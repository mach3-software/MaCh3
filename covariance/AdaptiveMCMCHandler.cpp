#include "covariance/AdaptiveMCMCHandler.h"

namespace adaptive_mcmc{

// ********************************************
AdaptiveMCMCHandler::AdaptiveMCMCHandler() {
// ********************************************
  start_adaptive_throw  = 0;
  start_adaptive_update = 0;
  end_adaptive_update   = 1;
  adaptive_update_step  = 1000;
  total_steps = 0;

  par_means = {};
  adaptive_covariance = nullptr;
}

// ********************************************
AdaptiveMCMCHandler::~AdaptiveMCMCHandler() {
// ********************************************
  if(adaptive_covariance != nullptr) {
    delete adaptive_covariance;
  }
}

// ********************************************
bool AdaptiveMCMCHandler::InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str, const int Npars) {
// ********************************************
  /*
   * HW: Idea is that adaption can simply read the YAML config
   * Options :
   *         External Info:
   * UseExternalMatrix [bool]     :    Use an external matrix
   * ExternalMatrixFileName [str] :    Name of file containing external info
   * ExternalMatrixName [str]     :    Name of external Matrix
   * ExternalMeansName [str]      :    Name of external means vector [for updates]
   *
   *         General Info:
   * DoAdaption [bool]            :    Do we want to do adaption?
   * AdaptionStartThrow [int]     :    Step we start throwing adaptive matrix from
   * AdaptionEndUpdate [int]      :    Step we stop updating adaptive matrix
   * AdaptionStartUpdate [int]    :    Do we skip the first N steps?
   * AdaptionUpdateStep [int]     :    Number of steps between matrix updates
   * Adaption blocks [vector<vector<int>>] : Splits the throw matrix into several block matrices
   */

  // setAdaptionDefaults();
  if(!adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]) {
    MACH3LOG_WARN("Adaptive Settings not found for {}, this is fine if you don't want adaptive MCMC", matrix_name_str);
    return false;
  }

  // We"re going to grab this info from the YAML manager
  if(!GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["DoAdaption"], false)) {
    MACH3LOG_WARN("Not using adaption for {}", matrix_name_str);
    return false;
  }

  start_adaptive_throw  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionStartThrow"], 10);
  start_adaptive_update = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionStartUpdate"], 0);
  end_adaptive_update   = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionEndUpdate"], 10000);
  adaptive_update_step  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionUpdateStep"], 100);

  // We also want to check for "blocks" by default all parameters "know" about each other
  // but we can split the matrix into independent block matrices

  // We"ll set a dummy variable here
  auto matrix_blocks = GetFromManager<std::vector<std::vector<int>>>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["MatrixBlocks"], {{}});

  SetAdaptiveBlocks(matrix_blocks, Npars);
  return true;
}

// ********************************************
void AdaptiveMCMCHandler::CreateNewAdaptiveCovariance(const int Npars) {
// ********************************************
  adaptive_covariance = new TMatrixDSym(Npars);
  adaptive_covariance->Zero();
  par_means = std::vector<double>(Npars, 0);
}

// ********************************************
void AdaptiveMCMCHandler::SetAdaptiveBlocks(std::vector<std::vector<int>> block_indices, const int Npars) {
// ********************************************
  /*
   *   In order to adapt efficient we want to setup our throw matrix to be a serious of block-diagonal (ish) matrices
   *
   *   To do this we set sub-block in the config by parameter index. For example having
   *   [[0,4],[4, 6]] in your config will set up two blocks one with all indices 0<=i<4 and the other with 4<=i<6
   */
  // Set up block regions
  adapt_block_matrix_indices = std::vector<int>(Npars, 0);

  // Should also make a matrix of block sizes
  adapt_block_sizes = std::vector<int>(block_indices.size()+1, 0);
  adapt_block_sizes[0] = Npars;

  if(block_indices.size()==0 || block_indices[0].size()==0) return;

  // Now we loop over our blocks
  for(int iblock=0; iblock<int(block_indices.size()); iblock++){
    // Loop over blocks in the block
    for(int isubblock=0; isubblock<int(block_indices[iblock].size())-1; isubblock+=2){
      int block_lb = block_indices[iblock][isubblock];
      int block_ub = block_indices[iblock][isubblock+1];

      if(block_lb > Npars || block_ub > Npars){
        MACH3LOG_ERROR("Cannot set matrix block with edges {}, {} for matrix of size {}",
                       block_lb, block_ub, Npars);
        throw MaCh3Exception(__FILE__, __LINE__);;
      }
      for(int ipar = block_lb; ipar < block_ub; ipar++){
        adapt_block_matrix_indices[ipar] = iblock+1;
        adapt_block_sizes[iblock+1] += 1;
        adapt_block_sizes[0] -= 1;
      }
    }
  }
}

// ********************************************
//HW: Truly adaptive MCMC!
void AdaptiveMCMCHandler::SaveAdaptiveToFile(const TString& outFileName, const TString& systematicName){
// ********************************************
  TFile* outFile = new TFile(outFileName, "UPDATE");
  if(outFile->IsZombie()){
    MACH3LOG_ERROR("Couldn't find {}", outFileName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  TVectorD* outMeanVec = new TVectorD(int(par_means.size()));
  for(int i = 0; i < int(par_means.size()); i++){
    (*outMeanVec)(i) = par_means[i];
  }
  outFile->cd();
  adaptive_covariance->Write(systematicName+"_postfit_matrix");
  outMeanVec->Write(systematicName+"_mean_vec");
  outFile->Close();
  delete outMeanVec;
  delete outFile;
}

// ********************************************
// HW : I would like this to be less painful to use!
// First things first we need setters
void AdaptiveMCMCHandler::SetThrowMatrixFromFile(const std::string& matrix_file_name,
                                                 const std::string& matrix_name,
                                                 const std::string& means_name,
                                                 bool& use_adaptive,
                                                 const int Npars) {
// ********************************************
  // Lets you set the throw matrix externally
  // Open file
  std::unique_ptr<TFile>matrix_file(new TFile(matrix_file_name.c_str()));
  use_adaptive = true;

  if(matrix_file->IsZombie()){
    MACH3LOG_ERROR("Couldn't find {}", matrix_file_name);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Next we grab our matrix
  adaptive_covariance = static_cast<TMatrixDSym*>(matrix_file->Get(matrix_name.c_str()));
  if(!adaptive_covariance){
    MACH3LOG_ERROR("Couldn't find {} in {}", matrix_name, matrix_file_name);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Finally we grab the means vector
  TVectorD* means_vector = static_cast<TVectorD*>(matrix_file->Get(means_name.c_str()));

  // This is fine to not exist!
  if(means_vector){
    // Yay our vector exists! Let's loop and fill it
    // Should check this is done
    if(means_vector->GetNrows()){
      MACH3LOG_ERROR("External means vec size ({}) != matrix size ({})", means_vector->GetNrows(), Npars);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    par_means = std::vector<double>(Npars);
    for(int i = 0; i < Npars; i++){
      par_means[i] = (*means_vector)(i);
    }
    MACH3LOG_INFO("Found Means in External File, Will be able to adapt");
  }
  // Totally fine if it doesn't exist, we just can't do adaption
  else{
    // We don't need a means vector, set the adaption=false
    MACH3LOG_WARN("Cannot find means vector in {}, therefore I will not be able to adapt!", matrix_file_name);
    use_adaptive = false;
  }

  matrix_file->Close();
  MACH3LOG_INFO("Set up matrix from external file");
}

// ********************************************
void AdaptiveMCMCHandler::UpdateAdaptiveCovariance(const std::vector<double>& _fCurrVal, const int Npars) {
// ********************************************
  std::vector<double> par_means_prev = par_means;

  int steps_post_burn = total_steps - start_adaptive_update;

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int iRow = 0; iRow < Npars; iRow++) {
    par_means[iRow] = (_fCurrVal[iRow]+par_means[iRow]*steps_post_burn)/(steps_post_burn+1);
  }

  //Now we update the covariances using cov(x,y)=E(xy)-E(x)E(y)
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int irow = 0; irow < Npars; irow++){
    int block = adapt_block_matrix_indices[irow];
    // int scale_factor = 5.76/double(adapt_block_sizes[block]);
    for(int icol = 0; icol <= irow; icol++){
      double cov_val=0;
      // Not in the same blocks
      if(adapt_block_matrix_indices[icol] == block){
        // Calculate Covariance for block
        // https://projecteuclid.org/journals/bernoulli/volume-7/issue-2/An-adaptive-Metropolis-algorithm/bj/1080222083.full
        cov_val = (*adaptive_covariance)(irow, icol)*Npars/5.6644;
        cov_val += par_means_prev[irow]*par_means_prev[icol]; //First we remove the current means
        cov_val = (cov_val*steps_post_burn+_fCurrVal[irow]*_fCurrVal[icol])/(steps_post_burn+1); //Now get mean(iRow*iCol)
        cov_val -= par_means[icol]*par_means[irow];
        cov_val*=5.6644/Npars;
      }
      (*adaptive_covariance)(icol, irow) = cov_val;
      (*adaptive_covariance)(irow, icol) = cov_val;
    }
  }
}

// ********************************************
bool AdaptiveMCMCHandler::IndivStepScaleAdapt() {
// ********************************************
  if(total_steps == start_adaptive_throw) return true;
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::UpdateMatrixAdapt() {
// ********************************************
  if(total_steps >= start_adaptive_throw &&
    (total_steps - start_adaptive_throw)%adaptive_update_step == 0) return true;
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::SkipAdaption() {
// ********************************************
  if(total_steps > end_adaptive_update ||
    total_steps< start_adaptive_update) return true;
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::AdaptionUpdate() {
// ********************************************
  if(total_steps <= start_adaptive_throw) return true;
  else return false;
}

// ********************************************
void AdaptiveMCMCHandler::Print() {
// ********************************************
  MACH3LOG_INFO("Adaptive MCMC Info:");
  MACH3LOG_INFO("Throwing from New Matrix from Step : {}", start_adaptive_throw);
  MACH3LOG_INFO("Adaption Matrix Start Update       : {}", start_adaptive_update);
  MACH3LOG_INFO("Adaption Matrix Ending Updates     : {}", end_adaptive_update);
  MACH3LOG_INFO("Steps Between Updates              : {}", adaptive_update_step);
}

} //end adaptive_mcmc
