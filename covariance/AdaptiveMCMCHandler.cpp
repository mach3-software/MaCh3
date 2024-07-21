#include "covariance/AdaptiveMCMCHandler.h"

namespace adaptive_mcmc{

// ********************************************
AdaptiveMCMCHandler::AdaptiveMCMCHandler() {
// ********************************************
  start_adaptive_throw  = 0;
  start_adaptive_update = 0;
  end_adaptive_update   = 1;
  adaptive_update_step  = 1000;

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
void AdaptiveMCMCHandler::InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str, const int Npars) {
// ********************************************
  //  setAdaptionDefaults();
  if(GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str], "")==""){
    MACH3LOG_WARN("Adaptive Settings not found for {}, this is fine if you don't want adaptive MCMC", matrix_name_str);
    return;
  }

  // We"re going to grab this info from the YAML manager
  if(GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["DoAdaption"], false)) {
    MACH3LOG_WARN("Not using adaption for {}", matrix_name_str);
    return;
  }

  start_adaptive_throw  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionStartThrow"], 10);
  start_adaptive_update = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionStartUpdate"], 0);
  end_adaptive_update   = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionEndUpdate"], 10000);
  adaptive_update_step  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AdaptionUpdateStep"], 100);


  // We also want to check for "blocks" by default all parameters "know" about each other
  // but we can split the matrix into independent block matrices

  // We"ll set a dummy variable here
  auto matrix_blocks = GetFromManager<std::vector<std::vector<int>>>(adapt_manager["AdaptionOptions"]["Settings"][matrix_name_str]["AdaptionUpdateStep"], {{}});

  SetAdaptiveBlocks(matrix_blocks, Npars);
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
  adapt_block_sizes = std::vector<int>((int)block_indices.size()+1, 0);
  adapt_block_sizes[0] = Npars;

  if(block_indices.size()==0 || block_indices[0].size()==0) return;

  // Now we loop over our blocks
  for(int iblock=0; iblock<(int)block_indices.size(); iblock++){
    // Loop over blocks in the block
    for(int isubblock=0; isubblock<(int)block_indices[iblock].size()-1; isubblock+=2){
      int block_lb = block_indices[iblock][isubblock];
      int block_ub = block_indices[iblock][isubblock+1];

      //std::cout<<block_lb<<" "<<block_ub<<std::endl;
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
    throw;
  }
  TVectorD* outMeanVec = new TVectorD((int)par_means.size());
  for(int i = 0; i < (int)par_means.size(); i++){
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
void AdaptiveMCMCHandler::Print() {
// ********************************************
  MACH3LOG_INFO("Adaptive MCMC Info:");
  MACH3LOG_INFO("Throwing from New Matrix from Step : {}", start_adaptive_throw);
  MACH3LOG_INFO("Adaption Matrix Start Update       : {}", start_adaptive_update);
  MACH3LOG_INFO("Adaption Matrix Ending Updates     : {}", end_adaptive_update);
  MACH3LOG_INFO("Steps Between Updates              : {}", adaptive_update_step);
}


} //end adaptive_mcmc
