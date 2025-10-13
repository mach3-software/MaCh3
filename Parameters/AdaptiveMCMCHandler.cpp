#include "Parameters/AdaptiveMCMCHandler.h"

namespace adaptive_mcmc{

/*
AMCMC Mode
*/
AdaptiveMode::AdaptiveMode(const std::vector<double>& point)
{
  n_steps = 0;
  means = point;
  variances.resize(point.size(), 0.0);

}

void AdaptiveMode::Update(const std::vector<double>& point){
  if(point.size() != means.size()){
    MACH3LOG_ERROR("Point and means vectors must be the same size");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  n_steps++;
  for(size_t i = 0; i < means.size(); ++i){
    double delta = point[i] - means[i];
    means[i] += delta / static_cast<double>(n_steps);
    variances[i] += delta * (point[i] - means[i]);
  }
}

double AdaptiveMode::NSigmaFromMode(std::vector<double> point) const{
  if(point.size() != means.size()){
    MACH3LOG_ERROR("Point and means vectors must be the same size");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  /// Need an "okay" estimate for low variance
  /// With low statistics, use a simple Euclidean distance with adaptive scaling
  if(n_steps < 10){
    double euclidean_dist = 0.0;
    for(size_t i = 0; i < means.size(); ++i){
      double delta = point[i] - means[i];
      euclidean_dist += delta * delta;
    }
    // Return normalized Euclidean distance scaled by dimension
    return std::sqrt(euclidean_dist / static_cast<double>(means.size()));
  }

  // With sufficient statistics, use proper Mahalanobis-like distance
  double nSigma = 0.0;
  for(size_t i = 0; i < means.size(); ++i){
    double var = variances[i] / (n_steps - 1);
    if(var > 1e-10){  // Numerical safety
      nSigma += ((point[i] - means[i]) * (point[i] - means[i])) / var;
    } else {
      // If variance is essentially zero, use the squared deviation directly
      double delta = point[i] - means[i];
      nSigma += delta * delta * 1e10;  // Penalize deviation from a precise mean
    }
  }
  return std::sqrt(nSigma / static_cast<double>(means.size()));  // Normalize by dimension
}


/*
Adaptive MCMC code
*/

// ********************************************
AdaptiveMCMCHandler::AdaptiveMCMCHandler() {
// ********************************************
  start_adaptive_throw  = 0;
  start_adaptive_update = 0;
  end_adaptive_update   = 1;
  adaptive_update_step  = 1000;
  total_steps = 0;

  par_means = {};
  mode_information.clear();  // Use clear() instead of assignment for vector of unique_ptr
  adaptive_covariance = nullptr;
  
  track_modes = false;  // Disabled by default
}

// ********************************************
AdaptiveMCMCHandler::~AdaptiveMCMCHandler() {
// ********************************************
  if(adaptive_covariance != nullptr) {
    delete adaptive_covariance;
  }
}

// ********************************************
bool AdaptiveMCMCHandler::InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str,
                                        const std::vector<double>* parameters, const std::vector<double>* fixed) {
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
   * AdaptionSaveNIterations [int]:    You don't have to save every adaptive stage so decide how often you want to save
   * Adaption blocks [vector<vector<int>>] : Splits the throw matrix into several block matrices
   * OuputFileName [std::string]  :    Name of the file that the adaptive matrices will be saved into
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

  if(!CheckNodeExists(adapt_manager, "AdaptionOptions", "Settings", "OutputFileName")) {
    MACH3LOG_ERROR("No OutputFileName specified in AdaptionOptions::Settings into your config file");
    MACH3LOG_ERROR("This is required if you are using adaptive MCMC");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  start_adaptive_throw  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["StartThrow"], 10);
  start_adaptive_update = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["StartUpdate"], 0);
  end_adaptive_update   = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["EndUpdate"], 10000);
  adaptive_update_step  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["UpdateStep"], 100);
  adaptive_save_n_iterations  = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["SaveNIterations"], -1);
  output_file_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Settings"]["OutputFileName"], "");

  #ifdef MPIENABLED
  /// Means we can do AMCMC with MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  output_file_name = output_file_name.substr(0, output_file_name.find_last_of('.')) + "_rank" + std::to_string(mpi_rank) + output_file_name.substr(output_file_name.find_last_of('.'));
  #endif

  // We also want to check for "blocks" by default all parameters "know" about each other
  // but we can split the matrix into independent block matrices
  SetParams(parameters);
  SetFixed(fixed);

  /// HW: This is technically wrong, should be across all systematics but will be addressed in a later PR
  adaption_scale = 2.38*2.38/GetNumParams(); 
  
  // Check for mode tracking configuration
  track_modes = GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Settings"]["TrackModes"], true);
  
  // We"ll set a dummy variable here
  auto matrix_blocks = GetFromManager<std::vector<std::vector<int>>>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["MatrixBlocks"], {{}});

  SetAdaptiveBlocks(matrix_blocks);
  return true;
}

// ********************************************
void AdaptiveMCMCHandler::CreateNewAdaptiveCovariance() {
// ********************************************
  adaptive_covariance = new TMatrixDSym(GetNumParams());
  adaptive_covariance->Zero();
  par_means = std::vector<double>(GetNumParams(), 0);
}

// ********************************************
void AdaptiveMCMCHandler::SetAdaptiveBlocks(const std::vector<std::vector<int>>& block_indices) {
// ********************************************
  /*
   *   In order to adapt efficient we want to setup our throw matrix to be a serious of block-diagonal (ish) matrices
   *
   *   To do this we set sub-block in the config by parameter index. For example having
   *   [[0,4],[4, 6]] in your config will set up two blocks one with all indices 0<=i<4 and the other with 4<=i<6
   */
  // Set up block regions
  adapt_block_matrix_indices = std::vector<int>(GetNumParams(), 0);

  // Should also make a matrix of block sizes
  adapt_block_sizes = std::vector<int>(block_indices.size()+1, 0);
  adapt_block_sizes[0] = GetNumParams();

  if(block_indices.size()==0 || block_indices[0].size()==0) return;

  int block_size = static_cast<int>(block_indices.size());
  // Now we loop over our blocks
  for(int iblock=0; iblock < block_size; iblock++){
    // Loop over blocks in the block
    int sub_block_size = static_cast<int>(block_indices.size()-1);
    for(int isubblock=0; isubblock < sub_block_size ; isubblock+=2){
      int block_lb = block_indices[iblock][isubblock];
      int block_ub = block_indices[iblock][isubblock+1];

      if(block_lb > GetNumParams() || block_ub > GetNumParams()){
        MACH3LOG_ERROR("Cannot set matrix block with edges {}, {} for matrix of size {}",
                       block_lb, block_ub, GetNumParams());
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
void AdaptiveMCMCHandler::SaveAdaptiveToFile(const std::string &outFileName,
                                             const std::string &systematicName, bool is_final) {
// ********************************************
  // Skip saving if adaptive_save_n_iterations is negative,
  // unless this is the final iteration (is_final overrides the condition)
  if (adaptive_save_n_iterations < 0 && !is_final) return;
  if (is_final ||
    (adaptive_save_n_iterations > 0 &&
    ((total_steps - start_adaptive_throw) / adaptive_update_step) % adaptive_save_n_iterations == 0)) {

    TFile *outFile = new TFile(outFileName.c_str(), "UPDATE");
    if (outFile->IsZombie()) {
      MACH3LOG_ERROR("Couldn't find {}", outFileName);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    TVectorD *outMeanVec = new TVectorD(int(par_means.size()));
    for (int i = 0; i < int(par_means.size()); i++) {
      (*outMeanVec)(i) = par_means[i];
    }

    std::string adaptive_cov_name = systematicName + "_posfit_matrix";
    std::string mean_vec_name = systematicName + "_mean_vec";
    if (!is_final) {
      // Some string to make the name of the saved adaptive matrix clear
      std::string total_steps_str =
          std::to_string(total_steps);
      std::string syst_name_str = systematicName;
      adaptive_cov_name =
          total_steps_str + '_' + syst_name_str + std::string("_throw_matrix");
      mean_vec_name =
          total_steps_str + '_' + syst_name_str + std::string("_mean_vec");
    }

    outFile->cd();
    adaptive_covariance->Write(adaptive_cov_name.c_str());
    outMeanVec->Write(mean_vec_name.c_str());
    outFile->Close();
    delete outMeanVec;
    delete outFile;
  }
}

// ********************************************
// HW : I would like this to be less painful to use!
// First things first we need setters
void AdaptiveMCMCHandler::SetThrowMatrixFromFile(const std::string& matrix_file_name,
                                                 const std::string& matrix_name,
                                                 const std::string& means_name,
                                                 bool& use_adaptive) {
// ********************************************
  // Lets you set the throw matrix externally
  // Open file
  auto matrix_file = std::make_unique<TFile>(matrix_file_name.c_str());
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
      MACH3LOG_ERROR("External means vec size ({}) != matrix size ({})", means_vector->GetNrows(), GetNumParams());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    par_means = std::vector<double>(GetNumParams());
    for(int i = 0; i < GetNumParams(); i++){
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
void AdaptiveMCMCHandler::UpdateAdaptiveCovariance() {
// ********************************************
  std::vector<double> par_means_prev = par_means;
  int steps_post_burn = total_steps - start_adaptive_update;

  // Update mode information if tracking is enabled
  if(track_modes) {
    UpdateModes(*_fCurrVal);
  }

  // Mode correction: "revert" the current point to the primary mode
  // This helps build a unimodal covariance even with mode flips
  std::vector<double> corrected_point(GetNumParams());
  bool use_mode_correction = track_modes && !mode_information.empty() && mode_information.size() > 1;
  
  if(use_mode_correction) {
    corrected_point = MoveToInitMode(*_fCurrVal);
  } else {
    corrected_point = *_fCurrVal;
  }

  // Step 1: Update means and compute deviations
  for (int i = 0; i < GetNumParams(); ++i) {
    if (IsFixed(i)) continue;

    // Use corrected point for mean calculation
    double point_val = corrected_point[i];
    par_means[i] = (point_val + par_means_prev[i]*steps_post_burn)/(steps_post_burn+1);
  }

  // Step 2: Update covariance
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < GetNumParams(); ++i) {
    if (IsFixed(i)) {
      (*adaptive_covariance)(i, i) = 1.0;
      continue;
    }

    int block_i = adapt_block_matrix_indices[i];

    for (int j = 0; j <= i; ++j) {
      if (IsFixed(j) || adapt_block_matrix_indices[j] != block_i) {
        (*adaptive_covariance)(i, j) = 0.0;
        (*adaptive_covariance)(j, i) = 0.0;
        continue;
      }

      double cov_prev = (*adaptive_covariance)(i, j);
      double cov_updated = 0.0;

      if (steps_post_burn > 0) {
        // Haario-style update with mode-corrected points
        double cov_t = cov_prev * (steps_post_burn - 1) / steps_post_burn;
        double prev_means_t = steps_post_burn * par_means_prev[i] * par_means_prev[j];
        double curr_means_t = (steps_post_burn + 1) * par_means[i] * par_means[j];
        double curr_step_t = corrected_point[i] * corrected_point[j];

        cov_updated = cov_t + adaption_scale * (prev_means_t - curr_means_t + curr_step_t) / steps_post_burn;
      }

      (*adaptive_covariance)(i, j) = cov_updated;
      (*adaptive_covariance)(j, i) = cov_updated;
    }
  }
}

// ********************************************
bool AdaptiveMCMCHandler::IndivStepScaleAdapt() const {
// ********************************************
  if(total_steps == start_adaptive_throw) return true;
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::UpdateMatrixAdapt() {
// ********************************************
  if(total_steps >= start_adaptive_throw &&
    // Check whether the number of steps is divisible by the adaptive update step
    // e.g. if adaptive_update_step = 1000 and (total_step - start_adpative_throw) is 5000 then this is true
    (total_steps - start_adaptive_throw)%adaptive_update_step == 0) {
    return true;
  } 
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::SkipAdaption() const {
// ********************************************
  if(total_steps > end_adaptive_update ||
    total_steps< start_adaptive_update) return true;
  else return false;
}

// ********************************************
bool AdaptiveMCMCHandler::AdaptionUpdate() const {
// ********************************************
  if(total_steps <= start_adaptive_throw) return true;
  else return false;
}

// ********************************************
void AdaptiveMCMCHandler::Print() const {
// ********************************************
  MACH3LOG_INFO("Adaptive MCMC Info:");
  MACH3LOG_INFO("Throwing from New Matrix from Step : {}", start_adaptive_throw);
  MACH3LOG_INFO("Adaption Matrix Start Update       : {}", start_adaptive_update);
  MACH3LOG_INFO("Adaption Matrix Ending Updates     : {}", end_adaptive_update);
  MACH3LOG_INFO("Steps Between Updates              : {}", adaptive_update_step);
  MACH3LOG_INFO("Saving matrices to file            : {}", output_file_name);
  MACH3LOG_INFO("Will only save every {} iterations", adaptive_save_n_iterations);
  MACH3LOG_INFO("Mode tracking enabled              : {}", track_modes);
  if(track_modes) {
    MACH3LOG_INFO("Number of detected modes           : {}", mode_information.size());
  }
}

// ********************************************
void AdaptiveMCMCHandler::PrintModeInfo() const {
// ********************************************
  if(!track_modes) {
    MACH3LOG_INFO("Mode tracking is disabled");
    return;
  }
  
  MACH3LOG_INFO("=== Mode Information ===");
  MACH3LOG_INFO("Total number of modes detected: {}", mode_information.size());
  
  for(size_t i = 0; i < mode_information.size(); ++i) {
    MACH3LOG_INFO("Mode {}: {} steps accumulated", i, mode_information[i]->GetNSteps());
    
    const auto& means = mode_information[i]->GetMeans();
    const auto& vars = mode_information[i]->GetVariances();
    
    // Print first few parameters for each mode
    size_t n_print = std::min(size_t(5), means.size());
    MACH3LOG_INFO("  First {} parameters:", n_print);
    for(size_t j = 0; j < n_print; ++j) {
      double std_dev = (mode_information[i]->GetNSteps() > 1) ? 
                       std::sqrt(vars[j] / (mode_information[i]->GetNSteps() - 1)) : 0.0;
      MACH3LOG_INFO("    Par {}: mean = {:.4f}, std = {:.4f}", j, means[j], std_dev);
    }
  }
  MACH3LOG_INFO("========================");
}

double AdaptiveMCMCHandler::CurrVal(const int par_index) const {
  /// HW Implemented as its own method to allow for
  /// different behaviour in the future
  return (*_fCurrVal)[par_index];
}

void AdaptiveMCMCHandler::UpdateModes(const std::vector<double>& point) {
  // Adaptive threshold: start large and decrease as we gain statistics
  // With low stats, we're more lenient about creating new modes
  double mode_threshold = 3.0;  // N-sigma threshold
  
  // Make threshold more stringent as we accumulate data
  if(total_steps > 1000) {
    mode_threshold = 2.5;
  }
  if(total_steps > 5000) {
    mode_threshold = 2.0;
  }
  
  // First point: create initial mode
  if(mode_information.empty()){
    std::unique_ptr<AdaptiveMode> new_mode = std::make_unique<AdaptiveMode>(point);
    mode_information.push_back(std::move(new_mode));
    MACH3LOG_DEBUG("Created initial mode 0");
    return;
  }

  // Find closest existing mode
  int closest_mode_idx = GetClosestMode(point);
  double distance_to_closest = mode_information[closest_mode_idx]->NSigmaFromMode(point);
  
  MACH3LOG_DEBUG("Distance to closest mode {}: {} sigma", closest_mode_idx, distance_to_closest);
  
  // If close enough to existing mode, update it
  if(distance_to_closest < mode_threshold) {
    mode_information[closest_mode_idx]->Update(point);
    MACH3LOG_DEBUG("Updated existing mode {}", closest_mode_idx);
  } 
  // If far from all existing modes, create new mode
  else {
    // Limit the number of modes to prevent over-segmentation with noise
    const size_t max_modes = 20;  // Configurable limit
    
    if(mode_information.size() < max_modes) {
      std::unique_ptr<AdaptiveMode> new_mode = std::make_unique<AdaptiveMode>(point);
      mode_information.push_back(std::move(new_mode));
      MACH3LOG_INFO("Created new mode {} at step {} ({} total modes)", 
                    mode_information.size()-1, total_steps, mode_information.size());
    } else {
      // If at max modes, force update of closest mode
      mode_information[closest_mode_idx]->Update(point);
      MACH3LOG_DEBUG("At max modes ({}), forcing update of mode {}", max_modes, closest_mode_idx);
    }
  }
  
  // Periodically merge modes that have become too similar
  if(total_steps % 1000 == 0 && mode_information.size() > 1) {
    MergeSimilarModes();
  }
}

int AdaptiveMCMCHandler::GetClosestMode(const std::vector<double>& point) const {
  if(mode_information.empty()){
    MACH3LOG_ERROR("Cannot get closest mode from empty mode list");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  int closest_idx = 0;
  double min_distance = mode_information[0]->NSigmaFromMode(point);
  
  for(size_t i = 1; i < mode_information.size(); ++i) {
    double distance = mode_information[i]->NSigmaFromMode(point);
    if(distance < min_distance) {
      min_distance = distance;
      closest_idx = static_cast<int>(i);
    }
  }
  
  return closest_idx;
}

void AdaptiveMCMCHandler::MergeSimilarModes() {
  // Check all pairs of modes and merge if they're too similar
  const double merge_threshold = 3.0;  // Modes within 3 sigma get merged
  
  for(size_t i = 0; i < mode_information.size(); ++i) {
    for(size_t j = i+1; j < mode_information.size(); ++j) {
      // Get means from both modes
      std::vector<double> means_i = mode_information[i]->GetMeans();
      double distance = mode_information[j]->NSigmaFromMode(means_i);
      
      if(distance < merge_threshold) {
        MACH3LOG_INFO("Merging mode {} into mode {} (distance: {} sigma)", j, i, distance);
        // Simply remove the redundant mode (the other will capture future points)
        mode_information.erase(mode_information.begin() + j);
        // After erasing, adjust indices
        --j;
      }
    }
  }
}

std::vector<double> AdaptiveMCMCHandler::MoveToInitMode(const std::vector<double>& point) const {
  // Find which mode this point belongs to
  int mode_idx = GetClosestMode(point);
  
  // Get the displacement from the mode center
  std::vector<double> mode_means = mode_information[mode_idx]->GetMeans();
  std::vector<double> corrected_point = point;
  
  if(mode_information.empty() || mode_idx != 0) {
    // If not in the primary (first) mode, reflect the point
    std::vector<double> primary_means = mode_information[0]->GetMeans();
    
    for(size_t i = 0; i < point.size(); ++i) {
      // Calculate displacement from current mode
      double displacement = point[i] - mode_means[i];
      // Apply same displacement from primary mode
      corrected_point[i] = primary_means[i] + displacement;
    }
    
    MACH3LOG_DEBUG("Corrected point from mode {} to mode 0", mode_idx);
  }
  
  return corrected_point;
}

} //end adaptive_mcmc
