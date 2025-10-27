#include "covariance/AdaptiveMCMCHandler.h"

namespace adaptive_mcmc
{

  // ********************************************
  AdaptiveMCMCHandler::AdaptiveMCMCHandler()
  {
    // ********************************************
    start_adaptive_throw = 0;
    start_adaptive_update = 0;
    end_adaptive_update = 1;
    adaptive_update_step = 1000;
    total_steps = 0;
    par_means = {};
    adaptive_covariance = nullptr;

    /// Robbins-Monro adaption
    use_robbins_monro = false;
    total_rm_restarts = 0;

    target_acceptance = 0.234;
    acceptance_rate_batch_size = 10000;
  }

  // ********************************************
  AdaptiveMCMCHandler::~AdaptiveMCMCHandler()
  {
    // ********************************************
    if (adaptive_covariance != nullptr)
    {
      delete adaptive_covariance;
    }
  }

  // ********************************************
  bool AdaptiveMCMCHandler::InitFromConfig(const YAML::Node &adapt_manager, const std::string &matrix_name_str,
                                           const std::vector<double> *parameters, const std::vector<double> *fixed)
  {
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
    if (!adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str])
    {
      MACH3LOG_WARN("Adaptive Settings not found for {}, this is fine if you don't want adaptive MCMC", matrix_name_str);
      return false;
    }

    // We"re going to grab this info from the YAML manager
    if (!GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["DoAdaption"], false))
    {
      MACH3LOG_WARN("Not using adaption for {}", matrix_name_str);
      return false;
    }

    if (!CheckNodeExists(adapt_manager, "AdaptionOptions", "Settings", "OutputFileName"))
    {
      MACH3LOG_ERROR("No OutputFileName specified in AdaptionOptions::Settings into your config file");
      MACH3LOG_ERROR("This is required if you are using adaptive MCMC");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    start_adaptive_throw = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["StartThrow"], 10);
    start_adaptive_update = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["StartUpdate"], 0);
    end_adaptive_update = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["EndUpdate"], 10000);
    adaptive_update_step = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["UpdateStep"], 100);
    adaptive_save_n_iterations = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["SaveNIterations"], -1);
    output_file_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Settings"]["OutputFileName"], "");

    // Check for Robbins-Monro adaption
    use_robbins_monro = GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Settings"]["UseRobbinsMonro"], false);
    target_acceptance = GetFromManager<double>(adapt_manager["AdaptionOptions"]["Settings"]["TargetAcceptance"], 0.234);
    if (target_acceptance <= 0 || target_acceptance >= 1)
    {
      MACH3LOG_ERROR("Target acceptance must be in (0,1), got {}", target_acceptance);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    acceptance_rate_batch_size = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["AcceptanceRateBatchSize"], 10000);
    total_rm_restarts = GetFromManager<int>(adapt_manager["AdaptionOptions"]["Settings"]["TotalRobbinsMonroRestarts"], 0);
    n_rm_restarts = 0;

    prev_step_accepted = false;

    // We also want to check for "blocks" by default all parameters "know" about each other
    // but we can split the matrix into independent block matrices
    SetParams(parameters);
    SetFixed(fixed);

    /// HW: This is technically wrong, should be across all systematics but will be addressed in a later PR
    if (use_robbins_monro)
    {
      MACH3LOG_INFO("Using Robbins-Monro for adaptive step size tuning with target acceptance {}", target_acceptance);
      // Now we need some thing
      CalculateRobbinsMonroStepLength();
    }

    // We'll use the same start scale for Robbins-Monro as standard adaption
    adaption_scale = 2.38 * 2.38 / GetNumParams();

    // We"ll set a dummy variable here
    auto matrix_blocks = GetFromManager<std::vector<std::vector<int>>>(adapt_manager["AdaptionOptions"]["Covariance"][matrix_name_str]["MatrixBlocks"], {{}});

    SetAdaptiveBlocks(matrix_blocks);
    return true;
  }

  // ********************************************
  void AdaptiveMCMCHandler::CreateNewAdaptiveCovariance()
  {
    // ********************************************
    adaptive_covariance = new TMatrixDSym(GetNumParams());
    adaptive_covariance->Zero();
    par_means = std::vector<double>(GetNumParams(), 0);
  }

  // ********************************************
  void AdaptiveMCMCHandler::SetAdaptiveBlocks(const std::vector<std::vector<int>> &block_indices)
  {
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
    adapt_block_sizes = std::vector<int>(block_indices.size() + 1, 0);
    adapt_block_sizes[0] = GetNumParams();

    if (block_indices.size() == 0 || block_indices[0].size() == 0)
      return;

    int block_size = static_cast<int>(block_indices.size());
    // Now we loop over our blocks
    for (int iblock = 0; iblock < block_size; iblock++)
    {
      // Loop over blocks in the block
      int sub_block_size = static_cast<int>(block_indices.size() - 1);
      for (int isubblock = 0; isubblock < sub_block_size; isubblock += 2)
      {
        int block_lb = block_indices[iblock][isubblock];
        int block_ub = block_indices[iblock][isubblock + 1];

        if (block_lb > GetNumParams() || block_ub > GetNumParams())
        {
          MACH3LOG_ERROR("Cannot set matrix block with edges {}, {} for matrix of size {}",
                         block_lb, block_ub, GetNumParams());
          throw MaCh3Exception(__FILE__, __LINE__);
          ;
        }
        for (int ipar = block_lb; ipar < block_ub; ipar++)
        {
          adapt_block_matrix_indices[ipar] = iblock + 1;
          adapt_block_sizes[iblock + 1] += 1;
          adapt_block_sizes[0] -= 1;
        }
      }
    }
  }

  // ********************************************
  // HW: Truly adaptive MCMC!
  void AdaptiveMCMCHandler::SaveAdaptiveToFile(const std::string &outFileName,
                                               const std::string &systematicName, bool is_final)
  {
    // ********************************************
    // Skip saving if adaptive_save_n_iterations is negative,
    // unless this is the final iteration (is_final overrides the condition)
    if (adaptive_save_n_iterations < 0 && !is_final)
      return;
    if (is_final ||
        (adaptive_save_n_iterations > 0 &&
         ((total_steps - start_adaptive_throw) / adaptive_update_step) % adaptive_save_n_iterations == 0))
    {

      TFile *outFile = new TFile(outFileName.c_str(), "UPDATE");
      if (outFile->IsZombie())
      {
        MACH3LOG_ERROR("Couldn't find {}", outFileName);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      TVectorD *outMeanVec = new TVectorD(int(par_means.size()));
      for (int i = 0; i < int(par_means.size()); i++)
      {
        (*outMeanVec)(i) = par_means[i];
      }

      std::string adaptive_cov_name = systematicName + "_posfit_matrix";
      std::string mean_vec_name = systematicName + "_mean_vec";
      if (!is_final)
      {
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
      auto adaptive_to_save = std::make_unique<TMatrixDSym>(*adaptive_covariance);
      // Multiply by scale
      // (*adaptive_to_save) *= GetAdaptionScale()/2.0;
      adaptive_to_save->Write(adaptive_cov_name.c_str());
      outMeanVec->Write(mean_vec_name.c_str());
      outFile->Close();
      delete outMeanVec;
      delete outFile;
    }
  }

  // ********************************************
  // HW : I would like this to be less painful to use!
  // First things first we need setters
  void AdaptiveMCMCHandler::SetThrowMatrixFromFile(const std::string &matrix_file_name,
                                                   const std::string &matrix_name,
                                                   const std::string &means_name,
                                                   bool &use_adaptive)
  {
    // ********************************************
    // Lets you set the throw matrix externally
    // Open file
    auto matrix_file = std::make_unique<TFile>(matrix_file_name.c_str());
    use_adaptive = true;

    if (matrix_file->IsZombie())
    {
      MACH3LOG_ERROR("Couldn't find {}", matrix_file_name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Next we grab our matrix
    adaptive_covariance = static_cast<TMatrixDSym *>(matrix_file->Get(matrix_name.c_str()));
    if (!adaptive_covariance)
    {
      MACH3LOG_ERROR("Couldn't find {} in {}", matrix_name, matrix_file_name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Finally we grab the means vector
    TVectorD *means_vector = static_cast<TVectorD *>(matrix_file->Get(means_name.c_str()));

    // This is fine to not exist!
    if (means_vector)
    {
      // Yay our vector exists! Let's loop and fill it
      // Should check this is done
      if (means_vector->GetNrows())
      {
        MACH3LOG_ERROR("External means vec size ({}) != matrix size ({})", means_vector->GetNrows(), GetNumParams());
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      par_means = std::vector<double>(GetNumParams());
      for (int i = 0; i < GetNumParams(); i++)
      {
        par_means[i] = (*means_vector)(i);
      }
      MACH3LOG_INFO("Found Means in External File, Will be able to adapt");
    }
    // Totally fine if it doesn't exist, we just can't do adaption
    else
    {
      // We don't need a means vector, set the adaption=false
      MACH3LOG_WARN("Cannot find means vector in {}, therefore I will not be able to adapt!", matrix_file_name);
      use_adaptive = false;
    }

    matrix_file->Close();
    MACH3LOG_INFO("Set up matrix from external file");
  }

  // ********************************************
  void AdaptiveMCMCHandler::UpdateAdaptiveCovariance()
  {
    // ********************************************
    std::vector<double> par_means_prev = par_means;
    int steps_post_burn = total_steps - start_adaptive_update;

    prev_step_accepted = false;

    // Step 1: Update means and compute deviations
    for (int i = 0; i < GetNumParams(); ++i)
    {
      if (IsFixed(i))
        continue;

      par_means[i] = (CurrVal(i) + par_means_prev[i] * steps_post_burn) / (steps_post_burn + 1);
      // Left over from cyclic means
    }

// Step 2: Update covariance
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < GetNumParams(); ++i)
    {
      if (IsFixed(i))
      {
        (*adaptive_covariance)(i, i) = 1.0;
        continue;
      }

      int block_i = adapt_block_matrix_indices[i];

      for (int j = 0; j <= i; ++j)
      {
        if (IsFixed(j) || adapt_block_matrix_indices[j] != block_i)
        {
          (*adaptive_covariance)(i, j) = 0.0;
          (*adaptive_covariance)(j, i) = 0.0;
          continue;
        }

        double cov_prev = (*adaptive_covariance)(i, j);
        double cov_updated = 0.0;

        if (steps_post_burn > 0)
        {
          // Haario-style update
          double cov_t = cov_prev * (steps_post_burn - 1) / steps_post_burn;
          double prev_means_t = steps_post_burn * par_means_prev[i] * par_means_prev[j];
          double curr_means_t = (steps_post_burn + 1) * par_means[i] * par_means[j];
          double curr_step_t = CurrVal(i) * CurrVal(j);

          cov_updated = cov_t + (prev_means_t - curr_means_t + curr_step_t) / steps_post_burn;
        }

        (*adaptive_covariance)(i, j) = cov_updated;
        (*adaptive_covariance)(j, i) = cov_updated;
      }
    }
  }

  // ********************************************
  bool AdaptiveMCMCHandler::IndivStepScaleAdapt() const
  {
    // ********************************************
    if (total_steps == start_adaptive_throw)
      return true;
    else
      return false;
  }

  // ********************************************
  bool AdaptiveMCMCHandler::UpdateMatrixAdapt()
  {
    // ********************************************
    if (total_steps >= start_adaptive_throw &&
        // Check whether the number of steps is divisible by the adaptive update step
        // e.g. if adaptive_update_step = 1000 and (total_step - start_adpative_throw) is 5000 then this is true
        (total_steps - start_adaptive_throw) % adaptive_update_step == 0)
    {
      return true;
    }
    else
      return false;
  }

  // ********************************************
  bool AdaptiveMCMCHandler::SkipAdaption() const
  {
    // ********************************************
    if (total_steps > end_adaptive_update ||
        total_steps < start_adaptive_update)
      return true;
    else
      return false;
  }

  // ********************************************
  bool AdaptiveMCMCHandler::AdaptionUpdate() const
  {
    // ********************************************
    if (total_steps <= start_adaptive_throw)
      return true;
    else
      return false;
  }

  // ********************************************
  void AdaptiveMCMCHandler::Print() const
  {
    // ********************************************
    MACH3LOG_INFO("Adaptive MCMC Info:");
    MACH3LOG_INFO("Throwing from New Matrix from Step : {}", start_adaptive_throw);
    MACH3LOG_INFO("Adaption Matrix Start Update       : {}", start_adaptive_update);
    MACH3LOG_INFO("Adaption Matrix Ending Updates     : {}", end_adaptive_update);
    MACH3LOG_INFO("Steps Between Updates              : {}", adaptive_update_step);
    MACH3LOG_INFO("Saving matrices to file            : {}", output_file_name);
    MACH3LOG_INFO("Will only save every {} iterations", adaptive_save_n_iterations);
    if (use_robbins_monro)
    {
      MACH3LOG_INFO("Using Robbins-Monro for adaptive step size tuning with target acceptance {}", target_acceptance);
    }
  }

  // ********************************************
  void AdaptiveMCMCHandler::CheckMatrixValidityForAdaption(const TMatrixDSym *covMatrix) const
  {
    // ********************************************
    int n = covMatrix->GetNrows();
    std::vector<std::vector<double>> corrMatrix(n, std::vector<double>(n));

    // Extract standard deviations from the diagonal
    std::vector<double> stdDevs(n);
    for (int i = 0; i < n; ++i)
    {
      stdDevs[i] = std::sqrt((*covMatrix)(i, i));
    }

    // Compute correlation for each element
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        double cov = (*covMatrix)(i, j);
        double corr = cov / (stdDevs[i] * stdDevs[j]);
        corrMatrix[i][j] = corr;
      }
    }

    // KS: These numbers are completely arbitrary
    constexpr double NumberOfParametersThreshold = 5;
    constexpr double CorrelationThreshold = 0.5;

    bool FlagMessage = false;
    // For each parameter, count how many others it is correlated with
    for (int i = 0; i < n; ++i)
    {
      if (IsFixed(i))
      {
        continue;
      }
      int block_i = adapt_block_matrix_indices[i];
      int correlatedCount = 0;
      std::vector<int> correlatedWith;

      for (int j = 0; j < n; ++j)
      {
        if (i == j || IsFixed(j) || adapt_block_matrix_indices[j] != block_i)
        {
          continue;
        }
        if (corrMatrix[i][j] > CorrelationThreshold)
        {
          correlatedCount++;
          correlatedWith.push_back(j);
        }
      }

      if (correlatedCount > NumberOfParametersThreshold)
      {
        MACH3LOG_WARN("Parameter {} is highly correlated with {} other parameters (above threshold {}).", i, correlatedCount, CorrelationThreshold);
        MACH3LOG_WARN("Correlated with parameters: {}.", fmt::join(correlatedWith, ", "));
        FlagMessage = true;
      }
    }
    if (FlagMessage)
    {
      MACH3LOG_WARN("Found highly correlated block, adaptive may not work as intended, proceed with caution");
      MACH3LOG_WARN("You can consider defying so called Adaption Block to not update this part of matrix");
    }
  }

  // ********************************************
  double AdaptiveMCMCHandler::CurrVal(const int par_index) const
  {
    // ********************************************
    /// HW Implemented as its own method to allow for
    /// different behaviour in the future
    return (*_fCurrVal)[par_index];
  }

  // ********************************************
  void AdaptiveMCMCHandler::CalculateRobbinsMonroStepLength()
  {
    // ********************************************
    /*
      Obtains the constant "step scale" for Robbins-Monro adaption
      for simplicity and because it has the degrees of freedom "baked in"
      it will be same across ALL blocks.
    */

    /// Firstly we need to calculate the alpha value, this is in some sense "optimal"
    double alpha = -ROOT::Math::normal_quantile(target_acceptance / 2, 1.0);

    int non_fixed_pars = GetNumParams() - GetNFixed();

    /// Now we can calculate the scale factor
    c_robbins_monro = (1 - 1 / non_fixed_pars) * std::sqrt(TMath::Pi() * 2 * std::exp(alpha * alpha));
    c_robbins_monro += 1 / (non_fixed_pars * target_acceptance * (1 - target_acceptance));

    MACH3LOG_INFO("Robbins-Monro scale factor set to {}", c_robbins_monro);
  }

  // ********************************************
  void AdaptiveMCMCHandler::UpdateRobbinsMonroScale()
  {
    // ********************************************
    /*
    Update the scale factor using Robbins-Monro.
    TLDR: If acceptance rate is too high, scale factor goes up, if too low goes down
          will pretty rapidly converge to the right value.
    */

    if ((total_steps > 0) && (total_steps % acceptance_rate_batch_size == 0))
    {
      n_rm_restarts++;
    }

    /// This allows to move towards a STABLE number of steps in the batch
    int batch_step;
    if (n_rm_restarts > total_rm_restarts)
    {
      batch_step = total_steps;
    }
    else
    {
      batch_step = total_steps % acceptance_rate_batch_size;
    }

    int non_fixed_pars = GetNumParams() - GetNFixed();

    double root_scale = std::sqrt(adaption_scale);

    /// Update the scale factor for Robbins-Monro adaption [200 is arbitrary for early stability but motivated by paper]
    double scale_factor = root_scale * c_robbins_monro / std::max(200.0, static_cast<double>(batch_step) / non_fixed_pars);

    /// Now we either increase or decrease the scale
    if (prev_step_accepted)
    {
      root_scale += (1 - target_acceptance) * scale_factor;
    }
    else
    {
      root_scale -= target_acceptance * scale_factor;
    }
    adaption_scale = root_scale * root_scale;
  }

} // end adaptive_mcmc
