#pragma once

// MaCh3 Includes
#include "manager/MaCh3Logger.h"
#include "manager/manager.h"
#include "covariance/CovarianceUtils.h"

namespace adaptive_mcmc{

/// @brief Contains information about adaptive covariance matrix
/// @see An adaptive Metropolis algorithm, H.Haario et al., 2001 for more info!
///@details struct encapsulating all adaptive MCMC information
class AdaptiveMCMCHandler{
 public:

  /// @brief Constructor
  AdaptiveMCMCHandler();

  /// @brief Destructor
  ~AdaptiveMCMCHandler();

  /// @brief Print all class members
  void Print();

  /// @brief Read initial values from config file
  /// @param adapt_manager Config file from which we update matrix
  void InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str, const int Npars);

  /// @brief If we don't have a covariance matrix to start from for adaptive tune we need to make one!
  void CreateNewAdaptiveCovariance(const int Npars);

  /// @brief HW: sets adaptive block matrix
  /// @param block_indices Values for sub-matrix blocks
  void SetAdaptiveBlocks(std::vector<std::vector<int>> block_indices, const int Npars);

  /// @brief HW: Save adaptive throw matrix to file
  void SaveAdaptiveToFile(const TString& outFileName, const TString& systematicName);

  /// @brief sets throw matrix from a file
  /// @param matrix_file_name name of file matrix lives in
  /// @param matrix_name name of matrix in file
  /// @param means_name name of means vec in file
  void SetThrowMatrixFromFile(const std::string& matrix_file_name,
                              const std::string& matrix_name,
                              const std::string& means_name,
                              bool& use_adaptive,
                              const int Npars);




  /// @brief Method to update adaptive MCMC
  /// @see https://projecteuclid.org/journals/bernoulli/volume-7/issue-2/An-adaptive-Metropolis-algorithm/bj/1080222083.full
  /// @param _fCurrVal Value of each parameter necessary for updating throw matrix
  void UpdateAdaptiveCovariance(const std::vector<double>& _fCurrVal, const int steps_post_burn, const int Npars);

  /// Meta variables related to adaption run time
  /// When do we start throwing
  int start_adaptive_throw;

  /// When do we stop update the adaptive matrix
  int start_adaptive_update;

  /// Steps between changing throw matrix
  int end_adaptive_update;

  /// Steps between changing throw matrix
  int adaptive_update_step;

  /// Indices for block-matrix adaption
  std::vector<int> adapt_block_matrix_indices;

  /// Size of blocks for adaption
  std::vector<int> adapt_block_sizes;

  // Variables directedly linked to adaption
  /// Mean values for all parameters
  std::vector<double> par_means;

  /// Full adaptive covariance matrix
  TMatrixDSym* adaptive_covariance;
};

} // adaptive_mcmc namespace
