#pragma once

// MaCh3 Includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"

namespace adaptive_mcmc{

/// @brief Contains information about adaptive covariance matrix
/// @cite haario2001adaptive
/// @author Henry Wallace
/// @details struct encapsulating all adaptive MCMC information
class AdaptiveMCMCHandler{
 public:

  /// @brief Constructor
  AdaptiveMCMCHandler();

  /// @brief Destructor
  virtual ~AdaptiveMCMCHandler();

  /// @brief Print all class members
  void Print();

  /// @brief Read initial values from config file
  /// @param adapt_manager Config file from which we update matrix
  bool InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str,
                      const std::vector<double>* parameters, const std::vector<double>* fixed);

  /// @brief If we don't have a covariance matrix to start from for adaptive tune we need to make one!
  void CreateNewAdaptiveCovariance();

  /// @brief HW: sets adaptive block matrix
  /// @param block_indices Values for sub-matrix blocks
  void SetAdaptiveBlocks(const std::vector<std::vector<int>>& block_indices);

  /// @brief HW: Save adaptive throw matrix to file
  void SaveAdaptiveToFile(const std::string& outFileName, const std::string& systematicName, const bool is_final=false);

  /// @brief sets throw matrix from a file
  /// @param matrix_file_name name of file matrix lives in
  /// @param matrix_name name of matrix in file
  /// @param means_name name of means vec in file
  void SetThrowMatrixFromFile(const std::string& matrix_file_name,
                              const std::string& matrix_name,
                              const std::string& means_name,
                              bool& use_adaptive);

  /// @brief Method to update adaptive MCMC
  /// @cite haario2001adaptive
  /// @param _fCurrVal Value of each parameter necessary for updating throw matrix
  void UpdateAdaptiveCovariance();

  /// @brief Tell whether we want reset step scale or not
  bool IndivStepScaleAdapt();

  /// @brief Tell whether matrix should be updated
  bool UpdateMatrixAdapt();

  /// @brief To be fair not a clue...
  bool AdaptionUpdate();

  /// @brief Tell if we are Skipping Adaption
  bool SkipAdaption();

  /// @brief Set the current values of the parameters
  void SetParams(const std::vector<double>* params){
    _fCurrVal = params;
  }

  /// @brief Set the fixed parameters
  void SetFixed(const std::vector<double>* fix){
    _fFixedPars = fix;
  }

  /// @brief Get the current values of the parameters
  int GetNumParams() const {
    return static_cast<int>(_fCurrVal->size());
  }

  /// @brief Check if a parameter is fixed
  bool IsFixed(const int ipar) const {
    if(!_fFixedPars){
      return false;
    }
    return ((*_fFixedPars)[ipar] < 0);
  }

  /// @brief Get Current value of parameter
  double CurrVal(const int par_index);

  /// Meta variables related to adaption run time
  /// When do we start throwing
  int start_adaptive_throw;

  /// When do we stop update the adaptive matrix
  int start_adaptive_update;

  /// Steps between changing throw matrix
  int end_adaptive_update;

  /// Steps between changing throw matrix
  int adaptive_update_step;

  /// If you don't want to save every adaption then
  /// you can specify this here
  int adaptive_save_n_iterations;

  /// Name of the file to save the adaptive matrices into
  std::string output_file_name;

  /// Indices for block-matrix adaption
  std::vector<int> adapt_block_matrix_indices;

  /// Size of blocks for adaption
  std::vector<int> adapt_block_sizes;

  // Variables directedly linked to adaption
  /// Mean values for all parameters
  std::vector<double> par_means;

  /// Full adaptive covariance matrix
  TMatrixDSym* adaptive_covariance;

  /// Total number of MCMC steps
  int total_steps;

  /// Scaling factor
  double adaption_scale;

  /// Vector of fixed parameters
  const std::vector<double>* _fFixedPars;

  /// Current values of parameters
  const std::vector<double>* _fCurrVal;
};

} // adaptive_mcmc namespace
