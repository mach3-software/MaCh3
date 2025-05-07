#pragma once

// MaCh3 Includes
#include "manager/manager.h"
#include "covariance/CovarianceUtils.h"
#include "unordered_set"

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
  bool InitFromConfig(const YAML::Node& adapt_manager, const std::string& matrix_name_str, std::vector<double>* parameters, std::vector<double>* fixed);

  /// @brief If we don't have a covariance matrix to start from for adaptive tune we need to make one!
  void CreateNewAdaptiveCovariance();

  /// @brief HW: sets adaptive block matrix
  /// @param block_indices Values for sub-matrix blocks
  void SetAdaptiveBlocks(std::vector<std::vector<int>> block_indices);

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
  void UpdateAdaptiveCovariance();

  /// @brief Tell whether we want reset step scale or not
  bool IndivStepScaleAdapt();

  /// @brief Tell whether matrix should be updated
  bool UpdateMatrixAdapt();

  /// @brief To be fair not a clue...
  bool AdaptionUpdate();

  /// @brief Tell if we are Skipping Adaption
  bool SkipAdaption();

  void SetParams(std::vector<double>* params){
    _fCurrVal = params;
  }

  void SetFixed(std::vector<double>* fix){
    errors = fix;
  }

  int GetNPars(){
    return static_cast<int>(_fCurrVal->size());
  }

  double CurrVal(int par_index);

  bool IsFixed(int ipar){

    if(!errors){
      return false;
    }

    return ((*errors)[ipar] < 0);
  }

  void SetScale(double scale){
    adaption_scale = scale;
  }


  double ModeAdjusted(int par_index);
  void SetBiModal(int par_index, double midpoint);

  /// Meta variables related to adaption run time
  /// When do we start throwing
  int start_adaptive_throw;

  /// When do we stop update the adaptive matrix
  int start_adaptive_update;

  /// Steps between changing throw matrix
  int end_adaptive_update;

  /// Steps between changing throw matrix
  int adaptive_update_step;

  /// If you don't want to save every adpation then
  /// you can specify this here
  int adaptive_save_n_iterations;

  /// Name of the file to save the adpative matrices into
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

  double adaption_scale;

  std::vector<double>* errors;
  std::vector<double>* _fCurrVal;
  std::unordered_map<int, double> bimodal_pars;

};

} // adaptive_mcmc namespace
