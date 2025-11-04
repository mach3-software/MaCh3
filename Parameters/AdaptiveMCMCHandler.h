#pragma once

// MaCh3 Includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"

_MaCh3_Safe_Include_Start_ //{
#include "Math/DistFunc.h"
_MaCh3_Safe_Include_End_ //}

namespace adaptive_mcmc
{
  /// @brief Contains information about adaptive covariance matrix
  /// @cite haario2001adaptive
  /// @cite Garthwaite01092016
  /// @author Henry Wallace
  /// @details struct encapsulating all adaptive MCMC information. Adaption can be done in the pure "Haario/Roberts-Rosenthal"
  /// style using a constant scale of 2.38^2/num_params or via a Robbins-Monro style adaption where the scale is updated.
  /// The latter approach allows you to target a specific acceptance rate resulting in potentially better mixing.
  class AdaptiveMCMCHandler
  {
  public:
    /// @brief Constructor
    AdaptiveMCMCHandler();

    /// @brief Destructor
    virtual ~AdaptiveMCMCHandler();

    /// @brief Print all class members
    void Print() const;

    /// @brief Read initial values from config file
    ///
    /// @param adapt_manager YAML node containing the configuration (`AdaptionOptions`).
    /// @param matrix_name_str Name of the covariance matrix block to configure.
    /// @param parameters Pointer to a vector of parameter values (nominal values).
    /// @param fixed Pointer to a vector of fixed parameter values.
    /// @return True if adaptive MCMC configuration was successfully initialized, false otherwise.
    bool InitFromConfig(const YAML::Node &adapt_manager, const std::string &matrix_name_str,
                        const std::vector<double> *parameters, const std::vector<double> *fixed);

    /// @brief If we don't have a covariance matrix to start from for adaptive tune we need to make one!
    void CreateNewAdaptiveCovariance();

    /// @brief HW: sets adaptive block matrix
    /// @param block_indices Values for sub-matrix blocks
    void SetAdaptiveBlocks(const std::vector<std::vector<int>> &block_indices);

    /// @brief HW: Save adaptive throw matrix to file
    void SaveAdaptiveToFile(const std::string &outFileName, const std::string &systematicName, const bool is_final = false);

    /// @brief sets throw matrix from a file
    /// @param matrix_file_name name of file matrix lives in
    /// @param matrix_name name of matrix in file
    /// @param means_name name of means vec in file
    void SetThrowMatrixFromFile(const std::string &matrix_file_name,
                                const std::string &matrix_name,
                                const std::string &means_name,
                                bool &use_adaptive);

    /// @brief Check if there are structures in matrix that could result in failure of adaption fits.
    /// this mostly cover cases of highly correlated block which seem to results in adaption failure
    /// @param covMatrix covariance matrix which we convert into correlation matrix
    void CheckMatrixValidityForAdaption(const TMatrixDSym *covMatrix) const;

    /// @brief Method to update adaptive MCMC
    /// @cite haario2001adaptive
    /// @param _fCurrVal Value of each parameter necessary for updating throw matrix
    void UpdateAdaptiveCovariance();

    /// @brief Tell whether we want reset step scale or not
    bool IndivStepScaleAdapt() const;

    /// @brief Tell whether matrix should be updated
    bool UpdateMatrixAdapt();

    /// @brief To be fair not a clue...
    bool AdaptionUpdate() const;

    /// @brief Tell if we are Skipping Adaption
    bool SkipAdaption() const;

    /// @brief Set the current values of the parameters
    void SetParams(const std::vector<double> *params)
    {
      _fCurrVal = params;
    }

    /// @brief Set the fixed parameters
    void SetFixed(const std::vector<double> *fix)
    {
      _fFixedPars = fix;
    }

    /// @brief Get the number of fixed parameters
    /// @ingroup ParameterHandlerGetters
    int GetNFixed() const
    {
      if (!_fFixedPars)
      {
        return 0;
      }
      int n_fixed = 0;
      for (int i = 0; i < static_cast<int>(_fFixedPars->size()); i++)
      {
        if ((*_fFixedPars)[i] < 0)
        {
          n_fixed++;
        }
      }
      return n_fixed;
    }

    /// @brief Get the current values of the parameters
    /// @ingroup ParameterHandlerGetters
    int GetNumParams() const
    {
      return static_cast<int>(_fCurrVal->size());
    }

    /// @brief Check if a parameter is fixed
    bool IsFixed(const int ipar) const
    {
      if (!_fFixedPars)
      {
        return false;
      }
      return ((*_fFixedPars)[ipar] < 0);
    }

    /// @brief Get Current value of parameter
    double CurrVal(const int par_index) const;

    /// @brief Get Total Number of Steps
    /// @ingroup ParameterHandlerGetters
    int GetTotalSteps() const
    {
      return total_steps;
    }

    /// @brief Change Total Number of Steps to new value
    void SetTotalSteps(const int nsteps)
    {
      total_steps = nsteps;
    }

    /// @brief Increase by one number of total steps
    void IncrementNSteps()
    {
      total_steps++;
    }

    void IncrementAcceptedSteps()
    {
      prev_step_accepted = true;
    }

    /// @brief Increase by one number of total steps
    /// @ingroup ParameterHandlerGetters
    TMatrixDSym *GetAdaptiveCovariance() const
    {
      return adaptive_covariance;
    }

    /// @brief Get the parameter means used in the adaptive handler
    /// @ingroup ParameterHandlerGetters
    std::vector<double> GetParameterMeans() const
    {
      return par_means;
    }

    /// @brief Get Name of Output File
    /// @ingroup ParameterHandlerGetters
    std::string GetOutFileName() const
    {
      return output_file_name;
    }

    /// @brief Get the current adaption scale
    double GetAdaptionScale()
    {
      if (use_robbins_monro)
      {
        UpdateRobbinsMonroScale();
      }
      return adaption_scale;
    }

    /// Use Robbins-Monro approach?
    bool GetUseRobbinsMonro() const
    {
      return use_robbins_monro;
    }

  private:
    /// Update the scale factor for Robbins-Monro adaption
    void UpdateRobbinsMonroScale();
    /// Calculate the constant step length for Robbins-Monro adaption
    void CalculateRobbinsMonroStepLength();

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
    TMatrixDSym *adaptive_covariance;

    /// Total number of MCMC steps
    int total_steps;

    /// Scaling factor
    double adaption_scale;

    /// Vector of fixed parameters
    const std::vector<double> *_fFixedPars;

    /// Current values of parameters
    const std::vector<double> *_fCurrVal;

    /// Acceptance rate in the current batch
    int acceptance_rate_batch_size;

    /// Use Robbins Monro https://arxiv.org/pdf/1006.3690
    bool use_robbins_monro;

    /// Target acceptance rate for Robbins Monro
    double target_acceptance;

    /// Constant "step scaling" factor for Robbins-Monro
    double c_robbins_monro;

    /// Number of restarts for Robbins Monro (so far)
    int n_rm_restarts;

    /// Total number of restarts ALLOWED for Robbins Monro
    int total_rm_restarts;

    /// Need to keep track of whether previous step was accepted for RM
    bool prev_step_accepted;
  };
} // adaptive_mcmc namespace
