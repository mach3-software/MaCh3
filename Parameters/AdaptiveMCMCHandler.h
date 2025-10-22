#pragma once

// MaCh3 Includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"
#include <deque> // NEW: For sliding window

namespace adaptive_mcmc
{

  // A rough approximation of a mode
  class AdaptiveMode
  {
  public:
    AdaptiveMode(const std::vector<double> &point);
    void Update(const std::vector<double> &point);
    double NSigmaFromMode(std::vector<double> point) const;

    /// @brief Get the mean values for this mode
    std::vector<double> GetMeans() const { return means; }

    /// @brief Get the variance values for this mode
    std::vector<double> GetVariances() const { return variances; }

    /// @brief Get the number of steps accumulated in this mode
    int GetNSteps() const { return n_steps; }

  protected:
    std::vector<double> means;
    std::vector<double> variances;
    int n_steps;
  };

  /// @brief Contains information about adaptive covariance matrix
  /// @cite haario2001adaptive
  /// @author Henry Wallace
  /// @details struct encapsulating all adaptive MCMC information
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
    bool InitFromConfig(const YAML::Node &adapt_manager, const std::string &matrix_name_str,
                        const std::vector<double> *parameters, const std::vector<double> *fixed);

    /// @brief If we don't have a covariance matrix to start from for adaptive tune we need to make one!
    void CreateNewAdaptiveCovariance();

    /// @brief HW: sets adaptive block matrix
    void SetAdaptiveBlocks(const std::vector<std::vector<int>> &block_indices);

    /// @brief HW: Save adaptive throw matrix to file
    void SaveAdaptiveToFile(const std::string &outFileName, const std::string &systematicName, const bool is_final = false);

    /// @brief sets throw matrix from a file
    void SetThrowMatrixFromFile(const std::string &matrix_file_name,
                                const std::string &matrix_name,
                                const std::string &means_name,
                                bool &use_adaptive);

    /// @brief Method to update adaptive MCMC
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

    /// @brief Get the current values of the parameters
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

    /// @brief Get adaptive covariance matrix
    TMatrixDSym *GetAdaptiveCovariance() const
    {
      return adaptive_covariance;
    }

    /// @brief Get the parameter means used in the adaptive handler
    std::vector<double> GetParameterMeans() const
    {
      return par_means;
    }

    /// @brief Get Name of Output File
    std::string GetOutFileName() const
    {
      return output_file_name;
    }

    /// @brief Enable or disable mode tracking
    void SetModeTracking(bool enable)
    {
      track_modes = enable;
    }

    /// @brief Check if mode tracking is enabled
    bool IsModeTrackingEnabled() const
    {
      return track_modes;
    }

    /// @brief Get the number of detected modes
    size_t GetNModes() const
    {
      return mode_information.size();
    }

    /// @brief Print information about detected modes
    void PrintModeInfo() const;

    /// @brief Set the size of the sliding window for covariance estimation
    void SetWindowSize(int size)
    {
      window_size = size;
    }

    /// @brief Get the current window size
    int GetWindowSize() const
    {
      return window_size;
    }

    /// @brief Set minimum steps required for stability before adapting
    void SetMinStableSteps(int steps)
    {
      min_stable_steps = steps;
    }

    /// @brief Set the mode jump detection threshold (in sigma)
    void SetModeJumpThreshold(double threshold)
    {
      mode_jump_threshold = threshold;
    }

    /// @brief Set the covariance shrink factor applied after mode jumps
    void SetPostJumpShrinkFactor(double factor)
    {
      post_jump_shrink_factor = factor;
    }

  private:
    /// Meta variables related to adaption run time
    int start_adaptive_throw;
    int start_adaptive_update;
    int end_adaptive_update;
    int adaptive_update_step;
    int adaptive_save_n_iterations;
    std::string output_file_name;

    /// Indices for block-matrix adaption
    std::vector<int> adapt_block_matrix_indices;
    std::vector<int> adapt_block_sizes;

    // Variables directly linked to adaption
    std::vector<double> par_means;
    TMatrixDSym *adaptive_covariance;
    int total_steps;
    double adaption_scale;

    /// Fixed and current parameter values
    const std::vector<double> *_fFixedPars;
    const std::vector<double> *_fCurrVal;

    /// Mode tracking
    bool track_modes;
    std::vector<std::unique_ptr<AdaptiveMode>> mode_information;
    void UpdateModes(const std::vector<double> &point);
    int GetClosestMode(const std::vector<double> &point) const;
    void MergeSimilarModes();
    std::vector<double> MoveToInitMode(const std::vector<double> &point) const;

    /// NEW: Windowed adaptation for multimodal distributions
    std::deque<std::vector<double>> recent_points; // Sliding window of recent MCMC points
    std::vector<double> prev_point;                // Previous MCMC point
    int steps_since_mode_change;                   // Steps since last mode jump
    int current_mode_idx;                          // Index of current mode

    // Configurable parameters for windowed adaptation
    int window_size;                    // Size of sliding window (default: 500)
    int min_stable_steps;               // Min steps in mode before adapting (default: 50)
    double mode_jump_threshold;         // N-sigma threshold for mode jump (default: 3.0)
    double post_jump_shrink_factor;     // Covariance shrink after jump (default: 0.5)
    double exponential_smoothing_alpha; // Weight for new covariance (default: 0.8)
  };

} // adaptive_mcmc namespace