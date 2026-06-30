/// @file SampleHandlerNuDockBase.h
/// @brief Client-side sample handler that delegates reweighting and likelihood
///        evaluation to a remote NuDock server.
/// @author Hank Hua

#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"
#include "NuDockFactory.h"
_MaCh3_Safe_Include_Start_ //{
#include "nudock.hpp"
_MaCh3_Safe_Include_End_ //}

/// @brief A SampleHandlerInterface class that acts as a client to an external NuDock server.
///
/// Instead of performing local event reweighting and likelihood calculation,
/// this sample forwards the current parameter values to an external NuDock server via
/// JSON requests and retrieves the computed likelihood. 
///
/// @author Hank Hua
class SampleHandlerNuDockBase : public SampleHandlerInterface {
public:
  /// @brief Construct the NuDock client sample handler.
  ///
  /// @details Reads the NuDockClient configuration from `configFile`, initialises the
  /// NuDock communication object, and caches the parameter indices that must
  /// be forwarded to the server on each reweight call.
  ///
  /// @param configFile Path to a YAML configuration file containing the
  ///                   NuDockClient block.
  /// @param xsec_cov   Pointer to the cross-section parameter handler from
  ///                   which current parameter values are read.
  SampleHandlerNuDockBase(std::string configFile, ParameterHandlerGeneric* xsec_cov);

  /// @brief Destructor.
  virtual ~SampleHandlerNuDockBase();

  /// @brief Send current parameter values to the NuDock server.
  ///
  /// @details Collects oscillation and systematic parameter values from the
  /// ParameterHandlerGeneric, converts oscillation parameters from MaCh3 to
  /// NuDock convention, and sends a "/set_parameters" request.
  void Reweight() override;

  /// @brief Retrieve the log-likelihood from the NuDock server.
  ///
  /// @details Sends a "/log_likelihood" request and converts the returned 2NLL value
  /// to MaCh3's NLL convention by dividing by 2.
  ///
  /// @return The log-likelihood value from the remote server.
  double GetLikelihood() const override;

  /// @copydoc SampleHandlerInterface::GetSampleTitle
  std::string GetSampleTitle([[maybe_unused]] const int Sample) const override { return "NuDockSample"; };

  /// @copydoc SampleHandlerInterface::GetName
  std::string GetName() const override { return "NuDockSample"; };

  /// @brief Get the likelihood for a specific sub-sample.
  /// @param isample Sub-sample index (unused -- delegates to GetLikelihood()).
  /// @return The total log-likelihood from the remote server.
  double GetSampleLikelihood([[maybe_unused]] const int isample) const override { return GetLikelihood(); };

  /// @brief Print event rates (no-op for the NuDock client).
  /// @param DataOnly Whether to print data-only rates (unused).
  void PrintRates([[maybe_unused]] const bool DataOnly = false) override { MACH3LOG_INFO("No rates to print for NuDock sample handler"); };

  /// @brief Get the number of oscillation channels for a sample.
  /// @param iSample Sample index (unused).
  /// @return Always returns 0 - oscillation channels are managed server-side.
  int GetNOscChannels([[maybe_unused]] const int iSample) const override { return 0; };

  /// @copydoc SampleHandlerInterface::GetKinVarName
  /// @note Functions for posterior predictive - left unimplemented in the base class since they are not needed for likelihood evaluation and may require experiment-specific handling
  /// @copydoc SampleHan
  std::string GetKinVarName([[maybe_unused]] const int iSample, [[maybe_unused]] const int Dimension) const override { return ""; };

  /// @copydoc SampleHandlerInterface::GetDataHist
  const TH1* GetDataHist([[maybe_unused]] const int Sample) override { return nullptr; };
  /// @copydoc SampleHandlerInterface::GetMCHist
  const TH1* GetMCHist([[maybe_unused]] const int Sample) override { return nullptr; };
  /// @copydoc SampleHandlerInterface::GetW2Hist
  const TH1* GetW2Hist([[maybe_unused]] const int Sample) override { return nullptr; };

  /// @copydoc SampleHandlerInterface::GetNDim
  int GetNDim([[maybe_unused]] const int Sample) const override { return 0; };
  /// @copydoc SampleHandlerInterface::GetFlavourName
  std::string GetFlavourName([[maybe_unused]] const int iSample, [[maybe_unused]] const int iChannel) const override { return ""; };

  /// @copydoc SampleHandlerInterface::ReturnKinematicParameterBinning
  std::vector<double> ReturnKinematicParameterBinning([[maybe_unused]] const int Sample, [[maybe_unused]] const std::string &KinematicParameter) const override { return {}; };
  /// @copydoc SampleHandlerInterface::Get1DVarHistByModeAndChannel
  std::unique_ptr<TH1> Get1DVarHistByModeAndChannel([[maybe_unused]] const int iSample, [[maybe_unused]] const std::string& ProjectionVar_Str,
                                            [[maybe_unused]] int kModeToFill = -1, [[maybe_unused]] int kChannelToFill = -1,
                                            [[maybe_unused]] int WeightStyle = 0) override { return nullptr; };
                                            /// @copydoc SampleHandlerInterface::Get2DVarHistByModeAndChannel
  std::unique_ptr<TH2> Get2DVarHistByModeAndChannel([[maybe_unused]] const int iSample, [[maybe_unused]] const std::string& ProjectionVar_StrX,
                                            [[maybe_unused]] const std::string& ProjectionVar_StrY, [[maybe_unused]] int kModeToFill = -1,
                                            [[maybe_unused]] int kChannelToFill = -1, [[maybe_unused]] int WeightStyle = 0) override { return nullptr;};
  /// @copydoc SampleHandlerInterface::Get1DVarHist
  std::unique_ptr<TH1> Get1DVarHist([[maybe_unused]] const int iSample, [[maybe_unused]] const std::string &ProjectionVar,
                           [[maybe_unused]] const std::vector<KinematicCut> &EventSelectionVec = {}, [[maybe_unused]] int WeightStyle = 0,
                           [[maybe_unused]] const std::vector<KinematicCut> &SubEventSelectionVec = {}) override {return nullptr;};
  /// @copydoc SampleHandlerInterface::Get2DVarHist
  std::unique_ptr<TH2> Get2DVarHist([[maybe_unused]] const int iSample, [[maybe_unused]] const std::string& ProjectionVarX,
                                    [[maybe_unused]] const std::string& ProjectionVarY, [[maybe_unused]] const std::vector< KinematicCut >& EventSelectionVec = {},
                                    [[maybe_unused]] int WeightStyle = 0, [[maybe_unused]] const std::vector< KinematicCut >& SubEventSelectionVec = {}) override { return nullptr; };

protected:
  /// @brief No-op -- memory cleanup is handled by the NuDock server.
  void CleanMemoryBeforeFit() override {};

  /// @brief Initialise the NuDock client connection and cache parameter indices.
  ///
  /// @details Calls InitialiseNuDockObj() and then gathers system parameter indices
  /// (kFunc, kNorm, kSpline) tagged with sample "NuDock" from the ParameterHandler.
  void Init();

  /// @brief Pointer to the NuDock client communication object.
  std::unique_ptr<NuDock> nudock_ptr;

  /// @brief Manager owning the NuDockClient configuration.
  std::unique_ptr<Manager> SampleManager;

  /// @brief Verbose logging flag, read from the NuDockClient config block.
  bool verbose;

  /// @brief Cached indices into the ParameterHandler for parameters sent to the server.
  std::vector<int> nudockParamInds;

  /// @brief Non-owning pointer to the cross-section ParameterHandler.
  ParameterHandlerGeneric* ParHandler;
};
