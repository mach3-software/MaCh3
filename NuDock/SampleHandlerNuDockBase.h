/// @file SampleHandlerNuDockBase.h
/// @brief Client-side sample handler that delegates reweighting and likelihood
///        evaluation to a remote NuDock server.
/// @author Hank Hua

#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerBase.h"
#include "NuDockFactory.h"
_MaCh3_Safe_Include_Start_ //{
#include "nudock.hpp"
_MaCh3_Safe_Include_End_ //}

/// @brief A SampleHandlerBase class that acts as a client to an external NuDock server.
///
/// Instead of performing local event reweighting and likelihood calculation,
/// this sample forwards the current parameter values to an external NuDock server via
/// JSON requests and retrieves the computed likelihood. 
///
/// @author Hank Hua
class SampleHandlerNuDockBase : public SampleHandlerBase {
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
  virtual void Reweight() override;

  /// @brief Retrieve the log-likelihood from the NuDock server.
  ///
  /// @details Sends a "/log_likelihood" request and converts the returned 2NLL value
  /// to MaCh3's NLL convention by dividing by 2.
  ///
  /// @return The log-likelihood value from the remote server.
  virtual double GetLikelihood() const override;

  /// @brief Get the title string for a given sample index.
  /// @param Sample Sample index (unused - always returns "NuDockSample").
  /// @return The fixed string "NuDockSample".
  virtual std::string GetSampleTitle(const int Sample) const override { (void)Sample; return "NuDockSample"; };

  /// @brief Get the name of this sample handler.
  /// @return The fixed string "NuDockSample".
  virtual std::string GetName() const override { return "NuDockSample"; };

  /// @brief Get the likelihood for a specific sub-sample.
  /// @param isample Sub-sample index (unused -- delegates to GetLikelihood()).
  /// @return The total log-likelihood from the remote server.
  virtual double GetSampleLikelihood(const int isample) const override { (void)isample; return GetLikelihood(); };

  /// @brief Print event rates (no-op for the NuDock client).
  /// @param DataOnly Whether to print data-only rates (unused).
  virtual void PrintRates(const bool DataOnly = false) override { (void)DataOnly; MACH3LOG_INFO("No rates to print for NuDock sample handler"); };

  /// @brief Get the number of oscillation channels for a sample.
  /// @param iSample Sample index (unused).
  /// @return Always returns 0 - oscillation channels are managed server-side.
  virtual int GetNOscChannels(const int iSample) const override { (void)iSample; return 0; };

  // Functions for posterior predictive - left unimplemented in the base class since they are not needed for likelihood evaluation and may require experiment-specific handling
  /// @brief Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
  /// @param iSample Sample index
  /// @param Dimension Dimension index
  virtual std::string GetKinVarName(const int iSample, const int Dimension) const override { (void)iSample; (void)Dimension; return ""; };

  /// @brief Get Data histogram
  virtual TH1* GetDataHist(const int Sample) override { (void)Sample; return nullptr; };
  /// @brief Get MC histogram
  virtual TH1* GetMCHist(const int Sample) override { (void)Sample; return nullptr; };
  /// @brief Get W2 histogram
  virtual TH1* GetW2Hist(const int Sample) override { (void)Sample; return nullptr; };

  /// @brief DB Function to differentiate 1D or 2D binning
  virtual int GetNDim(const int Sample) const override { (void)Sample; return 0; };
  virtual std::string GetFlavourName(const int iSample, const int iChannel) const override { (void)iSample; (void)iChannel; return ""; };

  /// @brief Return the binning used to draw a kinematic parameter
  virtual std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const override { (void)Sample; (void)KinematicParameter; return {}; };

  virtual TH1* Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                            int kModeToFill = -1, int kChannelToFill = -1, int WeightStyle = 0, TAxis* Axis = nullptr) override { 
    (void)iSample; (void)ProjectionVar_Str; (void)kModeToFill; (void)kChannelToFill; (void)WeightStyle; (void)Axis;
    return nullptr;
  };

  virtual TH2* Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                            const std::string& ProjectionVar_StrY, int kModeToFill = -1,
                                            int kChannelToFill = -1, int WeightStyle = 0,
                                            TAxis* AxisX = nullptr, TAxis* AxisY = nullptr) override { 
    (void)iSample; (void)ProjectionVar_StrX; (void)ProjectionVar_StrY; (void)kModeToFill; (void)kChannelToFill; (void)WeightStyle; (void)AxisX; (void)AxisY;
    return nullptr;
  };

  virtual TH1 *Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                           const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                            TAxis *Axis = nullptr, const std::vector<KinematicCut> &SubEventSelectionVec = {}) override {
    (void)iSample; (void)ProjectionVar; (void)EventSelectionVec; (void)WeightStyle; (void)Axis; (void)SubEventSelectionVec;
    return nullptr;
  };

  virtual TH2* Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                            const std::vector< KinematicCut >& EventSelectionVec = {},
                            int WeightStyle = 0, TAxis* AxisX = nullptr, TAxis* AxisY = nullptr,
                            const std::vector< KinematicCut >& SubEventSelectionVec = {}) override {
    (void)iSample; (void)ProjectionVarX; (void)ProjectionVarY; (void)EventSelectionVec; (void)WeightStyle; (void)AxisX; (void)AxisY; (void)SubEventSelectionVec;
    return nullptr;
  };


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
