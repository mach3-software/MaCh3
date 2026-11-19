#pragma once

//MaCh3 includes
#include "Samples/SampleInfo.h"
#include "Samples/EventInfo.h"

#include "Splines/BinnedSplineHandler.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"
#include "Samples/OscillationHandler.h"
#include "Samples/LikelihoodGradient.h"
#include "Samples/ReweightOptions.h"
#include "Samples/FastEventAccess.h"
#include "Samples/OscillationGradient.h"
#include "Samples/WeightParameterSnapshot.h"
#include "Samples/BinningHandler.h"
#include "Samples/SampleHandlerFunctional.h"

_MaCh3_Safe_Include_Start_ //{
#include "THStack.h"
#include "TLegend.h"
_MaCh3_Safe_Include_End_ //}

enum SamplePlotType {
  kModePlot = 0,
  kOscChannelPlot = 1
};

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
///
/// @details This class manages samples used in the analysis. It supports event-by-event reweighting and allows
/// kinematic shifts to be applied at the event level, which may cause migration between analysis bins.
///
/// The likelihood calculation is performed in binned space.
///
/// @author Dan Barrow
/// @author Ed Atkin
///
/// @ingroup SamplesAndParameters
class SampleHandlerBase :  public SampleHandlerInterface
{
 public:
  //######################################### Functions #########################################
  /// @brief Constructor
  /// @param ConfigFileName Name of config to initialise the sample object
  SampleHandlerBase(std::string ConfigFileName, ParameterHandlerGeneric* _ParHandler,
                  const std::shared_ptr<OscillationHandler>& OscillatorObj_ = nullptr);
  /// @brief destructor
  virtual ~SampleHandlerBase();

  /// @copydoc SampleHandlerInterface::GetNDim
  int GetNDim(const int Sample) const final { return SampleDetails[Sample].nDimensions; }
  /// @copydoc SampleHandlerInterface::GetName
  std::string GetName() const final;
  /// @copydoc SampleHandlerInterface::GetSampleTitle
  std::string GetSampleTitle(const int iSample) const final {return SampleDetails[iSample].SampleTitle;}
  /// @brief Sample name tag used only for getting relevant uncertainties
  /// @param iSample Index of the sample.
  std::string GetSampleName(const int iSample) const {return SampleDetails[iSample].SampleName;}
  /// @copydoc SampleHandlerInterface::GetKinVarName
  std::string GetKinVarName(const int iSample, const int Dimension) const final;

  /// @brief Computes and prints the integral breakdown of all modes and oscillation channels for a given sample.
  /// @param iSample Index of the sample.
  void PrintIntegral(const int iSample, const TString& OutputName="/dev/null", const int WeightStyle=0, const TString& OutputCSVName="/dev/null");
  /// @brief Print number of events for a given sample with specified cuts
  double GetSampleIntegral(const int iSample, const std::vector<KinematicCut> &EventSelectionVec, int WeightStyle);
  /// @brief Return sample integral after for given sample with specified restriction like only single mode
  std::vector<double> ReturnSampleIntegralBySelection(const int iSample, const SamplePlotType Selection1,
                                                      const int Selection2, const int WeightStyle);

  /// @brief DB: Add data for a given sample from a ROOT histogram.
  /// @param Sample Index of the sample.
  /// @param Data Pointer to a TH1 containing the data to be stored.
  void AddData(const int Sample, TH1* Data);

  /// @brief ETA:  Add data for a given sample from a raw array.
  /// @param Sample Index of the sample.
  /// @param Data_Array Vector containing the data values.
  void AddData(const int Sample, const std::vector<double>& Data_Array);

  /// @copydoc SampleHandlerInterface::PrintRates
  void PrintRates(const bool DataOnly = false) final;
  /// @copydoc SampleHandlerInterface::GetLikelihood
  double GetLikelihood() const override;
  /// @copydoc SampleHandlerInterface::SetComputeGradients
  void SetComputeGradients(const bool Enable) override;
  /// @copydoc SampleHandlerInterface::SetGradientsPaused
  void SetGradientsPaused(const bool Paused) override { GradientsPaused = Paused; }
  /// @copydoc SampleHandlerInterface::GetGradLikelihood
  std::vector<double> GetGradLikelihood() const override;
  /// @copydoc SampleHandlerInterface::GetGradParHandler
  const ParameterHandlerBase* GetGradParHandler() const override;
  /// @copydoc SampleHandlerInterface::GetGradientCoverage
  std::vector<bool> GetGradientCoverage() const override;
  /// @copydoc SampleHandlerInterface::EnableSurrogate
  bool EnableSurrogate(const int Stride) override;
  /// @copydoc SampleHandlerInterface::AnchorSurrogate
  void AnchorSurrogate() override;
  /// @copydoc SampleHandlerInterface::ReweightSurrogate
  void ReweightSurrogate() override;
  /// @copydoc SampleHandlerInterface::ReverseSurrogateLikelihood
  double ReverseSurrogateLikelihood() const override;
  /// @brief Compare the reweight backend in use with the plain event loop
  /// @details Moves the parameters through NPoints points, cycling between moving only the
  /// functional (shift) parameters, only the others, and both, which takes the backend through each
  /// of its passes. At every point it reweights with the backend and then with the plain loop, and
  /// compares the MC prediction and W2 bin by bin, the likelihood and, if gradients are on, the
  /// gradient. The parameters are put back and the handler reweighted before returning.
  /// @param NPoints Number of points
  /// @param Tolerance Differences above this are logged as errors
  /// @param Width Parameters are thrown this many prior widths around their prefit values
  /// @param Seed For the points
  /// @return The largest relative difference found, 0 if no backend is in use
  double CheckReweightBackend(int NPoints = 3, double Tolerance = 1e-9, double Width = 1., unsigned int Seed = 2026);
  /// @brief Name of the reweight backend in use, "legacy" for the plain event loop
  std::string GetReweightBackendName() const { return Backend ? Backend->Name() : "legacy"; }
  /// @brief How closely the reweight backend in use reproduces the plain event loop, see ReweightBackend::Precision
  double GetReweightBackendPrecision() const { return Backend ? Backend->Precision() : 0.; }
  /// @brief Bring every event's shifted kinematic variables up to date with the last reweight
  /// @details Only does anything with the opt-in linear shifts (Reweight: LinearShifts), which work
  /// out the shifted variables without writing them to the events. Core's own plotting and integral
  /// functions call this first; experiment code that reads shifted variables itself after a
  /// reweight has to call it too
  void SyncShiftedKinematics();
  /// @brief Forget everything the reweight backend reuses between reweights
  /// @details For an experiment that changes how events are shifted, cut or weighted outside the
  /// parameters (its own cached kinematics, say), so the next reweight starts from scratch
  void InvalidateReweightCache();

  /// @copydoc SampleHandlerInterface::GetSampleLikelihood
  double GetSampleLikelihood(const int iSample) const override;
  //===============================================================================

  /// @brief Get index of sample based on name
  /// @param SampleTitle The title of the sample to search for.
  int GetSampleIndex(const std::string& SampleTitle) const;

  /// @copydoc SampleHandlerInterface::GetDataHist
  const TH1* GetDataHist(const int Sample) final;
  /// @brief Get Data histogram by sample name
  /// @param Sample Sample title
  const TH1* GetDataHist(const std::string& Sample);

  /// @copydoc SampleHandlerInterface::GetMCHist
  const TH1* GetMCHist(const int Sample) final;
  /// @brief Get MC histogram by sample title
  /// @param Sample Sample name
  const TH1* GetMCHist(const std::string& Sample);

  /// @copydoc SampleHandlerInterface::GetW2Hist
  const TH1* GetW2Hist(const int Sample) final;
  /// @brief Get W2 histogram by sample name
  /// @param Sample Sample title
  const TH1* GetW2Hist(const std::string& Sample);
  /// @copydoc SampleHandlerInterface::Reweight
  void Reweight() override;
  /// @brief Computes the total event weight for a given entry.
  /// @param iEvent Event enumerator
  M3::float_t GetEventWeight(const int iEvent);

  /// @copydoc SampleHandlerInterface::GetNOscChannels
  int GetNOscChannels(const int iSample) const final {return static_cast<int>(SampleDetails[iSample].OscChannels.size());};
  /// @copydoc SampleHandlerInterface::GetFlavourName
  std::string GetFlavourName(const int iSample, const int iChannel) const final {
    if (iChannel < 0 || iChannel > GetNOscChannels(iSample)) {
      MACH3LOG_ERROR("Invalid Channel Requested: {}", iChannel);
      throw MaCh3Exception(__FILE__ , __LINE__);
    }
    return SampleDetails[iSample].OscChannels[iChannel].flavourName;
  }
  /// @copydoc SampleHandlerInterface::Get1DVarHist
  std::unique_ptr<TH1> Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                                    const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                                    const std::vector<KinematicCut> &SubEventSelectionVec = {}) final;
  /// @copydoc SampleHandlerInterface::Get2DVarHist
  std::unique_ptr<TH2> Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                                    const std::vector< KinematicCut >& EventSelectionVec = {},
                                    int WeightStyle = 0, const std::vector< KinematicCut >& SubEventSelectionVec = {}) final;
  /// @brief Construct vector of kinematic cuts that will be applied, on top of default cuts include stuff like cut on mode etc.
  /// @param iSample Index of the sample.
  std::vector<KinematicCut> BuildModeChannelSelection(const int iSample, const int kModeToFill, const int kChannelToFill) const;
  /// @brief Fill projection histogram by looping over all events, and skipping one which doesn't pass specified condition
  void Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {},
                          int WeightStyle=0);
  /// @brief Fill projection histogram by looping over all events, and skipping one which doesn't pass specified condition
  void Fill2DSubEventHist(const int iSample, TH2* _h2DVar, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {}, int WeightStyle = 0);
  /// @copydoc SampleHandlerInterface::Get1DVarHistByModeAndChannel
  std::unique_ptr<TH1> Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                                    const int kModeToFill = -1, const int kChannelToFill = -1,
                                                    const int WeightStyle = 0) final;
  /// @copydoc SampleHandlerInterface::Get2DVarHistByModeAndChannel
  std::unique_ptr<TH2> Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                                    const std::string& ProjectionVar_StrY, const int kModeToFill = -1,
                                                    const int kChannelToFill = -1, const int WeightStyle = 0) final;

  /// @brief KS: Return range for plot type, for example number of modes, osc channels etc
  /// @param TypeEnum Plot type enumerator see @ref SamplePlotType
  /// @param iSample Sample enumerator
  int GetRangeForPlotType(const SamplePlotType TypeEnum, const int iSample) const;

  std::vector<std::unique_ptr<TH1>> ReturnHistsBySelection1D(const int iSample, const std::string& KinematicProjection,
                                                             const SamplePlotType Selection1, const int Selection2 = -1,
                                                             const int WeightStyle = 0);
  std::vector<std::unique_ptr<TH2>> ReturnHistsBySelection2D(const int iSample, const std::string& KinematicProjectionX,
                                                             const std::string& KinematicProjectionY,
                                                             const SamplePlotType Selection1, const int Selection2 = -1,
                                                             const int WeightStyle=0);
  std::unique_ptr<THStack> ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection,
                                          const SamplePlotType Selection1, const int Selection2 = -1, const int WeightStyle = 0);
  /// @brief Return the legend used for stacked histograms with sample info
  const TLegend* ReturnStackHistLegend() const {return THStackLeg;}

  /// @brief ETA function to generically convert a string from param handler to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from param handler to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

  /// @copydoc SampleHandlerInterface::SaveAdditionalInfo
  void SaveAdditionalInfo(TDirectory* Dir) final;

  /// @brief JM: Convert a kinematic vector name to its corresponding integer ID.
  int ReturnKinematicVectorFromString(const std::string& KinematicStr) const;
  /// @brief JM: Convert a kinematic vector integer ID to its corresponding name as a string.
  std::string ReturnStringFromKinematicVector(const int KinematicVariable) const;
  /// @brief JM: Check if a kinematic parameter string corresponds to a subevent-level variable
  bool IsSubEventVarString(const std::string& VarStr) const;

  /// @brief Return total number of events
  unsigned int GetNEvents() const {return nEvents;}

  /// @brief Return array storing data entries for every bin
  auto GetDataArray() const {
    return SampleHandler_data;
  }
  /// @brief Return array storing MC entries for every bin
  auto GetMCArray() const {
    return SampleHandler_array;
  }
  /// @brief Return array storing W2 entries for every bin
  auto GetW2Array() const {
    return SampleHandler_array_w2;
  }
  /// @brief Return a sub-array for a given sample.
  std::vector<double> GetArrayForSample(const int Sample, std::vector<double> const & array) const;

  /// @brief Return array storing data entries for every bin
  /// @param Sample Sample index
  std::vector<double> GetDataArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_data);
  }
  /// @brief Return array storing MC entries for every bin
  /// @param Sample Sample index
  std::vector<double> GetMCArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_array);
  }
  /// @brief Return array storing W2 entries for single sample
  /// @param Sample Sample index
  std::vector<double> GetW2Array(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandler_array_w2);
  }
  /// @brief Loop over bins and checks if there are any which have 0 entries
  void CheckEmptyBins() const;

 protected:
  /// @brief including Dan's magic NuOscillator
  void InitialiseNuOscillatorObjects();
  /// @brief Get pointer to NuOscillator weight for a given event
  /// @param iEvent Event enumerator
  const M3::float_t* GetNuOscillatorPointers(const int iEvent) const;
  /// @brief Initialise pointer to oscillation weight to NuOscillator object
  void SetupNuOscillatorPointers();
  /// @brief Load information about sample handler and corresponding samples from config file
  void ReadConfig();
  /// @brief Initialise single sample from config file
  void LoadSingleSample(const int iSample, const YAML::Node& Settings);

  /// @brief DB Function to determine which weights apply to which types of samples
  virtual void AddAdditionalWeightPointers() = 0;
  /// @brief Promise that every weight AddAdditionalWeightPointers() adds is constant for the whole
  /// fit, and that CalcWeightFunc() does nothing
  /// @details The fast reweight backends need this: they multiply each event's constant weights out
  /// once, and sum events that share their varying weights. Off unless a handler opts in, in which
  /// case the plain event loop is used
  virtual bool AdditionalWeightsAreConstant() const { return false; }
  /// @brief Whether a functional shift (or ResetShifts / FinaliseShifts) can change this kinematic variable
  /// @details A cut on a variable no shift touches gives the same answer whatever the shifts do, so
  /// events failing one can never be selected and the fast backends leave them out entirely.
  /// Every variable counts as shifted unless a handler says otherwise
  virtual bool IsKinematicParameterShifted([[maybe_unused]] const int KinematicVariable) const { return true; }

  /// @brief Ensure Kinematic Map is setup and make sure it is initialised correctly
  void SetupKinematicMap();

  /// @todo abstract the spline initialisation completely to core
  /// @brief initialise your splineXX object and then use InitialiseSplineObject to conveniently setup everything up
  virtual void SetupSplines() = 0;

  //DB Require all objects to have a function which reads in the MC
  /// @brief Initialise any variables that your experiment specific SampleHandler needs
  virtual void Init() = 0;

  /// @brief Experiment specific setup, returns the number of events which were loaded
  virtual int SetupExperimentMC() = 0;

  /// @brief Function which translates experiment struct into core struct
  virtual void SetupMC() = 0;
  /// @brief Function responsible for loading data from file or loading from file
  virtual void InititialiseData() = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();

  //===============================================================================
  /// @brief Set pointers for each event to appropriate weights, for unbinned based on event number
  /// while for binned based on other kinematical properties
  void SetSplinePointers();
  /// @brief Retrieve the spline bin indices associated with a given event.
  /// @warning ThrowCrititcal argument will be eventually removed
  std::vector< SplineIndex > GetSplineBins(int Event, BinnedSplineHandler* BinnedSpline, bool& ThrowCrititcal) const;

  /// @brief Functions which find the nominal bin and bin edges
  void FindNominalBinAndEdges();

  /// @brief set the binning for 2D sample used for the likelihood calculation
  void SetBinning();

  /// @brief Initialise data, MC and W2 histograms
  void SetupReweightArrays();
  //===============================================================================

  // ----- start Functional Parameters -----
  /// @brief HH - a experiment-specific function where the maps to actual functions are set up
  virtual void RegisterFunctionalParameters(){};
  /// @brief Update the functional parameter values to the latest proposed values. Needs to be called before every new reweight so is called in FillArray
  virtual void PrepFunctionalParameters(){};

  /// @brief LP - Registration template function for multi-dimensional functional shifts
  /// @details This method is templated over the EventType and the response
  ///          function type because the response function type is a function
  ///          of the event type. The implementation checks via static assert
  ///          that the passed function object has the correct call signature.
  ///          The experiment class event vector is passed so that the function
  ///          objects that are used to apply the shifts are able to access the
  ///          full event information from the experiment class, rather than
  ///          only being able to respond to event properties defined in
  ///          Samples/EventInfo.h
  ///          A vector of input parameter names that should be consumed by the
  ///          functional shift is also passed. For a 1D version of this method,
  ///          see below.
  ///          This function might be used like below:
  ///
  ///          RegisterIndividualFunctionalParameter(ExptEventInfo, {"par1", "par2"},
  ///            [](std::vector<double> const & par_vals, ExptEventType & ev){
  ///              constexpr static const double a = 1.2345;
  ///              ev.property += par_vals[0] * a + par_vals[1];
  ///          });
  template <typename EventType, typename SFType>
  void RegisterIndividualFunctionalParameter(
      std::vector<EventType> &ExptEvents,
      std::vector<std::string> const &par_names,
      SFType shift_func);

  /// @brief LP - Registration template function for one-dimensional functional shifts
  /// @details See above for more details on the call signature and usage.
  template <typename EventType, typename SFType>
  void RegisterIndividualFunctionalParameter(std::vector<EventType> &ExptEvents,
                                             std::string const &par_name,
                                             SFType shift_func);
  /// @brief ETA - generic function applying shifts
  /// @note It is virtual so we can perform unorthodox shifts, ideally we should de-virtualise once we ensure we can support everything in core
  virtual void ApplyShifts(const int iEvent);
  /// @brief Copy the shift parameters' current values into the flat shift table ApplyShifts uses
  /// @details The table holds every shift's plain function and value, and every event's shifts as
  /// one flat list, so an event's shifts take one typed call rather than a std::function call each.
  /// Only used when every shift was registered as a plain function of one parameter
  void RefreshShiftTable();
  /// @brief HH - reset the shifted values to the original values
  virtual void ResetShifts([[maybe_unused]] const int iEvent) {};
  /// @brief LP - Optionally calculate derived observables after all shifts have been applied
  /// @details For example, have shifts that varied lepton energy and hadron energy separately
  ///          in a subclass implementation of this method you may add the shifted quantities
  ///          together to build a shifted neutrino energy estimator
  virtual void FinaliseShifts([[maybe_unused]] const int iEvent) {};
  // ----- end Functional Parameters -----

  /// @brief DB Function which determines if an event is selected based on @ref KinematicCut
  bool IsEventSelected(const int iSample, const int iEvent) _noexcept_;
  /// @brief JM Function which determines if a subevent is selected
  bool IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, unsigned const int iSubEvent, size_t nsubevents);

  /// @brief Check whether a normalisation systematic affects an event or not
  /// @param norm_parameters indexed [sample][param] describe norm params and associated kinematic cuts etc.
  void CalcNormsBins(std::vector<std::vector<NormParameter>>& norm_parameters, std::vector< std::vector< int > >& norms_bins);
  template <typename ParT> bool PassesSelection(const ParT& Par, std::size_t iEvent);
  /// @brief Calculate the total weight weight for a given event
  M3::float_t CalcWeightTotal(const EventInfo* _restrict_ MCEvent) const _noexcept_;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc([[maybe_unused]] const int iEvent) {return;};

  /// @brief Return the value of an associated kinematic parameter for an event
  double ReturnKinematicParameter(const std::string& KinematicParameter, int iEvent) const {
    return ReturnKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  /// @brief Return the value of an associated kinematic parameter for an event
  virtual double ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const = 0;

  // === JM declare the same functions for kinematic vectors ===
  std::vector<double> ReturnKinematicVector(const std::string& KinematicParameter, const int iEvent) const {
    return ReturnKinematicVector(ReturnKinematicVectorFromString(KinematicParameter), iEvent);
  }
  virtual std::vector<double> ReturnKinematicVector([[maybe_unused]] const int KinematicVariable,
                                                    [[maybe_unused]] const int iEvent) const {return {};};
  // ===========================================================

  /// @copydoc SampleHandlerInterface::ReturnKinematicParameterBinning
  std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const final;

  const double* GetPointerToKinematicParameter(const std::string& KinematicParameter, int iEvent) const {
    return GetPointerToKinematicParameter(ReturnKinematicParameterFromString(KinematicParameter), iEvent);
  }
  virtual const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const = 0;

  /// @brief Get pointer to oscillation channel associated with given event. Osc channel is const
  const double* GetPointerToOscChannel(const int iEvent) const;
  /// @brief Setup the norm parameters by assigning each event with bin
  void SetupNormParameters();
  /// @brief Setup the osc parameters
  void SetupOscParameters();
  //===============================================================================
  /// @brief Fill a histogram with the event-level information used in the fit
  /// @details
  /// DB Functions required for reweighting functions
  /// DB Replace previous implementation with reading bin contents from SampleHandler_array
  void FillHist(const int Sample, TH1* Hist, std::vector<double> &Array);

  /// @brief Function which does the core reweighting. This assumes that oscillation weights have
  /// already been calculated.
  void FillArray();

  /// @brief Helper function to reset histograms
  void ResetHistograms() _noexcept_;
  /// @brief Setup spline handler (both binned or unbinned)
  void InitialiseSplineObject();

  /// Contains all your splines (binned or unbinned) and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<SplineBase> SplineHandler;

  /// Contains oscillator handling calculating oscillation probabilities
  std::shared_ptr<OscillationHandler> Oscillator;

  //===============================================================================
  //DB Variables required for GetLikelihood
  /// KS: This stores binning information, in future could be come vector to store binning for every used sample
  std::unique_ptr<BinningHandler> Binning;
  /// DB Array to be filled after reweighting
  std::vector<double> SampleHandler_array;
  /// KS Array used for MC stat
  std::vector<double> SampleHandler_array_w2;
  /// DB Array to be filled in AddData
  std::vector<double> SampleHandler_data;
  //===============================================================================

  //===============================================================================
  /// Stores information about every MC event
  std::vector<EventInfo> MCEvents;
  /// Stores info about currently initialised sample
  std::vector<SampleInfo> SampleDetails;
  //===============================================================================

  //===============================================================================
  //DB Covariance Objects
  /// ETA - All experiments will need an xsec, det and osc cov
  ParameterHandlerGeneric *ParHandler = nullptr;

  /// HW The likelihood gradient, null unless SetComputeGradients(true). The reweight records each
  /// event's bin and weight in it, and GetGradLikelihood() runs its pass over the events
  std::unique_ptr<LikelihoodGradient> Gradient;
  /// HW Whether reweights currently leave the gradient out, see SetGradientsPaused()
  bool GradientsPaused = false;
  /// HW Whether Gradient holds the events of the last reweight, false after a paused one
  bool GradientCurrent = false;
  /// HW Per parameter, whether GetGradLikelihood() carries its derivative, see GetGradientCoverage()
  std::vector<bool> GradientCoverage;
  /// HW Derivatives of the oscillation weights, null unless gradients are on and the oscillation
  /// calculation is one OscillationGradient can differentiate
  std::unique_ptr<OscillationGradient> OscGradient;
  /// @brief Set up OscGradient, if the oscillation calculation allows it
  void SetupOscillationGradient();

  /// @name Reweight backend
  /// HW See ReweightBackend. Null means the plain event loop in FillArray, which is also the
  /// reference CheckReweightBackend() compares against
  /// @{
  /// Lets a backend shift, cut and bin this handler's events
  class BackendEventAccess final : public ReweightEventAccess {
   public:
    /// @brief Serve the given handler
    explicit BackendEventAccess(SampleHandlerBase& Handler_) : Handler(Handler_) {}
    /// @copydoc ReweightEventAccess::ShiftEvent
    void ShiftEvent(const int iEvent) override { Handler.ApplyShifts(iEvent); }
    /// @copydoc ReweightEventAccess::PassesCuts
    bool PassesCuts(const int iEvent) override {
      return Handler.IsEventSelected(Handler.MCEvents[static_cast<std::size_t>(iEvent)].NominalSample, iEvent);
    }
    /// @copydoc ReweightEventAccess::GlobalBinOf
    int GlobalBinOf(const int iEvent) override {
      const EventInfo& Event = Handler.MCEvents[static_cast<std::size_t>(iEvent)];
      return Handler.Binning->FindGlobalBin(Event.NominalSample, Event.KinVar, Event.NomBin);
    }
   private:
    /// The handler served
    SampleHandlerBase& Handler;
  };
  /// The backend in use, null for the plain event loop
  std::unique_ptr<ReweightBackend> Backend;
  /// @brief The surrogate for moves of the shifts alone, see SampleHandlerInterface::ReweightSurrogate
  struct SurrogateState {
    /// Whether EnableSurrogate succeeded
    bool On = false;
    /// Whether the anchor's W2 moves with the events (UpdateW2), or stays as it was filled
    bool MoveW2 = false;
    /// The prediction at the anchor
    std::vector<double> AnchorMC;
    /// W2 at the anchor
    std::vector<double> AnchorW2;
    /// The sampled migrations of the last ReweightSurrogate()
    std::vector<double> DeltaMC;
    /// Their squared weights
    std::vector<double> DeltaW2;
  } Surrogate;
  /// What Backend calls back into
  std::unique_ptr<BackendEventAccess> BackendAccess;
  /// The faster access Backend uses when it could be set up and checked, see FastEventAccess
  std::unique_ptr<FastEventAccess> FastAccess;
  /// @brief Set up FastAccess, using each of its shortcuts only if it agrees with the generic access
  /// for every live event
  void SetupFastEventAccess(const std::vector<char>& NeverSelected);
  /// @brief Log how much memory the reweight's own structures hold
  void ReportReweightMemory() const;
  /// @brief Measure the linear shift table and check it against the experiment's shifts, see LinearShiftTable
  void SetupLinearShifts(const std::vector<char>& NeverSelected);
  /// Whether the events' shifted variables may be behind the last reweight, see SyncShiftedKinematics
  bool KinematicsStale = false;
  /// Set while CheckReweightBackend() runs the plain loop for comparison
  bool UseLegacyLoop = false;
  /// Per event, where the weights AddAdditionalWeightPointers() added start and end in total_weight_pointers
  std::vector<std::uint16_t> ConstWeightStart;
  /// @brief Per event, how many weight pointers it has so far
  std::vector<std::uint16_t> WeightPointerCounts() const;
  /// See ConstWeightStart
  std::vector<std::uint16_t> ConstWeightEnd;
  /// The parameters the oscillation and spline weights were last evaluated at, see EvaluateWeights()
  WeightParameterSnapshot WeightSnapshot;
  /// @brief Choose and build the backend, as the sample config's Reweight section asks
  void SetupReweightBackend();
  /// @brief Compare the new backend with the plain loop at a few points: throw if it differs and was
  /// asked for by name, fall back to the plain loop if it was chosen automatically
  void SelfCheckReweightBackend(const ReweightOptions& Options);
  /// @brief Per event, non-zero if a cut on a variable no shift changes rejects it
  std::vector<char> FindNeverSelected();
  /// @brief The current values of every functional parameter, as the shifts see them
  std::vector<double> CurrentShiftKey() const;
  /// @brief Evaluate the oscillation and spline weights, unless with a backend in use no parameter
  /// they can depend on has moved since they were last evaluated
  void EvaluateWeights();
  /// @brief The plain event loop
  void FillArrayLegacy();
  /// @}

  //===============================================================================
  /// @brief Identifier of this Sample Handler, mostly used for fancy printing in FitterBase
  std::string SampleHandlerName;

  //===========================================================================
  /// @brief DB Vectors to store which kinematic cuts we apply.
  /// Gets used in IsEventSelected
  /// Read in from sample yaml file
  std::vector< std::vector< KinematicCut > > StoredSelection;
  /// @brief a way to store selection cuts which you may push back in the get1DVar functions
  /// most of the time this is just the same as StoredSelection
  std::vector< std::vector< KinematicCut > > Selection;
   //===========================================================================

  /// Mapping between string and kinematic enum
  const std::unordered_map<std::string, int>* KinematicParameters;
  /// Mapping between kinematic enum and string
  const std::unordered_map<int, std::string>* ReversedKinematicParameters;

  // === JM mapping between string and kinematic vector enum ===
  const std::unordered_map<std::string, int>* KinematicVectors;
  const std::unordered_map<int, std::string>* ReversedKinematicVectors;
  // ===========================================================

  /// The manager object used to read the sample yaml file
  std::unique_ptr<Manager> SampleManager;

  std::unordered_map<std::string, double> _modeNomWeightMap;

  //===============================================================================
  /// DB: Legend associated with stacked histograms produced by this class.
  TLegend* THStackLeg = nullptr;

  /// KS:Super hacky to update W2 or not
  bool FirstTimeW2;
  /// KS:Super hacky to update W2 or not
  bool UpdateW2;

  /// @brief Retrieve the initial neutrino PDG code associated with a given input file name.
  NuPDG GetInitPDGFromFileName(const std::string& FileName) const {return FileToInitPDGMap.at(FileName);}
  /// @brief Retrieve the final neutrino PDG code associated with a given input file name.
  NuPDG GetFinalPDGFromFileName(const std::string& FileName) const {return FileToFinalPDGMap.at(FileName);}

 private:
   /// @brief Temporarily extend Selection for a given sample with additional cuts.
   /// Returns the original Selection so the caller can restore it later.
   std::vector<std::vector<KinematicCut>> ApplyTemporarySelection(const int iSample,
                                                                  const std::vector<KinematicCut>& ExtraCuts);
  /// Mapping from input file names to initial neutrino PDG codes.
  std::unordered_map<std::string, NuPDG> FileToInitPDGMap;
  /// Mapping from input file names to final neutrino PDG codes.
  std::unordered_map<std::string, NuPDG> FileToFinalPDGMap;

  /// @brief Helper object for storing/updating information related to functional shift parameters
  M3::detail::Functional functional;

  /// Number of MC events are there
  unsigned int nEvents;
};

//template implementation details live here
#include "Samples/SampleHandlerBase.txx"
