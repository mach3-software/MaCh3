#pragma once

//MaCh3 includes
#include "Samples/SampleInfo.h"
#include "Samples/EventInfo.h"

#include "Splines/BinnedSplineHandler.h"
#include "Samples/DeviceReweight.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerInterface.h"
#include "Samples/OscillationHandler.h"
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
  /// @copydoc SampleHandlerInterface::GetGradLikelihood
  std::vector<double> GetGradLikelihood() const override;
  /// @copydoc SampleHandlerInterface::SetComputeGradients
  void SetComputeGradients(const bool Enable) override;
  /// @copydoc SampleHandlerInterface::SetGradientsPaused
  void SetGradientsPaused(const bool Paused) override { GradientsPaused = Paused; }
  /// @copydoc SampleHandlerInterface::GetGradParHandler
  const ParameterHandlerBase* GetGradParHandler() const override { return ParHandler; }
  /// @copydoc SampleHandlerInterface::GetSampleLikelihood
  double GetSampleLikelihood(const int iSample) const override;
  //===============================================================================

  /// @brief Route every reweight through the plain reference loop instead of the optimised paths
  /// @details The reference loop is the event loop as it was before any of the reweight
  /// optimisations: every event is shifted, selected, binned and weighted on every reweight, the
  /// oscillator and splines are always evaluated, and the gradient is filled with atomics. It
  /// reads none of the caches and leaves them as they were, so the optimised paths carry on
  /// where they left off once it is switched off again. Much slower, for validation only
  void SetReferencePath(const bool Enable) { ReferencePath = Enable; }
  /// @brief Check the optimised reweight against the reference loop, see SetReferencePath
  /// @details Moves the parameters through NPoints points, cycling between moving only the
  /// functional (shift) parameters, only the others, and both, which is how the optimised code
  /// goes through each of its paths: a full pass with event weights reused, a cached pass over
  /// grouped events, and a full pass with everything recomputed. Flat-prior parameters are
  /// thrown uniformly in their bounds, the rest by Gaussian about the prefit value with
  /// Width times their prior error, clipped to their bounds. At every point it reweights through
  /// the optimised path and then the reference one and compares the MC prediction and the W2
  /// array bin by bin, relative to each bin, the gradient array relative to its largest entry,
  /// if gradients are on, and the likelihood. The parameters are put back and the handler
  /// reweighted at them before returning
  /// @return The largest discrepancy found, in the relative terms above. Anything above
  /// Tolerance is also logged as an error with the point, bin and values
  double CheckAgainstReference(const int NPoints = 9, const double Tolerance = 1e-9,
                               const double Width = 1., const unsigned int Seed = 2026);

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
  /// @brief Return flat array of d(bin content)/d(parameter) for every bin
  /// @details Indexed as [iBin * GetNGradPars() + iPar], empty unless SetComputeGradients(true)
  const std::vector<double>& GetGradArray() const {
    return SampleHandler_gradients;
  }
  /// @brief Return the number of parameters GetGradArray() is indexed by, 0 when gradients are off
  int GetNGradPars() const {
    return nGradPars;
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
  /// @brief Promise that every weight AddAdditionalWeightPointers registers is constant for the
  /// whole fit, and that CalcWeightFunc does nothing
  /// @details With that promise, events which share a bin and every other weight pointer can be
  /// summed into one before the event loop, see SetupEventGrouping. Off unless a handler opts in
  virtual bool AdditionalWeightsAreConstant() const { return false; }
  /// @brief Whether any functional shift (or ResetShifts / FinaliseShifts) can change this
  /// kinematic variable
  /// @details A selection cut on a variable no shift touches gives the same answer whatever the
  /// shifts do, so events failing one can never be selected and the full pass skips them outright,
  /// shifts included. Every variable counts as shifted unless a handler says otherwise
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

  /// @brief Accumulate one event's weight derivatives into the flat gradient array
  /// @param MCEvent Event whose weight derivatives we want
  /// @param TotalWeight The event's total weight, as returned by CalcWeightTotal
  /// @param GlobalBin Global bin the event falls into
  /// @param GradArray Flat NBins * nGradPars array to accumulate into, see SampleHandler_gradients
  /// @details For \f$W = \prod_k w_k\f$ the chain rule gives
  /// \f$\partial W/\partial x_p = W \sum_k w_k^{-1} \partial w_k/\partial x_p\f$. Only call this with
  /// TotalWeight > 0, which guarantees no \f$w_k\f$ is zero and so the division is safe.
  void FillGradientArray(const EventInfo* _restrict_ MCEvent, const double TotalWeight,
                         const int GlobalBin, double* _restrict_ GradArray) const _noexcept_;
  /// @brief As FillGradientArray, into one bin's row, for a caller which is the only writer of it
  /// @param MCEvent Event whose weight derivatives we want
  /// @param TotalWeight The event's total weight, as returned by CalcWeightTotal, > 0
  /// @param GradRow The nGradPars entries of SampleHandler_gradients for the event's bin
  void FillGradientRow(const EventInfo* _restrict_ MCEvent, const double TotalWeight,
                       double* _restrict_ GradRow) const _noexcept_;

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
  /// HW Flat NBins * nGradPars array of d(bin content)/d(parameter), indexed as
  /// [iBin * nGradPars + iPar] so the parameters of one bin are contiguous.
  /// Empty unless SetComputeGradients(true) was called
  std::vector<double> SampleHandler_gradients;
  /// HW Number of parameters SampleHandler_gradients is indexed by, 0 when gradients are off
  int nGradPars = 0;
  /// HW Whether to fill SampleHandler_gradients during the reweight, see SetComputeGradients
  bool ComputeGradients = false;
  /// HW Whether reweights currently leave the gradient out, see SetGradientsPaused
  bool GradientsPaused = false;
  /// HW Whether reweights go through the plain reference loop, see SetReferencePath
  bool ReferencePath = false;
  /// @brief The reference loop FillArray runs instead of the optimised paths with ReferencePath set
  void FillArrayReference();
  /// HW Whether SampleHandler_gradients belongs to the last reweight, false after a paused one
  bool GradientsCurrent = false;

  /// @name Kinematic cache
  /// HW Which events pass the selection, and which bin each lands in, depend only on the shifted
  /// kinematics, and those depend only on the functional parameters. While they stay put (every
  /// leapfrog of a NUTS trajectory, any fit holding them fixed) both answers are the same on every
  /// reweight, so FillArray works them out once and reuses them, see KinematicCacheMatches
  /// @warning Assumes a shift depends on nothing but its own functional parameter values and the
  /// nominal event. A derived PrepFunctionalParameters which feeds the shifts anything else breaks
  /// that, and has to call InvalidateKinematicCache whenever what it feeds them changes
  /// @{
  /// Value of every functional parameter, in registration order, when EventBinCache was filled
  std::vector<double> KinCacheKey;
  /// Whether EventBinCache holds a complete pass of FillArray at KinCacheKey
  bool KinCacheValid = false;
  /// Global bin of every event at KinCacheKey, kKinCacheNotSelected if the selection rejects it
  std::vector<int> EventBinCache;
  /// Events which pass the selection and land inside the binning, grouped by global bin, built
  /// from EventBinCache the first time it is reused so a reweight which never repeats pays nothing
  std::vector<int> CachedEvents;
  /// Where each global bin's events start in CachedEvents, one entry per bin plus one at the end
  std::vector<int> CachedBinStart;
  /// Global bins with at least one cached event, most events first, the order threads take them in
  std::vector<int> CachedBinOrder;
  /// Whether CachedEvents, CachedBinStart and CachedBinOrder match EventBinCache
  bool CachedListReady = false;
  /// Marker in EventBinCache for an event the selection rejects; UnderOverFlowBin means selected
  /// but outside the binning
  static constexpr int kKinCacheNotSelected = M3::UnderOverFlowBin - 1;
  /// @brief Whether the cache was filled at the current functional parameter values
  /// @details Compares exactly, there is no tolerance. On a mismatch it records the current
  /// values as the new key and marks the cache as invalid until the next full pass fills it
  bool KinematicCacheMatches();
  /// @brief Build CachedEvents, CachedBinStart and CachedBinOrder from EventBinCache
  void BuildCachedEventList();
  /// @brief Forget the cache, so the next FillArray does a full pass
  void InvalidateKinematicCache() { KinCacheValid = false; CachedListReady = false; }

  /// HW Event grouping. Inside one bin, events whose varying weights are the very same pointers
  /// differ only by their constant weights, so their sum is (sum of constant weights) times one
  /// product of the shared weights. The cached loop then runs over groups rather than events,
  /// typically an order of magnitude fewer. Only with AdditionalWeightsAreConstant()
  /// Where each event's AddAdditionalWeightPointers entries start and end in total_weight_pointers
  std::vector<std::uint16_t> AdditionalWeightStart;
  /// See AdditionalWeightStart
  std::vector<std::uint16_t> AdditionalWeightEnd;
  /// Whether SetupEventGrouping has run since the weight lists last changed
  bool GroupingChecked = false;
  /// Whether the cached loop runs over groups
  bool GroupingEnabled = false;
  /// Identifier of each event's list of varying weights, the same for events whose lists match
  std::vector<int> EventSignature;
  /// Product of each event's constant weights
  std::vector<double> EventConstWeight;
  /// Every event, sorted by signature, so the bins CachedEvents is built into come out sorted too
  std::vector<int> EventsBySignature;
  /// One event standing for each group, whose varying weights the group shares
  std::vector<int> GroupRep;
  /// Signature of each group, indexing the flat lists below
  std::vector<int> GroupSig;
  /// Per signature, its shared weight pointers, flat: SigWeightStart[s] to SigWeightStart[s+1]-1
  std::vector<const M3::float_t*> SigWeightPtr;
  /// See SigWeightPtr
  std::vector<int> SigWeightStart;
  /// Per signature, its weights with a derivative: parameter, derivative and weight pointers, flat
  struct SigDeriv { int Par; const M3::float_t* Deriv; const M3::float_t* Weight; };
  /// See SigDeriv
  std::vector<SigDeriv> SigDerivs;
  /// See SigWeightStart
  std::vector<int> SigDerivStart;
  /// Sum over a group of the constant weights, and of their squares for the W2 array
  std::vector<double> GroupConst;
  /// See GroupConst
  std::vector<double> GroupConst2;
  /// Where each global bin's groups start in GroupRep, one entry per bin plus one at the end
  std::vector<int> BinGroupStart;
  /// @brief Work out signatures once, and whether grouping can be used at all
  void SetupEventGrouping();
  /// @name Bookkeeping the device reweight relies on, kept in every build so the class is the same
  /// @{
  /// Bumped whenever SetupEventGrouping rebuilds the signatures, and BuildEventGroups the groups,
  /// so a device copy knows when it is stale
  unsigned int SignatureVersion = 0;
  /// See SignatureVersion
  unsigned int GroupsVersion = 0;
  /// Bumped by AddData, see SignatureVersion
  unsigned int DataVersion = 0;
  /// Whether the last reweight's gradient came from the device, see DeviceFill
  bool GradFromDevice = false;
  /// Whether ResetHistograms left the per bin gradient array for the CPU pass to wipe
  bool GradNeedsZero = false;
  /// @brief Whether ResetHistograms should leave the gradient array alone, see GradNeedsZero
  bool GradResetDeferred() const {
    #ifdef MaCh3_CUDA
    return DeviceUsable;
    #else
    return false;
    #endif
  }
  /// @brief Wipe the per bin gradient array if ResetHistograms left it, before a CPU pass fills it
  void ZeroGradientsIfDeferred() {
    if (!GradNeedsZero) return;
    std::fill(SampleHandler_gradients.begin(), SampleHandler_gradients.end(), 0.0);
    GradNeedsZero = false;
  }
  /// @}

  #ifdef MaCh3_CUDA
  /// @name GPU device reweight, see Samples/DeviceReweight.h
  /// Prototype: only in a CUDA build, and only with MACH3_DEVICE_REWEIGHT=1 in the environment
  /// (M3::DeviceReweightRequested). Each part falls back to the CPU path when it doesn't fit a handler
  /// @{
  /// The device copy, built on the first cached reweight once grouping is set up
  std::unique_ptr<DeviceReweight> Device;
  /// Whether SetupDevice has run for the current signatures, and whether it succeeded
  bool DeviceChecked = false;
  /// See DeviceChecked
  bool DeviceUsable = false;
  /// SignatureVersion and GroupsVersion the device copy was built from
  unsigned int DeviceSignatureVersion = 0;
  /// See DeviceSignatureVersion
  unsigned int DeviceGroupsVersion = 0;
  /// See DeviceSignatureVersion
  unsigned int DeviceDataVersion = 0;
  /// Weights the host gathers into the pool each reweight (oscillation and anything else)
  std::vector<const M3::float_t*> DeviceGatheredPtrs;
  /// Staging for DeviceGatheredPtrs' values, the parameter values and the gradient
  std::vector<double> DeviceGatheredVals;
  /// See DeviceGatheredVals
  std::vector<double> DeviceParVals;
  /// See DeviceGatheredVals
  std::vector<double> DeviceGrad;
  /// @brief Build the device copy from the signatures, see DeviceReweight
  void SetupDevice();
  /// @brief Bring the device copy up to date (grouping, setup, data). Returns whether it can be used
  bool DevicePrepare();
  /// @brief Copy the parameter values, gathered weights and shift values into the staging vectors, and
  /// find the spline segments. Returns the segments, nullptr without splines
  const short int* StageDeviceInputs();
  /// @brief The cached reweight on the device, if it is usable. Returns whether it ran
  bool DeviceFill(const bool Grads);
  /// Whether full passes (the random walk sweeps) run on the device too, see SetupDeviceSweep
  bool DeviceSweepUsable = false;
  /// Whether EventBinCache's live entries are on the device, from the last device full pass, see PullEventBins
  bool EventBinsOnDevice = false;
  /// Whether the device has built its groups from the last device full pass, see DeviceFill
  bool DeviceGroupsFromFull = false;
  /// Live event indices of the device full pass, to put its bins back into EventBinCache
  std::vector<int> DeviceLiveEvents;
  /// Shift parameter values, staged for the device
  std::vector<double> DeviceShiftVals;
  /// @brief Turn every shift into per event coefficients, check them, and upload. Leaves
  /// DeviceSweepUsable false (full passes stay on the CPU) if anything doesn't fit, see DeviceSweepStatic
  void SetupDeviceSweep();
  /// Most variables (binning variable plus shifted cuts) an event can have on the device, see DeviceReweight.cu
  static constexpr int kDeviceMaxVariables = 8;
  /// @brief Fill the per sample parts of Sweep (shifted cuts, bin edges) and, per sample, the indices
  /// into StoredSelection of the cuts on shifted variables
  void DescribeSweepSamples(DeviceSweepStatic& Sweep, std::vector<std::vector<int>>& ShiftedCuts);
  /// @brief Host pointers to an event's device variables: the binning variable, then its shifted cuts' variables
  void GetSweepVariables(const int iEvent, const std::vector<std::vector<int>>& ShiftedCuts, const int NVariables,
                         std::vector<double*>& Variables) const;
  /// @brief Find each live event's nominal variables and the coefficient of every shift on them
  void ProbeShiftCoefficients(DeviceSweepStatic& Sweep, const std::vector<std::vector<int>>& ShiftedCuts);
  /// @brief Compare the coefficients' prediction with ApplyShifts at random shift values. Returns the
  /// worst relative mismatch and counts the values compared in NChecked
  double CheckShiftLinearity(const DeviceSweepStatic& Sweep, const std::vector<std::vector<int>>& ShiftedCuts, long& NChecked);
  /// @brief A full pass on the device, if it is usable. Returns whether it ran
  bool DeviceFullPass(const bool Grads);
  /// @brief Copy the live events' bins of the last device full pass into EventBinCache
  void PullEventBins();
  /// @}
  #endif

  /// HW Faster full pass, for the reweights where the shifts have moved (random walk sweeps).
  /// Whether ApplyShifts can use functional.apply_all, see UpdateFastShifts
  bool FastShifts = false;
  /// @brief Decide FastShifts: every registered shift has a plain function on one event vector
  void UpdateFastShifts();
  /// Per event, pointers to the variables of its sample's selection cuts, in cut order, so the
  /// full pass reads them directly rather than through ReturnKinematicParameter
  std::vector<const double*> SelectionPtr;
  /// Where each event's entries start in SelectionPtr, one per event plus one at the end
  std::vector<int> SelectionPtrStart;
  /// Whether SetupSelectionPointers has run since the cuts last changed, and whether it succeeded
  bool SelectionPtrChecked = false;
  /// See SelectionPtrChecked
  bool SelectionPtrReady = false;
  /// Events which fail a cut on a variable no shift touches, so can never be selected. Built with
  /// SelectionPtr, empty when nothing can be ruled out
  std::vector<char> NeverSelected;
  /// The other events, which the full pass loops over when NeverSelected is in use
  std::vector<int> LiveEvents;
  /// EventBinCache's storage when the never selected events were last marked in it
  const int* DeadBinsWritten = nullptr;
  /// @brief Ask the experiment once per event and cut for a pointer to the cut variable, and check
  /// each one reads the same value ReturnKinematicParameter returns
  void SetupSelectionPointers();
  /// @brief IsEventSelected through SelectionPtr, for the stored cuts
  bool IsEventSelectedFast(const int iSample, const int iEvent) const _noexcept_ {
    const auto& Cuts = StoredSelection[iSample];
    const double* const* Ptr = SelectionPtr.data() + SelectionPtrStart[iEvent];
    const int NCuts = static_cast<int>(Cuts.size());
    for (int iCut = 0; iCut < NCuts; ++iCut) {
      const double Val = *Ptr[iCut];
      if ((Val < Cuts[iCut].LowerBound) || (Val >= Cuts[iCut].UpperBound)) return false;
    }
    return true;
  }
  /// Weight of every event, kept across reweights in which no weight bearing parameter moved, so
  /// the random walk sweeps, which only move shifts, reuse them. Only with AdditionalWeightsAreConstant()
  std::vector<M3::float_t> EventWeightCache;
  /// WeightCacheVersion at which each EventWeightCache entry was filled
  std::vector<unsigned int> EventWeightStamp;
  /// Bumped whenever a weight bearing parameter moves, invalidating every entry at once
  unsigned int WeightCacheVersion = 1;
  /// Values of the weight bearing parameters, every non functional one, at WeightCacheVersion
  std::vector<double> WeightCacheKey;
  /// Indices of the weight bearing parameters in ParHandler
  std::vector<int> WeightCachePars;
  /// Whether the full pass reads and fills EventWeightCache
  bool WeightCacheEnabled = false;
  /// WeightCacheVersion at which the oscillator and splines were last evaluated, 0 for never
  unsigned int WeightsEvaluatedVersion = 0;
  /// @brief Compare the weight bearing parameters with WeightCacheKey, bump the version if any moved
  void UpdateWeightCache();

  /// HW Bin lookup for the full pass without the per-event vectors FindGlobalBin goes through.
  /// Only when every sample has uniform 1D binning; the algorithm is FindBin's, unchanged
  struct FastAxis {
    /// Bin edges, NBins + 1 of them
    const double* Edges = nullptr;
    /// FindBin's neighbour lookup, one per bin
    const BinShiftLookup* Lookup = nullptr;
    /// Number of bins
    int NBins = 0;
    /// Stride and global offset, as FindGlobalBin applies them
    int Stride = 1;
    /// See Stride
    int Offset = 0;
  };
  /// One per sample
  std::vector<FastAxis> FastAxes;
  /// HW Every event's NominalSample, in a plain array so the full pass needn't touch MCEvents
  /// for events whose weight is cached
  std::vector<int> EventSampleIdx;
  /// Per event, the binning variable FindGlobalBin would read
  std::vector<const double*> FastBinVar;
  /// Per event, its nominal bin
  std::vector<int> FastNomBin;
  /// Whether SetupFastBinning has run, and whether it can be used
  bool FastBinChecked = false;
  /// See FastBinChecked
  bool FastBinReady = false;
  /// @brief Build the flat lookup and check it against FindGlobalBin for every event
  void SetupFastBinning();
  /// @brief FindGlobalBin through the flat lookup
  int FindGlobalBinFast(const int iSample, const int iEvent) const _noexcept_ {
    const FastAxis& Axis = FastAxes[iSample];
    const double Var = *FastBinVar[iEvent];
    if (Var < Axis.Edges[0] || Var >= Axis.Edges[Axis.NBins]) return M3::UnderOverFlowBin;
    const int NomBin = FastNomBin[iEvent];
    int Bin;
    if (NomBin > M3::UnderOverFlowBin && Var < Axis.Lookup[NomBin].upper_binedge && Var >= Axis.Lookup[NomBin].lower_binedge) {
      Bin = NomBin;
    } else if (NomBin > M3::UnderOverFlowBin && Var < Axis.Lookup[NomBin].lower_binedge && Var >= Axis.Lookup[NomBin].lower_lower_binedge) {
      Bin = NomBin - 1;
    } else if (NomBin > M3::UnderOverFlowBin && Var < Axis.Lookup[NomBin].upper_upper_binedge && Var >= Axis.Lookup[NomBin].upper_binedge) {
      Bin = NomBin + 1;
    } else {
      Bin = static_cast<int>(std::upper_bound(Axis.Edges, Axis.Edges + Axis.NBins + 1, Var) - Axis.Edges) - 1;
    }
    return Bin * Axis.Stride + Axis.Offset;
  }
  /// @brief Merge the events of each bin of CachedEvents into groups
  void BuildEventGroups();
  /// @}
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
