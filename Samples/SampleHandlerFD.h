#pragma once

//MaCh3 includes
#include "Splines/BinnedSplineHandler.h"
#include "Parameters/ParameterHandlerGeneric.h"
#include "Samples/SampleHandlerBase.h"
#include "Samples/OscillationHandler.h"
#include "Samples/FarDetectorCoreInfoStruct.h"
#include "Samples/BinningHandler.h"

_MaCh3_Safe_Include_Start_ //{
#include "THStack.h"
#include "TLegend.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Class responsible for handling implementation of samples used in analysis, reweighting and returning LLH
/// @author Dan Barrow
/// @author Ed Atkin
///
/// @ingroup SamplesAndParameters
class SampleHandlerFD :  public SampleHandlerBase
{
 public:
  //######################################### Functions #########################################
  /// @brief Constructor
  /// @param ConfigFileName Name of config to initialise the sample object
  SampleHandlerFD(std::string ConfigFileName, ParameterHandlerGeneric* xsec_cov,
                  const std::shared_ptr<OscillationHandler>& OscillatorObj_ = nullptr);
  /// @brief destructor
  virtual ~SampleHandlerFD();

  /// @brief DB Function to differentiate 1D or 2D binning
  int GetNDim(const int Sample) const override { return SampleDetails[Sample].nDimensions; }
  std::string GetName() const override;
  std::string GetSampleTitle(const int Sample) const override {return SampleDetails[Sample].SampleTitle;}

  /// @brief Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
  /// @param iSample Sample index
  /// @param Dimension Dimension index
  std::string GetKinVarName(const int iSample, const int Dimension) const override;

  std::string GetXBinVarName(const int Sample) const {return GetKinVarName(Sample, 0);}
  std::string GetYBinVarName(const int Sample) const {return GetKinVarName(Sample, 1);}
  /// @brief Get pointer to binning handler
  const BinningHandler* GetBinningHandler() const {return Binning.get();}

  void PrintIntegral(const int iSample, const TString& OutputName="/dev/null", const int WeightStyle=0, const TString& OutputCSVName="/dev/null");

  //===============================================================================
  // DB Reweighting and Likelihood functions

  //ETA - abstract these to SampleHandlerFDBase
  //DB Require these four functions to allow conversion from TH1(2)D to array for multi-threaded GetLikelihood
  void AddData(const int Sample, TH1* Data);
  void AddData(const int Sample, const std::vector<double>& Data_Array);

  /// @brief Helper function to print rates for the samples with LLH
  /// @param DataOnly whether to print data only rates
  void PrintRates(const bool DataOnly = false) override;
  /// @brief DB Multi-threaded GetLikelihood
  double GetLikelihood() const override;
  /// @brief Get likelihood for single sample
  double GetSampleLikelihood(const int isample) const override;
  //===============================================================================

  /// @brief Get index of sample based on name
  /// @param SampleTitle The title of the sample to search for.
  int GetSampleIndex(const std::string& SampleTitle) const;

  /// @brief Get Data histogram
  TH1* GetDataHist(const int Sample) override;
  TH1* GetDataHist(const std::string& Sample);

  /// @brief Get MC histogram
  TH1* GetMCHist(const int Sample) override;
  TH1* GetMCHist(const std::string& Sample);

  /// @brief Get W2 histogram
  TH1* GetW2Hist(const int Sample) override;
  TH1* GetW2Hist(const std::string& Sample);

  void Reweight() override;
  M3::float_t GetEventWeight(const int iEntry);

  /// @brief including Dan's magic NuOscillator
  void InitialiseNuOscillatorObjects();
  void SetupNuOscillatorPointers();
  const M3::float_t* GetNuOscillatorPointers(const int iEvent) const;

  void ReadConfig();
  void LoadSingleSample(const int iSample, const YAML::Node& Settings);

  int GetNOscChannels(const int iSample) const override {return static_cast<int>(SampleDetails[iSample].OscChannels.size());};

  std::string GetFlavourName(const int iSample, const int iChannel) const override {
    if (iChannel < 0 || iChannel > GetNOscChannels(iSample)) {
      MACH3LOG_ERROR("Invalid Channel Requested: {}", iChannel);
      throw MaCh3Exception(__FILE__ , __LINE__);
    }
    return SampleDetails[iSample].OscChannels[iChannel].flavourName;
  }
  /// @brief Temporarily extend Selection for a given sample with additional cuts.
  /// Returns the original Selection so the caller can restore it later.
  std::vector<std::vector<KinematicCut>> ApplyTemporarySelection(const int iSample,
                                                                 const std::vector<KinematicCut>& ExtraCuts);
  TH1 *Get1DVarHist(const int iSample, const std::string &ProjectionVar,
                    const std::vector<KinematicCut> &EventSelectionVec = {}, int WeightStyle = 0,
                    TAxis *Axis = nullptr, const std::vector<KinematicCut> &SubEventSelectionVec = {}) override;
  TH2* Get2DVarHist(const int iSample, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                    const std::vector< KinematicCut >& EventSelectionVec = {},
                    int WeightStyle = 0, TAxis* AxisX = nullptr, TAxis* AxisY = nullptr,
                    const std::vector< KinematicCut >& SubEventSelectionVec = {}) override;
  std::vector<KinematicCut> BuildModeChannelSelection(const int iSample, const int kModeToFill, const int kChannelToFill) const;

  void Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {},
                          int WeightStyle=0);
  void Fill2DSubEventHist(const int iSample, TH2D* _h2DVar, const std::string& ProjectionVarX, const std::string& ProjectionVarY,
                          const std::vector< KinematicCut >& SubEventSelectionVec = {}, int WeightStyle = 0);

  TH1* Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                    int kModeToFill = -1, int kChannelToFill = -1, int WeightStyle = 0, TAxis* Axis = nullptr) override;
  TH2* Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                    const std::string& ProjectionVar_StrY, int kModeToFill = -1,
                                    int kChannelToFill = -1, int WeightStyle = 0,
                                    TAxis* AxisX = nullptr, TAxis* AxisY = nullptr) override;

  TH1 *GetModeHist1D(const int iSample, int s, int m, int style = 0) {
    return Get1DVarHistByModeAndChannel(iSample, GetXBinVarName(iSample), m, s, style);
  }
  TH2 *GetModeHist2D(const int iSample, int s, int m, int style = 0) {
    return Get2DVarHistByModeAndChannel(iSample, GetXBinVarName(iSample),GetYBinVarName(iSample), m, s, style);
  }

  std::vector<TH1*> ReturnHistsBySelection1D(const int iSample, const std::string& KinematicProjection,
                                             int Selection1, int Selection2 = -1,
                                             int WeightStyle = 0, TAxis* Axis = nullptr);
  std::vector<TH2*> ReturnHistsBySelection2D(const int iSample, const std::string& KinematicProjectionX,
                                             const std::string& KinematicProjectionY,
                                             int Selection1, int Selection2=-1, int WeightStyle=0,
                                             TAxis* XAxis = nullptr, TAxis* YAxis = nullptr);
  THStack* ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection,
                                          int Selection1, int Selection2 = -1, int WeightStyle = 0, TAxis* Axis = nullptr);
  TLegend* ReturnStackHistLegend() {return THStackLeg;}

  /// @brief ETA function to generically convert a string from xsec cov to a kinematic type
  int ReturnKinematicParameterFromString(const std::string& KinematicStr) const;
  /// @brief ETA function to generically convert a kinematic type from xsec cov to a string
  std::string ReturnStringFromKinematicParameter(const int KinematicVariable) const;

  /// @brief Store additional info in a chan
  void SaveAdditionalInfo(TDirectory* Dir) override;

  // === JM declare the same functions for kinematic vectors ===
  int ReturnKinematicVectorFromString(const std::string& KinematicStr) const;
  std::string ReturnStringFromKinematicVector(const int KinematicVariable) const;
  // ===========================================================
  /// @brief JM Check if a kinematic parameter string corresponds to a subevent-level variable
  bool IsSubEventVarString(const std::string& VarStr);

  /// @brief Return array storing data entries for every bin
  auto GetDataArray() const {
    return SampleHandlerFD_data;
  }
  /// @brief Return array storing MC entries for every bin
  auto GetMCArray() const {
    return SampleHandlerFD_array;
  }
  /// @brief Return array storing W2 entries for every bin
  auto GetW2Array() const {
    return SampleHandlerFD_array_w2;
  }
  /// @brief Return a sub-array for a given sample.
  std::vector<double> GetArrayForSample(const int Sample, std::vector<double> const & array) const;

  /// @brief Return array storing data entries for every bin
  std::vector<double> GetDataArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandlerFD_data);
  }
  /// @brief Return array storing MC entries for every bin
  std::vector<double> GetMCArray(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandlerFD_array);
  }
  /// @brief Return array storing W2 entries for single sample
  std::vector<double> GetW2Array(const int Sample) const {
    return GetArrayForSample(Sample, SampleHandlerFD_array_w2);
  }

 protected:
  /// @brief DB Function to determine which weights apply to which types of samples
  virtual void AddAdditionalWeightPointers() = 0;

  /// @brief Ensure Kinematic Map is setup and make sure it is initialised correctly
  void SetupKinematicMap();

  /// @todo abstract the spline initialisation completely to core
  /// @brief initialise your splineXX object and then use InitialiseSplineObject to conviently setup everything up
  virtual void SetupSplines() = 0;

  //DB Require all objects to have a function which reads in the MC
  /// @brief Initialise any variables that your experiment specific SampleHandler needs
  virtual void Init() = 0;

  /// @brief Experiment specific setup, returns the number of events which were loaded
  virtual int SetupExperimentMC() = 0;

  /// @brief Function which translates experiment struct into core struct
  virtual void SetupFDMC() = 0;

  /// @brief Function which does a lot of the lifting regarding the workflow in creating different MC objects
  void Initialise();

  /// @brief Contains all your splines (binned or unbinned) and handles the setup and the returning of weights from spline evaluations
  std::unique_ptr<SplineBase> SplineHandler;

  /// @brief Contains oscillator handling calculating oscillation probabilities
  std::shared_ptr<OscillationHandler> Oscillator;
  //===============================================================================
  /// @brief Set pointers for each event to appropriate weights, for unbinned based on event number
  /// while for binned based on other kinematical properties
  void SetSplinePointers();

  //Functions which find the nominal bin and bin edges
  void FindNominalBinAndEdges();

  /// @brief set the binning for 2D sample used for the likelihood calculation
  void SetBinning();

  /// @brief Initialise data, MC and W2 histograms
  void SetupReweightArrays();
  //===============================================================================

  /// @brief DB Function which determines if an event is selected, where Selection double looks like {{ND280KinematicTypes Var1, douuble LowBound}
  bool IsEventSelected(const int iSample, const int iEvent) _noexcept_;
  /// @brief JM Function which determines if a subevent is selected
  bool IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, unsigned const int iSubEvent, size_t nsubevents);

  // ----- start Functional Parameters -----
  /// @brief HH - a experiment-specific function where the maps to actual functions are set up
  virtual void RegisterFunctionalParameters() = 0;

  std::vector<std::function<void(const int)>> event_shift_functions;
  std::vector<std::vector<int>> events_and_shifts;

  /// @brief HH - a helper function for RegisterFunctionalParameter
  template <typename EventType, typename SFType>
  void RegisterIndividualFunctionalParameter(
      std::vector<EventType> &ExptEvents,
      std::vector<std::string> const &par_names,
      SFType shift_func) {

    static_assert(
        std::is_same_v<std::function<void(std::vector<double const *> const &,
                                          EventType &)>,
                       decltype(std::function(shift_func))>,
        "Function call signature for single parameter Functional shift must be "
        "void(std::vector<double const *> const &, EventType &). -- note the "
        "two consts on the parameter value vector");

    auto sample_func_pars =
        ParHandler->GetFunctionalParametersFromSampleName(SampleHandlerName);

    std::stringstream ss_pars, ss_miss;
    std::vector<FunctionalParameter const *> matched_pars;
    for (auto const &par_name : par_names) {
      ss_pars << par_name << " ";
      bool found = false;
      for (auto const &fp : sample_func_pars) {
        if (fp.name == par_name) {
          matched_pars.push_back(&fp);
          found = true;
        }
      }
      if (!found) {
        ss_miss << par_name;
      }
    }

    // allows experiments to effectively disable functional parameters by not
    // supplying the YAML defining the parameters
    if (!matched_pars.size()) {
      MACH3LOG_INFO(
          "Functional shift consuming parameters: [ {}], doesn't apply to "
          "sample handler: {}",
          ss_pars.str(), SampleHandlerName);
      return;
    } else if (matched_pars.size() !=
               par_names.size()) { // not well defined how to procede with
                                   // partially defined parameter sets, so don't
      MACH3LOG_ERROR(
          "Functional shift consuming parameters: [ {}], only partially "
          "applys to sample handler: {}, missed parameters: [ {}]",
          ss_pars.str(), SampleHandlerName, ss_miss.str());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (!events_and_shifts.size()) {
      events_and_shifts.resize(ExptEvents.size());
    } else if (events_and_shifts.size() != ExptEvents.size()) {
      MACH3LOG_ERROR("When registering functional shift consuming parameters: "
                     "[ {}], SampleHandler: {} has an allocated event map of "
                     "size: {}, but passed a vector of experiment events size: "
                     "{}. SampleHandler must have a unique set of event "
                     "indices so this indicates something has gone wrong.",
                     ss_pars.str(), SampleHandlerName,
                     events_and_shifts.size(), ExptEvents.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (ExptEvents.size() != MCSamples.size()) {
      MACH3LOG_ERROR("When registering functional shift consuming parameters: "
                     "[ {}], SampleHandler: {} knows about {} MCEvents, but "
                     "passed a vector of experiment events size: "
                     "{}. SampleHandler must have a unique set of event "
                     "indices so this indicates something has gone wrong.",
                     ss_pars.str(), SampleHandlerName, MCSamples.size(),
                     ExptEvents.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // This functional shift is correctly configured for this SampleHandler
    std::vector<double const *> par_vals;
    for (auto const &fp : sample_func_pars) {
      par_vals.push_back(ParHandler->RetPointer(fp.index));
    }

    MACH3LOG_INFO("Registered functional shift consuming parameters: "
                  "[ {}] for SampleHandler: {}, with {} par vals.",
                  ss_pars.str(), SampleHandlerName, par_vals.size());

    int iShift = int(event_shift_functions.size());
    event_shift_functions.push_back(
      [shift_func, par_vals, &ExptEvents](int iEvent) {
        shift_func(par_vals, ExptEvents[iEvent]);
      });

    // For each event, make a vector of pointers to the functional parameters
    int NEvents = GetNEvents();
    for (int iEvent = 0; iEvent < NEvents; ++iEvent) {
      // Now loop over the functional parameters and get a vector of enums
      // corresponding to the functional parameters
      int nmatch = 0;
      for (auto const &par : matched_pars) {
        if (!MatchCondition(par->modes, static_cast<int>(std::round(
                                            *(MCSamples[iEvent].mode))))) {
          MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent,
                         *(MCSamples[iEvent].mode), par->name);
          break;
        }
        if (!PassesSelection((*par), iEvent)) {
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}",
                         iEvent, par->name);
          break;
        }
        nmatch++;
      }
      if (!nmatch) {
        continue;
      }
      if (nmatch != matched_pars.size()) {
        MACH3LOG_ERROR(
            "When determining wether functional shift consuming parameters: "
            "[ {}], in SampleHandler: {} should apply to event: {}, only {}/{} "
            "parameters matched. Partially applied shifts are ill-defined.",
            ss_pars.str(), SampleHandlerName, iEvent, nmatch,
            matched_pars.size());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      events_and_shifts[iEvent].push_back(iShift);
    }
  }

  template <typename EventType, typename SFType>
  void RegisterIndividualFunctionalParameter(std::vector<EventType> &ExptEvents,
                                             std::string const &par_name,
                                             SFType shift_func) {

    static_assert(
        std::is_same_v<std::function<void(double const *, EventType &)>,
                       decltype(std::function(shift_func))>,
        "Function call signature for single parameter Functional shift must be "
        "void(double const *, EventType&).");

    RegisterIndividualFunctionalParameter(
        ExptEvents,
        std::vector<std::string>{
            par_name,
        },
        [=](std::vector<double const *> const &par_vals, EventType &ev) {
          shift_func(par_vals[0], ev);
        });
  }
  /// @brief ETA - generic function applying shifts
  void ApplyShifts(const int iEvent);
  /// @brief HH - reset the shifted values to the original values
  virtual void ResetShifts(const int iEvent) {(void)iEvent;};
  /// @brief LP - Optionally calculate derived observables after all shifts have been applied
  /// @details LP - For example, have shifts that varied lepton energy and hadron energy separately
  ///               in a subclass implementation of this method you may add the shifted quantities
  ///               together to build a shifted neutrino energy estimator
  virtual void FinaliseShifts(const int iEvent) {(void)iEvent;};
  // ----- end Functional Parameters -----

  /// @brief Check whether a normalisation systematic affects an event or not
  void CalcNormsBins(std::vector<NormParameter>& norm_parameters, std::vector< std::vector< int > >& norms_bins);
  template <typename ParT>
  bool PassesSelection(const ParT &Par, std::size_t iEvent) {
    // ***************************************************************************
    bool IsSelected = true;
    if (Par.hasKinBounds) {
      const auto &kinVars = Par.KinematicVarStr;
      const auto &selection = Par.Selection;

      for (std::size_t iKinPar = 0; iKinPar < kinVars.size(); ++iKinPar) {
        const double kinVal = ReturnKinematicParameter(
            kinVars[iKinPar], static_cast<int>(iEvent));

        bool passedAnyBound = false;
        const auto &boundsList = selection[iKinPar];

        for (const auto &bounds : boundsList) {
          if (kinVal > bounds[0] && kinVal <= bounds[1]) {
            passedAnyBound = true;
            break;
          }
        }

        if (!passedAnyBound) {
          MACH3LOG_TRACE("Event {}, missed kinematic check ({}) for dial {}",
                         iEvent, kinVars[iKinPar], Par.name);
          IsSelected = false;
          break;
        }
      }
    }
    return IsSelected;
  }

  /// @brief Calculate the total weight weight for a given event
  M3::float_t CalcWeightTotal(const EventInfo* _restrict_ MCEvent) const;

  /// @brief Calculate weights for function parameters
  ///
  /// First you need to setup additional pointers in you experiment code in SetupWeightPointers
  /// Then in this function you can calculate whatever fancy function you want by filling weight to which you have pointer
  /// This way func weight shall be used in GetEventWeight
  virtual void CalcWeightFunc(int iEvent) {return; (void)iEvent;};

  /// @brief Return the value of an associated kinematic parameter for an event
  virtual double ReturnKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual double ReturnKinematicParameter(int KinematicVariable, int iEvent) = 0;

  // === JM declare the same functions for kinematic vectors ===
  virtual std::vector<double> ReturnKinematicVector(std::string KinematicParameter, int iEvent) {return {}; (void)KinematicParameter; (void)iEvent;};
  virtual std::vector<double> ReturnKinematicVector(int KinematicVariable, int iEvent) {return {}; (void)KinematicVariable; (void)iEvent;};
  // ===========================================================

  /// @brief Return the binning used to draw a kinematic parameter
  std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const override;

  virtual const double* GetPointerToKinematicParameter(std::string KinematicParamter, int iEvent) = 0;
  virtual const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) = 0;

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
  /// DB Replace previous implementation with reading bin contents from SampleHandlerFD_array
  void FillHist(const int Sample, TH1* Hist, std::vector<double> &Array);

  /// @brief DB Nice new multi-threaded function which calculates the event weights and fills the relevant bins of an array
#ifdef MULTITHREAD
  /// @brief Function which does the core reweighting, fills the @ref SampleHandlerFD::SampleHandlerFD_array
  /// vector with the weight calculated from reweighting but multithreaded
  void FillArray_MP();
#endif
  /// @brief Function which does the core reweighting, fills the @ref SampleHandlerFD::SampleHandlerFD_array
  /// vector with the weight calculated from reweighting
  void FillArray();

  /// @brief Helper function to reset histograms
  void ResetHistograms();

  //===============================================================================
  //DB Variables required for GetLikelihood
  /// KS: This stores binning information, in future could be come vector to store binning for every used sample
  std::unique_ptr<BinningHandler> Binning;
  /// DB Array to be filled after reweighting
  std::vector<double> SampleHandlerFD_array;
  /// KS Array used for MC stat
  std::vector<double> SampleHandlerFD_array_w2;
  /// DB Array to be filled in AddData
  std::vector<double> SampleHandlerFD_data;
  //===============================================================================

  //===============================================================================
  /// Stores information about every MC event
  std::vector<EventInfo> MCSamples;
  /// Stores info about currently initialised sample
  std::vector<SampleInfo> SampleDetails;
  //===============================================================================

  //===============================================================================
  //DB Covariance Objects
  /// ETA - All experiments will need an xsec, det and osc cov
  ParameterHandlerGeneric *ParHandler = nullptr;

  //===============================================================================
  /// @brief A unique ID for each sample based on which we can define what systematic should be applied
  std::string SampleHandlerName;

  //===========================================================================
  //DB Vectors to store which kinematic cuts we apply
  //like in XsecNorms but for events in sample. Read in from sample yaml file
  //What gets used in IsEventSelected, which gets set equal to user input plus
  //all the vectors in StoreSelection

  /// @brief What gets pulled from config options, these are constant after loading in
  /// this is of length 3: 0th index is the value, 1st is lower bound, 2nd is upper bound
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
  void InitialiseSplineObject();

  std::unordered_map<std::string, double> _modeNomWeightMap;

  //===============================================================================
  /// DB Miscellaneous Variables
  TLegend* THStackLeg = nullptr;
  //===============================================================================

  /// KS:Super hacky to update W2 or not
  bool FirstTimeW2;
  /// KS:Super hacky to update W2 or not
  bool UpdateW2;

  /// @brief Retrieve the initial neutrino PDG code associated with a given input file name.
  NuPDG GetInitPDGFromFileName(const std::string& FileName) const {return FileToInitPDGMap.at(FileName);}
  /// @brief Retrieve the final neutrino PDG code associated with a given input file name.
  NuPDG GetFinalPDGFromFileName(const std::string& FileName) const {return FileToFinalPDGMap.at(FileName);}

 private:
  std::unordered_map<std::string, NuPDG> FileToInitPDGMap;
  std::unordered_map<std::string, NuPDG> FileToFinalPDGMap;

  enum FDPlotType {
    kModePlot = 0,
    kOscChannelPlot = 1
  };
};
