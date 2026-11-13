#include "SampleHandlerBase.h"
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Splines/UnbinnedSplineHandler.h"

#include "TRandom3.h"

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <memory>
#include <numeric>

// ************************************************
SampleHandlerBase::SampleHandlerBase(std::string ConfigFileName, ParameterHandlerGeneric* _ParHandler,
                                     const std::shared_ptr<OscillationHandler>& OscillatorObj_) : SampleHandlerInterface() {
// ************************************************
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Creating SampleHandlerBase object");

  //ETA - safety feature so you can't pass a NULL _ParHandler
  if(!_ParHandler) {
    MACH3LOG_WARN("You've passed me a nullptr ParameterHandler so I will not use any model parameter");
  }
  ParHandler = _ParHandler;
  nEvents = 0;
  nSamples = 0;

  if (OscillatorObj_ != nullptr) {
    MACH3LOG_WARN("You have passed an Oscillator object through the constructor of a SampleHandlerBase object - this will be used for all oscillation channels");
    Oscillator = OscillatorObj_;
    if(!ParHandler) {
      MACH3LOG_CRITICAL("You've passed me a nullptr to ParamHandler while non null to Oscillator");
      MACH3LOG_CRITICAL("Make up you mind");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  KinematicParameters = nullptr;
  ReversedKinematicParameters = nullptr;
  KinematicVectors = nullptr;
  ReversedKinematicVectors = nullptr;

  SampleHandlerName = "";
  SampleManager = std::make_unique<Manager>(ConfigFileName.c_str());
  Binning = std::make_unique<BinningHandler>();
  // Variables related to MC stat
  FirstTimeW2 = true;
  UpdateW2 = false;
}

SampleHandlerBase::~SampleHandlerBase() {
  MACH3LOG_DEBUG("I'm deleting SampleHandlerBase");
  if(THStackLeg != nullptr) delete THStackLeg;
}

// ************************************************
void SampleHandlerBase::ReadConfig() {
// ************************************************
  auto ModeName = Get<std::string>(SampleManager->raw()["MaCh3ModeConfig"], __FILE__ , __LINE__);
  Modes = std::make_unique<MaCh3Modes>(getenv("MACH3")+std::string("/") + ModeName);
  //SampleName has to be provided in the sample yaml otherwise this will throw an exception
  SampleHandlerName = Get<std::string>(SampleManager->raw()["SampleHandlerName"], __FILE__ , __LINE__);

  fTestStatistic = static_cast<TestStatistic>(SampleManager->GetMCStatLLH());
  if (CheckNodeExists(SampleManager->raw(), "LikelihoodOptions")) {
    UpdateW2 = GetFromManager<bool>(SampleManager->raw()["LikelihoodOptions"]["UpdateW2"], false, __FILE__ , __LINE__);
  }

  if (!CheckNodeExists(SampleManager->raw(), "BinningFile")){
    MACH3LOG_ERROR("BinningFile not given in for sample handler {}, ReturnKinematicParameterBinning will not work", SampleHandlerName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto EnabledSasmples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  // Get number of samples and resize relevant objects
  nSamples = static_cast<M3::int_t>(EnabledSasmples.size());
  if (nSamples == 0){
    MACH3LOG_ERROR("No samples for Sample Handler {}, please double check sample config", GetName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SampleDetails.resize(nSamples);
  StoredSelection.resize(nSamples);
  for (int iSample = 0; iSample < nSamples; iSample++)
  {
    auto SampleSettings = SampleManager->raw()[EnabledSasmples[iSample]];
    LoadSingleSample(iSample, SampleSettings);
  } // end loop over enabling samples

  // EM: initialise the mode weight map
  for( int iMode=0; iMode < Modes->GetNModes(); iMode++ ) {
    _modeNomWeightMap[Modes->GetMaCh3ModeName(iMode)] = 1.0;
  }

  // EM: multiply by the nominal weight specified in the sample config file
  if ( SampleManager->raw()["NominalWeights"] ) {
    for( int iMode=0; iMode<Modes->GetNModes(); iMode++ ) {
      std::string modeStr = Modes->GetMaCh3ModeName(iMode);
      if( SampleManager->raw()["NominalWeights"][modeStr] ) {
        double modeWeight = SampleManager->raw()["NominalWeights"][modeStr].as<double>();
        _modeNomWeightMap[Modes->GetMaCh3ModeName(iMode)] *= modeWeight;
      }
    }
  }

  // EM: print em out
  MACH3LOG_INFO("  Nominal mode weights to apply: ");
  for(int iMode=0; iMode<Modes->GetNModes(); iMode++ ) {
    std::string modeStr = Modes->GetMaCh3ModeName(iMode);
    MACH3LOG_INFO("    - {}: {}", modeStr, _modeNomWeightMap.at(modeStr));
  }
}


// ************************************************
void SampleHandlerBase::LoadSingleSample(const int iSample, const YAML::Node& SampleSettings) {
// ************************************************
  SampleInfo SingleSample;
  //SampleTitle has to be provided in the sample yaml otherwise this will throw an exception
  SingleSample.SampleTitle = Get<std::string>(SampleSettings["SampleTitle"], __FILE__ , __LINE__);
  // Sample name used for defying syst uncertainties, if not specified use SampleHandler name
  SingleSample.SampleName = GetFromManager<std::string>(SampleSettings["SampleName"], GetName(), __FILE__ , __LINE__);

  Binning->SetupSampleBinning(SampleSettings["Binning"], SingleSample);

  auto MCFilePrefix  = Get<std::string>(SampleSettings["InputFiles"]["mtupleprefix"], __FILE__, __LINE__);
  auto MCFileSuffix  = Get<std::string>(SampleSettings["InputFiles"]["mtuplesuffix"], __FILE__, __LINE__);
  auto SplinePrefix  = Get<std::string>(SampleSettings["InputFiles"]["splineprefix"], __FILE__, __LINE__);
  auto SplineSuffix  = Get<std::string>(SampleSettings["InputFiles"]["splinesuffix"], __FILE__, __LINE__);

  int NChannels = static_cast<M3::int_t>(SampleSettings["OscChannels"].size());
  SingleSample.OscChannels.reserve(NChannels);

  YAML::Node OscChannelsConfig;
  // KS: We first check whether OscChannel are defined individually for this sample or taken from list
  if(SampleSettings["OscChannels"].IsScalar()) {
    auto PredeterminedChannelsName = Get<std::string>(SampleSettings["OscChannels"], __FILE__, __LINE__);
    if(!SampleManager->raw()["OscChannels"]) {
      MACH3LOG_ERROR("Trying to use Predetermined OscChannels however such field doesn't exist in config for SampleHandler: {}", GetName());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if(!SampleManager->raw()["OscChannels"][PredeterminedChannelsName]) {
      MACH3LOG_ERROR("I didn't find PredeterminedChannelsName called: {}", PredeterminedChannelsName);
      MACH3LOG_ERROR("However I have PredeterminedChannelsName known as:");
      for (const auto& item : SampleManager->raw()["OscChannels"]) {
        MACH3LOG_ERROR("{}", item.first.as<std::string>());
      }
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    OscChannelsConfig = SampleManager->raw()["OscChannels"][PredeterminedChannelsName];
  } else {
    OscChannelsConfig = SampleSettings["OscChannels"];
  }
  int OscChannelCounter = 0;
  for (auto const &osc_channel : OscChannelsConfig) {
    OscChannelInfo OscInfo;
    OscInfo.flavourName       = Get<std::string>(osc_channel["Name"], __FILE__ , __LINE__);
    OscInfo.flavourName_Latex = Get<std::string>(osc_channel["LatexName"], __FILE__ , __LINE__);
    OscInfo.InitPDG           = GetFromManager(osc_channel["nutype"], 0, __FILE__,__LINE__);
    OscInfo.FinalPDG          = GetFromManager(osc_channel["oscnutype"], 0, __FILE__,__LINE__);
    OscInfo.ChannelIndex      = OscChannelCounter++;

    for (const auto& Existing : SingleSample.OscChannels) {
      if (Existing.InitPDG == OscInfo.InitPDG && Existing.FinalPDG == OscInfo.FinalPDG) {
        MACH3LOG_ERROR("Duplicate oscillation channel detected! InitPDG = {}, FinalPDG = {}"
                       "already defined in channel {} for sample {}",
                       OscInfo.InitPDG, OscInfo.FinalPDG, Existing.ChannelIndex, SingleSample.SampleTitle);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    auto MCFileNames = Get<std::vector<std::string>>(osc_channel["mtuplefile"], __FILE__ , __LINE__);
    for(size_t iFile = 0; iFile < MCFileNames.size(); iFile++){
      std::string FileName = MCFilePrefix + MCFileNames[iFile] + MCFileSuffix;
      MCFileNames[iFile] = FileName;
      FileToInitPDGMap[FileName] = NuPDG(OscInfo.InitPDG);
      FileToFinalPDGMap[FileName] = NuPDG(OscInfo.FinalPDG);
    }

    SingleSample.OscChannels.push_back(std::move(OscInfo));
    SingleSample.mc_files.push_back(MCFileNames);
    SingleSample.spline_files.push_back(SplinePrefix+osc_channel["splinefile"].as<std::string>()+SplineSuffix);
  }
  //Now grab the selection cuts from the manager
  for ( auto const &SelectionCuts : SampleSettings["SelectionCuts"]) {
    auto TempBoundsVec = GetBounds(SelectionCuts["Bounds"]);
    KinematicCut CutObj;
    CutObj.LowerBound = TempBoundsVec[0];
    CutObj.UpperBound = TempBoundsVec[1];
    CutObj.ParamToCutOnIt = ReturnKinematicParameterFromString(SelectionCuts["KinematicStr"].as<std::string>());
    MACH3LOG_INFO("Adding cut on {} with bounds {} to {}", SelectionCuts["KinematicStr"].as<std::string>(), TempBoundsVec[0], TempBoundsVec[1]);
    StoredSelection[iSample].emplace_back(CutObj);
  }
  Selection = StoredSelection;
  // The cache records which events these cuts select
  InvalidateKinematicCache();
  SelectionPtrChecked = false;
  /// Add new sample
  SampleDetails[iSample] = std::move(SingleSample);
}

// ************************************************
void SampleHandlerBase::Initialise() {
// ************************************************
  TStopwatch clock;
  clock.Start();

  //First grab all the information from your sample config via your manager
  ReadConfig();

  //Now initialise all the variables you will need
  Init();

  nEvents = SetupExperimentMC();
  MCEvents.resize(nEvents);
  SetupMC();

  MACH3LOG_INFO("=============================================");
  MACH3LOG_INFO("Total number of events is: {}", GetNEvents());
  SetupOscParameters();
  MACH3LOG_INFO("Setting up Sample Binning..");
  SetBinning();
  MACH3LOG_INFO("Setting up Splines..");
  SetupSplines();
  MACH3LOG_INFO("Setting up Normalisation Pointers..");
  SetupNormParameters();
  MACH3LOG_INFO("Setting up Functional Pointers..");
  RegisterFunctionalParameters();
  MACH3LOG_INFO("Setting up Additional Weight Pointers..");
  AdditionalWeightStart.resize(GetNEvents());
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    AdditionalWeightStart[iEvent] = static_cast<std::uint16_t>(MCEvents[iEvent].total_weight_pointers.size());
  }
  AddAdditionalWeightPointers();
  AdditionalWeightEnd.resize(GetNEvents());
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    AdditionalWeightEnd[iEvent] = static_cast<std::uint16_t>(MCEvents[iEvent].total_weight_pointers.size());
  }
  GroupingChecked = false;
  MACH3LOG_INFO("Setting up Kinematic Map..");
  SetupKinematicMap();
  clock.Stop();
  MACH3LOG_INFO("Finished loading MC for {}, it took {:.2f}s to finish", GetName(), clock.RealTime());
  MACH3LOG_INFO("Initialising Data");
  InititialiseData();
  MACH3LOG_INFO("=======================================================");
  CheckEmptyBins();
}

// ************************************************
void SampleHandlerBase::SetupKinematicMap() {
// ************************************************
  if(KinematicParameters == nullptr || ReversedKinematicParameters == nullptr) {
    MACH3LOG_ERROR("Map KinematicParameters or ReversedKinematicParameters hasn't been initialised");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  // KS: Ensure maps exist correctly
  for (const auto& pair : *KinematicParameters) {
    const auto& key = pair.first;
    const auto& value = pair.second;

    auto it = ReversedKinematicParameters->find(value);
    if (it == ReversedKinematicParameters->end() || it->second != key) {
      MACH3LOG_ERROR("Mismatch found: {} -> {} but {} -> {}",
                     key, value, value, (it != ReversedKinematicParameters->end() ? it->second : "NOT FOUND"));
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  // KS: Ensure some MaCh3 specific variables are defined
  std::vector<std::string> Vars = {"Mode", "OscillationChannel", "TargetNucleus"};
  for(size_t iVar = 0; iVar < Vars.size(); iVar++) {
    try {
      ReturnKinematicParameterFromString(Vars[iVar]);
    } catch (const MaCh3Exception&) {
      MACH3LOG_ERROR("MaCh3 expected variable: {} not found in KinematicParameters.", Vars[iVar]);
      MACH3LOG_ERROR("All keys in KinematicParameters:");
      for (const auto& pair : *KinematicParameters) {
        MACH3LOG_ERROR("Key: {}", pair.first);
      }
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// ************************************************
void SampleHandlerBase::FillHist(const int Sample, TH1* Hist, std::vector<double> &Array) {
// ************************************************
  int Dimension = GetNDim(Sample);
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  if(Dimension == 1) {
    Hist->Reset();
    for (int xBin = 0; xBin < Binning->GetNAxisBins(Sample, 0); ++xBin) {
      const int idx = Binning->GetGlobalBinSafe(Sample, {xBin});
      Hist->SetBinContent(xBin + 1, Array[idx]);
    }
  } else if (Dimension == 2) {
    Hist->Reset();
    if(Binning->IsUniform(Sample)) {
      for (int yBin = 0; yBin < Binning->GetNAxisBins(Sample, 1); ++yBin) {
        for (int xBin = 0; xBin < Binning->GetNAxisBins(Sample, 0); ++xBin) {
          const int idx = Binning->GetGlobalBinSafe(Sample, {xBin, yBin});
          Hist->SetBinContent(xBin + 1, yBin + 1, Array[idx]);
        }
      }
    } else {
      for (int iBin = 0; iBin < Binning->GetNBins(Sample); ++iBin) {
        const int idx = iBin + Binning->GetSampleStartBin(Sample);
        //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
        Hist->SetBinContent(iBin + 1, Array[idx]);
      }
    }
  } else {
    for (int iBin = 0; iBin < Binning->GetNBins(Sample); ++iBin) {
      const int idx = iBin + Binning->GetSampleStartBin(Sample);
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      Hist->SetBinContent(iBin + 1, Array[idx]);
    }
  }
}
#pragma GCC diagnostic pop


// ************************************************
bool SampleHandlerBase::IsEventSelected(const int iSample, const int iEvent) _noexcept_ {
// ************************************************
  const auto& SampleSelection = Selection[iSample];
  const int SelectionSize = static_cast<int>(SampleSelection.size());
  for (int iSelection = 0; iSelection < SelectionSize; ++iSelection) {
    const auto& Cut = SampleSelection[iSelection];
    const double Val = ReturnKinematicParameter(Cut.ParamToCutOnIt, iEvent);
    if ((Val < Cut.LowerBound) || (Val >= Cut.UpperBound)) {
      return false;
    }
  }
  //DB To avoid unnecessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

// === JM Define function to check if sub-event is selected ===
bool SampleHandlerBase::IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, const unsigned int iSubEvent, size_t nsubevents) {
  for (unsigned int iSelection=0;iSelection < SubEventCuts.size() ;iSelection++) {
    std::vector<double> Vec = ReturnKinematicVector(SubEventCuts[iSelection].ParamToCutOnIt, iEvent);
    if (nsubevents != Vec.size()) {
      MACH3LOG_ERROR("Cannot apply kinematic cut on {} as it is of different size to plotting variable");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const double Val = Vec[iSubEvent];
    if ((Val < SubEventCuts[iSelection].LowerBound) || (Val >= SubEventCuts[iSelection].UpperBound)) {
      return false;
    }
  }
  //DB To avoid unnecessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}
// ===========================================================

//************************************************
// Reweight function
void SampleHandlerBase::Reweight() {
//************************************************
  //KS: Reset the histograms before reweight
  ResetHistograms();

  //update the list of parameter values passed to shift functionals
  functional.update_vals();
  // Update the functional parameter values to the latest proposed values
  PrepFunctionalParameters();

  //HW The device reweight evaluates the splines itself, so when it is going to take this pass (the
  //shifts haven't moved) the CPU needn't
  bool OnDevice = false;
  #ifdef MaCh3_CUDA
  //HW Only while the device copy matches the grouping: SetComputeGradients and the like mark the
  //grouping stale, and the device is rebuilt from it in the pass itself
  OnDevice = DeviceUsable && !ReferencePath && GroupingChecked && (DeviceSignatureVersion == SignatureVersion) &&
             (DeviceSweepUsable || KinematicCacheMatches());
  #endif

  //HW Oscillation and spline weights depend only on non functional parameters. If none of those
  //moved since they were last evaluated (every random walk sweep), they are already right.
  //The reference loop always evaluates them and leaves the cache's bookkeeping alone
  if (!ReferencePath) UpdateWeightCache();
  const bool WeightsCurrent = !ReferencePath && WeightCacheEnabled && (WeightsEvaluatedVersion == WeightCacheVersion);

  //You only need to do these things if Oscillator has been initialised
  //if not then you're not considering oscillations
  if (Oscillator && !WeightsCurrent) Oscillator->Evaluate();

  // Calculate weight coming from all splines if we initialised handler
  if(SplineHandler && !WeightsCurrent && !OnDevice) SplineHandler->Evaluate();
  //HW After the reference loop the oscillator and splines may hold weights for parameters the cache
  //has no record of, and after a device pass the CPU spline weights weren't evaluated at all. Either
  //way the next CPU pass has to evaluate them again
  WeightsEvaluatedVersion = (WeightCacheEnabled && !ReferencePath && !OnDevice) ? WeightCacheVersion : 0;

  //KS: If using CPU this does nothing, if on GPU need to make sure we finished copying memory from
  if(SplineHandler) SplineHandler->SynchroniseMemTransfer();

  // Here we fill weights to MC predictions weights for splines and osc already have been filled
  FillArray();

  //KS: If you want to not update W2 wights then uncomment this line
  if(!UpdateW2) FirstTimeW2 = false;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"
// ************************************************
void SampleHandlerBase::FillArray() {
// ************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;
  if (ReferencePath) {
    FillArrayReference();
    return;
  }

  // NOTE comment below is left for historical reasons
  //DB - Brain dump of speedup ideas
  //
  //Those relevant to reweighting
  // 1. Don't bother storing and calculating NC signal events - Implemented and saves marginal s/step
  // 2. Loop over spline event weight calculation in the following event loop - Currently done in splineSKBase->calcWeight() where multi-threading won't be optimised - Implemented and saves 0.3s/step
  // 3. Inline getDiscVar or somehow include that calculation inside the multi-threading - Implemented and saves about 0.01s/step
  // 4. Include isCC inside SKMCStruct so don't have to have several 'if' statements determine if oscillation weight needs to be set to 1.0 for NC events - Implemented and saves marginal s/step
  // 5. Do explicit check on adjacent bins when finding event XBin instead of looping over all BinEdge indices - Implemented but doesn't significantly affect s/step
  //
  //Other aspects
  // 1. Order minituples in Y-axis variable as this will *hopefully* reduce cache misses inside SampleHandler_array_class[yBin][xBin]
  //
  // We will hit <0.1 s/step eventually! :D
  [[maybe_unused]] const auto TotalBins = Binning->GetNBins();
  const unsigned int NumberOfEvents = GetNEvents();

  double* _restrict_ MC_Array_for_reduction = SampleHandler_array.data();
  double* _restrict_ W2_array_for_reduction = SampleHandler_array_w2.data();
  //HW Null unless SetComputeGradients(true). Kept out of the reduction below, per-thread copies of
  //a NBins * nGradPars array would dwarf the MC prediction, so FillGradientArray uses atomics
  double* _restrict_ Grad_Array = (ComputeGradients && !GradientsPaused) ? SampleHandler_gradients.data() : nullptr;
  GradientsCurrent = (Grad_Array != nullptr);

  //HW The functional parameters haven't moved since the last full pass, so every event's shifted
  //kinematics, selection and bin are what that pass found. Skip all three and only redo the weights,
  //and only for the events which can contribute. Same sums as the full pass, just fewer of them
  GradFromDevice = false;
  if (KinematicCacheMatches()) {
    #ifdef MaCh3_CUDA
    if (DeviceFill(Grad_Array != nullptr)) return;
    #endif
    if (!CachedListReady) BuildCachedEventList();
    ZeroGradientsIfDeferred();
    const int NBinsUsed = static_cast<int>(CachedBinOrder.size());
    const int* _restrict_ Events = CachedEvents.data();
    const int* _restrict_ BinStart = CachedBinStart.data();
    const int* _restrict_ BinOrder = CachedBinOrder.data();
    const size_t RowLength = static_cast<size_t>(nGradPars);

    //HW One bin per iteration, so exactly one thread writes each bin of the MC and W2 arrays and
    //each row of the gradient: no reduction copies and, above all, no atomics on the gradient,
    //which is what kept the full pass from scaling with threads. Dynamic, since bins differ a lot
    //in how many events they hold, and the busiest are handed out first
    #ifdef MULTITHREAD
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (int iOrder = 0; iOrder < NBinsUsed; ++iOrder) {
      const int GlobalBin = BinOrder[iOrder];
      double* _restrict_ GradRow = (Grad_Array != nullptr) ? Grad_Array + static_cast<size_t>(GlobalBin) * RowLength : nullptr;
      double BinSum = 0.;
      double BinSumW2 = 0.;
      if (GroupingEnabled) {
        const M3::float_t* const* _restrict_ WPtr = SigWeightPtr.data();
        const SigDeriv* _restrict_ DList = SigDerivs.data();
        for (int iGroup = BinGroupStart[GlobalBin]; iGroup < BinGroupStart[GlobalBin + 1]; ++iGroup) {
          const int Sig = GroupSig[iGroup];
          double Shared = 1.;
          for (int k = SigWeightStart[Sig]; k < SigWeightStart[Sig + 1]; ++k) Shared *= *WPtr[k];
          //HW Every event of a group has constant weights of the same sign, so this skips exactly the
          //events the per event test below would
          const double GroupWeight = GroupConst[iGroup] * Shared;
          if (GroupWeight <= 0.) continue;
          BinSum += GroupWeight;
          if (FirstTimeW2) BinSumW2 += GroupConst2[iGroup] * Shared * Shared;
          if (GradRow != nullptr) {
            // Constant weights have no derivative, so only the shared ones contribute
            for (int k = SigDerivStart[Sig]; k < SigDerivStart[Sig + 1]; ++k) {
              GradRow[DList[k].Par] += GroupWeight * (*DList[k].Deriv) / *DList[k].Weight;
            }
          }
        }
        MC_Array_for_reduction[GlobalBin] += BinSum;
        if (FirstTimeW2) W2_array_for_reduction[GlobalBin] += BinSumW2;
        continue;
      }
      for (int iCached = BinStart[GlobalBin]; iCached < BinStart[GlobalBin + 1]; ++iCached) {
        const int iEvent = Events[iCached];
        const EventInfo* _restrict_ MCEvent = &MCEvents[iEvent];

        // Virtual by default does nothing, has to happen before CalcWeightTotal
        CalcWeightFunc(iEvent);

        const M3::float_t totalweight = CalcWeightTotal(MCEvent);
        if (totalweight <= 0.){
          continue;
        }
        BinSum += totalweight;
        if (FirstTimeW2) BinSumW2 += totalweight*totalweight;
        if (GradRow != nullptr) FillGradientRow(MCEvent, totalweight, GradRow);
      }
      MC_Array_for_reduction[GlobalBin] += BinSum;
      if (FirstTimeW2) W2_array_for_reduction[GlobalBin] += BinSumW2;
    }
    return;
  }

  #ifdef MaCh3_CUDA
  if (DeviceFullPass(Grad_Array != nullptr)) {
    // EventBinCache now holds every event at KinCacheKey, though the live ones are still on the device
    KinCacheValid = true;
    return;
  }
  EventBinsOnDevice = false;
  #endif
  ZeroGradientsIfDeferred();
  if (EventBinCache.size() != NumberOfEvents) EventBinCache.resize(NumberOfEvents);
  int* _restrict_ EventBin = EventBinCache.data();

  UpdateFastShifts();
  if (!SelectionPtrChecked) SetupSelectionPointers();
  UpdateWeightCache();
  const bool FastSelection = SelectionPtrReady;
  if (!FastBinChecked) SetupFastBinning();
  const bool FastBins = FastBinReady;
  if (EventSampleIdx.size() != NumberOfEvents) {
    EventSampleIdx.resize(NumberOfEvents);
    for (unsigned int i = 0; i < NumberOfEvents; ++i) EventSampleIdx[i] = MCEvents[i].NominalSample;
  }
  const int* _restrict_ EventSample = EventSampleIdx.data();
  const bool UseWeightCache = WeightCacheEnabled;
  M3::float_t* _restrict_ WeightCache = UseWeightCache ? EventWeightCache.data() : nullptr;
  unsigned int* _restrict_ WeightStamp = UseWeightCache ? EventWeightStamp.data() : nullptr;
  const unsigned int Version = WeightCacheVersion;

  //HW Events failing a cut no shift can change are never selected, so they needn't be visited at
  //all, beyond recording that once. Only with the pointer based cuts, which built the list
  const bool SkipDead = FastSelection && !LiveEvents.empty();
  if (SkipDead && DeadBinsWritten != EventBin) {
    for (unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) if (NeverSelected[iEvent]) EventBin[iEvent] = kKinCacheNotSelected;
    DeadBinsWritten = EventBin;
  }
  const unsigned int NLoop = SkipDead ? static_cast<unsigned int>(LiveEvents.size()) : NumberOfEvents;
  const int* _restrict_ Live = SkipDead ? LiveEvents.data() : nullptr;

  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:MC_Array_for_reduction[:TotalBins], W2_array_for_reduction[:TotalBins])
  #endif
  for (unsigned int iLoop = 0; iLoop < NLoop; ++iLoop) {
    const unsigned int iEvent = SkipDead ? static_cast<unsigned int>(Live[iLoop]) : iLoop;
    //ETA - generic functions to apply shifts to kinematic variables
    // Apply this before IsEventSelected is called.
    ApplyShifts(iEvent);

    const EventInfo* _restrict_ MCEvent = &MCEvents[iEvent];
    const int Sample = EventSample[iEvent];
    //ETA - generic functions to apply shifts to kinematic variable
    const bool Selected = FastSelection ? IsEventSelectedFast(Sample, static_cast<int>(iEvent))
                                        : IsEventSelected(Sample, iEvent);
    if(!Selected){
      EventBin[iEvent] = kKinCacheNotSelected;
      continue;
    }

    //DB Find the relevant bin in the PDF for each event
    //HW Ahead of the weight check, since the bin doesn't depend on the weight and the cache needs
    //it for events whose weight is zero now but needn't be on the next reweight
    const int GlobalBin = FastBins ? FindGlobalBinFast(Sample, static_cast<int>(iEvent))
                                   : Binning->FindGlobalBin(Sample, MCEvent->KinVar, MCEvent->NomBin);
    EventBin[iEvent] = GlobalBin;

    M3::float_t totalweight;
    if (UseWeightCache && WeightStamp[iEvent] == Version) {
      totalweight = WeightCache[iEvent];
    } else {
      // Virtual by default does nothing, has to happen before CalcWeightTotal
      CalcWeightFunc(iEvent);
      totalweight = CalcWeightTotal(MCEvent);
      if (UseWeightCache) {
        WeightCache[iEvent] = totalweight;
        WeightStamp[iEvent] = Version;
      }
    }
    //DB Catch negative total weights and skip any event with a negative weight. Previously we would set weight to zero and continue but that is inefficient
    if (totalweight <= 0.){
      continue;
    }

    //ETA - we can probably remove this final if check on the -1?
    //Maybe we can add an overflow bin to the array and assign any events to this bin?
    //Might save us an extra if call?
    //DB Fill relevant part of thread array
    if (GlobalBin > M3::UnderOverFlowBin) {
      MC_Array_for_reduction[GlobalBin] += totalweight;
      //HW totalweight > 0 here, so FillGradientArray can safely divide by the individual weights
      if (Grad_Array != nullptr) FillGradientArray(MCEvent, totalweight, GlobalBin, Grad_Array);
      if (FirstTimeW2) W2_array_for_reduction[GlobalBin] += totalweight*totalweight;
    }
  }
  // EventBinCache now holds every event at KinCacheKey, which KinematicCacheMatches set above
  KinCacheValid = true;
}
#pragma GCC diagnostic pop

// ************************************************
void SampleHandlerBase::FillArrayReference() {
// ************************************************
  //HW The event loop as it stood before the optimisations, kept as the reference they are checked
  //against: no caches, no grouping, no skipped events, generic shifts, cuts and bin lookup
  [[maybe_unused]] const auto TotalBins = Binning->GetNBins();
  const unsigned int NumberOfEvents = GetNEvents();
  double* _restrict_ MC_Array_for_reduction = SampleHandler_array.data();
  double* _restrict_ W2_array_for_reduction = SampleHandler_array_w2.data();
  double* _restrict_ Grad_Array = (ComputeGradients && !GradientsPaused) ? SampleHandler_gradients.data() : nullptr;
  GradientsCurrent = (Grad_Array != nullptr);
  GradFromDevice = false;
  ZeroGradientsIfDeferred();

  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:MC_Array_for_reduction[:TotalBins], W2_array_for_reduction[:TotalBins])
  #endif
  for (unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    ApplyShifts(iEvent);
    const EventInfo* _restrict_ MCEvent = &MCEvents[iEvent];
    if (!IsEventSelected(MCEvent->NominalSample, iEvent)) continue;

    CalcWeightFunc(iEvent);
    const M3::float_t totalweight = CalcWeightTotal(MCEvent);
    if (totalweight <= 0.) continue;

    const int GlobalBin = Binning->FindGlobalBin(MCEvent->NominalSample, MCEvent->KinVar, MCEvent->NomBin);
    if (GlobalBin > M3::UnderOverFlowBin) {
      MC_Array_for_reduction[GlobalBin] += totalweight;
      if (Grad_Array != nullptr) FillGradientArray(MCEvent, totalweight, GlobalBin, Grad_Array);
      if (FirstTimeW2) W2_array_for_reduction[GlobalBin] += totalweight*totalweight;
    }
  }
}

#ifdef MaCh3_CUDA
namespace {

/// Works out where in the device pool (see DeviceReweight.h) the value behind each weight or derivative
/// pointer lives: a parameter value, a spline weight or derivative, one of the constants 1 and 0, or
/// otherwise a value the host gathers each reweight (oscillation weights and the like). The constants
/// sit after the gathered values, whose number is only known once every pointer has been seen, so
/// they get placeholder indices until ResolveConstants
class PoolIndexer {
 public:
  PoolIndexer(const M3::float_t* ParBase, const int NPars, const M3::float_t* SplineWeights,
              const M3::float_t* SplineDerivs, const int NSplines)
      : ParBase(ParBase), NPars(NPars), SplineWeights(SplineWeights), SplineDerivs(SplineDerivs), NSplines(NSplines) {}

  /// The pool index of the value P points at
  int Index(const M3::float_t* P) {
    if (P >= ParBase && P < ParBase + NPars) return static_cast<int>(P - ParBase);
    if (NSplines > 0) {
      if (P >= SplineWeights && P < SplineWeights + NSplines) return NPars + static_cast<int>(P - SplineWeights);
      if (P >= SplineDerivs && P < SplineDerivs + NSplines) return NPars + NSplines + static_cast<int>(P - SplineDerivs);
    }
    if (P == &M3::Unity) return kUnityPlaceholder;
    if (P == &M3::Zero) return kZeroPlaceholder;

    // Anything else is gathered: the same pointer always gets the same slot
    auto Found = GatheredSlot.find(P);
    if (Found != GatheredSlot.end()) return Found->second;
    const int Slot = NPars + 2 * NSplines + static_cast<int>(Gathered.size());
    GatheredSlot.emplace(P, Slot);
    Gathered.push_back(P);
    return Slot;
  }

  /// Swap the placeholders for the constants' real indices, the last two in the pool
  void ResolveConstants(std::vector<int>& Indices) const {
    const int PoolSize = NPars + 2 * NSplines + static_cast<int>(Gathered.size()) + 2;
    for (int& k : Indices) {
      if (k == kUnityPlaceholder) k = PoolSize - 2;
      else if (k == kZeroPlaceholder) k = PoolSize - 1;
    }
  }

  /// The pointers whose values the host gathers each reweight, in pool order
  const std::vector<const M3::float_t*>& GatheredPointers() const { return Gathered; }

 private:
  static constexpr int kUnityPlaceholder = -1;
  static constexpr int kZeroPlaceholder = -2;
  const M3::float_t* ParBase;
  int NPars;
  const M3::float_t* SplineWeights;
  const M3::float_t* SplineDerivs;
  int NSplines;
  std::unordered_map<const M3::float_t*, int> GatheredSlot;
  std::vector<const M3::float_t*> Gathered;
};

/// Copy the spline monolith into the device inputs (the coefficients are copied, as the device keeps its own)
void FillDeviceSplineInputs(DeviceReweightStatic& Inputs, const BinnedSplineHandler::ArrayView& View) {
  Inputs.NSplines = static_cast<int>(View.NSplines);
  Inputs.NSplineParams = View.NParams;
  Inputs.MaxKnots = View.MaxKnots;
  Inputs.KnotX.assign(View.KnotX, View.KnotX + View.NKnotX);
  Inputs.SegmentCoeffs.assign(View.SegmentCoeffs, View.SegmentCoeffs + View.NSegmentCoeffs);
  Inputs.SplineCoeffStart = *View.SplineCoeffStart;
  Inputs.SplineParam = *View.SplineParam;
  Inputs.SplineParamIndex = View.SplineParamIndex;
}

}  // namespace

// ************************************************
void SampleHandlerBase::SetupDevice() {
// ************************************************
  DeviceChecked = true;
  DeviceUsable = false;
  Device.reset();

  // 1. Can this handler use the device at all?
  if (!M3::DeviceReweightRequested()) return;
  if (!GroupingEnabled) {
    MACH3LOG_WARN("{}: MACH3_DEVICE_REWEIGHT is set, but the device reweight needs event grouping "
                  "(AdditionalWeightsAreConstant), so this handler stays on the CPU", GetName());
    return;
  }
  auto* Binned = dynamic_cast<BinnedSplineHandler*>(SplineHandler.get());
  if (SplineHandler && Binned == nullptr) {
    MACH3LOG_WARN("{}: the device reweight only knows binned splines, so this handler stays on the CPU", GetName());
    return;
  }

  // 2. Sizes, data and splines
  DeviceReweightStatic Inputs;
  Inputs.NBins = Binning->GetNBins();
  Inputs.Data = SampleHandler_data;
  Inputs.NPars = ParHandler->GetNumParams();
  Inputs.LowMCBound = M3::_LOW_MC_BOUND_;
  BinnedSplineHandler::ArrayView View;
  if (Binned != nullptr) {
    View = Binned->GetArrayView();
    FillDeviceSplineInputs(Inputs, View);
  }

  // 3. Every signature's weights and derivatives as pool indices instead of pointers
  PoolIndexer Pool(ParHandler->RetPointer(0), Inputs.NPars, View.Weights, View.Grads, Inputs.NSplines);
  Inputs.SigWeightStart = SigWeightStart;
  Inputs.SigWeightIdx.reserve(SigWeightPtr.size());
  for (const M3::float_t* P : SigWeightPtr) Inputs.SigWeightIdx.push_back(Pool.Index(P));
  Inputs.SigDerivStart = SigDerivStart;
  for (const auto& D : SigDerivs) {
    Inputs.SigDerivPar.push_back(D.Par);
    Inputs.SigDerivIdx.push_back(Pool.Index(D.Deriv));
    Inputs.SigDerivWeightIdx.push_back(Pool.Index(D.Weight));
  }
  Pool.ResolveConstants(Inputs.SigWeightIdx);
  Pool.ResolveConstants(Inputs.SigDerivIdx);
  Pool.ResolveConstants(Inputs.SigDerivWeightIdx);
  DeviceGatheredPtrs = Pool.GatheredPointers();
  Inputs.NGathered = static_cast<int>(DeviceGatheredPtrs.size());

  // 4. The device copy. If the GPU can't take it, this handler simply stays on the CPU
  try {
    Device = std::make_unique<DeviceReweight>(Inputs);
  } catch (const std::exception& E) {
    MACH3LOG_ERROR("{}: couldn't set up the device reweight, so this handler stays on the CPU: {}", GetName(), E.what());
    Device.reset();
    return;
  }

  // 5. Staging vectors and bookkeeping
  DeviceGatheredVals.assign(DeviceGatheredPtrs.size(), 0.);
  DeviceParVals.assign(static_cast<size_t>(Inputs.NPars), 0.);
  DeviceGrad.assign(static_cast<size_t>(Inputs.NPars), 0.);
  DeviceSignatureVersion = SignatureVersion;
  DeviceGroupsVersion = 0;
  DeviceDataVersion = DataVersion;
  DeviceUsable = true;

  // 6. Full passes on the device too, if the shifts allow it
  SetupDeviceSweep();

  MACH3LOG_INFO("{}: device reweight on, {} parameters, {} splines, {} gathered weights, {} signatures, {:.1f} MB on the device",
                GetName(), Inputs.NPars, Inputs.NSplines, Inputs.NGathered, SigWeightStart.size() - 1,
                static_cast<double>(Device->DeviceBytes()) / 1048576.);
}

// ************************************************
bool SampleHandlerBase::DevicePrepare() {
// ************************************************
  if (ReferencePath) return false;
  //HW Rebuild the grouping first if something (SetComputeGradients) marked it stale, so the device
  //copy is rebuilt from the current signatures rather than kept from old ones
  if (!GroupingChecked) SetupEventGrouping();
  if (!GroupingEnabled) return false;
  if (!DeviceChecked || DeviceSignatureVersion != SignatureVersion) SetupDevice();
  if (!DeviceUsable) return false;
  if (DeviceDataVersion != DataVersion) { Device->SetData(SampleHandler_data); DeviceDataVersion = DataVersion; }
  return true;
}

// ************************************************
const short int* SampleHandlerBase::StageDeviceInputs() {
// ************************************************
  const auto& Prop = ParHandler->GetParPropVec();
  for (size_t i = 0; i < DeviceParVals.size(); ++i) DeviceParVals[i] = static_cast<double>(Prop[i]);
  for (size_t i = 0; i < DeviceGatheredPtrs.size(); ++i) DeviceGatheredVals[i] = static_cast<double>(*DeviceGatheredPtrs[i]);
  for (size_t k = 0; k < DeviceShiftVals.size(); ++k) DeviceShiftVals[k] = functional.shifts[k].par_vals[0];
  if (auto* Binned = dynamic_cast<BinnedSplineHandler*>(SplineHandler.get())) return Binned->FindSegments();
  return nullptr;
}

// ************************************************
bool SampleHandlerBase::DeviceFill(const bool Grads) {
// ************************************************
  if (!DevicePrepare()) return false;
  //HW If the last full pass ran on the device its event bins are still there, and the device builds the
  //groups itself. Otherwise the kinematic cache came from a CPU pass, and the host's groups go up
  if (EventBinsOnDevice && DeviceSweepUsable) {
    if (!DeviceGroupsFromFull) {
      Device->BuildGroupsFromFullPass();
      DeviceGroupsFromFull = true;
      DeviceGroupsVersion = 0;
    }
  } else {
    if (!CachedListReady) BuildCachedEventList();
    if (DeviceGroupsVersion != GroupsVersion) {
      DeviceReweightGroups G;
      G.Sig = GroupSig;
      G.Const = GroupConst;
      G.Const2 = GroupConst2;
      G.BinStart = BinGroupStart;
      Device->SetGroups(G, SampleHandler_data);
      DeviceGroupsVersion = GroupsVersion;
    }
  }
  const short int* Segments = StageDeviceInputs();
  Device->Run(DeviceParVals.data(), Segments, DeviceGatheredVals.data(), FirstTimeW2, Grads,
              SampleHandler_array.data(), SampleHandler_array_w2.data(), DeviceGrad.data());
  GradFromDevice = Grads;
  return true;
}

// ************************************************
void SampleHandlerBase::SetupDeviceSweep() {
// ************************************************
  //HW A full pass on the device needs every shifted variable to move linearly and additively with the
  //shift parameters, x = nominal + sum_k shift_k * coefficient_k. The shifts themselves are arbitrary
  //C++ functions, so the coefficients are found by probing them (ProbeShiftCoefficients) and the
  //linear form is then checked against the CPU's own ApplyShifts (CheckShiftLinearity)
  DeviceSweepUsable = false;
  EventBinsOnDevice = false;

  // 1. Requirements on the cuts, the binning and the shifts
  if (!SelectionPtrChecked) SetupSelectionPointers();
  if (!FastBinChecked) SetupFastBinning();
  if (!SelectionPtrReady || !FastBinReady) {
    MACH3LOG_INFO("{}: full passes stay on the CPU, the device needs pointer based cuts and uniform 1D binning", GetName());
    return;
  }
  for (const auto& Shift : functional.shifts) {
    if (Shift.par_vals.size() != 1) {
      MACH3LOG_INFO("{}: full passes stay on the CPU, a shift takes more than one parameter", GetName());
      return;
    }
  }

  // 2. Per sample: the cuts that can move and the binning
  DeviceSweepStatic Sweep;
  Sweep.NShifts = static_cast<int>(functional.shifts.size());
  Sweep.UnderOverFlowBin = M3::UnderOverFlowBin;
  Sweep.NotSelected = kKinCacheNotSelected;
  std::vector<std::vector<int>> ShiftedCuts;
  DescribeSweepSamples(Sweep, ShiftedCuts);
  if (Sweep.NVariables > kDeviceMaxVariables) {
    MACH3LOG_INFO("{}: full passes stay on the CPU, more than {} shifted cuts per sample", GetName(), kDeviceMaxVariables - 1);
    return;
  }

  // 3. The live events: those not already ruled out by a cut no shift can change
  DeviceLiveEvents = LiveEvents;
  if (DeviceLiveEvents.empty()) {
    DeviceLiveEvents.resize(GetNEvents());
    std::iota(DeviceLiveEvents.begin(), DeviceLiveEvents.end(), 0);
  }
  Sweep.NLive = static_cast<int>(DeviceLiveEvents.size());

  // 4. Every live event's nominal variables and shift coefficients, then the check
  ProbeShiftCoefficients(Sweep, ShiftedCuts);
  long NChecked = 0;
  const double Worst = CheckShiftLinearity(Sweep, ShiftedCuts, NChecked);
  if (!(Worst < 1e-12)) {
    MACH3LOG_WARN("{}: full passes stay on the CPU, the shifts aren't additive and linear (worst mismatch {:.3g})", GetName(), Worst);
    return;
  }

  // 5. Upload
  try {
    Device->SetSweep(Sweep);
  } catch (const std::exception& E) {
    MACH3LOG_ERROR("{}: couldn't set up device full passes, they stay on the CPU: {}", GetName(), E.what());
    return;
  }
  DeviceShiftVals.assign(static_cast<size_t>(Sweep.NShifts), 0.);
  DeviceSweepUsable = true;
  MACH3LOG_INFO("{}: full passes on the device too, {} live events, {} shifted variables per event, {} coefficients, "
                "linearity checked to {:.2g} over {} values, {:.1f} MB on the device in all", GetName(), Sweep.NLive,
                Sweep.NVariables, Sweep.Coef.size(), Worst, NChecked, static_cast<double>(Device->DeviceBytes()) / 1048576.);
}

// ************************************************
void SampleHandlerBase::DescribeSweepSamples(DeviceSweepStatic& Sweep, std::vector<std::vector<int>>& ShiftedCuts) {
// ************************************************
  //HW An event's variables, in the order the device reads them: 0 is the binning variable, 1, 2, ...
  //are the sample's cuts on shifted variables. Cuts on unshifted variables were settled for good when
  //the never selectable events were removed (LiveEvents), so the device doesn't need them
  const int NSamples = GetNSamples();
  ShiftedCuts.assign(static_cast<size_t>(NSamples), {});
  Sweep.NVariables = 1;
  Sweep.CutStart.assign(1, 0);
  Sweep.EdgeStart.assign(1, 0);
  for (int iSample = 0; iSample < NSamples; ++iSample) {
    std::vector<int>& Shifted = ShiftedCuts[static_cast<size_t>(iSample)];
    const auto& Cuts = StoredSelection[static_cast<size_t>(iSample)];
    for (size_t iCut = 0; iCut < Cuts.size(); ++iCut) {
      if (!IsKinematicParameterShifted(Cuts[iCut].ParamToCutOnIt)) continue;
      Shifted.push_back(static_cast<int>(iCut));
      Sweep.CutVariable.push_back(static_cast<int>(Shifted.size()));  // the variable slot this cut reads
      Sweep.CutLo.push_back(Cuts[iCut].LowerBound);
      Sweep.CutHi.push_back(Cuts[iCut].UpperBound);
    }
    Sweep.CutStart.push_back(static_cast<int>(Sweep.CutVariable.size()));
    Sweep.NVariables = std::max(Sweep.NVariables, 1 + static_cast<int>(Shifted.size()));

    const FastAxis& Axis = FastAxes[static_cast<size_t>(iSample)];
    Sweep.Edges.insert(Sweep.Edges.end(), Axis.Edges, Axis.Edges + Axis.NBins + 1);
    Sweep.EdgeStart.push_back(static_cast<int>(Sweep.Edges.size()));
    Sweep.Stride.push_back(Axis.Stride);
    Sweep.Offset.push_back(Axis.Offset);
  }
}

// ************************************************
void SampleHandlerBase::GetSweepVariables(const int iEvent, const std::vector<std::vector<int>>& ShiftedCuts,
                                          const int NVariables, std::vector<double*>& Variables) const {
// ************************************************
  // Where each of an event's device variables lives on the host (nullptr for slots its sample doesn't use)
  Variables.assign(static_cast<size_t>(NVariables), nullptr);
  Variables[0] = const_cast<double*>(FastBinVar[static_cast<size_t>(iEvent)]);
  const std::vector<int>& Shifted = ShiftedCuts[static_cast<size_t>(MCEvents[iEvent].NominalSample)];
  const size_t First = static_cast<size_t>(SelectionPtrStart[iEvent]);
  for (size_t j = 0; j < Shifted.size(); ++j) {
    Variables[1 + j] = const_cast<double*>(SelectionPtr[First + static_cast<size_t>(Shifted[j])]);
  }
}

// ************************************************
void SampleHandlerBase::ProbeShiftCoefficients(DeviceSweepStatic& Sweep, const std::vector<std::vector<int>>& ShiftedCuts) {
// ************************************************
  //HW With an event's variables set to zero, applying one shift with its parameter at 1 leaves exactly
  //that shift's coefficient in each variable it moves. Repeat for each shift the event has
  const size_t NVar = static_cast<size_t>(Sweep.NVariables);
  const size_t NLive = static_cast<size_t>(Sweep.NLive);
  const std::vector<double> AtOne{1.};
  Sweep.Sample.resize(NLive);
  Sweep.Sig.resize(NLive);
  Sweep.ConstWeight.resize(NLive);
  Sweep.Nominal.assign(NLive * NVar, 0.);
  Sweep.CoefStart.assign(NLive * NVar + 1, 0);
  Sweep.CoefShift.clear();
  Sweep.Coef.clear();

  std::vector<double*> Variables;
  // Found coefficients, (shift, coefficient) per variable of the current event
  std::vector<std::vector<std::pair<int, double>>> Found(NVar);
  for (size_t i = 0; i < NLive; ++i) {
    const int iEvent = DeviceLiveEvents[i];
    Sweep.Sample[i] = MCEvents[iEvent].NominalSample;
    Sweep.Sig[i] = EventSignature[static_cast<size_t>(iEvent)];
    Sweep.ConstWeight[i] = EventConstWeight[static_cast<size_t>(iEvent)];
    GetSweepVariables(iEvent, ShiftedCuts, Sweep.NVariables, Variables);

    static const std::vector<int> NoShifts;
    const std::vector<int>& Shifts = functional.event_shifts.empty() ? NoShifts : functional.event_shifts[static_cast<size_t>(iEvent)];

    // Nominal values, with the shifts undone
    if (!Shifts.empty()) ResetShifts(iEvent);
    for (size_t j = 0; j < NVar; ++j) {
      Sweep.Nominal[i * NVar + j] = (Variables[j] != nullptr) ? *Variables[j] : 0.;
      Found[j].clear();
    }

    // One shift at a time on zeroed variables
    for (const int k : Shifts) {
      ResetShifts(iEvent);
      for (double* V : Variables) if (V != nullptr) *V = 0.;
      functional.shifts[static_cast<size_t>(k)].apply(AtOne, iEvent);
      for (size_t j = 0; j < NVar; ++j) {
        if (Variables[j] == nullptr) continue;
        const double C = *Variables[j];
        if (C != 0.) Found[j].push_back({k, C});
      }
    }

    // Put back the nominal values, whether or not ResetShifts covers every variable that was zeroed
    if (!Shifts.empty()) ResetShifts(iEvent);
    for (size_t j = 0; j < NVar; ++j) {
      if (Variables[j] != nullptr) *Variables[j] = Sweep.Nominal[i * NVar + j];
    }

    for (size_t j = 0; j < NVar; ++j) {
      for (const auto& [k, C] : Found[j]) {
        Sweep.CoefShift.push_back(k);
        Sweep.Coef.push_back(C);
      }
      Sweep.CoefStart[i * NVar + j + 1] = static_cast<int>(Sweep.Coef.size());
    }
  }
}

// ************************************************
double SampleHandlerBase::CheckShiftLinearity(const DeviceSweepStatic& Sweep, const std::vector<std::vector<int>>& ShiftedCuts,
                                              long& NChecked) {
// ************************************************
  //HW At random shift values the CPU's own ApplyShifts has to give what the coefficients say, for every
  //live event and variable. A shift that is nonlinear, doesn't add, or is touched by FinaliseShifts fails.
  //Returns the worst relative mismatch
  const int NShifts = Sweep.NShifts;
  const size_t NVar = static_cast<size_t>(Sweep.NVariables);
  std::vector<double> Saved(static_cast<size_t>(NShifts));
  for (int k = 0; k < NShifts; ++k) Saved[static_cast<size_t>(k)] = functional.shifts[static_cast<size_t>(k)].par_vals[0];
  const bool WasReference = ReferencePath;
  ReferencePath = true;  // ApplyShifts then reads par_vals rather than the fast table

  TRandom3 Rng(77);
  double Worst = 0.;
  NChecked = 0;
  std::vector<double*> Variables;
  std::vector<double> Vals(static_cast<size_t>(NShifts));
  for (int Trial = 0; Trial < 3; ++Trial) {
    for (int k = 0; k < NShifts; ++k) {
      Vals[static_cast<size_t>(k)] = 0.1 * Rng.Gaus();
      functional.shifts[static_cast<size_t>(k)].par_vals[0] = Vals[static_cast<size_t>(k)];
    }
    for (size_t i = 0; i < static_cast<size_t>(Sweep.NLive); ++i) {
      const int iEvent = DeviceLiveEvents[i];
      ApplyShifts(iEvent);
      GetSweepVariables(iEvent, ShiftedCuts, Sweep.NVariables, Variables);
      for (size_t j = 0; j < NVar; ++j) {
        if (Variables[j] == nullptr) continue;
        const size_t Cell = i * NVar + j;
        double Predicted = Sweep.Nominal[Cell];
        for (int c = Sweep.CoefStart[Cell]; c < Sweep.CoefStart[Cell + 1]; ++c) {
          Predicted = Predicted + Vals[static_cast<size_t>(Sweep.CoefShift[static_cast<size_t>(c)])] * Sweep.Coef[static_cast<size_t>(c)];
        }
        const double Actual = *Variables[j];
        Worst = std::max(Worst, std::fabs(Predicted - Actual) / std::max(1., std::fabs(Actual)));
        ++NChecked;
      }
    }
  }

  // Leave everything as it was
  for (int k = 0; k < NShifts; ++k) functional.shifts[static_cast<size_t>(k)].par_vals[0] = Saved[static_cast<size_t>(k)];
  if (!functional.event_shifts.empty()) {
    for (const int iEvent : DeviceLiveEvents) {
      if (!functional.event_shifts[static_cast<size_t>(iEvent)].empty()) ResetShifts(iEvent);
    }
  }
  ReferencePath = WasReference;
  InvalidateKinematicCache();
  return Worst;
}

// ************************************************
bool SampleHandlerBase::DeviceFullPass(const bool Grads) {
// ************************************************
  if (!DevicePrepare() || !DeviceSweepUsable) return false;
  const short int* Segments = StageDeviceInputs();
  Device->RunFull(DeviceParVals.data(), Segments, DeviceGatheredVals.data(), DeviceShiftVals.data(), FirstTimeW2, Grads,
                  SampleHandler_array.data(), SampleHandler_array_w2.data(), DeviceGrad.data());
  GradFromDevice = Grads;
  EventBinsOnDevice = true;
  DeviceGroupsFromFull = false;
  return true;
}

// ************************************************
void SampleHandlerBase::PullEventBins() {
// ************************************************
  const unsigned int NumberOfEvents = GetNEvents();
  if (EventBinCache.size() != NumberOfEvents) EventBinCache.resize(NumberOfEvents);
  // Everything not live failed a cut no shift can change
  if (DeviceLiveEvents.size() != NumberOfEvents) std::fill(EventBinCache.begin(), EventBinCache.end(), kKinCacheNotSelected);
  std::vector<int> Bins(DeviceLiveEvents.size());
  Device->GetEventBins(Bins.data());
  for (size_t i = 0; i < Bins.size(); ++i) EventBinCache[static_cast<size_t>(DeviceLiveEvents[i])] = Bins[i];
  // The CPU full pass skips marking the dead events again while this pointer matches, so make it look again
  DeadBinsWritten = nullptr;
  EventBinsOnDevice = false;
}
#endif

// ************************************************
double SampleHandlerBase::CheckAgainstReference(const int NPoints, const double Tolerance,
                                                const double Width, const unsigned int Seed) {
// ************************************************
  const int NPar = ParHandler->GetNumParams();
  std::vector<double> Saved(static_cast<size_t>(NPar));
  for (int i = 0; i < NPar; ++i) Saved[static_cast<size_t>(i)] = ParHandler->GetParProp(i);
  const bool WasReference = ReferencePath;
  ReferencePath = false;

  TRandom3 Rng(Seed);
  // 0: only the shifts move, 1: only the rest, 2: both
  auto Throw = [&](const int What) {
    for (int i = 0; i < NPar; ++i) {
      if (ParHandler->IsParameterFixed(i)) continue;
      const bool Shift = (ParHandler->GetParamType(i) == SystType::kFunc);
      if ((What == 0 && !Shift) || (What == 1 && Shift)) continue;
      const double Lo = ParHandler->GetLowerBound(i), Hi = ParHandler->GetUpperBound(i);
      double Val;
      if (ParHandler->GetFlatPrior(i)) Val = Rng.Uniform(std::max(Lo, -10.), std::min(Hi, 10.));
      else Val = ParHandler->GetParPreFit(i) + Width * ParHandler->GetDiagonalError(i) * Rng.Gaus();
      ParHandler->SetParProp(i, std::min(Hi, std::max(Lo, Val)));
    }
  };
  const char* Moves[] = {"shifts only", "non-shifts only", "everything"};

  const bool Grads = ComputeGradients && !GradientsPaused;
  const bool W2 = FirstTimeW2;
  double Worst = 0.;
  for (int iPoint = 0; iPoint < NPoints; ++iPoint) {
    // Start with a move of everything, so the first optimised reweight is a full pass from scratch
    const int What = (iPoint + 2) % 3;
    Throw(What);

    Reweight();
    const std::vector<double> MC = SampleHandler_array, MCW2 = SampleHandler_array_w2, Grad = SampleHandler_gradients;
    const double LLH = GetLikelihood();
    //HW A device reweight sums d(-lnL)/d(parameter) directly and never fills the per bin array, so
    //compare what the sampler actually uses
    const bool OptFromDevice = GradFromDevice;
    const bool CompareLLHGrad = Grads && fTestStatistic == kPoisson;
    const std::vector<double> LLHGrad = CompareLLHGrad ? GetGradLikelihood() : std::vector<double>();
    ReferencePath = true;
    Reweight();
    ReferencePath = false;

    double PointWorst = 0.;
    auto Note = [&](const double Diff, const char* What2, const int Index, const double Opt, const double Ref) {
      PointWorst = std::max(PointWorst, Diff);
      if (Diff > Tolerance) {
        MACH3LOG_ERROR("{}: point {} ({}), {} {}: optimised {:.17g}, reference {:.17g}, relative difference {:.3g}",
                       GetName(), iPoint, Moves[What], What2, Index, Opt, Ref, Diff);
      }
    };
    for (size_t b = 0; b < MC.size(); ++b) {
      const double Ref = SampleHandler_array[b];
      Note(std::fabs(MC[b] - Ref) / std::max(std::fabs(Ref), 1e-300), "bin", static_cast<int>(b), MC[b], Ref);
      if (W2) {
        const double Ref2 = SampleHandler_array_w2[b];
        Note(std::fabs(MCW2[b] - Ref2) / std::max(std::fabs(Ref2), 1e-300), "W2 bin", static_cast<int>(b), MCW2[b], Ref2);
      }
    }
    if (CompareLLHGrad) {
      const std::vector<double> RefGrad = GetGradLikelihood();
      double Scale = 0.;
      for (const double G : RefGrad) Scale = std::max(Scale, std::fabs(G));
      for (size_t k = 0; k < RefGrad.size(); ++k) {
        Note(std::fabs(LLHGrad[k] - RefGrad[k]) / std::max(Scale, 1e-300), "-lnL gradient entry",
             static_cast<int>(k), LLHGrad[k], RefGrad[k]);
      }
    }
    if (Grads && !OptFromDevice) {
      double Scale = 0.;
      for (const double G : SampleHandler_gradients) Scale = std::max(Scale, std::fabs(G));
      for (size_t k = 0; k < Grad.size(); ++k) {
        Note(std::fabs(Grad[k] - SampleHandler_gradients[k]) / std::max(Scale, 1e-300), "gradient entry",
             static_cast<int>(k), Grad[k], SampleHandler_gradients[k]);
      }
    }
    const double RefLLH = GetLikelihood();
    Note(std::fabs(LLH - RefLLH) / std::max(std::fabs(RefLLH), 1.), "-lnL", 0, LLH, RefLLH);
    MACH3LOG_INFO("{}: point {} moving {}, -lnL {:.10f}, largest relative difference from the reference loop {:.3g}",
                  GetName(), iPoint, Moves[What], RefLLH, PointWorst);
    Worst = std::max(Worst, PointWorst);
  }

  for (int i = 0; i < NPar; ++i) ParHandler->SetParProp(i, Saved[static_cast<size_t>(i)]);
  ReferencePath = WasReference;
  Reweight();
  if (Worst > Tolerance) {
    MACH3LOG_ERROR("{}: the optimised reweight differs from the reference loop by up to {:.3g}, tolerance {:.3g}",
                   GetName(), Worst, Tolerance);
  } else {
    MACH3LOG_INFO("{}: the optimised reweight agrees with the reference loop to {:.3g} over {} points",
                  GetName(), Worst, NPoints);
  }
  return Worst;
}

// ************************************************
bool SampleHandlerBase::KinematicCacheMatches() {
// ************************************************
  size_t NVals = 0;
  for (const auto& Shift : functional.shifts) NVals += Shift.par_vals.size();

  bool Same = KinCacheValid && (KinCacheKey.size() == NVals);
  size_t iVal = 0;
  for (const auto& Shift : functional.shifts) {
    if (!Same) break;
    for (const double Val : Shift.par_vals) {
      //HW Exact on purpose: a parameter a sampler holds fixed comes back bit for bit identical, and
      //anything else has moved. A NaN never matches, so it always gets a full pass
      if (KinCacheKey[iVal++] != Val) { Same = false; break; }
    }
  }
  if (Same) return true;

  KinCacheKey.clear();
  KinCacheKey.reserve(NVals);
  for (const auto& Shift : functional.shifts) {
    KinCacheKey.insert(KinCacheKey.end(), Shift.par_vals.begin(), Shift.par_vals.end());
  }
  InvalidateKinematicCache();
  return false;
}

// ************************************************
void SampleHandlerBase::UpdateFastShifts() {
// ************************************************
  FastShifts = (functional.apply_all != nullptr) && functional.apply_all_possible && !functional.shifts.empty();
  for (const auto& Shift : functional.shifts) {
    if (Shift.raw == nullptr) { FastShifts = false; break; }
  }
  if (!FastShifts) return;

  // The event lists only change at registration, so they are flattened once
  const size_t NEvents = functional.event_shifts.size();
  if (functional.flat_start.size() != NEvents + 1) {
    functional.flat_start.assign(NEvents + 1, 0);
    functional.flat_index.clear();
    for (size_t iEvent = 0; iEvent < NEvents; ++iEvent) {
      functional.flat_index.insert(functional.flat_index.end(), functional.event_shifts[iEvent].begin(), functional.event_shifts[iEvent].end());
      functional.flat_start[iEvent + 1] = static_cast<int>(functional.flat_index.size());
    }
  }
  // The values change every reweight; update_vals has already run
  functional.raw_table.resize(functional.shifts.size());
  for (size_t iShift = 0; iShift < functional.shifts.size(); ++iShift) {
    functional.raw_table[iShift].raw = functional.shifts[iShift].raw;
    functional.raw_table[iShift].value = functional.shifts[iShift].par_vals[0];
  }
}

// ************************************************
void SampleHandlerBase::SetupSelectionPointers() {
// ************************************************
  SelectionPtrChecked = true;
  SelectionPtrReady = false;
  const int NumberOfEvents = static_cast<int>(GetNEvents());
  SelectionPtrStart.assign(static_cast<size_t>(NumberOfEvents) + 1, 0);
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    SelectionPtrStart[iEvent + 1] = SelectionPtrStart[iEvent] + static_cast<int>(StoredSelection[MCEvents[iEvent].NominalSample].size());
  }
  SelectionPtr.assign(static_cast<size_t>(SelectionPtrStart[NumberOfEvents]), nullptr);
  try {
    for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
      const auto& Cuts = StoredSelection[MCEvents[iEvent].NominalSample];
      for (size_t iCut = 0; iCut < Cuts.size(); ++iCut) {
        const double* Ptr = GetPointerToKinematicParameter(Cuts[iCut].ParamToCutOnIt, iEvent);
        const double Val = ReturnKinematicParameter(Cuts[iCut].ParamToCutOnIt, iEvent);
        //HW A variable which only exists as a computed value can't be read through a pointer
        if (Ptr == nullptr || !(*Ptr == Val || (std::isnan(*Ptr) && std::isnan(Val)))) {
          MACH3LOG_INFO("{}: cut variable {} of event {} can't be read through a pointer, selection cuts "
                        "go through ReturnKinematicParameter", GetName(), Cuts[iCut].ParamToCutOnIt, iEvent);
          return;
        }
        SelectionPtr[static_cast<size_t>(SelectionPtrStart[iEvent]) + iCut] = Ptr;
      }
    }
  } catch (const std::exception& Err) {
    MACH3LOG_INFO("{}: no pointer to a selection cut variable ({}), selection cuts go through "
                  "ReturnKinematicParameter", GetName(), Err.what());
    return;
  }
  SelectionPtrReady = true;
  MACH3LOG_INFO("{}: full event passes read selection cut variables through pointers", GetName());

  // Cuts on variables no shift touches decide once and for all
  NeverSelected.assign(static_cast<size_t>(NumberOfEvents), 0);
  LiveEvents.clear();
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    const auto& Cuts = StoredSelection[MCEvents[iEvent].NominalSample];
    for (size_t iCut = 0; iCut < Cuts.size(); ++iCut) {
      if (IsKinematicParameterShifted(Cuts[iCut].ParamToCutOnIt)) continue;
      const double Val = *SelectionPtr[static_cast<size_t>(SelectionPtrStart[iEvent]) + iCut];
      if (Val < Cuts[iCut].LowerBound || Val >= Cuts[iCut].UpperBound) { NeverSelected[iEvent] = 1; break; }
    }
    if (!NeverSelected[iEvent]) LiveEvents.push_back(iEvent);
  }
  if (static_cast<int>(LiveEvents.size()) == NumberOfEvents) {
    NeverSelected.clear();
    LiveEvents.clear();
  } else {
    MACH3LOG_INFO("{}: {} of {} events fail a cut no shift can change, full event passes skip them",
                  GetName(), NumberOfEvents - static_cast<int>(LiveEvents.size()), NumberOfEvents);
  }
}

// ************************************************
void SampleHandlerBase::SetupFastBinning() {
// ************************************************
  FastBinChecked = true;
  FastBinReady = false;
  const int NSamples = GetNSamples();
  FastAxes.assign(static_cast<size_t>(NSamples), FastAxis{});
  for (int iSample = 0; iSample < NSamples; ++iSample) {
    const auto& Info = Binning->GetSampleBinningInfo(iSample);
    if (Info.BinEdges.size() != 1 || !Info.Uniform || Info.BinLookup.size() != 1 || Info.Strides.size() != 1) return;
    FastAxis& Axis = FastAxes[static_cast<size_t>(iSample)];
    Axis.Edges = Info.BinEdges[0].data();
    Axis.Lookup = Info.BinLookup[0].data();
    Axis.NBins = Info.AxisNBins[0];
    Axis.Stride = Info.Strides[0];
    Axis.Offset = Info.GlobalOffset;
    if (static_cast<int>(Info.BinEdges[0].size()) != Axis.NBins + 1 || static_cast<int>(Info.BinLookup[0].size()) != Axis.NBins) return;
  }
  const int NumberOfEvents = static_cast<int>(GetNEvents());
  FastBinVar.resize(static_cast<size_t>(NumberOfEvents));
  FastNomBin.resize(static_cast<size_t>(NumberOfEvents));
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    const EventInfo& Event = MCEvents[iEvent];
    if (Event.KinVar.size() != 1 || Event.NomBin.size() != 1) return;
    FastBinVar[iEvent] = Event.KinVar[0];
    FastNomBin[iEvent] = Event.NomBin[0];
    //HW Same answer as FindGlobalBin for every event at the current values, or the fast lookup isn't used
    if (FindGlobalBinFast(Event.NominalSample, iEvent) != Binning->FindGlobalBin(Event.NominalSample, Event.KinVar, Event.NomBin)) {
      MACH3LOG_WARN("{}: fast bin lookup disagrees with FindGlobalBin for event {}, not using it", GetName(), iEvent);
      return;
    }
  }
  FastBinReady = true;
  MACH3LOG_INFO("{}: full event passes use the flat 1D bin lookup", GetName());
}

// ************************************************
void SampleHandlerBase::UpdateWeightCache() {
// ************************************************
  WeightCacheEnabled = false;
  if (!AdditionalWeightsAreConstant() || ParHandler == nullptr) return;
  if (WeightCachePars.empty()) {
    for (int i = 0; i < ParHandler->GetNumParams(); ++i) {
      if (ParHandler->GetParamType(i) != SystType::kFunc) WeightCachePars.push_back(i);
    }
  }
  const auto& Values = ParHandler->GetParPropVec();
  bool Same = (WeightCacheKey.size() == WeightCachePars.size());
  for (size_t k = 0; Same && k < WeightCachePars.size(); ++k) {
    if (WeightCacheKey[k] != Values[static_cast<size_t>(WeightCachePars[k])]) Same = false;
  }
  if (!Same) {
    WeightCacheKey.resize(WeightCachePars.size());
    for (size_t k = 0; k < WeightCachePars.size(); ++k) WeightCacheKey[k] = Values[static_cast<size_t>(WeightCachePars[k])];
    // Every stamp from before is now stale. Wrapping round would make the oldest ones look fresh
    if (++WeightCacheVersion == 0) {
      std::fill(EventWeightStamp.begin(), EventWeightStamp.end(), 0u);
      WeightCacheVersion = 1;
    }
  }
  if (EventWeightCache.size() != GetNEvents()) {
    EventWeightCache.assign(GetNEvents(), 0.);
    EventWeightStamp.assign(GetNEvents(), 0u);
  }
  WeightCacheEnabled = true;
}

// ************************************************
void SampleHandlerBase::SetupEventGrouping() {
// ************************************************
  GroupingChecked = true;
  GroupingEnabled = false;
  const int NumberOfEvents = static_cast<int>(GetNEvents());
  if (!AdditionalWeightsAreConstant()) return;
  if (static_cast<int>(AdditionalWeightStart.size()) != NumberOfEvents || static_cast<int>(AdditionalWeightEnd.size()) != NumberOfEvents) return;

  struct PtrListHash {
    size_t operator()(const std::vector<const M3::float_t*>& V) const noexcept {
      size_t h = 1469598103934665603ULL;
      for (const auto* P : V) { h ^= reinterpret_cast<std::uintptr_t>(P); h *= 1099511628211ULL; }
      return h;
    }
  };
  std::unordered_map<std::vector<const M3::float_t*>, int, PtrListHash> Signatures;
  std::vector<int> SignatureFirstEvent;
  EventSignature.resize(static_cast<size_t>(NumberOfEvents));
  EventConstWeight.resize(static_cast<size_t>(NumberOfEvents));
  std::vector<const M3::float_t*> Shared;
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    const auto& Pointers = MCEvents[iEvent].total_weight_pointers;
    const int ConstStart = AdditionalWeightStart[iEvent];
    const int ConstEnd = AdditionalWeightEnd[iEvent];
    Shared.assign(Pointers.begin(), Pointers.begin() + ConstStart);
    Shared.insert(Shared.end(), Pointers.begin() + ConstEnd, Pointers.end());
    double Const = 1.;
    for (int iWeight = ConstStart; iWeight < ConstEnd; ++iWeight) Const *= *(Pointers[iWeight]);
    EventConstWeight[iEvent] = Const;

    const auto Found = Signatures.emplace(Shared, static_cast<int>(Signatures.size()));
    EventSignature[iEvent] = Found.first->second;
    if (Found.second) { SignatureFirstEvent.push_back(iEvent); continue; }

    //HW Same weights should mean the same derivatives, but the group sum leans on it, so check
    const int iFirst = SignatureFirstEvent[static_cast<size_t>(Found.first->second)];
    const auto& Derivs = MCEvents[iEvent].total_weight_gradients;
    const auto& FirstDerivs = MCEvents[iFirst].total_weight_gradients;
    const bool Aligned = (Derivs.size() == Pointers.size()) && (FirstDerivs.size() == MCEvents[iFirst].total_weight_pointers.size());
    if (!Aligned) continue;
    const int FirstStart = AdditionalWeightStart[iFirst];
    const int FirstEnd = AdditionalWeightEnd[iFirst];
    std::vector<std::pair<int, const M3::float_t*>> A(Derivs.begin(), Derivs.begin() + ConstStart), B(FirstDerivs.begin(), FirstDerivs.begin() + FirstStart);
    A.insert(A.end(), Derivs.begin() + ConstEnd, Derivs.end());
    B.insert(B.end(), FirstDerivs.begin() + FirstEnd, FirstDerivs.end());
    if (A != B) {
      MACH3LOG_WARN("{}: events {} and {} share every weight but not every derivative, so they can't "
                    "be summed together. Event grouping is off", GetName(), iFirst, iEvent);
      return;
    }
  }

  // Counting sort on the signature
  const int NSignatures = static_cast<int>(Signatures.size());
  std::vector<int> Start(static_cast<size_t>(NSignatures) + 1, 0);
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) ++Start[static_cast<size_t>(EventSignature[iEvent]) + 1];
  for (int iSig = 0; iSig < NSignatures; ++iSig) Start[iSig + 1] += Start[iSig];
  EventsBySignature.resize(static_cast<size_t>(NumberOfEvents));
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) EventsBySignature[static_cast<size_t>(Start[EventSignature[iEvent]]++)] = iEvent;

  // Flat copies of every signature's shared weights and derivatives, from its first event
  SigWeightPtr.clear(); SigDerivs.clear();
  SigWeightStart.assign(static_cast<size_t>(NSignatures) + 1, 0);
  SigDerivStart.assign(static_cast<size_t>(NSignatures) + 1, 0);
  for (int iSig = 0; iSig < NSignatures; ++iSig) {
    const int iEvent = SignatureFirstEvent[static_cast<size_t>(iSig)];
    const auto& Pointers = MCEvents[iEvent].total_weight_pointers;
    const auto& Derivs = MCEvents[iEvent].total_weight_gradients;
    const bool HaveDerivs = (Derivs.size() == Pointers.size());
    const int ConstStart = AdditionalWeightStart[iEvent];
    const int ConstEnd = AdditionalWeightEnd[iEvent];
    for (int iWeight = 0; iWeight < static_cast<int>(Pointers.size()); ++iWeight) {
      if (iWeight >= ConstStart && iWeight < ConstEnd) continue;
      SigWeightPtr.push_back(Pointers[iWeight]);
      if (HaveDerivs && Derivs[iWeight].second != &M3::Zero) {
        SigDerivs.push_back({Derivs[iWeight].first, Derivs[iWeight].second, Pointers[iWeight]});
      }
    }
    SigWeightStart[iSig + 1] = static_cast<int>(SigWeightPtr.size());
    SigDerivStart[iSig + 1] = static_cast<int>(SigDerivs.size());
  }

  GroupingEnabled = true;
  ++SignatureVersion;
  MACH3LOG_INFO("{}: {} events share {} distinct sets of varying weights, the kinematic cache will sum "
                "events with the same set and bin before the event loop", GetName(), NumberOfEvents, NSignatures);
}

// ************************************************
void SampleHandlerBase::BuildEventGroups() {
// ************************************************
  const int NBins = Binning->GetNBins();
  GroupRep.clear();
  GroupSig.clear();
  GroupConst.clear();
  GroupConst2.clear();
  BinGroupStart.assign(static_cast<size_t>(NBins) + 1, 0);
  for (int iBin = 0; iBin < NBins; ++iBin) {
    BinGroupStart[iBin] = static_cast<int>(GroupRep.size());
    int iCached = CachedBinStart[iBin];
    const int End = CachedBinStart[iBin + 1];
    while (iCached < End) {
      // One run of a signature, split by the sign of the constant weights
      const int Sig = EventSignature[CachedEvents[iCached]];
      int PosRep = -1, NegRep = -1;
      double PosC = 0., PosC2 = 0., NegC = 0., NegC2 = 0.;
      for (; iCached < End && EventSignature[CachedEvents[iCached]] == Sig; ++iCached) {
        const int iEvent = CachedEvents[iCached];
        const double C = EventConstWeight[iEvent];
        // A zero constant weight makes every total weight zero, which the event loop skips anyway
        if (C > 0.) { if (PosRep < 0) PosRep = iEvent; PosC += C; PosC2 += C * C; }
        else if (C < 0.) { if (NegRep < 0) NegRep = iEvent; NegC += C; NegC2 += C * C; }
      }
      if (PosRep >= 0) { GroupRep.push_back(PosRep); GroupSig.push_back(Sig); GroupConst.push_back(PosC); GroupConst2.push_back(PosC2); }
      if (NegRep >= 0) { GroupRep.push_back(NegRep); GroupSig.push_back(Sig); GroupConst.push_back(NegC); GroupConst2.push_back(NegC2); }
    }
  }
  BinGroupStart[NBins] = static_cast<int>(GroupRep.size());
  ++GroupsVersion;
}

// ************************************************
void SampleHandlerBase::BuildCachedEventList() {
// ************************************************
  #ifdef MaCh3_CUDA
  if (EventBinsOnDevice) PullEventBins();
  #endif
  //HW Counting sort on the bin: one pass to count, one to place, events keep their order in a bin
  const int NBins = Binning->GetNBins();
  const int NumberOfEvents = static_cast<int>(EventBinCache.size());
  CachedBinStart.assign(static_cast<size_t>(NBins) + 1, 0);
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    // Events outside the binning (UnderOverFlowBin) never fill anything, so they are left out too
    const int Bin = EventBinCache[iEvent];
    if (Bin > M3::UnderOverFlowBin) ++CachedBinStart[static_cast<size_t>(Bin) + 1];
  }
  for (int iBin = 0; iBin < NBins; ++iBin) CachedBinStart[iBin + 1] += CachedBinStart[iBin];

  CachedEvents.resize(static_cast<size_t>(CachedBinStart[NBins]));
  std::vector<int> Next(CachedBinStart.begin(), CachedBinStart.end() - 1);
  if (!GroupingChecked) SetupEventGrouping();
  //HW Placing events in signature order leaves each bin sorted by signature, so a group is a run
  for (int iOrder = 0; iOrder < NumberOfEvents; ++iOrder) {
    const int iEvent = GroupingEnabled ? EventsBySignature[iOrder] : iOrder;
    const int Bin = EventBinCache[iEvent];
    if (Bin > M3::UnderOverFlowBin) CachedEvents[static_cast<size_t>(Next[Bin]++)] = iEvent;
  }
  if (GroupingEnabled) BuildEventGroups();

  CachedBinOrder.clear();
  for (int iBin = 0; iBin < NBins; ++iBin) {
    if (CachedBinStart[iBin + 1] > CachedBinStart[iBin]) CachedBinOrder.push_back(iBin);
  }
  std::stable_sort(CachedBinOrder.begin(), CachedBinOrder.end(), [this](const int a, const int b) {
    return (CachedBinStart[a + 1] - CachedBinStart[a]) > (CachedBinStart[b + 1] - CachedBinStart[b]);
  });
  CachedListReady = true;
}

// **************************************************
// Helper function to reset the data and MC histograms
void SampleHandlerBase::ResetHistograms() _noexcept_ {
// **************************************************
  // DB Reset values stored in PDF array to 0.
  // Don't openMP this; no significant gain
  const int nBins = Binning->GetNBins();
  std::fill_n(SampleHandler_array.begin(), nBins, 0.0);
  if (FirstTimeW2) {
    std::fill_n(SampleHandler_array_w2.begin(), nBins, 0.0);
  }
  //HW The gradients accumulate the same way the MC prediction does, so they need wiping too. With
  //the device reweight on most passes never touch the array, so the CPU passes wipe it themselves
  //(ZeroGradientsIfDeferred) rather than every reweight paying for it
  if (ComputeGradients && !GradientsPaused) {
    if (GradResetDeferred()) GradNeedsZero = true;
    else std::fill(SampleHandler_gradients.begin(), SampleHandler_gradients.end(), 0.0);
  }
} // end function

// ***************************************************************************
void SampleHandlerBase::ApplyShifts(const int iEvent) {
// ***************************************************************************
  // KS: If there are no shifts then there is no point in resetting which can be costly.
  if(!functional.event_shifts.size() || !functional.event_shifts[iEvent].size()) {
    return;
  }

  // Given a sample and event, apply the shifts to the event based on the vector of functional parameter enums
  // First reset shifted array back to nominal values
  ResetShifts(iEvent);

  if (FastShifts && !ReferencePath) {
    functional.apply_all(functional, functional.apply_all_events, iEvent);
  } else {
    for (auto const &iShift : functional.event_shifts[iEvent]) {
      auto & shift = functional.shifts[iShift];
      shift.apply(shift.par_vals, iEvent);
    }
  }

  FinaliseShifts(iEvent);
}

// ***************************************************************************
// Calculate the spline weight for one event
M3::float_t SampleHandlerBase::CalcWeightTotal(const EventInfo* _restrict_ MCEvent) const _noexcept_ {
// ***************************************************************************
  M3::float_t TotalWeight = 1.0;
  const int TotalWeights = static_cast<int>(MCEvent->total_weight_pointers.size());
  //DB Loop over stored pointers
  #ifdef MULTITHREAD
  #pragma omp simd reduction(*:TotalWeight)
  #endif
  for (int iWeight = 0; iWeight < TotalWeights; ++iWeight) {
    TotalWeight *= *(MCEvent->total_weight_pointers[iWeight]);
  }

  return TotalWeight;
}

// ***************************************************************************
// Accumulate the derivative of one event's total weight into the flat gradient array
void SampleHandlerBase::FillGradientArray(const EventInfo* _restrict_ MCEvent, const double TotalWeight,
                                          const int GlobalBin, double* _restrict_ GradArray) const _noexcept_ {
// ***************************************************************************
  //HW Flat array, so this bin's nGradPars entries start here and sit next to each other
  const size_t BinOffset = static_cast<size_t>(GlobalBin) * static_cast<size_t>(nGradPars);
  const size_t TotalWeights = MCEvent->total_weight_pointers.size();

  //DB Loop over stored pointers
  //HW No omp simd, several weights can land on the same parameter so the scatter isn't vectorisable
  for (size_t iWeight = 0; iWeight < TotalWeights; ++iWeight) {
    const auto& [iPar, DerivPointer] = MCEvent->total_weight_gradients[iWeight];
    //HW Weights we don't track a derivative for (oscillation, unbinned splines) point at M3::Zero
    if (DerivPointer == &M3::Zero) continue;

    const double Weight = *(MCEvent->total_weight_pointers[iWeight]);
    const double Deriv = *DerivPointer;
    //HW dW/dx_p = W * sum_k (1/w_k) dw_k/dx_p
    const double Contribution = TotalWeight * Deriv / Weight;

    #ifdef MULTITHREAD
    #pragma omp atomic
    #endif
    GradArray[BinOffset + static_cast<size_t>(iPar)] += Contribution;
  }
}

// ***************************************************************************
// Same sum as FillGradientArray, for a caller that owns the row so needs no atomics
void SampleHandlerBase::FillGradientRow(const EventInfo* _restrict_ MCEvent, const double TotalWeight,
                                        double* _restrict_ GradRow) const _noexcept_ {
// ***************************************************************************
  const size_t TotalWeights = MCEvent->total_weight_pointers.size();
  for (size_t iWeight = 0; iWeight < TotalWeights; ++iWeight) {
    const auto& [iPar, DerivPointer] = MCEvent->total_weight_gradients[iWeight];
    if (DerivPointer == &M3::Zero) continue;
    GradRow[iPar] += TotalWeight * (*DerivPointer) / *(MCEvent->total_weight_pointers[iWeight]);
  }
}

// ***************************************************************************
// Setup the osc parameters
void SampleHandlerBase::SetupOscParameters() {
// ***************************************************************************
  // KS: Only make sense to setup osc if you have ParHandler
  if(ParHandler == nullptr ) return;

  auto OscParams = ParHandler->GetOscParsFromSampleName(GetSampleName(0));
  if (OscParams.size() > 0) {
    MACH3LOG_INFO("Setting up NuOscillator..");
    if (Oscillator != nullptr) {
      MACH3LOG_INFO("You have passed an OscillatorBase object through the constructor of a SampleHandlerFD object - this will be used for all oscillation channels");
      if(Oscillator->isEqualBinningPerOscChannel() != true) {
        MACH3LOG_ERROR("Trying to run shared NuOscillator without EqualBinningPerOscChannel, this will not work");
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      if(OscParams.size() != Oscillator->GetOscParamsSize()){
        MACH3LOG_ERROR("SampleHandler {} has {} osc params, while shared NuOsc has {} osc params", GetName(),
                        OscParams.size(), Oscillator->GetOscParamsSize());
        MACH3LOG_ERROR("This indicate misconfiguration in your Osc yaml");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } else {
      InitialiseNuOscillatorObjects();
    }
    SetupNuOscillatorPointers();
  } else{
    MACH3LOG_WARN("Didn't find any oscillation params, thus will not enable oscillations");
    if(CheckNodeExists(SampleManager->raw(), "NuOsc")){
      MACH3LOG_ERROR("However config for SampleHandler {} has 'NuOsc' field", GetName());
      MACH3LOG_ERROR("This may indicate misconfiguration");
      MACH3LOG_ERROR("Either remove 'NuOsc' field from SampleHandler config or check your model.yaml and include oscillation for sample");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}


// ***************************************************************************
// Setup the norm parameters
void SampleHandlerBase::SetupNormParameters() {
// ***************************************************************************
  if(ParHandler == nullptr) return;
  std::vector< std::vector< int > > norms_bins(GetNEvents());

  std::vector< std::vector<NormParameter>> norm_parameters(GetNSamples());

  for (int iSample = 0; iSample < GetNSamples(); ++iSample) {
   norm_parameters[iSample] = ParHandler->GetNormParsFromSampleName(GetSampleName(iSample));
  }
  if(!ParHandler) {
    MACH3LOG_ERROR("ParHandler is not setup!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Assign norm bins in MCEvents tree
  CalcNormsBins(norm_parameters, norms_bins);

  //DB Attempt at reducing impact of SystematicHandlerGeneric::calcReweight()
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    size_t counter = 0;
    const size_t offset = MCEvents[iEvent].total_weight_pointers.size();
    const size_t addSize = norms_bins[iEvent].size();
    MCEvents[iEvent].total_weight_pointers.resize(offset + addSize);
    MCEvents[iEvent].total_weight_gradients.resize(offset + addSize);

    for(auto const & norm_bin: norms_bins[iEvent]) {
      MCEvents[iEvent].total_weight_pointers[offset + counter] = ParHandler->RetPointer(norm_bin);
      //HW A norm weight is just w = x, so dw/dx is unity. FillGradientArray applies the 1/w
      MCEvents[iEvent].total_weight_gradients[offset + counter] = {norm_bin, &M3::Unity};
      counter += 1;
    }
  }
}

// ************************************************
//A way to check whether a normalisation parameter applies to an event or not
void SampleHandlerBase::CalcNormsBins(std::vector <std::vector<NormParameter>>& norm_parameters, std::vector< std::vector< int > >& norms_bins) {
// ************************************************
  for(unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent){
    std::vector< int > NormBins = {};
    if (ParHandler) {
      const auto SampleId = MCEvents[iEvent].NominalSample;
      auto& NormParam = norm_parameters[SampleId];
      // Skip oscillated NC events
      // Not strictly needed, but these events don't get included in oscillated predictions, so
      // no need to waste our time calculating and storing information about xsec/flux etc. parameters
      // that will never be used.
      if (MCEvents[iEvent].isNC && (MCEvents[iEvent].nupdg != MCEvents[iEvent].nupdgUnosc) ) {
        MACH3LOG_TRACE("Event {}, missed NC/signal check", iEvent);
        continue;
      } //DB Abstract check on MaCh3Modes to determine which apply to neutral current
      for (std::vector<NormParameter>::iterator it = NormParam.begin(); it != NormParam.end(); ++it) {
        //Now check that the target of an interaction matches with the normalisation parameters
        const int Target = static_cast<int>(std::round(ReturnKinematicParameter("TargetNucleus", iEvent)));
        bool TargetMatch = MatchCondition(it->targets, Target);
        if (!TargetMatch) {
          MACH3LOG_TRACE("Event {}, missed target check ({}) for dial {}", iEvent, Target, it->name);
          continue;
        }

        //Now check that the neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourMatch = MatchCondition(it->pdgs, MCEvents[iEvent].nupdg);
        if (!FlavourMatch) {
          MACH3LOG_TRACE("Event {}, missed PDG check ({}) for dial {}", iEvent,MCEvents[iEvent].nupdg, it->name);
          continue;
        }

        //Now check that the unoscillated neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourUnoscMatch = MatchCondition(it->preoscpdgs, MCEvents[iEvent].nupdgUnosc);
        if (!FlavourUnoscMatch){
          MACH3LOG_TRACE("Event {}, missed FlavourUnosc check ({}) for dial {}", iEvent,MCEvents[iEvent].nupdgUnosc, it->name);
          continue;
        }

        //Now check that the mode of an interaction matches with the normalisation parameters
        const int Mode = static_cast<int>(std::round(ReturnKinematicParameter("Mode", iEvent)));
        bool ModeMatch = MatchCondition(it->modes, Mode);
        if (!ModeMatch) {
          MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, Mode, it->name);
          continue;
        }

        //Now check whether the norm has kinematic bounds
        //i.e. does it only apply to events in a particular kinematic region?
        // Now check whether within kinematic bounds
        bool IsSelected = PassesSelection((*it), iEvent);
        // Need to then break the event loop
        if(!IsSelected){
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", iEvent, it->name);
          continue;
        }
        // Now set 'index bin' for each normalisation parameter
        // All normalisations are just 1 bin for 2015, so bin = index (where index is just the bin for that normalisation)
        int bin = it->index;

        NormBins.push_back(bin);
        MACH3LOG_TRACE("Event {}, will be affected by dial {}", iEvent, it->name);
      } // end iteration over norm_parameters
    } // end if (ParHandler)
    norms_bins[iEvent] = NormBins;
  }//end loop over events
}

// ************************************************
void SampleHandlerBase::SetupReweightArrays() {
// ************************************************
  SampleHandler_array = std::vector<double>(Binning->GetNBins(),0);
  SampleHandler_array_w2 = std::vector<double>(Binning->GetNBins(),0);
  SampleHandler_data = std::vector<double>(Binning->GetNBins(),0);
  //HW Keep the gradient array in step with the binning if gradients were switched on before setup
  if (ComputeGradients) SetComputeGradients(true);
}

// ************************************************
void SampleHandlerBase::SetComputeGradients(const bool Enable) {
// ************************************************
  ComputeGradients = Enable;
  if (SplineHandler) SplineHandler->SetComputeDerivatives(Enable);
  // Switching gradients on pads the derivative lists, which the grouping check reads
  GroupingChecked = false;
  InvalidateKinematicCache();
  if (!Enable) {
    nGradPars = 0;
    SampleHandler_gradients.clear();
    SampleHandler_gradients.shrink_to_fit();
    return;
  }
  if (ParHandler == nullptr) {
    MACH3LOG_ERROR("Can't track gradients for {} without a parameter handler", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //HW No binning yet, SetupReweightArrays will come back to us once there is
  if (Binning == nullptr) return;

  /// @todo HW: UnbinnedSplineHandler doesn't expose spline derivatives, so gradient based fitters
  /// can't be used with it yet. See SetSplinePointers
  if (dynamic_cast<UnbinnedSplineHandler*>(SplineHandler.get()) != nullptr) {
    MACH3LOG_WARN("{} uses unbinned splines, whose derivatives are not propagated. The gradient will "
                  "be wrong for every spline parameter", GetName());
  }

  nGradPars = ParHandler->GetNumParams();
  if (nGradPars <= 0) {
    MACH3LOG_ERROR("Parameter handler {} has no parameters to take gradients with respect to", ParHandler->GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //HW Experiments push their own functional weights into total_weight_pointers in SetupWeightPointers,
  //which core never sees, so pad those with a zero derivative to keep the two arrays parallel
  size_t NMissingDerivs = 0;
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    auto& Pointers = MCEvents[iEvent].total_weight_pointers;
    auto& Derivs = MCEvents[iEvent].total_weight_gradients;
    if (Derivs.size() == Pointers.size()) continue;
    if (Derivs.size() > Pointers.size()) {
      MACH3LOG_ERROR("Event {} of {} has {} weight derivatives but only {} weights, so they can't be "
                     "matched up", iEvent, GetName(), Derivs.size(), Pointers.size());
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    NMissingDerivs += Pointers.size() - Derivs.size();
    Derivs.resize(Pointers.size(), {0, &M3::Zero});
  }
  if (NMissingDerivs > 0) {
    MACH3LOG_WARN("{} has {} event weights with no registered derivative, most likely functional "
                  "parameters. Their gradient is dropped", GetName(), NMissingDerivs);
  }

  //HW Bin major, so one bin's parameters stay contiguous
  SampleHandler_gradients.assign(static_cast<size_t>(Binning->GetNBins()) * static_cast<size_t>(nGradPars), 0.);
  MACH3LOG_INFO("{} is now tracking gradients for {} parameters over {} bins",
                GetName(), nGradPars, Binning->GetNBins());
}

// ************************************************
void SampleHandlerBase::SetBinning() {
// ************************************************
  for(int iSample = 0; iSample < GetNSamples(); iSample++)
  {
    int Dimension = GetNDim(iSample);
    std::string HistTitle = GetSampleTitle(iSample);

    auto* SamDet = &SampleDetails[iSample];
    if(Dimension == 1) {
      auto XVec = Binning->GetBinEdges(iSample, 0);
      SamDet->DataHist = new TH1D(("d" + HistTitle).c_str(), HistTitle.c_str(), static_cast<int>(XVec.size()-1), XVec.data());
      SamDet->MCHist   = new TH1D(("h" + HistTitle).c_str(), HistTitle.c_str(), static_cast<int>(XVec.size()-1), XVec.data());
      SamDet->W2Hist   = new TH1D(("w" + HistTitle).c_str(), HistTitle.c_str(), static_cast<int>(XVec.size()-1), XVec.data());

      // Set all titles so most of projections don't have empty titles...
      SamDet->DataHist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->DataHist->GetYaxis()->SetTitle("Events");
      SamDet->MCHist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->MCHist->GetYaxis()->SetTitle("Events");
      SamDet->W2Hist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->W2Hist->GetYaxis()->SetTitle("Events");
    } else if (Dimension == 2){
      if(Binning->IsUniform(iSample)) {
        auto XVec = Binning->GetBinEdges(iSample, 0);
        auto YVec = Binning->GetBinEdges(iSample, 1);
        int nX = static_cast<int>(XVec.size() - 1);
        int nY = static_cast<int>(YVec.size() - 1);

        SamDet->DataHist = new TH2D(("d" + HistTitle).c_str(), HistTitle.c_str(), nX, XVec.data(), nY, YVec.data());
        SamDet->MCHist   = new TH2D(("h" + HistTitle).c_str(), HistTitle.c_str(), nX, XVec.data(), nY, YVec.data());
        SamDet->W2Hist   = new TH2D(("w" + HistTitle).c_str(), HistTitle.c_str(), nX, XVec.data(), nY, YVec.data());
      } else {
        auto AddBinsToTH2Poly = [](TH2Poly* hist, const std::vector<BinInfo>& bins) {
          for (const auto& bin : bins) {
            double xLow  = bin.Extent[0][0];
            double xHigh = bin.Extent[0][1];
            double yLow  = bin.Extent[1][0];
            double yHigh = bin.Extent[1][1];

            double x[4] = {xLow, xHigh, xHigh, xLow};
            double y[4] = {yLow, yLow, yHigh, yHigh};

            hist->AddBin(4, x, y);
          }
        };
        // Create all three histograms
        SamDet->DataHist = new TH2Poly();
        SamDet->DataHist->SetName(("d" + HistTitle).c_str());
        SamDet->DataHist->SetTitle(HistTitle.c_str());

        SamDet->MCHist = new TH2Poly();
        SamDet->MCHist->SetName(("h" + HistTitle).c_str());
        SamDet->MCHist->SetTitle(HistTitle.c_str());

        SamDet->W2Hist = new TH2Poly();
        SamDet->W2Hist->SetName(("w" + HistTitle).c_str());
        SamDet->W2Hist->SetTitle(HistTitle.c_str());

        // Add bins to each
        AddBinsToTH2Poly(static_cast<TH2Poly*>(SamDet->DataHist), Binning->GetNonUniformBins(iSample));
        AddBinsToTH2Poly(static_cast<TH2Poly*>(SamDet->MCHist),   Binning->GetNonUniformBins(iSample));
        AddBinsToTH2Poly(static_cast<TH2Poly*>(SamDet->W2Hist),   Binning->GetNonUniformBins(iSample));
      }

      // Set all titles so most of projections don't have empty titles...
      SamDet->DataHist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->DataHist->GetYaxis()->SetTitle(SamDet->VarStr[1].c_str());
      SamDet->MCHist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->MCHist->GetYaxis()->SetTitle(SamDet->VarStr[1].c_str());
      SamDet->W2Hist->GetXaxis()->SetTitle(SamDet->VarStr[0].c_str());
      SamDet->W2Hist->GetYaxis()->SetTitle(SamDet->VarStr[1].c_str());
    } else {
      int nbins = Binning->GetNBins(iSample);
      SamDet->DataHist = new TH1D(("d" + HistTitle).c_str(), HistTitle.c_str(), nbins, 0, nbins);
      SamDet->MCHist   = new TH1D(("h" + HistTitle).c_str(), HistTitle.c_str(), nbins, 0, nbins);
      SamDet->W2Hist   = new TH1D(("w" + HistTitle).c_str(), HistTitle.c_str(), nbins, 0, nbins);

      for(int iBin = 0; iBin < nbins; iBin++) {
        auto BinName = Binning->GetBinName(iSample, iBin);
        SamDet->DataHist->GetXaxis()->SetBinLabel(iBin+1, BinName.c_str());
        SamDet->MCHist->GetXaxis()->SetBinLabel(iBin+1, BinName.c_str());
        SamDet->W2Hist->GetXaxis()->SetBinLabel(iBin+1, BinName.c_str());
      }

      // Set all titles so most of projections don't have empty titles...
      SamDet->DataHist->GetYaxis()->SetTitle("Events");
      SamDet->MCHist->GetYaxis()->SetTitle("Events");
      SamDet->W2Hist->GetYaxis()->SetTitle("Events");
    }

    SamDet->DataHist->SetDirectory(nullptr);
    SamDet->MCHist->SetDirectory(nullptr);
    SamDet->W2Hist->SetDirectory(nullptr);
  }

  //Set the number of X and Y bins now
  SetupReweightArrays();
  FindNominalBinAndEdges();
}

// ************************************************
void SampleHandlerBase::FindNominalBinAndEdges() {
// ************************************************
  for (unsigned int event_i = 0; event_i < GetNEvents(); event_i++) {
    int Sample = MCEvents[event_i].NominalSample;
    const int dim = GetNDim(Sample);
    MCEvents[event_i].KinVar.resize(dim);
    MCEvents[event_i].NomBin.resize(dim);

    auto SetNominalBin = [&](int bin, int max_bins, int& out_bin) {
      if (bin >= 0 && bin < max_bins) {
        out_bin = bin;
      } else {
        out_bin = M3::UnderOverFlowBin;  // Out of bounds
      }
    };

    // Find nominal bin for each dimension
    for(int iDim = 0; iDim < dim; iDim++) {
      MCEvents[event_i].KinVar[iDim] = GetPointerToKinematicParameter(GetKinVarName(Sample, iDim), event_i);
      if (std::isnan(*MCEvents[event_i].KinVar[iDim]) || std::isinf(*MCEvents[event_i].KinVar[iDim])) {
        MACH3LOG_ERROR("Variable {} for sample {} and dimension {} is ill-defined and equal to {}",
                       GetKinVarName(Sample, iDim), GetSampleTitle(Sample), dim, *MCEvents[event_i].KinVar[iDim]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      const int bin = Binning->FindNominalBin(Sample, iDim, *MCEvents[event_i].KinVar[iDim]);
      int NBins_i = static_cast<int>(Binning->GetBinEdges(Sample, iDim).size() - 1);
      SetNominalBin(bin, NBins_i, MCEvents[event_i].NomBin[iDim]);
    }
  }
}

// ************************************************
int SampleHandlerBase::GetSampleIndex(const std::string& SampleTitle) const {
// ************************************************
  for (M3::int_t iSample = 0; iSample < nSamples; ++iSample) {
    if (SampleTitle == GetSampleTitle(iSample)) {
      return iSample;
    }
  }
  MACH3LOG_ERROR("Sample name not found: {}", SampleTitle);
  throw MaCh3Exception(__FILE__, __LINE__);
}

// ************************************************
const TH1* SampleHandlerBase::GetW2Hist(const int Sample) {
// ************************************************
  FillHist(Sample, SampleDetails[Sample].W2Hist, SampleHandler_array_w2);
  if(SampleDetails[Sample].W2Hist == nullptr) {
    MACH3LOG_ERROR("Can't access {} for {}Dimensions", __func__, GetNDim(Sample));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleDetails[Sample].W2Hist;
}

// ************************************************
const TH1* SampleHandlerBase::GetW2Hist(const std::string& Sample) {
// ************************************************
  const int Index = GetSampleIndex(Sample);
  return GetW2Hist(Index);
}

// ************************************************
const TH1* SampleHandlerBase::GetMCHist(const int Sample) {
// ************************************************
  FillHist(Sample, SampleDetails[Sample].MCHist, SampleHandler_array);
  if(SampleDetails[Sample].MCHist == nullptr) {
    MACH3LOG_ERROR("Can't access {} for {}Dimensions", __func__, GetNDim(Sample));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleDetails[Sample].MCHist;
}

// ************************************************
const TH1* SampleHandlerBase::GetMCHist(const std::string& Sample) {
// ************************************************
  const int Index = GetSampleIndex(Sample);
  return GetMCHist(Index);
}

// ************************************************
const TH1* SampleHandlerBase::GetDataHist(const int Sample) {
// ************************************************
  if(SampleDetails[Sample].DataHist == nullptr) {
    MACH3LOG_ERROR("Can't access {} for {}Dimensions", __func__, GetNDim(Sample));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleDetails[Sample].DataHist;
}

// ************************************************
const TH1* SampleHandlerBase::GetDataHist(const std::string& Sample) {
// ************************************************
  int Index = GetSampleIndex(Sample);
  return GetDataHist(Index);
}

// ************************************************
void SampleHandlerBase::AddData(const int Sample, TH1* Data) {
// ************************************************
  ++DataVersion;
  int Dim = GetNDim(Sample);
  MACH3LOG_INFO("Adding {}D data histogram: {} with {:.2f} events", Dim, Data->GetTitle(), Data->Integral());
  // delete old histogram
  delete SampleDetails[Sample].DataHist;
  SampleDetails[Sample].DataHist = static_cast<TH1*>(Data->Clone());

  if(!SampleHandler_data.size()) {
    MACH3LOG_ERROR("SampleHandler_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto ChecHistType = [&](const std::string& Type, const int Dimen, const TH1* Hist,
                          const std::string& file, const int line) {
    if (std::string(Hist->ClassName()) != Type) {
      MACH3LOG_ERROR("Expected {} for {}D sample, got {}", Type, Dimen, Hist->ClassName());
      throw MaCh3Exception(file, line);
    }
  };

  if (Dim == 1) {
    // Ensure we really have a TH1D
    ChecHistType("TH1D", Dim, SampleDetails[Sample].DataHist, __FILE__, __LINE__);
    M3::CheckBinningMatch(static_cast<TH1D*>(SampleDetails[Sample].DataHist),
                          static_cast<TH1D*>(SampleDetails[Sample].MCHist), __FILE__, __LINE__);
    for (int xBin = 0; xBin < Binning->GetNAxisBins(Sample, 0); ++xBin) {
      const int idx = Binning->GetGlobalBinSafe(Sample, {xBin});
      // ROOT histograms are 1-based, so bin index + 1
      SampleHandler_data[idx] = SampleDetails[Sample].DataHist->GetBinContent(xBin + 1);
    }
    SampleDetails[Sample].DataHist->GetXaxis()->SetTitle(GetKinVarName(Sample, 0).c_str());
    SampleDetails[Sample].DataHist->GetYaxis()->SetTitle("Number of Events");
  } else if (Dim == 2) {
    if(Binning->IsUniform(Sample)) {
      ChecHistType("TH2D", Dim, SampleDetails[Sample].DataHist, __FILE__, __LINE__);
      M3::CheckBinningMatch(static_cast<TH2D*>(SampleDetails[Sample].DataHist),
                            static_cast<TH2D*>(SampleDetails[Sample].MCHist), __FILE__, __LINE__);
      for (int yBin = 0; yBin < Binning->GetNAxisBins(Sample, 1); ++yBin) {
        for (int xBin = 0; xBin < Binning->GetNAxisBins(Sample, 0); ++xBin) {
          const int idx = Binning->GetGlobalBinSafe(Sample, {xBin, yBin});
          //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
          SampleHandler_data[idx] = SampleDetails[Sample].DataHist->GetBinContent(xBin + 1, yBin + 1);
        }
      }
    } else{
      ChecHistType("TH2Poly", Dim, SampleDetails[Sample].DataHist, __FILE__, __LINE__);
      M3::CheckBinningMatch(static_cast<TH2Poly*>(SampleDetails[Sample].DataHist),
                            static_cast<TH2Poly*>(SampleDetails[Sample].MCHist), __FILE__, __LINE__);
      for (int iBin = 0; iBin < Binning->GetNBins(Sample); ++iBin) {
        const int idx = iBin + Binning->GetSampleStartBin(Sample);
        //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
        SampleHandler_data[idx] = SampleDetails[Sample].DataHist->GetBinContent(iBin + 1);
      }
    }

    SampleDetails[Sample].DataHist->GetXaxis()->SetTitle(GetKinVarName(Sample, 0).c_str());
    SampleDetails[Sample].DataHist->GetYaxis()->SetTitle(GetKinVarName(Sample, 1).c_str());
    SampleDetails[Sample].DataHist->GetZaxis()->SetTitle("Number of Events");
  } else {
    ChecHistType("TH1D", Dim, SampleDetails[Sample].DataHist, __FILE__, __LINE__);
    M3::CheckBinningMatch(static_cast<TH1D*>(SampleDetails[Sample].DataHist),
                          static_cast<TH1D*>(SampleDetails[Sample].MCHist), __FILE__, __LINE__);
    for (int iBin = 0; iBin < Binning->GetNBins(Sample); ++iBin) {
      const int idx = iBin + Binning->GetSampleStartBin(Sample);
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      SampleHandler_data[idx] = SampleDetails[Sample].DataHist->GetBinContent(iBin + 1);
    }
  }
}

// ************************************************
void SampleHandlerBase::AddData(const int Sample, const std::vector<double>& Data_Array) {
// ************************************************
  ++DataVersion;
  const int Start = Binning->GetSampleStartBin(Sample);
  const int End = Binning->GetSampleEndBin(Sample);
  const int ExpectedSize = End - Start;

  if (static_cast<int>(Data_Array.size()) != ExpectedSize) {
    MACH3LOG_ERROR("Data_Array size mismatch for sample '{}'.", GetSampleTitle(Sample));
    MACH3LOG_ERROR("Expected size: {}, received size: {}.", ExpectedSize, Data_Array.size());
    MACH3LOG_ERROR("Start bin: {}, End bin: {}.", Start, End);
    MACH3LOG_ERROR("This likely indicates a binning or sample slicing bug.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::copy_n(Data_Array.begin(), End-Start, SampleHandler_data.begin() + Start);

  FillHist(Sample, SampleDetails[Sample].DataHist, SampleHandler_data);
}

// ************************************************
void SampleHandlerBase::InitialiseNuOscillatorObjects() {
  // ************************************************
  auto NuOscillatorConfigFile = Get<std::string>(SampleManager->raw()["NuOsc"]["NuOscConfigFile"], __FILE__ , __LINE__);
  auto EqualBinningPerOscChannel = Get<bool>(SampleManager->raw()["NuOsc"]["EqualBinningPerOscChannel"], __FILE__ , __LINE__);

  // TN override the sample setting if not using binned oscillation
  if (EqualBinningPerOscChannel) {
    if (YAML::LoadFile(NuOscillatorConfigFile)["General"]["CalculationType"].as<std::string>() == "Unbinned") {
      MACH3LOG_WARN("Tried using EqualBinningPerOscChannel while using Unbinned oscillation calculation, changing EqualBinningPerOscChannel to false");
      EqualBinningPerOscChannel = false;
    }
  }

  std::vector<OscillationParameter> OscParams{};
  for (int iSample = 0; iSample < GetNSamples(); iSample++) {
    auto sample_osc_pars = ParHandler->GetOscParsFromSampleName(GetSampleName(iSample));
    if (!iSample) {
      if (sample_osc_pars.empty()) {
        MACH3LOG_ERROR("OscParams is empty for SampleHandler/Sample {}/{}.", GetName(), GetSampleName(iSample));
        MACH3LOG_ERROR("This indicates an error in your oscillation YAML configuration.");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      OscParams = sample_osc_pars;
    } else {
      if (sample_osc_pars.size() != OscParams.size()) {
        MACH3LOG_ERROR("SampleHandler/Sample {}/{} has {} osc params while Sample {} has {}",
                        GetName(), GetSampleName(iSample), sample_osc_pars.size(), 0, GetSampleName(0));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }

  std::vector<const M3::float_t*> OscParamsValues(OscParams.size());
  std::vector<std::string> NuOscName(OscParams.size());
  for(size_t ij = 0; ij < OscParams.size(); ij++){
    OscParamsValues[ij] = ParHandler->RetPointer(OscParams[ij].index);
    NuOscName[ij] = OscParams[ij].NuOscName;
  }

  if (EqualBinningPerOscChannel) {
    Oscillator = std::make_shared<OscillationHandler>(NuOscillatorConfigFile, EqualBinningPerOscChannel, OscParamsValues,
                                                      NuOscName, GetNOscChannels(0));
  } else {
    // aim is to build energy/cosz arrays per oscillation channel so that
    // NuOscillator can decide on the binning
    std::vector<std::map<std::pair<int, int>, std::vector<M3::float_t>>>
    channel_energy_array(GetNSamples()),
    channel_cosz_array(GetNSamples());
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
      if (MCEvents[iEvent].NominalSample >= GetNSamples()) {
        MACH3LOG_ERROR("Encountered Event with NominalSample: {}, but "
        "SampleHandler {} only has {} samples.",
        MCEvents[iEvent].NominalSample, GetName(),
                        GetNSamples());
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      if (MCEvents[iEvent].isNC) { // skip NC events
        continue;
      }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
      channel_energy_array[MCEvents[iEvent].NominalSample] // sample
      [std::make_pair(MCEvents[iEvent].nupdgUnosc,
                      MCEvents[iEvent].nupdg)] // osc chan
                      .push_back(M3::float_t(MCEvents[iEvent].enu_true)); // energy

      // DB Atmospheric only part, can only happen if truecz has been
      // initialised within the experiment specific code
      if (MCEvents[0].coszenith_true != M3::_BAD_DOUBLE_) {
        channel_cosz_array[MCEvents[iEvent].NominalSample] // sample
        [std::make_pair(MCEvents[iEvent].nupdgUnosc,
                        MCEvents[iEvent].nupdg)] // osc chan
                        .push_back(M3::float_t(MCEvents[iEvent].coszenith_true)); // cosz
      }
#pragma GCC diagnostic pop
    }

    for (int iSample = 0; iSample < GetNSamples(); iSample++) {
      if (!channel_energy_array[iSample].size()) { // sample had no events
        continue;
      }

      // loop through channels and check they were defined in config, or add
      // them if not
      for (auto const &[osc_flavs, energies] : channel_energy_array[iSample]) {
        bool found = false;
        for (size_t i = 0; i < SampleDetails[iSample].OscChannels.size(); ++i) {
          if (osc_flavs.first == SampleDetails[iSample].OscChannels[i].InitPDG &&
            osc_flavs.second == SampleDetails[iSample].OscChannels[i].FinalPDG) {
            found = true;
            break;
          }
        }
        if(found) { // already have this channel defined.
          continue;
        }

        static const std::map<int, std::pair<std::string, std::string>> spec_defnames = {
          {12, {"nue", "#nu_{e}"}},
          {-12, {"nueb", "#bar{#nu}_{e}"}},
          {14, {"numu", "#nu_{#mu}"}},
          {-14, {"nuemub", "#bar{#nu}_{#mu}"}},
          {16, {"nutau", "#nu_{#tau}"}},
          {-16, {"nutaub", "#bar{#nu}_{#tau}"}},
        };

        OscChannelInfo oci;
        oci.InitPDG = osc_flavs.first;
        oci.FinalPDG = osc_flavs.second;
        oci.flavourName = fmt::format("{}_x_{}", spec_defnames.at(osc_flavs.first).first,
                                      spec_defnames.at(osc_flavs.second).first);
        oci.flavourName_Latex =
        fmt::format("{} #rightarrow {}", spec_defnames.at(osc_flavs.first).second,
                    spec_defnames.at(osc_flavs.second).second);
        oci.ChannelIndex = double(SampleDetails[iSample].OscChannels.size());

        SampleDetails[iSample].OscChannels.push_back(oci);
      }

      // set up oscillator for this sample
      if (!Oscillator) {
        Oscillator = std::make_shared<OscillationHandler>(NuOscillatorConfigFile, EqualBinningPerOscChannel, OscParamsValues,
                                                          NuOscName, GetNOscChannels(iSample));
      } else {
        Oscillator->AddSample(NuOscillatorConfigFile, GetNOscChannels(iSample));
      }

      // set oscillator binning
      for (auto &[osc_flavs, energies] : channel_energy_array[iSample]) {
        auto iChannel = GetOscChannel(SampleDetails[iSample].OscChannels,
                                      osc_flavs.first, osc_flavs.second);

        std::sort(energies.begin(), energies.end());
        std::sort(channel_cosz_array[iSample][osc_flavs].begin(),
                  channel_cosz_array[iSample][osc_flavs].end());
        Oscillator->SetOscillatorBinning(iSample, iChannel, energies,
                                         channel_cosz_array[iSample][osc_flavs]);
      }
    } // end loop over samples
  } // end EqualBinningPerOscChannel
}

// ************************************************
void SampleHandlerBase::SetupNuOscillatorPointers() {
// ************************************************
  auto AddOscPointer = GetFromManager<bool>(SampleManager->raw()["NuOsc"]["AddOscPointer"], true, __FILE__ , __LINE__);
  // Sometimes we don't want to add osc pointer to allow for some experiment specific handling of osc weight, like for example being able to shift osc weigh
  if(!AddOscPointer) {
    return;
  }
  for (unsigned int iEvent=0;iEvent<GetNEvents();iEvent++) {
   const M3::float_t* osc_w_pointer = GetNuOscillatorPointers(iEvent);

   // KS: Do not add unity
   if (osc_w_pointer != &M3::Unity) {
     MCEvents[iEvent].total_weight_pointers.push_back(osc_w_pointer);
     //HW NuOscillator doesn't hand back derivatives, so oscillation parameter gradients are dropped
     MCEvents[iEvent].total_weight_gradients.push_back({0, &M3::Zero});
   }
  }
}

// ************************************************
const M3::float_t* SampleHandlerBase::GetNuOscillatorPointers(const int iEvent) const {
// ************************************************
  const M3::float_t* osc_w_pointer = &M3::Unity;
  if (MCEvents[iEvent].isNC) {
    if (MCEvents[iEvent].nupdg != MCEvents[iEvent].nupdgUnosc) {
      osc_w_pointer = &M3::Zero;
    } else {
      osc_w_pointer = &M3::Unity;
    }
  } else {
    int InitFlav = M3::_BAD_INT_;
    int FinalFlav = M3::_BAD_INT_;

    InitFlav =  M3::Utils::PDGToNuOscillatorFlavour(MCEvents[iEvent].nupdgUnosc);
    FinalFlav = M3::Utils::PDGToNuOscillatorFlavour(MCEvents[iEvent].nupdg);

    if (InitFlav == M3::_BAD_INT_ || FinalFlav == M3::_BAD_INT_) {
      MACH3LOG_ERROR("Something has gone wrong in the mapping between MCEvents.nutype and the enum used within NuOscillator");
      MACH3LOG_ERROR("MCEvents.nupdgUnosc: {}", MCEvents[iEvent].nupdgUnosc);
      MACH3LOG_ERROR("InitFlav: {}", InitFlav);
      MACH3LOG_ERROR("MCEvents.nupdg: {}", MCEvents[iEvent].nupdg);
      MACH3LOG_ERROR("FinalFlav: {}", FinalFlav);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const int Sample = MCEvents[iEvent].NominalSample;

    const int OscIndex = GetOscChannel(SampleDetails[Sample].OscChannels, MCEvents[iEvent].nupdgUnosc, MCEvents[iEvent].nupdg);
    //Can only happen if truecz has been initialised within the experiment specific code
    if (MCEvents[iEvent].coszenith_true != M3::_BAD_DOUBLE_) {
      //Atmospherics
      osc_w_pointer = Oscillator->GetNuOscillatorPointers(Sample, OscIndex, InitFlav, FinalFlav, FLOAT_T(MCEvents[iEvent].enu_true), FLOAT_T(MCEvents[iEvent].coszenith_true));
    } else {
      //Beam
      osc_w_pointer = Oscillator->GetNuOscillatorPointers(Sample, OscIndex, InitFlav, FinalFlav, FLOAT_T(MCEvents[iEvent].enu_true));
    }
  } // end if NC
  return osc_w_pointer;
}

// ************************************************
std::string SampleHandlerBase::GetName() const {
// ************************************************
  //ETA - extra safety to make sure SampleHandlerName is actually set
  // probably unnecessary due to the requirement for it to be in the yaml config
  if(SampleHandlerName.length() == 0) {
    MACH3LOG_ERROR("No sample name provided");
    MACH3LOG_ERROR("Please provide a SampleHandlerName in your configuration file: {}", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleHandlerName;
}

// ************************************************
M3::float_t SampleHandlerBase::GetEventWeight(const int iEntry) {
// ************************************************
  /// KS: @warning we have to here recalculate weight and cap because there is possibility weight wasn't calculated during FillArray because it didn't fulfil IsEventSelected

  // Virtual by default does nothing, has to happen before CalcWeightTotal
  CalcWeightFunc(iEntry);

  const EventInfo* _restrict_ MCEvent = &MCEvents[iEntry];
  M3::float_t totalweight = CalcWeightTotal(MCEvent);

  //DB Catch negative total weights and skip any event with a negative weight. Previously we would set weight to zero and continue but that is inefficient
  if (totalweight <= 0.){
    totalweight = 0.;
  }
  return totalweight;
}

// ************************************************
std::vector< SplineIndex > SampleHandlerBase::GetSplineBins(int Event, BinnedSplineHandler* BinnedSpline, bool& ThrowCrititcal) const {
// ************************************************
  const int SampleIndex = MCEvents[Event].NominalSample;
  const auto SampleTitle = GetSampleTitle(SampleIndex);
  bool NoOscChannels = false;
  if(Oscillator == nullptr && GetNOscChannels(SampleIndex) == 1) {
    MACH3LOG_DEBUG("Assuming there are no osc channels in {}", __func__);
    NoOscChannels = true;
  }
  const int OscIndex = NoOscChannels ? 0 : GetOscChannel(SampleDetails[SampleIndex].OscChannels,
                                                         MCEvents[Event].nupdgUnosc, MCEvents[Event].nupdg);
  const int Mode = static_cast<int>(std::round(ReturnKinematicParameter("Mode", Event)));
  const double Etrue = MCEvents[Event].enu_true;
  std::vector< SplineIndex > EventSplines;
  switch(GetNDim(SampleIndex)) {
    case 1:
      EventSplines = BinnedSpline->GetEventSplines(SampleTitle, OscIndex, Mode, {Etrue, *(MCEvents[Event].KinVar[0]), 0.});
      break;
    case 2:
      EventSplines = BinnedSpline->GetEventSplines(SampleTitle, OscIndex, Mode, {Etrue, *(MCEvents[Event].KinVar[0]), *(MCEvents[Event].KinVar[1])});
      break;
    default:
      if(ThrowCrititcal) {
        MACH3LOG_CRITICAL("MaCh3 doesn't support binned splines for more than 2D while you use {}", GetNDim(SampleIndex));
        MACH3LOG_CRITICAL("Will use 2D like approach");
        ThrowCrititcal = false;
      }
      EventSplines = BinnedSpline->GetEventSplines(SampleTitle, OscIndex, Mode, {Etrue, *(MCEvents[Event].KinVar[0]), *(MCEvents[Event].KinVar[1])});
      break;
  }
  return EventSplines;
}

// ************************************************
void SampleHandlerBase::SetSplinePointers() {
// ************************************************
  //Now loop over events and get the spline bin for each event
  if (auto BinnedSpline = dynamic_cast<BinnedSplineHandler*>(SplineHandler.get())) {
    bool ThrowCrititcal = true;

    std::vector< std::vector<SplineParameter> > SplineParsVec(GetNSamples());
    for (int iSample = 0; iSample < GetNSamples(); ++iSample) {
      SplineParsVec[iSample] = ParHandler->GetSplineParsFromSampleName(GetSampleName(iSample));
    }
    for (unsigned int j = 0; j < GetNEvents(); ++j) {
      auto EventSplines = GetSplineBins(j, BinnedSpline, ThrowCrititcal);
      const int NSplines = static_cast<int>(EventSplines.size());
      if(NSplines == 0) continue;
      auto& w_pointers = MCEvents[j].total_weight_pointers;
      auto& w_grad_pointers = MCEvents[j].total_weight_gradients;
      w_pointers.reserve(w_pointers.size() + NSplines);
      w_grad_pointers.reserve(w_grad_pointers.size() + static_cast<size_t>(NSplines));
      const auto SampleId = MCEvents[j].NominalSample;
      for(int spline = 0; spline < NSplines; spline++) {
        int SystIndex = EventSplines[spline].iSyst;

        bool IsSelected = PassesSelection(SplineParsVec[SampleId][SystIndex], j);
        // Need to then break the event loop
        if(!IsSelected){
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", j, SplineParsVec[SampleId][SystIndex].name);
          continue;
        }
        //Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
        w_pointers.push_back(BinnedSpline->RetPointer(EventSplines[spline]));

        //HW Same for the gradients, each spline only ever affects one systematic. iSyst indexes this
        //sample's spline params, so take .index off it for the global one the gradient array uses
        const int GlobalParIndex = SplineParsVec[static_cast<size_t>(SampleId)][static_cast<size_t>(SystIndex)].index;
        w_grad_pointers.push_back({GlobalParIndex, BinnedSpline->RetGradPointer(EventSplines[static_cast<size_t>(spline)])});
      } // end loop over splines
      w_pointers.shrink_to_fit();
      w_grad_pointers.shrink_to_fit();
    } // end loop over events
  } else if (auto UnbinnedSpline = dynamic_cast<UnbinnedSplineHandler*>(SplineHandler.get())) {
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
      MCEvents[iEvent].total_weight_pointers.push_back(UnbinnedSpline->RetPointer(iEvent));
      //HW Unbinned splines expose no derivative, so this weight's gradient is dropped
      MCEvents[iEvent].total_weight_gradients.push_back({0, &M3::Zero});
    }
  } else {
    MACH3LOG_ERROR("Not supported splines");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
double SampleHandlerBase::GetSampleLikelihood(const int isample) const {
// ************************************************
  const int Start = Binning->GetSampleStartBin(isample);
  const int End = Binning->GetSampleEndBin(isample);

  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:negLogL)
  #endif
  for (int idx = Start; idx < End; ++idx)
  {
    double const &DataVal = SampleHandler_data[idx];
    double const &MCPred = SampleHandler_array[idx];
    double const &w2 = SampleHandler_array_w2[idx];

    //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
    negLogL += GetTestStatLLH(DataVal, MCPred, w2);
  }
  return negLogL;
}

// ************************************************
double SampleHandlerBase::GetLikelihood() const {
// ************************************************
  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:negLogL)
  #endif
  for (int idx = 0; idx < Binning->GetNBins(); ++idx)
  {
    double const &DataVal = SampleHandler_data[idx];
    double const &MCPred = SampleHandler_array[idx];
    double const &w2 = SampleHandler_array_w2[idx];

    //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
    negLogL += GetTestStatLLH(DataVal, MCPred, w2);
  }
  return negLogL;
}

// ************************************************
std::vector<double> SampleHandlerBase::GetGradLikelihood() const {
// ************************************************
  if (!ComputeGradients) {
    MACH3LOG_ERROR("Asked {} for the likelihood gradient without calling SetComputeGradients(true) first", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if (!GradientsCurrent) {
    MACH3LOG_ERROR("Asked {} for the likelihood gradient, but the last reweight had gradients paused, "
                   "see SetGradientsPaused", GetName());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //HW The pre-factor below is the Poisson derivative, and the MC stat statistics would also pick up
  //a w2 dependence we don't track derivatives for
  if (fTestStatistic != kPoisson) {
    MACH3LOG_ERROR("Likelihood gradients are only implemented for the {} test statistic, you asked for {}",
                   TestStatistic_ToString(kPoisson), TestStatistic_ToString(fTestStatistic));
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  //HW The device reweight sums the gradient itself, with the same pre-factor and the same skipped bins
  #ifdef MaCh3_CUDA
  if (GradFromDevice) return DeviceGrad;
  #endif

  const int NBins = Binning->GetNBins();
  std::vector<double> gradLogL(static_cast<size_t>(nGradPars), 0.);
  double* _restrict_ gradOut = gradLogL.data();

  //HW Not multithreaded over bins, that would need a reduction over the whole gradLogL array
  for (int idx = 0; idx < NBins; ++idx)
  {
    const double DataVal = SampleHandler_data[static_cast<size_t>(idx)];
    const double MCPred = SampleHandler_array[static_cast<size_t>(idx)];
    //HW GetPoissonLLH clamps the prediction below this bound, so the test statistic is flat there
    if (MCPred < M3::_LOW_MC_BOUND_) continue;

    // -logL = MC - D + D*log(D/MC), so d(-logL)/d(MC) is this pre-factor
    const double PreFactor = 1. - (DataVal / MCPred);
    //HW Flat array, so this bin's nGradPars entries sit next to each other
    const double* _restrict_ Grad = SampleHandler_gradients.data() + static_cast<size_t>(idx) * static_cast<size_t>(nGradPars);

    for (int iPar = 0; iPar < nGradPars; ++iPar) {
      gradOut[iPar] += PreFactor * Grad[iPar];
    }
  }
  return gradLogL;
}

// ************************************************
void SampleHandlerBase::SaveAdditionalInfo(TDirectory* Dir) {
// ************************************************
  Dir->cd();

  YAML::Node Config = SampleManager->raw();
  TMacro ConfigSave = YAMLtoTMacro(Config, (std::string("Config_") + GetName()));
  ConfigSave.Write();

  for(int iSample = 0; iSample < GetNSamples(); iSample++)
  {
    std::unique_ptr<TH1> data_hist;

    if (GetNDim(iSample) == 1) {
      data_hist = M3::Clone<TH1D>(dynamic_cast<const TH1D*>(GetDataHist(iSample)), "data_" + GetSampleTitle(iSample));
      data_hist->GetXaxis()->SetTitle(GetKinVarName(iSample, 0).c_str());
      data_hist->GetYaxis()->SetTitle("Number of Events");
    } else if (GetNDim(iSample) == 2) {
      if(Binning->IsUniform(iSample)) {
        data_hist = M3::Clone<TH2D>(dynamic_cast<const TH2D*>(GetDataHist(iSample)), "data_" + GetSampleTitle(iSample));
      } else {
        data_hist = M3::Clone<TH2Poly>(dynamic_cast<const TH2Poly*>(GetDataHist(iSample)), "data_" + GetSampleTitle(iSample));
      }
      data_hist->GetXaxis()->SetTitle(GetKinVarName(iSample, 0).c_str());
      data_hist->GetYaxis()->SetTitle(GetKinVarName(iSample, 1).c_str());
      data_hist->GetZaxis()->SetTitle("Number of Events");
    } else {
      data_hist = M3::Clone<TH1D>(dynamic_cast<const TH1D*>(GetDataHist(iSample)), "data_" + GetSampleTitle(iSample));
      int nbins = Binning->GetNBins(iSample);
      for(int iBin = 0; iBin < nbins; iBin++) {
        auto BinName = Binning->GetBinName(iSample, iBin);
        data_hist->GetXaxis()->SetBinLabel(iBin+1, BinName.c_str());
      }
      data_hist->GetYaxis()->SetTitle("Number of Events");
    }

    if (!data_hist) {
      MACH3LOG_ERROR("nullptr data hist :(");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    data_hist->SetTitle(("data_" + GetSampleTitle(iSample)).c_str());
    data_hist->Write();
  }
}

// ************************************************
void SampleHandlerBase::InitialiseSplineObject() {
// ************************************************
  if(auto BinnedSplines = dynamic_cast<BinnedSplineHandler*>(SplineHandler.get())) {
    bool LoadSplineFile = GetFromManager<bool>(SampleManager->raw()["InputFiles"]["LoadSplineFile"], false, __FILE__, __LINE__);
    bool PrepSplineFile = GetFromManager<bool>(SampleManager->raw()["InputFiles"]["PrepSplineFile"], false, __FILE__, __LINE__);
    auto SplineFileName = GetFromManager<std::string>(SampleManager->raw()["InputFiles"]["SplineFileName"],
                                                      (SampleHandlerName + "_SplineFile.root"), __FILE__, __LINE__);
    if(!LoadSplineFile) {
      for(int iSample = 0; iSample < GetNSamples(); iSample++) {
        std::vector<std::string> spline_filepaths = SampleDetails[iSample].spline_files;

        //Keep a track of the spline variables
        std::vector<std::string> SplineVarNames = {"TrueNeutrinoEnergy"};
        if (GetNDim(iSample) == 1) {
          SplineVarNames.push_back(GetKinVarName(iSample, 0));
        } else if (GetNDim(iSample) == 2) {
          SplineVarNames.push_back(GetKinVarName(iSample, 0));
          SplineVarNames.push_back(GetKinVarName(iSample, 1));
        } else {
          MACH3LOG_CRITICAL("{} Not implemented for dimension {}, will use 2D", __func__, GetNDim(iSample));
          SplineVarNames.push_back(GetKinVarName(iSample, 0));
          SplineVarNames.push_back(GetKinVarName(iSample, 1));
        }
        BinnedSplines->AddSample(GetSampleName(iSample), GetSampleTitle(iSample), spline_filepaths, SplineVarNames);
      }
      BinnedSplines->CountNumberOfLoadedSplines(false, 1);
      BinnedSplines->TransferToMonolith();
      if(PrepSplineFile) BinnedSplines->PrepareSplineFile(SplineFileName);
    } else {
      // KS: Skip default spline loading and use flattened spline format allowing to read stuff much faster
      BinnedSplines->LoadSplineFile(SplineFileName);
    }
    MACH3LOG_INFO("--------------------------------");
    MACH3LOG_INFO("Setup Far Detector splines");

    SetSplinePointers();

    BinnedSplines->CleanUpMemory();
  } else if (dynamic_cast<UnbinnedSplineHandler*>(SplineHandler.get())) {
    SetSplinePointers();
  } else {
    MACH3LOG_ERROR("Unsupported spline type encountered.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
std::vector<std::vector<KinematicCut>> SampleHandlerBase::ApplyTemporarySelection(const int iSample,
                                         const std::vector<KinematicCut>& ExtraCuts) {
// ************************************************
  // Backup current selection
  auto originalSelection = Selection;

  //DB Add all the predefined selections to the selection vector which will be applied
  auto selectionToApply = Selection;

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (const auto& cut : ExtraCuts) {
    selectionToApply[iSample].emplace_back(cut);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = std::move(selectionToApply);

  return originalSelection;
}

// ************************************************
// === JM adjust GetNDVarHist functions to allow for subevent-level plotting ===
std::unique_ptr<TH1> SampleHandlerBase::Get1DVarHist(const int iSample, const std::string& ProjectionVar_Str,
                                                     const std::vector< KinematicCut >& EventSelectionVec,
                                                     const int WeightStyle,
                                                     const std::vector< KinematicCut >& SubEventSelectionVec) {
// ************************************************
  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  auto tmp_Selection = ApplyTemporarySelection(iSample, EventSelectionVec);

  //DB Define the histogram which will be returned
  std::vector<double> xBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_Str);
  auto _h1DVar = std::make_unique<TH1D>("", "", int(xBinEdges.size())-1, xBinEdges.data());
  _h1DVar->SetDirectory(nullptr);
  _h1DVar->GetXaxis()->SetTitle(ProjectionVar_Str.c_str());
  _h1DVar->GetYaxis()->SetTitle("Events");

  if (IsSubEventVarString(ProjectionVar_Str)) {
    Fill1DSubEventHist(iSample, _h1DVar.get(), ProjectionVar_Str, SubEventSelectionVec, WeightStyle);
  } else {
    //DB Grab the associated enum with the argument string
    int ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);

    //DB Loop over all events
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
      const int EventSample = MCEvents[iEvent].NominalSample;
      if(EventSample != iSample) continue;
      if (IsEventSelected(EventSample, iEvent)) {
        double Weight = GetEventWeight(iEvent);
        if (WeightStyle == 1) {
          Weight = 1.;
        }
        double Var = ReturnKinematicParameter(ProjectionVar_Int, iEvent);
        _h1DVar->Fill(Var, Weight);
      }
    }
  }
  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h1DVar;
}

// ************************************************
void SampleHandlerBase::Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar_Str,
                                           const std::vector< KinematicCut >& SubEventSelectionVec, const int WeightStyle) {
// ************************************************
  int ProjectionVar_Int = ReturnKinematicVectorFromString(ProjectionVar_Str);

  //JM Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    const int EventSample = MCEvents[iEvent].NominalSample;
    if(EventSample != iSample) continue;
    if (IsEventSelected(EventSample, iEvent)) {
      double Weight = GetEventWeight(iEvent);
      if (WeightStyle == 1) {
        Weight = 1.;
      }
      std::vector<double> Vec = ReturnKinematicVector(ProjectionVar_Int,iEvent);
      size_t nsubevents = Vec.size();
      //JM Loop over all subevents in event
      for (unsigned int iSubEvent = 0; iSubEvent < nsubevents; iSubEvent++) {
        if (IsSubEventSelected(SubEventSelectionVec, iEvent, iSubEvent, nsubevents)) {
          double Var = Vec[iSubEvent];
          _h1DVar->Fill(Var,Weight);
        }
      }
    }
  }
}

// ************************************************
std::unique_ptr<TH2> SampleHandlerBase::Get2DVarHist(const int iSample,
                                                     const std::string& ProjectionVar_StrX,
                                                     const std::string& ProjectionVar_StrY,
                                                     const std::vector< KinematicCut >& EventSelectionVec,
                                                     const int WeightStyle,
                                                     const std::vector< KinematicCut >& SubEventSelectionVec) {
// ************************************************
  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  auto tmp_Selection = ApplyTemporarySelection(iSample, EventSelectionVec);

  //DB Define the histogram which will be returned
  //KS: If we use 2D non uniform binning and wanting to plot fit variables use TH2Poly everything else uses TH2D with binning from config
  std::unique_ptr<TH2> _h2DVar;
  if(GetNDim(iSample) == 2 && !Binning->IsUniform(iSample) &&
     ProjectionVar_StrX == GetKinVarName(iSample, 0) && ProjectionVar_StrY == GetKinVarName(iSample, 1) ) {
    _h2DVar = std::unique_ptr<TH2>(static_cast<TH2*>(M3::Clone(GetMCHist(iSample)).release()));
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_StrX);
    std::vector<double> yBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_StrY);
    _h2DVar = std::make_unique<TH2D>("", "", int(xBinEdges.size())-1, xBinEdges.data(), int(yBinEdges.size())-1, yBinEdges.data());
  }
  _h2DVar->SetDirectory(nullptr);
  _h2DVar->GetXaxis()->SetTitle(ProjectionVar_StrX.c_str());
  _h2DVar->GetYaxis()->SetTitle(ProjectionVar_StrY.c_str());
  _h2DVar->GetZaxis()->SetTitle("Events");

  bool IsSubEventHist = IsSubEventVarString(ProjectionVar_StrX) || IsSubEventVarString(ProjectionVar_StrY);
  if (IsSubEventHist) Fill2DSubEventHist(iSample, _h2DVar.get(), ProjectionVar_StrX, ProjectionVar_StrY, SubEventSelectionVec, WeightStyle);
  else {
    //DB Grab the associated enum with the argument string
    int ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
    int ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

    //DB Loop over all events
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
      const int EventSample = MCEvents[iEvent].NominalSample;
      if(EventSample != iSample) continue;
      if (IsEventSelected(EventSample, iEvent)) {
        double Weight = GetEventWeight(iEvent);
        if (WeightStyle == 1) {
          Weight = 1.;
        }
        double VarX = ReturnKinematicParameter(ProjectionVar_IntX, iEvent);
        double VarY = ReturnKinematicParameter(ProjectionVar_IntY, iEvent);
        _h2DVar->Fill(VarX,VarY,Weight);
      }
    }
  }
  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h2DVar;
}

// ************************************************
void SampleHandlerBase::Fill2DSubEventHist(const int iSample, TH2* _h2DVar,
                                         const std::string& ProjectionVar_StrX,
                                         const std::string& ProjectionVar_StrY,
                                         const std::vector< KinematicCut >& SubEventSelectionVec,
                                         int WeightStyle) {
// ************************************************
  bool IsSubEventVarX = IsSubEventVarString(ProjectionVar_StrX);
  bool IsSubEventVarY = IsSubEventVarString(ProjectionVar_StrY);

  int ProjectionVar_IntX, ProjectionVar_IntY;
  if (IsSubEventVarX) ProjectionVar_IntX = ReturnKinematicVectorFromString(ProjectionVar_StrX);
  else ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  if (IsSubEventVarY) ProjectionVar_IntY = ReturnKinematicVectorFromString(ProjectionVar_StrY);
  else ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

  //JM Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    const int EventSample = MCEvents[iEvent].NominalSample;
    if(EventSample != iSample) continue;
    if (IsEventSelected(EventSample, iEvent)) {
      double Weight = GetEventWeight(iEvent);
      if (WeightStyle == 1) {
        Weight = 1.;
      }
      std::vector<double> VecX = {}, VecY = {};
      double VarX = M3::_BAD_DOUBLE_, VarY = M3::_BAD_DOUBLE_;
      size_t nsubevents = 0;
      // JM Three cases: subeventX vs eventY || eventX vs subeventY || subeventX vs subeventY
      if (IsSubEventVarX && !IsSubEventVarY) {
        VecX = ReturnKinematicVector(ProjectionVar_IntX, iEvent);
        VarY = ReturnKinematicParameter(ProjectionVar_IntY, iEvent);
        nsubevents = VecX.size();
      }
      else if (!IsSubEventVarX && IsSubEventVarY) {
        VecY = ReturnKinematicVector(ProjectionVar_IntY, iEvent);
        VarX = ReturnKinematicParameter(ProjectionVar_IntX, iEvent);
        nsubevents = VecY.size();
      }
      else {
        VecX = ReturnKinematicVector(ProjectionVar_IntX, iEvent);
        VecY = ReturnKinematicVector(ProjectionVar_IntY, iEvent);
        if (VecX.size() != VecY.size()) {
          MACH3LOG_ERROR("Cannot plot {} of size {} against {} of size {}", ProjectionVar_StrX, VecX.size(), ProjectionVar_StrY, VecY.size());
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        nsubevents = VecX.size();
      }
      //JM Loop over all subevents in event
      for (unsigned int iSubEvent = 0; iSubEvent < nsubevents; iSubEvent++) {
        if (IsSubEventSelected(SubEventSelectionVec, iEvent, iSubEvent, nsubevents)) {
          if (IsSubEventVarX) VarX = VecX[iSubEvent];
          if (IsSubEventVarY) VarY = VecY[iSubEvent];
          _h2DVar->Fill(VarX,VarY,Weight);
        }
      }
    }
  }
}

// ************************************************
int SampleHandlerBase::ReturnKinematicParameterFromString(const std::string& KinematicParameterStr) const {
// ************************************************
  auto it = KinematicParameters->find(KinematicParameterStr);
  if (it != KinematicParameters->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameterStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ************************************************
std::string SampleHandlerBase::ReturnStringFromKinematicParameter(const int KinematicParameter) const {
// ************************************************
  auto it = ReversedKinematicParameters->find(KinematicParameter);
  if (it != ReversedKinematicParameters->end()) {
    return it->second;
  }

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameter);
  throw MaCh3Exception(__FILE__, __LINE__);

  return "";
}

// ************************************************
// === JM define KinematicVector-to-string mapping functions  ===
int SampleHandlerBase::ReturnKinematicVectorFromString(const std::string& KinematicVectorStr) const {
// ************************************************
  auto it = KinematicVectors->find(KinematicVectorStr);
  if (it != KinematicVectors->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Vector: {}", KinematicVectorStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ************************************************
std::string SampleHandlerBase::ReturnStringFromKinematicVector(const int KinematicVector) const {
// ************************************************
  auto it = ReversedKinematicVectors->find(KinematicVector);
  if (it != ReversedKinematicVectors->end()) {
    return it->second;
  }

  MACH3LOG_ERROR("Did not recognise Kinematic Vector: {}", KinematicVector);
  throw MaCh3Exception(__FILE__, __LINE__);

  return "";
}

// ************************************************
std::vector<double> SampleHandlerBase::ReturnKinematicParameterBinning(const int Sample, const std::string& KinematicParameter) const {
// ************************************************
  // If any of fit based variables return them
  /// @todo might be useful to allow overwriting this
  if(Binning->IsUniform(Sample)) {
    for(int iDim = 0; iDim < GetNDim(Sample); iDim++) {
      if(KinematicParameter == GetKinVarName(Sample, iDim)) {
        return Binning->GetBinEdges(Sample, iDim);
      }
    } // end loop over Dimension
  } // end loop over IsUniform

  auto MakeBins = [](int nBins) {
    std::vector<double> bins(nBins + 1);
    for (int i = 0; i <= nBins; ++i)
      bins[i] = static_cast<double>(i) - 0.5;
    return bins;
  };

  if (KinematicParameter == "OscillationChannel") {
    return MakeBins(GetNOscChannels(Sample));
  } else if (KinematicParameter == "Mode") {
    return MakeBins(Modes->GetNModes());
  }

  std::vector<double> BinningVect;
  // We first check if binning for a sample has been specified
  auto BinningConfig = M3OpenConfig(SampleManager->raw()["BinningFile"].as<std::string>());
  bool found_range_specifier = false;
  if(BinningConfig[GetSampleTitle(Sample)] && BinningConfig[GetSampleTitle(Sample)][KinematicParameter]){
    BinningVect = BuildBinEdgesFromNode(BinningConfig[GetSampleTitle(Sample)][KinematicParameter], found_range_specifier);
  } else {
    BinningVect = BuildBinEdgesFromNode(BinningConfig[KinematicParameter], found_range_specifier);
  }

  // Ensure binning is increasing
  auto IsIncreasing = [](const std::vector<double>& vec) {
    for (size_t i = 1; i < vec.size(); ++i) {
      if (vec[i] <= vec[i-1]) {
        return false;
      }
    }
    return true;
  };

  if (!IsIncreasing(BinningVect)) {
    MACH3LOG_ERROR("Binning for {} is not increasing [{}]", KinematicParameter, fmt::join(BinningVect, ", "));
    if(found_range_specifier){
      MACH3LOG_ERROR("A bin range specifier was found. Please carefully check the number of square brackets used.");
    }
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return BinningVect;
}

// ************************************************
bool SampleHandlerBase::IsSubEventVarString(const std::string& VarStr) const {
// ************************************************
  if (KinematicVectors == nullptr) return false;

  if (KinematicVectors->count(VarStr)) {
    if (!KinematicParameters->count(VarStr)) return true;
    else {
      MACH3LOG_ERROR("Attempted to plot kinematic variable {}, but it appears in both KinematicVectors and KinematicParameters", VarStr);
      throw MaCh3Exception(__FILE__,__LINE__);
    }
  }
  return false;
}

// ************************************************
std::vector<KinematicCut> SampleHandlerBase::BuildModeChannelSelection(const int iSample, const int kModeToFill, const int kChannelToFill) const {
// ************************************************
  bool fChannel;
  bool fMode;

  if (kChannelToFill != -1) {
    if (kChannelToFill > GetNOscChannels(iSample)) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}", GetNOscChannels(iSample));
      MACH3LOG_ERROR("kChannelToFill given:{}", kChannelToFill);
      MACH3LOG_ERROR("Exiting.");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill != -1) {
    if (!(kModeToFill >= 0) && (kModeToFill < Modes->GetNModes())) {
      MACH3LOG_ERROR("Required mode is not available. kModeToFill should be between 0 and {}", Modes->GetNModes());
      MACH3LOG_ERROR("kModeToFill given:{}", kModeToFill);
      MACH3LOG_ERROR("Exiting..");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector< KinematicCut > SelectionVec;

  if (fMode) {
    KinematicCut SelecMode;
    SelecMode.ParamToCutOnIt = ReturnKinematicParameterFromString("Mode");
    SelecMode.LowerBound = kModeToFill;
    SelecMode.UpperBound = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    KinematicCut SelecChannel;
    SelecChannel.ParamToCutOnIt = ReturnKinematicParameterFromString("OscillationChannel");
    SelecChannel.LowerBound = kChannelToFill;
    SelecChannel.UpperBound = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return SelectionVec;
}

// ************************************************
std::unique_ptr<TH1> SampleHandlerBase::Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
                                                                     const int kModeToFill, const int kChannelToFill, const int WeightStyle) {
// ************************************************
  auto SelectionVec = BuildModeChannelSelection(iSample, kModeToFill, kChannelToFill);
  return Get1DVarHist(iSample, ProjectionVar_Str, SelectionVec, WeightStyle);
}

// ************************************************
std::unique_ptr<TH2> SampleHandlerBase::Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX,
                                                                     const std::string& ProjectionVar_StrY, const int kModeToFill,
                                                                     const int kChannelToFill, const int WeightStyle) {
// ************************************************
  auto SelectionVec = BuildModeChannelSelection(iSample, kModeToFill, kChannelToFill);
  return Get2DVarHist(iSample, ProjectionVar_StrX, ProjectionVar_StrY, SelectionVec, WeightStyle);
}

// ************************************************
void SampleHandlerBase::PrintIntegral(const int iSample, const TString& OutputFileName, const int WeightStyle, const TString& OutputCSVFileName) {
// ************************************************
  constexpr int space = 14;

  bool printToFile=false;
  if (OutputFileName.CompareTo("/dev/null")) {printToFile = true;}

  bool printToCSV=false;
  if(OutputCSVFileName.CompareTo("/dev/null")) printToCSV=true;

  std::ofstream outfile;
  if (printToFile) {
    outfile.open(OutputFileName.Data(), std::ios_base::app);
    outfile.precision(7);
  }

  std::ofstream outcsv;
  if(printToCSV){
    outcsv.open(OutputCSVFileName, std::ios_base::app); // Appened to CSV
    outcsv.precision(7);
  }

  double PDFIntegral = 0.;

  std::vector< std::vector< std::unique_ptr<TH1> > > IntegralList;
  IntegralList.resize(Modes->GetNModes());

  std::vector<double> ChannelIntegral;
  ChannelIntegral.resize(GetNOscChannels(iSample));
  for (unsigned int i=0;i<ChannelIntegral.size();i++) {ChannelIntegral[i] = 0.;}

  for (int i=0;i<Modes->GetNModes();i++) {
    if (GetNDim(iSample) == 1) {
      IntegralList[i] = ReturnHistsBySelection1D(iSample, GetKinVarName(iSample, 0), SamplePlotType::kOscChannelPlot, i, WeightStyle);
    } else {
      IntegralList[i] = CastVector<TH2, TH1>(ReturnHistsBySelection2D(iSample, GetKinVarName(iSample, 0), GetKinVarName(iSample, 1),
                                                                      SamplePlotType::kOscChannelPlot, i, WeightStyle));
    }
  }

  MACH3LOG_INFO("-------------------------------------------------");

  if (printToFile) {
    outfile << "\\begin{table}[ht]" << std::endl;
    outfile << "\\begin{center}" << std::endl;
    outfile << "\\caption{Integral breakdown for sample: " << GetSampleTitle(iSample) << "}" << std::endl;
    outfile << "\\label{" << GetSampleTitle(iSample) << "-EventRate}" << std::endl;

    TString nColumns;
    for (int i=0;i<GetNOscChannels(iSample);i++) {nColumns+="|c";}
    nColumns += "|c|";
    outfile << "\\begin{tabular}{|l" << nColumns.Data() << "}" << std::endl;
    outfile << "\\hline" << std::endl;
  }

  if(printToCSV){
    // HW Probably a better way but oh well, here I go making MaCh3 messy again
    outcsv<<"Integral Breakdown for sample :"<<GetSampleTitle(iSample)<<"\n";
  }

  MACH3LOG_INFO("Integral breakdown for sample: {}", GetSampleTitle(iSample));
  MACH3LOG_INFO("");

  if (printToFile) {outfile << std::setw(space) << "Mode:";}
  if(printToCSV) {outcsv<<"Mode,";}

  std::string table_headings = fmt::format("| {:<8} |", "Mode");
  std::string table_footline = "------------"; //Scalable table horizontal line
  for (int i = 0;i < GetNOscChannels(iSample); i++) {
    table_headings += fmt::format(" {:<17} |", GetFlavourName(iSample, i));
    table_footline += "--------------------";
    if (printToFile) {outfile << "&" << std::setw(space) << SampleDetails[iSample].OscChannels[i].flavourName_Latex << " ";}
    if (printToCSV)  {outcsv << GetFlavourName(iSample, i) << ",";}
  }
  if (printToFile) {outfile << "&" << std::setw(space) << "Total:" << "\\\\ \\hline" << std::endl;}
  if (printToCSV)  {outcsv <<"Total\n";}
  table_headings += fmt::format(" {:<10} |", "Total");
  table_footline += "-------------";

  MACH3LOG_INFO("{}", table_headings);
  MACH3LOG_INFO("{}", table_footline);

  for (unsigned int i=0;i<IntegralList.size();i++) {
    double ModeIntegral = 0;
    if (printToFile) {outfile << std::setw(space) << Modes->GetMaCh3ModeName(i);}
    if(printToCSV)   {outcsv << Modes->GetMaCh3ModeName(i) << ",";}

    table_headings = fmt::format("| {:<8} |", Modes->GetMaCh3ModeName(i)); //Start string with mode name

    for (unsigned int j=0;j<IntegralList[i].size();j++) {
      double Integral = IntegralList[i][j]->Integral();

      if (Integral < 1e-100) {Integral=0;}

      ModeIntegral += Integral;
      ChannelIntegral[j] += Integral;
      PDFIntegral += Integral;

      if (printToFile) {outfile << "&" << std::setw(space) << Form("%4.5f",Integral) << " ";}
      if (printToCSV)  {outcsv << Form("%4.5f", Integral) << ",";}

      table_headings += fmt::format(" {:<17.4f} |", Integral);
    }
    if (printToFile) {outfile << "&" << std::setw(space) << Form("%4.5f",ModeIntegral) <<  " \\\\ \\hline" << std::endl;}
    if (printToCSV)  {outcsv << Form("%4.5f", ModeIntegral) << "\n";}

    table_headings += fmt::format(" {:<10.4f} |", ModeIntegral);

    MACH3LOG_INFO("{}", table_headings);
  }

  if (printToFile) {outfile << std::setw(space) << "Total:";}
  if (printToCSV)  {outcsv << "Total,";}

  //Clear the table_headings to print last row of totals
  table_headings = fmt::format("| {:<8} |", "Total");
  for (unsigned int i=0;i<ChannelIntegral.size();i++) {
    if (printToFile) {outfile << "&" << std::setw(space) << Form("%4.5f",ChannelIntegral[i]) << " ";}
    if (printToCSV)  {outcsv << Form("%4.5f", ChannelIntegral[i]) << ",";}
    table_headings += fmt::format(" {:<17.4f} |", ChannelIntegral[i]);
  }
  if (printToFile) {outfile << "&" << std::setw(space) << Form("%4.5f",PDFIntegral) << " \\\\ \\hline" << std::endl;}
  if (printToCSV)  {outcsv << Form("%4.5f", PDFIntegral) << "\n\n\n\n";} // Let's have a few new lines!

  table_headings += fmt::format(" {:<10.4f} |", PDFIntegral);
  MACH3LOG_INFO("{}", table_headings);
  MACH3LOG_INFO("{}", table_footline);

  if (printToFile) {
    outfile << "\\end{tabular}" << std::endl;
    outfile << "\\end{center}" << std::endl;
    outfile << "\\end{table}" << std::endl;
  }

  MACH3LOG_INFO("");

  if (printToFile) {
    outfile << std::endl;
    outfile.close();
  }
}

// ************************************************
int SampleHandlerBase::GetRangeForPlotType(const SamplePlotType TypeEnum, const int iSample) const {
// ************************************************
  switch (TypeEnum) {
    case SamplePlotType::kModePlot:
      return Modes->GetNModes();
      break;
    case SamplePlotType::kOscChannelPlot:
      return GetNOscChannels(iSample);
      break;
    default:
      MACH3LOG_ERROR("You've passed me a SamplePlotType with value {} which was not implemented.", static_cast<int>(TypeEnum));
      throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
std::vector<std::unique_ptr<TH1>> SampleHandlerBase::ReturnHistsBySelection1D(const int iSample, const std::string& KinematicProjection,
                                                                              const SamplePlotType Selection1, const int Selection2, const int WeightStyle) {
// ************************************************
  std::vector<std::unique_ptr<TH1>> hHistList;
  std::string legendEntry;

  if (THStackLeg != nullptr) {delete THStackLeg;}
  THStackLeg = new TLegend(0.1,0.1,0.9,0.9);

  const int iMax = GetRangeForPlotType(Selection1, iSample);
  for (int i=0;i<iMax;i++) {
    if (Selection1 == SamplePlotType::kModePlot) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(iSample, KinematicProjection,i,Selection2,WeightStyle));
      THStackLeg->AddEntry(hHistList[i].get(), (Modes->GetMaCh3ModeName(i)+Form(" : (%4.2f)",hHistList[i]->Integral())).c_str(),"f");

      hHistList[i]->SetFillColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
      hHistList[i]->SetLineColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
    }
    if (Selection1 == SamplePlotType::kOscChannelPlot) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(iSample, KinematicProjection,Selection2,i,WeightStyle));
      THStackLeg->AddEntry(hHistList[i].get(),(GetFlavourName(iSample, i)+Form(" | %4.2f",hHistList[i]->Integral())).c_str(),"f");
    }
  }

  return hHistList;
}

// ************************************************
std::vector<std::unique_ptr<TH2>> SampleHandlerBase::ReturnHistsBySelection2D(const int iSample, const std::string& KinematicProjectionX,
                                                                              const std::string& KinematicProjectionY, const SamplePlotType Selection1,
                                                                              const int Selection2, const int WeightStyle) {
// ************************************************
  std::vector<std::unique_ptr<TH2>> hHistList;

  const int iMax = GetRangeForPlotType(Selection1, iSample);
  for (int i=0;i<iMax;i++) {
    if (Selection1 == SamplePlotType::kModePlot) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(iSample, KinematicProjectionX,KinematicProjectionY,i,Selection2,WeightStyle));
    }
    if (Selection1 == SamplePlotType::kOscChannelPlot) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(iSample, KinematicProjectionX,KinematicProjectionY,Selection2,i,WeightStyle));
    }
  }

  return hHistList;
}

// ************************************************
std::unique_ptr<THStack> SampleHandlerBase::ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection,
                                                                           const SamplePlotType Selection1, int Selection2, int WeightStyle) {
// ************************************************
  auto HistList = ReturnHistsBySelection1D(iSample, KinematicProjection, Selection1, Selection2, WeightStyle);
  auto StackHist = std::make_unique<THStack>((GetSampleTitle(iSample)+"_"+KinematicProjection+"_Stack").c_str(),"");
  // Note: we use .release() to transfer ownership of each TH1 to THStack.
  for (unsigned int i=0;i<HistList.size();i++) {
    StackHist->Add(HistList[i].release());
  }
  return StackHist;
}

// ************************************************
const double* SampleHandlerBase::GetPointerToOscChannel(const int iEvent) const {
// ************************************************
  auto& OscillationChannels = SampleDetails[MCEvents[iEvent].NominalSample].OscChannels;
  const int Channel = GetOscChannel(OscillationChannels, MCEvents[iEvent].nupdgUnosc, MCEvents[iEvent].nupdg);
  return &(OscillationChannels[Channel].ChannelIndex);
}

// ***************************************************************************
void SampleHandlerBase::CheckEmptyBins() const {
// ***************************************************************************
  const auto TotalBins = Binning->GetNBins();
  int iCounter = 0;
  for(int iBin = 0; iBin < TotalBins; iBin++) {
    if(SampleHandler_array[iBin] == 0) {
      iCounter++;
      MACH3LOG_DEBUG("Bin {}, for sample {}, has 0 entries",
                     Binning->GetBinName(iBin), GetSampleTitle(Binning->GetSampleIndex(iBin)));
    }
  }

  if(iCounter > 0){
    MACH3LOG_WARN("Found in total {} ({:.2f}%) empty bins for SampleHandler: {}",
                  iCounter, 100.0 * static_cast<double>(iCounter) / TotalBins, GetName());
  }
}

// ***************************************************************************
// Helper function to print rates for the samples with LLH
void SampleHandlerBase::PrintRates(const bool DataOnly) {
// ***************************************************************************
  if (!SampleHandler_data.size()) {
    MACH3LOG_ERROR("Data sample is empty!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Printing for {}", GetName());

  if (!DataOnly) {
    const std::string sep_full(81, '-');
    MACH3LOG_INFO("{}", sep_full);
    MACH3LOG_INFO("{:<40}{:<15}{:<15}{:<10}|", "Sample", "Data", "MC", "-LLH");
  } else {
    const std::string sep_data(56, '-');
    MACH3LOG_INFO("{}", sep_data);
    MACH3LOG_INFO("{:<40}{:<15}|", "Sample", "Data");
  }

  double sumData = 0.0;
  double sumMC = 0.0;
  double likelihood = 0.0;

  for (int iSample = 0; iSample < GetNSamples(); ++iSample) {
    std::string name = GetSampleTitle(iSample);
    std::vector<double> DataArray = GetDataArray(iSample);
    double dataIntegral = std::accumulate(DataArray.begin(), DataArray.end(), 0.0);
    sumData += dataIntegral;
    if (!DataOnly) {
      std::vector<double> MCArray = GetMCArray(iSample);
      double mcIntegral = std::accumulate(MCArray.begin(), MCArray.end(), 0.0);
      sumMC += mcIntegral;
      likelihood = GetSampleLikelihood(iSample);

      MACH3LOG_INFO("{:<40}{:<15.2f}{:<15.2f}{:<10.2f}|", name, dataIntegral, mcIntegral, likelihood);
    } else {
      MACH3LOG_INFO("{:<40}{:<15.2f}|", name, dataIntegral);
    }
  }
  if (!DataOnly) {
    likelihood = GetLikelihood();
    MACH3LOG_INFO("{:<40}{:<15.2f}{:<15.2f}{:<10.2f}|", "Total", sumData, sumMC, likelihood);
    const std::string sep_full(81, '-');
    MACH3LOG_INFO("{}", sep_full);
  } else {
    MACH3LOG_INFO("{:<40}{:<20.2f}|", "Total", sumData);
    const std::string sep_data(56, '-');
    MACH3LOG_INFO("{}", sep_data);
  }
}

// ***************************************************************************
// Return Kinematic Variable name for specified sample and dimension for example "Reconstructed_Neutrino_Energy"
std::string SampleHandlerBase::GetKinVarName(const int iSample, const int Dimension) const {
// ***************************************************************************
  if(Dimension > GetNDim(iSample)) {
    MACH3LOG_ERROR("Asking for dimension {}, while sample: {} only has {}", Dimension, GetSampleTitle(iSample), GetNDim(iSample));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleDetails[iSample].VarStr[Dimension];
}

// ***************************************************************************
std::vector<double> SampleHandlerBase::GetArrayForSample(const int Sample, std::vector<double> const & array) const {
// ***************************************************************************
  const int Start = Binning->GetSampleStartBin(Sample);
  const int End   = Binning->GetSampleEndBin(Sample);

  return std::vector<double>(array.begin() + Start, array.begin() + End);
}
