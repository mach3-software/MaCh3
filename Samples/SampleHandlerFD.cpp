#include "SampleHandlerFD.h"
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include <cstddef>

#include <algorithm>
#include <memory>

// ************************************************
SampleHandlerFD::SampleHandlerFD(std::string ConfigFileName, ParameterHandlerGeneric* xsec_cov,
                                 const std::shared_ptr<OscillationHandler>& OscillatorObj_) : SampleHandlerBase() {
// ************************************************
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Creating SampleHandlerFD object");

  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(!xsec_cov){
    MACH3LOG_ERROR("You've passed me a nullptr to a SystematicHandlerGeneric... I need this to setup splines!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  ParHandler = xsec_cov;

  nSamples = 1;

  if (OscillatorObj_ != nullptr) {
    MACH3LOG_WARN("You have passed an Oscillator object through the constructor of a SampleHandlerFD object - this will be used for all oscillation channels");
    Oscillator = OscillatorObj_;
  }

  KinematicParameters = nullptr;
  ReversedKinematicParameters = nullptr;
  KinematicVectors = nullptr;
  ReversedKinematicVectors = nullptr;

  SampleHandlerFD_array = nullptr;
  SampleHandlerFD_data = nullptr;
  SampleHandlerFD_array_w2 = nullptr;
  SampleHandlerName = "";
  SampleManager = std::make_unique<manager>(ConfigFileName.c_str());
  Binning = std::make_unique<BinningHandler>();
  // Variables related to MC stat
  FirstTimeW2 = true;
  UpdateW2 = false;
}

SampleHandlerFD::~SampleHandlerFD() {
  MACH3LOG_DEBUG("I'm deleting SampleHandlerFD");
  
  if (SampleHandlerFD_array != nullptr) delete[] SampleHandlerFD_array;
  if (SampleHandlerFD_array_w2 != nullptr) delete[] SampleHandlerFD_array_w2;
  //ETA - there is a chance that you haven't added any data...
  if (SampleHandlerFD_data != nullptr) delete[] SampleHandlerFD_data;

  if(THStackLeg != nullptr) delete THStackLeg;
}

void SampleHandlerFD::ReadConfig()
{
  auto ModeName = Get<std::string>(SampleManager->raw()["MaCh3ModeConfig"], __FILE__ , __LINE__);
  Modes = std::make_unique<MaCh3Modes>(ModeName);
  //SampleName has to be provided in the sample yaml otherwise this will throw an exception
  SampleHandlerName = Get<std::string>(SampleManager->raw()["SampleHandlerName"], __FILE__ , __LINE__);

  fTestStatistic = static_cast<TestStatistic>(SampleManager->GetMCStatLLH());
  if (CheckNodeExists(SampleManager->raw(), "LikelihoodOptions")) {
    UpdateW2 = GetFromManager<bool>(SampleManager->raw()["LikelihoodOptions"]["UpdateW2"], false);
  }

  if (!CheckNodeExists(SampleManager->raw(), "BinningFile")){
    MACH3LOG_ERROR("BinningFile not given in for sample handler {}, ReturnKinematicParameterBinning will not work", SampleHandlerName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto EnabledSasmples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  // Get number of samples and resize relevant objects
  nSamples = static_cast<M3::int_t>(EnabledSasmples.size());
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
void SampleHandlerFD::LoadSingleSample(const int iSample, const YAML::Node& SampleSettings) {
// ************************************************
  SampleInfo SingleSample;
  //SampleTitle has to be provided in the sample yaml otherwise this will throw an exception
  SingleSample.SampleTitle = Get<std::string>(SampleSettings["SampleTitle"], __FILE__ , __LINE__);


  Binning->SetupSampleBinning(SampleSettings["Binning"], SingleSample);

  auto mtupleprefix  = Get<std::string>(SampleSettings["InputFiles"]["mtupleprefix"], __FILE__, __LINE__);
  auto mtuplesuffix  = Get<std::string>(SampleSettings["InputFiles"]["mtuplesuffix"], __FILE__, __LINE__);
  auto splineprefix  = Get<std::string>(SampleSettings["InputFiles"]["splineprefix"], __FILE__, __LINE__);
  auto splinesuffix  = Get<std::string>(SampleSettings["InputFiles"]["splinesuffix"], __FILE__, __LINE__);

  int NChannels = static_cast<M3::int_t>(SampleSettings["SubSamples"].size());
  SingleSample.OscChannels.reserve(NChannels);

  int OscChannelCounter = 0;
  for (auto const &osc_channel : SampleSettings["SubSamples"]) {
    std::string MTupleFileName = mtupleprefix+osc_channel["mtuplefile"].as<std::string>()+mtuplesuffix;

    OscChannelInfo OscInfo;
    OscInfo.flavourName       = osc_channel["Name"].as<std::string>();
    OscInfo.flavourName_Latex = osc_channel["LatexName"].as<std::string>();
    OscInfo.InitPDG           = static_cast<NuPDG>(osc_channel["nutype"].as<int>());
    OscInfo.FinalPDG          = static_cast<NuPDG>(osc_channel["oscnutype"].as<int>());
    OscInfo.ChannelIndex      = OscChannelCounter;

    SingleSample.OscChannels.push_back(std::move(OscInfo));

    FileToInitPDGMap[MTupleFileName] = static_cast<NuPDG>(osc_channel["nutype"].as<int>());
    FileToFinalPDGMap[MTupleFileName] = static_cast<NuPDG>(osc_channel["oscnutype"].as<int>());

    SingleSample.mc_files.push_back(MTupleFileName);
    SingleSample.spline_files.push_back(splineprefix+osc_channel["splinefile"].as<std::string>()+splinesuffix);
    OscChannelCounter++;
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
  /// Add new sample
  // KS: Important to first call move and then InitialiseHistograms otheriwse histograms will be deleted
  // need safer solution but do this like this for now
  SampleDetails[iSample] = std::move(SingleSample);
  SampleDetails[iSample].InitialiseHistograms();
}

void SampleHandlerFD::Initialise() {
  TStopwatch clock;
  clock.Start();

  //First grab all the information from your sample config via your manager
  ReadConfig();

  //Now initialise all the variables you will need
  Init();

  nEvents = SetupExperimentMC();
  MCSamples.resize(nEvents);
  SetupFDMC();

  MACH3LOG_INFO("=============================================");
  MACH3LOG_INFO("Total number of events is: {}", GetNEvents());

  auto OscParams = ParHandler->GetOscParsFromSampleName(SampleHandlerName);
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
  MACH3LOG_INFO("Setting up Sample Binning..");
  SetBinning();
  MACH3LOG_INFO("Setting up Splines..");
  SetupSplines();
  MACH3LOG_INFO("Setting up Normalisation Pointers..");
  SetupNormParameters();
  MACH3LOG_INFO("Setting up Functional Pointers..");
  SetupFunctionalParameters();
  MACH3LOG_INFO("Setting up Additional Weight Pointers..");
  AddAdditionalWeightPointers();
  MACH3LOG_INFO("Setting up Kinematic Map..");
  SetupKinematicMap();
  clock.Stop();
  MACH3LOG_INFO("Finished loading MC for {}, it took {:.2f}s to finish", GetName(), clock.RealTime());
  MACH3LOG_INFO("=======================================================");
}

// ************************************************
void SampleHandlerFD::SetupKinematicMap() {
// ************************************************
  if(KinematicParameters == nullptr || ReversedKinematicParameters == nullptr) {
    MACH3LOG_INFO("Map KinematicParameters or ReversedKinematicParameters hasn't been initialised");
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
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// ************************************************
void SampleHandlerFD::FillMCHist(const int Sample, const int Dimension) {
// ************************************************
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  if(Dimension == 1){
    SampleDetails[Sample]._hPDF1D->Reset();
    for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
      for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
        const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
        SampleDetails[Sample]._hPDF1D->SetBinContent(xBin + 1, SampleHandlerFD_array[idx]);
      }
    }
  } else if (Dimension == 2) {
    SampleDetails[Sample]._hPDF2D->Reset();
    for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
      for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
        const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
        SampleDetails[Sample]._hPDF2D->SetBinContent(xBin + 1, yBin + 1, SampleHandlerFD_array[idx]);
      }
    }
  } else {
    MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, Dimension);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}
#pragma GCC diagnostic pop


// ************************************************
bool SampleHandlerFD::IsEventSelected(const int iSample, const int iEvent) {
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
bool SampleHandlerFD::IsSubEventSelected(const std::vector<KinematicCut> &SubEventCuts, const int iEvent, const unsigned int iSubEvent, size_t nsubevents) {
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
void SampleHandlerFD::Reweight() {
//************************************************
  //KS: Reset the histograms before reweight 
  ResetHistograms();
  
  //You only need to do these things if Oscillator has been initialised
  //if not then you're not considering oscillations
  if (Oscillator) Oscillator->Evaluate();

  // Calculate weight coming from all splines if we initialised handler
  if(SplineHandler) SplineHandler->Evaluate();

  #ifdef MULTITHREAD
  // Call entirely different routine if we're running with openMP
  FillArray_MP();
  #else
  FillArray();
  #endif

  //KS: If you want to not update W2 wights then uncomment this line
  if(!UpdateW2) FirstTimeW2 = false;
}

//************************************************
// Function which does the core reweighting. This assumes that oscillation weights have
// already been calculated and stored in SampleHandlerFD.osc_w[iEvent]. This
// function takes advantage of most of the things called in setupSKMC to reduce reweighting time.
// It also follows the ND code reweighting pretty closely. This function fills the SampleHandlerFD
// array array which is binned to match the sample binning, such that bin[1][1] is the
// equivalent of SampleDetails._hPDF2D->GetBinContent(2,2) {Noticing the offset}
void SampleHandlerFD::FillArray() {
//************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;
  
  PrepFunctionalParameters();

  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    ApplyShifts(iEvent);
    FarDetectorCoreInfo* MCEvent = &MCSamples[iEvent];

    if (!IsEventSelected(MCEvent->NominalSample, iEvent)) {
      continue;
    }

    // Virtual by default does nothing, has to happen before CalcWeightTotal
    CalcWeightFunc(iEvent);

    const M3::float_t totalweight = CalcWeightTotal(MCEvent);
    //DB Catch negative total weights and skip any event with a negative weight. Previously we would set weight to zero and continue but that is inefficient
    if (totalweight <= 0.){
      continue;
    }

    //DB Find the relevant bin in the PDF for each event
    const int GlobalBin = Binning->FindGlobalBin(MCEvent->NominalSample, MCEvent->KinVar, MCEvent->NomBin);

    //DB Fill relevant part of thread array
    if (GlobalBin > M3::UnderOverFlowBin) {
      SampleHandlerFD_array[GlobalBin] += totalweight;
      if (FirstTimeW2) SampleHandlerFD_array_w2[GlobalBin] += totalweight*totalweight;
    }
  }
}

#ifdef MULTITHREAD
// ************************************************ 
/// Multithreaded version of fillArray @see fillArray()
void SampleHandlerFD::FillArray_MP() {
// ************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;

  PrepFunctionalParameters();

  //This is stored as [y][x] due to shifts only occurring in the x variable (Erec/Lep mom) - I believe this will help reduce cache misses
  double* SampleHandlerFD_array_private = nullptr;
  double* SampleHandlerFD_array_private_w2 = nullptr;
  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
  #pragma omp parallel private(SampleHandlerFD_array_private, SampleHandlerFD_array_private_w2)
  {
    // private to each thread
    // ETA - maybe we can use parallel firstprivate to initialise these?
    SampleHandlerFD_array_private = new double[Binning->GetNBins()];
    SampleHandlerFD_array_private_w2 = new double[Binning->GetNBins()];

    std::fill_n(SampleHandlerFD_array_private, Binning->GetNBins(), 0.0);
    std::fill_n(SampleHandlerFD_array_private_w2, Binning->GetNBins(), 0.0);

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
    // 1. Order minituples in Y-axis variable as this will *hopefully* reduce cache misses inside SampleHandlerFD_array_class[yBin][xBin]
    //
    // We will hit <0.1 s/step eventually! :D

    const unsigned int NumberOfEvents = GetNEvents();
    #pragma omp for
    for (unsigned int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
      //ETA - generic functions to apply shifts to kinematic variables
      // Apply this before IsEventSelected is called.
      ApplyShifts(iEvent);

      FarDetectorCoreInfo* MCEvent = &MCSamples[iEvent];
      //ETA - generic functions to apply shifts to kinematic variable
      //this is going to be slow right now due to string comps under the hood.
      //Need to implement a more efficient version of event-by-event cut checks
      if(!IsEventSelected(MCEvent->NominalSample, iEvent)){
        continue;
      }

      // Virtual by default does nothing, has to happen before CalcWeightTotal
      CalcWeightFunc(iEvent);

      const M3::float_t totalweight = CalcWeightTotal(MCEvent);
      //DB Catch negative total weights and skip any event with a negative weight. Previously we would set weight to zero and continue but that is inefficient
      if (totalweight <= 0.){
        continue;
      }

      //DB Find the relevant bin in the PDF for each event
      const int GlobalBin = Binning->FindGlobalBin(MCEvent->NominalSample, MCEvent->KinVar, MCEvent->NomBin);

      //ETA - we can probably remove this final if check on the -1?
      //Maybe we can add an overflow bin to the array and assign any events to this bin?
      //Might save us an extra if call?
      //DB Fill relevant part of thread array
      if (GlobalBin > M3::UnderOverFlowBin) {
        SampleHandlerFD_array_private[GlobalBin] += totalweight;
        SampleHandlerFD_array_private_w2[GlobalBin] += totalweight*totalweight;
      }
    }
    //End of Calc Weights and fill Array
    //==================================================
    // DB Copy contents of 'SampleHandlerFD_array_private' into 'SampleHandlerFD_array' which can then be used in GetLikelihood
    for (int idx = 0; idx < Binning->GetNBins(); ++idx) {
      #pragma omp atomic
      SampleHandlerFD_array[idx] += SampleHandlerFD_array_private[idx];
      if (FirstTimeW2) {
        #pragma omp atomic
        SampleHandlerFD_array_w2[idx] += SampleHandlerFD_array_private_w2[idx];
      }
    }

    delete[] SampleHandlerFD_array_private;
    delete[] SampleHandlerFD_array_private_w2;
  } //end of parallel region
}
#endif

// **************************************************
// Helper function to reset the data and MC histograms
void SampleHandlerFD::ResetHistograms() {
// **************************************************  
  //DB Reset values stored in PDF array to 0.
  // Don't openMP this; no significant gain
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int i = 0; i < Binning->GetNBins(); ++i) {
    SampleHandlerFD_array[i] = 0.;
    if (FirstTimeW2) SampleHandlerFD_array_w2[i] = 0.;
  }
} // end function

void SampleHandlerFD::RegisterIndividualFunctionalParameter(const std::string& fpName, int fpEnum, FuncParFuncType fpFunc){
  // Add protections to not add the same functional parameter twice
  if (funcParsNamesMap.find(fpName) != funcParsNamesMap.end()) {
    MACH3LOG_ERROR("Functional parameter {} already registered in funcParsNamesMap with enum {}", fpName, funcParsNamesMap[fpName]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (std::find(funcParsNamesVec.begin(), funcParsNamesVec.end(), fpName) != funcParsNamesVec.end()) {
    MACH3LOG_ERROR("Functional parameter {} already in funcParsNamesVec", fpName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (funcParsFuncMap.find(fpEnum) != funcParsFuncMap.end()) {
    MACH3LOG_ERROR("Functional parameter enum {} already registered in funcParsFuncMap", fpEnum);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  funcParsNamesMap[fpName] = fpEnum;
  funcParsNamesVec.push_back(fpName);
  funcParsFuncMap[fpEnum] = fpFunc;
}

void SampleHandlerFD::SetupFunctionalParameters() {
  funcParsVec = ParHandler->GetFunctionalParametersFromSampleName(SampleHandlerName);
  // RegisterFunctionalParameters is implemented in experiment-specific code, 
  // which calls RegisterIndividualFuncPar to populate funcParsNamesMap, funcParsNamesVec, and funcParsFuncMap
  RegisterFunctionalParameters();
  funcParsMap.resize(funcParsNamesMap.size());
  funcParsGrid.resize(GetNEvents());

  // For every functional parameter in XsecCov that matches the name in funcParsNames, add it to the map
  for (FunctionalParameter & fp : funcParsVec) {
    for (std::string name : funcParsNamesVec) {
      if (fp.name == name) {
        MACH3LOG_INFO("Adding functional parameter: {} to funcParsMap with key: {}", fp.name, funcParsNamesMap[fp.name]);
        fp.funcPtr = &funcParsFuncMap[funcParsNamesMap[fp.name]];
        funcParsMap[static_cast<std::size_t>(funcParsNamesMap[fp.name])] = &fp;
        continue;
      }
    }
    // If we don't find a match, we need to throw an error
    if (funcParsMap[static_cast<std::size_t>(funcParsNamesMap[fp.name])] == nullptr) {
      MACH3LOG_ERROR("Functional parameter {} not found, did you define it in RegisterFunctionalParameters()?", fp.name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  // Mostly the same as CalcXsecNormsBins
  // For each event, make a vector of pointers to the functional parameters
  for (std::size_t iEvent = 0; iEvent < static_cast<std::size_t>(GetNEvents()); ++iEvent) {
    // Now loop over the functional parameters and get a vector of enums corresponding to the functional parameters
    for (std::vector<FunctionalParameter>::iterator it = funcParsVec.begin(); it != funcParsVec.end(); ++it) {
      // Check whether the interaction modes match
      bool ModeMatch = MatchCondition((*it).modes, static_cast<int>(std::round(*(MCSamples[iEvent].mode))));
      if (!ModeMatch) {
        MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, *(MCSamples[iEvent].mode), (*it).name);
        continue;
      }
      // Now check whether within kinematic bounds
      bool IsSelected = true;
      if ((*it).hasKinBounds) {
        const auto& kinVars = (*it).KinematicVarStr;
        const auto& selection = (*it).Selection;

        for (std::size_t iKinPar = 0; iKinPar < kinVars.size(); ++iKinPar) {
          const double kinVal = ReturnKinematicParameter(kinVars[iKinPar], static_cast<int>(iEvent));

          bool passedAnyBound = false;
          const auto& boundsList = selection[iKinPar];

          for (const auto& bounds : boundsList) {
            if (kinVal > bounds[0] && kinVal <= bounds[1]) {
              passedAnyBound = true;
              break;
            }
          }

          if (!passedAnyBound) {
            MACH3LOG_TRACE("Event {}, missed kinematic check ({}) for dial {}",
                           iEvent, kinVars[iKinPar], (*it).name);
            IsSelected = false;
            break;
          }
        }
      }
      // Need to then break the event loop
      if(!IsSelected){
        MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", iEvent, (*it).name);
        continue;
      }
      auto funcparenum = funcParsNamesMap[(*it).name];
      funcParsGrid.at(iEvent).push_back(funcparenum);
    }
  }
  MACH3LOG_INFO("Finished setting up functional parameters");
}

void SampleHandlerFD::ApplyShifts(int iEvent) {
  // Given a sample and event, apply the shifts to the event based on the vector of functional parameter enums
  // First reset shifted array back to nominal values
  resetShifts(iEvent);

  const std::vector<int>& fpEnums = funcParsGrid[iEvent];
  const std::size_t nShifts = fpEnums.size();

  for (std::size_t i = 0; i < nShifts; ++i) {
    const int fpEnum = fpEnums[i];
    FunctionalParameter *fp = funcParsMap[static_cast<std::size_t>(fpEnum)];
    // if (fp->funcPtr) {
    //   (*fp->funcPtr)(fp->valuePtr, iEvent);
    // } else {
    //   MACH3LOG_ERROR("Functional parameter function pointer for {} is null for event {} in sample {}", fp->name, iEvent);
    //   throw MaCh3Exception(__FILE__, __LINE__);
    // }
    (*fp->funcPtr)(fp->valuePtr, iEvent);
  }
}

// ***************************************************************************
// Calculate the spline weight for one event
M3::float_t SampleHandlerFD::CalcWeightTotal(const FarDetectorCoreInfo* _restrict_ MCEvent) const {
// ***************************************************************************
  M3::float_t TotalWeight = 1.0;
  const int nNorms = static_cast<int>(MCEvent->xsec_norm_pointers.size());
  //Loop over stored normalisation and function pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iParam = 0; iParam < nNorms; ++iParam)
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wuseless-cast"
    TotalWeight *= static_cast<M3::float_t>(*(MCEvent->xsec_norm_pointers[iParam]));
    #pragma GCC diagnostic pop
  }

  const int TotalWeights = static_cast<int>(MCEvent->total_weight_pointers.size());
  //DB Xsec syst
  //Loop over stored spline pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iWeight = 0; iWeight < TotalWeights; ++iWeight) {
    TotalWeight *= *(MCEvent->total_weight_pointers[iWeight]);
  }

  return TotalWeight;
}

// ***************************************************************************
// Setup the norm parameters
void SampleHandlerFD::SetupNormParameters() {
// ***************************************************************************
  std::vector< std::vector< int > > xsec_norms_bins(GetNEvents());

  std::vector<NormParameter> norm_parameters = ParHandler->GetNormParsFromSampleName(GetName());

  if(!ParHandler){
    MACH3LOG_ERROR("ParHandler is not setup!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Assign xsec norm bins in MCSamples tree
  CalcNormsBins(norm_parameters, xsec_norms_bins);

  //DB Attempt at reducing impact of SystematicHandlerGeneric::calcReweight()
  int counter;
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    counter = 0;

    MCSamples[iEvent].xsec_norm_pointers.resize(xsec_norms_bins[iEvent].size());
    for(auto const & norm_bin: xsec_norms_bins[iEvent]) {
      MCSamples[iEvent].xsec_norm_pointers[counter] = ParHandler->RetPointer(norm_bin);
      counter += 1;
    }
  }
}

// ************************************************
//A way to check whether a normalisation parameter applies to an event or not
void SampleHandlerFD::CalcNormsBins(std::vector<NormParameter>& norm_parameters, std::vector< std::vector< int > >& xsec_norms_bins) {
// ************************************************
  #ifdef DEBUG
  std::vector<int> VerboseCounter(norm_parameters.size(), 0);
  #endif
  for(unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent){
    std::vector< int > NormBins = {};
    if (ParHandler) {
      // Skip oscillated NC events
      // Not strictly needed, but these events don't get included in oscillated predictions, so
      // no need to waste our time calculating and storing information about xsec parameters
      // that will never be used.
      if (MCSamples[iEvent].isNC && (*MCSamples[iEvent].nupdg != *MCSamples[iEvent].nupdgUnosc) ) {
        MACH3LOG_TRACE("Event {}, missed NC/signal check", iEvent);
        continue;
      } //DB Abstract check on MaCh3Modes to determine which apply to neutral current
      for (std::vector<NormParameter>::iterator it = norm_parameters.begin(); it != norm_parameters.end(); ++it) {
        //Now check that the target of an interaction matches with the normalisation parameters
        bool TargetMatch = MatchCondition((*it).targets, *(MCSamples[iEvent].Target));
        if (!TargetMatch) {
          MACH3LOG_TRACE("Event {}, missed target check ({}) for dial {}", iEvent, *(MCSamples[iEvent].Target), (*it).name);
          continue;
        }

        //Now check that the neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourMatch = MatchCondition((*it).pdgs, *(MCSamples[iEvent].nupdg));
        if (!FlavourMatch) {
          MACH3LOG_TRACE("Event {}, missed PDG check ({}) for dial {}", iEvent,(*MCSamples[iEvent].nupdg), (*it).name);
          continue;
        }

        //Now check that the unoscillated neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourUnoscMatch = MatchCondition((*it).preoscpdgs, *(MCSamples[iEvent].nupdgUnosc));
        if (!FlavourUnoscMatch){
          MACH3LOG_TRACE("Event {}, missed FlavourUnosc check ({}) for dial {}", iEvent,(*MCSamples[iEvent].nupdgUnosc), (*it).name);
          continue;
        }

        //Now check that the mode of an interaction matches with the normalisation parameters
        bool ModeMatch = MatchCondition((*it).modes, static_cast<int>(std::round(*(MCSamples[iEvent].mode))));
        if (!ModeMatch) {
          MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, *(MCSamples[iEvent].mode), (*it).name);
          continue;
        }

        //Now check whether the norm has kinematic bounds
        //i.e. does it only apply to events in a particular kinematic region?
        // Now check whether within kinematic bounds
        bool IsSelected = true;
        if ((*it).hasKinBounds) {
          const auto& kinVars = (*it).KinematicVarStr;
          const auto& selection = (*it).Selection;

          for (std::size_t iKinPar = 0; iKinPar < kinVars.size(); ++iKinPar) {
            const double kinVal = ReturnKinematicParameter(kinVars[iKinPar], static_cast<int>(iEvent));

            bool passedAnyBound = false;
            const auto& boundsList = selection[iKinPar];

            for (const auto& bounds : boundsList) {
              if (kinVal > bounds[0] && kinVal <= bounds[1]) {
                passedAnyBound = true;
                break;
              }
            }

            if (!passedAnyBound) {
              MACH3LOG_TRACE("Event {}, missed kinematic check ({}) for dial {}",
                             iEvent, kinVars[iKinPar], (*it).name);
              IsSelected = false;
              break;
            }
          }
        }
        // Need to then break the event loop
        if(!IsSelected){
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", iEvent, (*it).name);
          continue;
        }
        // Now set 'index bin' for each normalisation parameter
        // All normalisations are just 1 bin for 2015, so bin = index (where index is just the bin for that normalisation)
        int bin = (*it).index;

        NormBins.push_back(bin);
        MACH3LOG_TRACE("Event {}, will be affected by dial {}", iEvent, (*it).name);
        #ifdef DEBUG
        VerboseCounter[std::distance(norm_parameters.begin(), it)]++;
        #endif
        //}
      } // end iteration over norm_parameters
    } // end if (ParHandler)
    xsec_norms_bins[iEvent] = NormBins;
  }//end loop over events
  #ifdef DEBUG
  MACH3LOG_DEBUG("┌──────────────────────────────────────────────────────────┐");
  for (std::size_t i = 0; i < norm_parameters.size(); ++i) {
    const auto& norm = norm_parameters[i];
    double eventRatio = static_cast<double>(VerboseCounter[i]) / static_cast<double>(GetNEvents());

    MACH3LOG_DEBUG("│ Param {:<15}, affects {:<8} events ({:>6.2f}%) │",
                  ParHandler->GetParFancyName(norm.index), VerboseCounter[i], eventRatio);
  }
  MACH3LOG_DEBUG("└──────────────────────────────────────────────────────────┘");
  #endif
}

// ************************************************
void SampleHandlerFD::SetupReweightArrays() {
// ************************************************
  SampleHandlerFD_array = new double[Binning->GetNBins()];
  SampleHandlerFD_array_w2 = new double[Binning->GetNBins()];
  SampleHandlerFD_data = new double[Binning->GetNBins()];

  for (int i = 0; i < Binning->GetNBins(); ++i) {
    SampleHandlerFD_array[i] = 0.0;
    SampleHandlerFD_array_w2[i] = 0.0;
    SampleHandlerFD_data[i] = 0.0;
  }
}

// ************************************************
void SampleHandlerFD::SetBinning() {
// ************************************************
  for(int iSample = 0; iSample < GetNsamples(); iSample++)
  {
    auto XVec = Binning->GetXBinEdges(iSample);
    auto YVec = Binning->GetYBinEdges(iSample);

    SampleDetails[iSample]._hPDF1D->Reset();
    SampleDetails[iSample]._hPDF1D->SetBins(static_cast<int>(XVec.size()-1), XVec.data());
    SampleDetails[iSample].dathist->SetBins(static_cast<int>(XVec.size()-1), XVec.data());

    SampleDetails[iSample]._hPDF2D->Reset();
    SampleDetails[iSample]._hPDF2D->SetBins(static_cast<int>(XVec.size()-1), XVec.data(),
                                static_cast<int>(YVec.size()-1), YVec.data());
    SampleDetails[iSample].dathist2d->SetBins(static_cast<int>(XVec.size()-1), XVec.data(),
                                  static_cast<int>(YVec.size()-1), YVec.data());

  }

  //Set the number of X and Y bins now
  SetupReweightArrays();
  FindNominalBinAndEdges();
}

// ************************************************
void SampleHandlerFD::FindNominalBinAndEdges() {
// ************************************************
  for (unsigned int event_i = 0; event_i < GetNEvents(); event_i++) {
    int Sample = MCSamples[event_i].NominalSample;
    const int dim = GetNDim(Sample);
    MCSamples[event_i].KinVar.resize(dim);
    MCSamples[event_i].NomBin.resize(dim);

    auto SetNominalBin = [&](int bin, int max_bins, int& out_bin) {
      if (bin >= 0 && bin < max_bins) {
        out_bin = bin;
      } else {
        out_bin = M3::UnderOverFlowBin;  // Out of bounds
      }
    };

    // Set and validate x_var (common to both cases)
    MCSamples[event_i].KinVar[0] = GetPointerToKinematicParameter(GetXBinVarName(Sample), event_i);
    if (std::isnan(*MCSamples[event_i].KinVar[0]) || std::isinf(*MCSamples[event_i].KinVar[0])) {
      MACH3LOG_ERROR("X var for event {} is ill-defined and equal to {}", event_i, *MCSamples[event_i].KinVar[0]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (dim == 1) {
      // 1D case
      int bin_x = Binning->FindNominalBin(Sample, 0, *MCSamples[event_i].KinVar[0]);
      SetNominalBin(bin_x, Binning->GetNXBins(Sample), MCSamples[event_i].NomBin[0]);
    }
    else if (dim == 2){
      // 2D case - set and validate y_var
      MCSamples[event_i].KinVar[1] = GetPointerToKinematicParameter(GetYBinVarName(Sample), event_i);
      if (std::isnan(*MCSamples[event_i].KinVar[1]) || std::isinf(*MCSamples[event_i].KinVar[1])) {
        MACH3LOG_ERROR("Y var for event {} is ill-defined and equal to {}", event_i, *MCSamples[event_i].KinVar[1]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      // Get global bin and convert to x/y bins
      int bin_x = Binning->FindNominalBin(Sample, 0, *MCSamples[event_i].KinVar[0]);
      int bin_y = Binning->FindNominalBin(Sample, 1, *MCSamples[event_i].KinVar[1]);

      SetNominalBin(bin_x, Binning->GetNXBins(Sample), MCSamples[event_i].NomBin[0]);
      SetNominalBin(bin_y, Binning->GetNYBins(Sample), MCSamples[event_i].NomBin[1]);
    } else {
      MACH3LOG_ERROR("Unsupported dimensionality: {}", dim);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// ************************************************
int SampleHandlerFD::GetSampleIndex(const std::string& SampleTitle) const {
// ************************************************
  for (M3::int_t iSample = 0; iSample < nSamples; ++iSample) {
    if (SampleTitle == GetSampleTitle(iSample)) {
      return iSample;
    }
  }
  MACH3LOG_ERROR("Sample name not found: {}", SampleTitle);
  throw MaCh3Exception(__FILE__, __LINE__);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// ************************************************
TH1* SampleHandlerFD::GetW2Hist(const int Sample, const int Dimension) {
// ************************************************
  if(Dimension == 1) {
    TH1D* W2Hist = dynamic_cast<TH1D*>(SampleDetails[Sample]._hPDF1D->Clone((SampleDetails[Sample]._hPDF1D->GetName() + std::string("_W2")).c_str()));
    if (!W2Hist) {
      MACH3LOG_ERROR("Failed to cast");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    W2Hist->Reset();
    for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
      for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
        const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
        W2Hist->SetBinContent(idx + 1, SampleHandlerFD_array_w2[idx]);
      }
    }
    return W2Hist;
  } else if(Dimension == 2) {
    TH2D* W2Hist = dynamic_cast<TH2D*>(SampleDetails[Sample]._hPDF2D->Clone((SampleDetails[Sample]._hPDF2D->GetName() + std::string("_W2")).c_str()));
    if (!W2Hist) {
      MACH3LOG_ERROR("Failed to cast");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    W2Hist->Reset();
    for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
      for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
        const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
        W2Hist->SetBinContent(xBin + 1, yBin + 1, SampleHandlerFD_array_w2[idx]);
      }
    }
    return W2Hist;
  } else{
    MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, Dimension);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}
#pragma GCC diagnostic pop

// ************************************************
TH1* SampleHandlerFD::GetW2Hist(const std::string& Sample, const int Dimension) {
// ************************************************
  const int Index = GetSampleIndex(Sample);
  return GetW2Hist(Index, Dimension);
}

// ************************************************
TH1* SampleHandlerFD::GetMCHist(const int Sample, const int Dimension) {
// ************************************************
  FillMCHist(Sample, Dimension);

  if(Dimension == 1) {
    return SampleDetails[Sample]._hPDF1D;
  } else if(Dimension == 2) {
    return SampleDetails[Sample]._hPDF2D;
  } else{
    MACH3LOG_ERROR("Asking for {} with N Dimension = {}. This is not implemented", __func__, Dimension);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
TH1* SampleHandlerFD::GetMCHist(const std::string& Sample, const int Dimension) {
// ************************************************
  const int Index = GetSampleIndex(Sample);
  return GetMCHist(Index, Dimension);
}

// ************************************************
TH1* SampleHandlerFD::GetDataHist(const int Sample, const int Dimension) {
// ************************************************
  if(Dimension == 1) {
    return SampleDetails[Sample].dathist;
  } else if(Dimension == 2) {
    return SampleDetails[Sample].dathist2d;
  } else{
    MACH3LOG_ERROR("Asdking for {} with N Dimension = {}. This is not implemented", __func__, Dimension);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
TH1* SampleHandlerFD::GetDataHist(const std::string& Sample, const int Dimension) {
// ************************************************
  int Index = GetSampleIndex(Sample);
  return GetDataHist(Index, Dimension);
}

void SampleHandlerFD::AddData(const int Sample, std::vector<double> &data) {
  if (SampleDetails[Sample].dathist == nullptr) {
    MACH3LOG_ERROR("Data hist hasn't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SampleDetails[Sample].dathist2d = nullptr;
  SampleDetails[Sample].dathist->Reset();

  if (GetNDim(Sample) != 1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram when the number of dimensions for this sample is {}", GetNDim(Sample));
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (auto const& data_point : data){
    SampleDetails[Sample].dathist->Fill(data_point);
  }

  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  // Assuming nBins == nXBins here, because you have 1D data
  for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
    for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
      const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
      // ROOT histograms are 1-based, so bin index + 1
      SampleHandlerFD_data[idx] = SampleDetails[Sample].dathist->GetBinContent(xBin + 1);
    }
  }
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
void SampleHandlerFD::AddData(const int Sample, std::vector< std::vector <double> > &data) {
  if (SampleDetails[Sample].dathist2d == nullptr) {
    MACH3LOG_ERROR("Data hist hasn't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SampleDetails[Sample].dathist = nullptr;
  SampleDetails[Sample].dathist2d->Reset();

  if (GetNDim(Sample) != 2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram when the number of dimensions for this sample is {}", GetNDim(Sample));
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //TODO: this assumes that std::vector is of length 2 and then both data.at(0) and data
  //ETA: I think this might just be wrong? We should probably just make this AddData(std::vector<double> data_x, std::vector<double> data_y)
  // or maybe something like AddData(std::vector<std::pair<double, double>> data)?
  for (int i = 0; i < int(data.size()); i++) {
    SampleDetails[Sample].dathist2d->Fill(data.at(0)[i],data.at(1)[i]);
  }

  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
    for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
      //Need to cast to an int (Int_t) for ROOT
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
      SampleHandlerFD_data[idx] = SampleDetails[Sample].dathist2d->GetBinContent(xBin + 1, yBin + 1);
    }
  }
}

void SampleHandlerFD::AddData(const int Sample, TH1D* Data) {
  MACH3LOG_INFO("Adding 1D data histogram: {} with {:.2f} events", Data->GetTitle(), Data->Integral());
  if (SampleDetails[Sample].dathist != nullptr) {
    delete SampleDetails[Sample].dathist;
  }
  SampleDetails[Sample].dathist2d = nullptr;
  SampleDetails[Sample].dathist = static_cast<TH1D*>(Data->Clone());

  if (GetNDim(Sample) != 1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram in a 2D sample - Quitting"); 
    MACH3LOG_ERROR("The number of dimensions for this sample is {}", GetNDim(Sample));
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
    
  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  SampleDetails[Sample].dathist->GetXaxis()->SetTitle(GetXBinVarName(Sample).c_str());
  SampleDetails[Sample].dathist->GetYaxis()->SetTitle("Number of Events");

  for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
    for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
      const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
      // ROOT histograms are 1-based, so bin index + 1
      SampleHandlerFD_data[idx] = SampleDetails[Sample].dathist->GetBinContent(xBin + 1);
    }
  }
}

void SampleHandlerFD::AddData(const int Sample, TH2D* Data) {
  MACH3LOG_INFO("Adding 2D data histogram: {} with {:.2f} events", Data->GetTitle(), Data->Integral());
  if (SampleDetails[Sample].dathist2d != nullptr) {
    delete SampleDetails[Sample].dathist2d;
  }
  SampleDetails[Sample].dathist2d = static_cast<TH2D*>(Data->Clone());
  SampleDetails[Sample].dathist = nullptr;

  if (GetNDim(Sample) != 2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram in a 1D sample - Quitting"); 
    throw MaCh3Exception(__FILE__ , __LINE__ );}
   
  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }


  SampleDetails[Sample].dathist2d->GetXaxis()->SetTitle(GetXBinVarName(Sample).c_str());
  SampleDetails[Sample].dathist2d->GetYaxis()->SetTitle(GetYBinVarName(Sample).c_str());
  SampleDetails[Sample].dathist2d->GetZaxis()->SetTitle("Number of Events");

  for (int yBin = 0; yBin < Binning->GetNYBins(Sample); ++yBin) {
    for (int xBin = 0; xBin < Binning->GetNXBins(Sample); ++xBin) {
      const int idx = Binning->GetGlobalBinSafe(Sample, xBin, yBin);
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      SampleHandlerFD_data[idx] = SampleDetails[Sample].dathist2d->GetBinContent(xBin + 1, yBin + 1);
    }
  }
}
#pragma GCC diagnostic pop

// ************************************************
void SampleHandlerFD::InitialiseNuOscillatorObjects() {
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
  std::vector<const double*> OscParams = ParHandler->GetOscParsFromSampleName(SampleHandlerName);
  if (OscParams.empty()) {
    MACH3LOG_ERROR("OscParams is empty for sample '{}'.", GetName());
    MACH3LOG_ERROR("This likely indicates an error in your oscillation YAML configuration.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  Oscillator = std::make_shared<OscillationHandler>(NuOscillatorConfigFile, EqualBinningPerOscChannel, OscParams, GetNOscChannels(0));
  // Add samples only if we don't use same binning
  if(!EqualBinningPerOscChannel) {
    // KS: Start from 1 because sample 0 already added
    for(int iSample = 1; iSample < GetNsamples(); iSample++) {
      Oscillator->AddSample(NuOscillatorConfigFile, GetNOscChannels(iSample));
    }
    for(int iSample = 0; iSample < GetNsamples(); iSample++) {
      for(int iChannel = 0; iChannel < GetNOscChannels(iSample); iChannel++) {
        std::vector<M3::float_t> EnergyArray;
        std::vector<M3::float_t> CosineZArray;

        for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
          if(MCSamples[iEvent].NominalSample != iSample) continue;
          // KS: This is bit weird but we basically loop over all events and push to vector only these which are part of a given OscChannel
          const int Channel = GetOscChannel(SampleDetails[MCSamples[iEvent].NominalSample].OscChannels, (*MCSamples[iEvent].nupdgUnosc), (*MCSamples[iEvent].nupdg));
          //DB Remove NC events from the arrays which are handed to the NuOscillator objects
          if (!MCSamples[iEvent].isNC && Channel == iChannel) {
            EnergyArray.push_back(M3::float_t(*(MCSamples[iEvent].rw_etru)));
          }
        }
        std::sort(EnergyArray.begin(),EnergyArray.end());

        //============================================================================
        //DB Atmospheric only part, can only happen if truecz has been initialised within the experiment specific code
        if (*(MCSamples[0].rw_truecz) != M3::_BAD_DOUBLE_) {
          for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
            if(MCSamples[iEvent].NominalSample != iSample) continue;
            // KS: This is bit weird but we basically loop over all events and push to vector only these which are part of a given OscChannel
            const int Channel = GetOscChannel(SampleDetails[MCSamples[iEvent].NominalSample].OscChannels, (*MCSamples[iEvent].nupdgUnosc), (*MCSamples[iEvent].nupdg));
            //DB Remove NC events from the arrays which are handed to the NuOscillator objects
            if (!MCSamples[iEvent].isNC && Channel == iChannel) {
              CosineZArray.push_back(M3::float_t(*(MCSamples[iEvent].rw_truecz)));
            }
          }
          std::sort(CosineZArray.begin(),CosineZArray.end());
        }
        Oscillator->SetOscillatorBinning(iSample, iChannel, EnergyArray, CosineZArray);
      } // end loop over channels
    } // end loop over samples
  }
}

void SampleHandlerFD::SetupNuOscillatorPointers() {
  auto AddOscPointer = GetFromManager<bool>(SampleManager->raw()["NuOsc"]["AddOscPointer"], true, __FILE__ , __LINE__);
  // Sometimes we don't want to add osc pointer to allow for some experiment specific handling of osc weight, like for example being able to shift osc weigh
  if(!AddOscPointer) {
    return;
  }
  for (unsigned int iEvent=0;iEvent<GetNEvents();iEvent++) {
   const M3::float_t* osc_w_pointer = GetNuOscillatorPointers(iEvent);

   // KS: Do not add unity
   if (osc_w_pointer != &M3::Unity) {
     MCSamples[iEvent].total_weight_pointers.push_back(osc_w_pointer);
   }
  }
}

const M3::float_t* SampleHandlerFD::GetNuOscillatorPointers(const int iEvent) const {
  const M3::float_t* osc_w_pointer = &M3::Unity;
  if (MCSamples[iEvent].isNC) {
    if (*MCSamples[iEvent].nupdg != *MCSamples[iEvent].nupdgUnosc) {
      osc_w_pointer = &M3::Zero;
    } else {
      osc_w_pointer = &M3::Unity;
    }
  } else {
    int InitFlav = M3::_BAD_INT_;
    int FinalFlav = M3::_BAD_INT_;

    InitFlav =  MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples[iEvent].nupdgUnosc));
    FinalFlav = MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples[iEvent].nupdg));

    if (InitFlav == M3::_BAD_INT_ || FinalFlav == M3::_BAD_INT_) {
      MACH3LOG_ERROR("Something has gone wrong in the mapping between MCSamples.nutype and the enum used within NuOscillator");
      MACH3LOG_ERROR("MCSamples.nupdgUnosc: {}", (*MCSamples[iEvent].nupdgUnosc));
      MACH3LOG_ERROR("InitFlav: {}", InitFlav);
      MACH3LOG_ERROR("MCSamples.nupdg: {}", (*MCSamples[iEvent].nupdg));
      MACH3LOG_ERROR("FinalFlav: {}", FinalFlav);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const int Sample = MCSamples[iEvent].NominalSample;

    const int OscIndex = GetOscChannel(SampleDetails[Sample].OscChannels, (*MCSamples[iEvent].nupdgUnosc), (*MCSamples[iEvent].nupdg));
    //Can only happen if truecz has been initialised within the experiment specific code
    if (*(MCSamples[iEvent].rw_truecz) != M3::_BAD_DOUBLE_) {
      //Atmospherics
      osc_w_pointer = Oscillator->GetNuOscillatorPointers(Sample, OscIndex, InitFlav, FinalFlav, FLOAT_T(*(MCSamples[iEvent].rw_etru)), FLOAT_T(*(MCSamples[iEvent].rw_truecz)));
    } else {
      //Beam
      osc_w_pointer = Oscillator->GetNuOscillatorPointers(Sample, OscIndex, InitFlav, FinalFlav, FLOAT_T(*(MCSamples[iEvent].rw_etru)));
    }
  } // end if NC
  return osc_w_pointer;
}

std::string SampleHandlerFD::GetName() const {
  //ETA - extra safety to make sure SampleHandlerName is actually set
  // probably unnecessary due to the requirement for it to be in the yaml config
  if(SampleHandlerName.length() == 0) {
    MACH3LOG_ERROR("No sample name provided");
    MACH3LOG_ERROR("Please provide a SampleHandlerName in your configuration file: {}", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleHandlerName;
}

M3::float_t SampleHandlerFD::GetEventWeight(const int iEntry) {
  // KS: WARNING we have to here recalculate weight and cap because there is possibility weight wasn't calculated during FillArray because it didn't fulfil IsEventSelected

  // Virtual by default does nothing, has to happen before CalcWeightTotal
  CalcWeightFunc(iEntry);

  const FarDetectorCoreInfo* MCEvent = &MCSamples[iEntry];
  M3::float_t totalweight = CalcWeightTotal(MCEvent);

  //DB Catch negative total weights and skip any event with a negative weight. Previously we would set weight to zero and continue but that is inefficient
  if (totalweight <= 0.){
    totalweight = 0.;
  }
  return totalweight;
}

// Finds the binned spline that an event should apply to and stored them in a
// a vector for easy evaluation in the fillArray() function.
void SampleHandlerFD::FillSplineBins() {
  //Now loop over events and get the spline bin for each event
  for (unsigned int j = 0; j < GetNEvents(); ++j) {
    const int SampleIndex = MCSamples[j].NominalSample;
    const int OscIndex = GetOscChannel(SampleDetails[SampleIndex].OscChannels, (*MCSamples[j].nupdgUnosc), (*MCSamples[j].nupdg));

    std::vector< std::vector<int> > EventSplines;
    switch(GetNDim(SampleIndex)){
      case 1:
        EventSplines = SplineHandler->GetEventSplines(GetSampleTitle(SampleIndex), OscIndex, int(*(MCSamples[j].mode)), *(MCSamples[j].rw_etru), *(MCSamples[j].KinVar[0]), 0.);
        break;
      case 2:
        EventSplines = SplineHandler->GetEventSplines(GetSampleTitle(SampleIndex), OscIndex, int(*(MCSamples[j].mode)), *(MCSamples[j].rw_etru), *(MCSamples[j].KinVar[0]), *(MCSamples[j].KinVar[1]));
        break;
      default:
        MACH3LOG_ERROR("Error in assigning spline bins because nDimensions = {}", GetNDim(SampleIndex));
        MACH3LOG_ERROR("MaCh3 only supports splines binned in Etrue + the sample binning");
        MACH3LOG_ERROR("Please check the sample binning you specified in your sample config ");
        throw MaCh3Exception(__FILE__, __LINE__);
        break;
    }
    const int NSplines = static_cast<int>(EventSplines.size());
    if(NSplines == 0) continue;
    const int PointersBefore = static_cast<int>(MCSamples[j].total_weight_pointers.size());
    MCSamples[j].total_weight_pointers.resize(PointersBefore + NSplines);

    for(int spline = 0; spline < NSplines; spline++) {
      //Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
      MCSamples[j].total_weight_pointers[PointersBefore+spline] = SplineHandler->retPointer(EventSplines[spline][0], EventSplines[spline][1],
                                                                    EventSplines[spline][2], EventSplines[spline][3],
                                                                    EventSplines[spline][4], EventSplines[spline][5],
                                                                    EventSplines[spline][6]);
    }
  }
}

// ************************************************
double SampleHandlerFD::GetSampleLikelihood(const int isample) const {
// ************************************************
  const int Start = Binning->GetSampleStartBin(isample);
  const int End = Binning->GetSampleEndBin(isample);

  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:negLogL)
  #endif
  for (int idx = Start; idx < End; ++idx)
  {
    const double DataVal = SampleHandlerFD_data[idx];
    const double MCPred = SampleHandlerFD_array[idx];
    const double w2 = SampleHandlerFD_array_w2[idx];

    //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
    negLogL += GetTestStatLLH(DataVal, MCPred, w2);
  }
  return negLogL;
}

// ************************************************
double SampleHandlerFD::GetLikelihood() const {
// ************************************************
  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:negLogL)
  #endif
  for (int idx = 0; idx < Binning->GetNBins(); ++idx)
  {
    const double DataVal = SampleHandlerFD_data[idx];
    const double MCPred = SampleHandlerFD_array[idx];
    const double w2 = SampleHandlerFD_array_w2[idx];

    //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
    negLogL += GetTestStatLLH(DataVal, MCPred, w2);
  }
  return negLogL;
}

// ************************************************
void SampleHandlerFD::SaveAdditionalInfo(TDirectory* Dir) {
// ************************************************
  Dir->cd();

  YAML::Node Config = SampleManager->raw();
  TMacro ConfigSave = YAMLtoTMacro(Config, (std::string("Config_") + GetName()));
  ConfigSave.Write();

  for(int iSample = 0; iSample < GetNsamples(); iSample++)
  {
    std::unique_ptr<TH1> data_hist;

    if (GetNDim(iSample) == 1) {
      data_hist = M3::Clone<TH1D>(dynamic_cast<TH1D*>(GetDataHist(iSample, 1)), "data_" + GetSampleTitle(iSample));
      data_hist->GetXaxis()->SetTitle(GetXBinVarName(iSample).c_str());
      data_hist->GetYaxis()->SetTitle("Number of Events");
    } else if (GetNDim(iSample) == 2) {
      data_hist = M3::Clone<TH2D>(dynamic_cast<TH2D*>(GetDataHist(iSample, 2)), "data_" + GetSampleTitle(iSample));
      data_hist->GetXaxis()->SetTitle(GetXBinVarName(iSample).c_str());
      data_hist->GetYaxis()->SetTitle(GetYBinVarName(iSample).c_str());
      data_hist->GetZaxis()->SetTitle("Number of Events");
    } else {
      MACH3LOG_ERROR("Not implemented for dimension {}", GetNDim(iSample));
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (!data_hist) {
      MACH3LOG_ERROR("nullptr data hist :(");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    data_hist->SetTitle(("data_" + GetSampleTitle(iSample)).c_str());
    data_hist->Write();
  }
}

void SampleHandlerFD::InitialiseSplineObject() {
  bool LoadSplineFile = GetFromManager<bool>(SampleManager->raw()["InputFiles"]["LoadSplineFile"], false, __FILE__, __LINE__);
  bool PrepSplineFile = GetFromManager<bool>(SampleManager->raw()["InputFiles"]["PrepSplineFile"], false, __FILE__, __LINE__);
  auto SplineFileName = GetFromManager<std::string>(SampleManager->raw()["InputFiles"]["SplineFileName"],
                                                    (SampleHandlerName + "_SplineFile.root"), __FILE__, __LINE__);
  if(!LoadSplineFile) {
    for(int iSample = 0; iSample < GetNsamples(); iSample++) {
      std::vector<std::string> spline_filepaths = SampleDetails[iSample].spline_files;

      //Keep a track of the spline variables
      std::vector<std::string> SplineVarNames = {"TrueNeutrinoEnergy"};
      if(GetXBinVarName(iSample).length() > 0){
        SplineVarNames.push_back(GetXBinVarName(iSample));
      }
      if(GetYBinVarName(iSample).length() > 0){
        SplineVarNames.push_back(GetYBinVarName(iSample));
      }
      SplineHandler->AddSample(SampleHandlerName, GetSampleTitle(iSample), spline_filepaths, SplineVarNames);
    }
    SplineHandler->CountNumberOfLoadedSplines(false, 1);
    SplineHandler->TransferToMonolith();
    if(PrepSplineFile) SplineHandler->PrepareSplineFile(SplineFileName);
  } else {
    // KS: Skip default spline loading and use flattened spline format allowing to read stuff much faster
    SplineHandler->LoadSplineFile(SplineFileName);
  }
  MACH3LOG_INFO("--------------------------------");
  MACH3LOG_INFO("Setup Far Detector splines");

  FillSplineBins();

  SplineHandler->cleanUpMemory();
}

// ************************************************
// === JM adjust GetNDVarHist functions to allow for subevent-level plotting ===
TH1* SampleHandlerFD::Get1DVarHist(const int iSample, const std::string& ProjectionVar_Str, const std::vector< KinematicCut >& EventSelectionVec,
    int WeightStyle, TAxis* Axis, const std::vector< KinematicCut >& SubEventSelectionVec) {
// ************************************************
  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector< KinematicCut > > tmp_Selection = Selection;
  std::vector< std::vector< KinematicCut > > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  SelectionVecToApply = tmp_Selection;

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<EventSelectionVec.size();iSelec++) {
    SelectionVecToApply[iSample].emplace_back(EventSelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH1D* _h1DVar = nullptr;;
  if (Axis) {
    _h1DVar = new TH1D("","",Axis->GetNbins(),Axis->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_Str);
    _h1DVar = new TH1D("", "", int(xBinEdges.size())-1, xBinEdges.data());
  }
  _h1DVar->GetXaxis()->SetTitle(ProjectionVar_Str.c_str());

  if (IsSubEventVarString(ProjectionVar_Str)) {
    Fill1DSubEventHist(iSample, _h1DVar, ProjectionVar_Str, SubEventSelectionVec, WeightStyle);
  } else {
    //DB Grab the associated enum with the argument string
    int ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);

    //DB Loop over all events
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
      const int EventSample = MCSamples[iEvent].NominalSample;
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

void SampleHandlerFD::Fill1DSubEventHist(const int iSample, TH1D* _h1DVar, const std::string& ProjectionVar_Str, const std::vector< KinematicCut >& SubEventSelectionVec, int WeightStyle) {
  int ProjectionVar_Int = ReturnKinematicVectorFromString(ProjectionVar_Str);

  //JM Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    const int EventSample = MCSamples[iEvent].NominalSample;
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
TH2* SampleHandlerFD::Get2DVarHist(const int iSample, const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY,
    const std::vector< KinematicCut >& EventSelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY, const std::vector< KinematicCut >& SubEventSelectionVec) {
// ************************************************
  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector< KinematicCut > > tmp_Selection = Selection;
  std::vector< std::vector< KinematicCut > > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  SelectionVecToApply = tmp_Selection;

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<EventSelectionVec.size();iSelec++) {
    SelectionVecToApply[iSample].emplace_back(EventSelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH2D* _h2DVar = nullptr;
  if (AxisX && AxisY) {
    _h2DVar = new TH2D("","",AxisX->GetNbins(),AxisX->GetXbins()->GetArray(),AxisY->GetNbins(),AxisY->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_StrX);
    std::vector<double> yBinEdges = ReturnKinematicParameterBinning(iSample, ProjectionVar_StrY);
    _h2DVar = new TH2D("", "", int(xBinEdges.size())-1, xBinEdges.data(), int(yBinEdges.size())-1, yBinEdges.data());
  }
  _h2DVar->GetXaxis()->SetTitle(ProjectionVar_StrX.c_str());
  _h2DVar->GetYaxis()->SetTitle(ProjectionVar_StrY.c_str());

  bool IsSubEventHist = IsSubEventVarString(ProjectionVar_StrX) || IsSubEventVarString(ProjectionVar_StrY);
  if (IsSubEventHist) Fill2DSubEventHist(iSample, _h2DVar, ProjectionVar_StrX, ProjectionVar_StrY, SubEventSelectionVec, WeightStyle);
  else {
    //DB Grab the associated enum with the argument string
    int ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
    int ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

    //DB Loop over all events
    for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
      const int EventSample = MCSamples[iEvent].NominalSample;
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

void SampleHandlerFD::Fill2DSubEventHist(const int iSample, TH2D* _h2DVar, const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY,
    const std::vector< KinematicCut >& SubEventSelectionVec, int WeightStyle) {
  bool IsSubEventVarX = IsSubEventVarString(ProjectionVar_StrX);
  bool IsSubEventVarY = IsSubEventVarString(ProjectionVar_StrY);   

  int ProjectionVar_IntX, ProjectionVar_IntY;
  if (IsSubEventVarX) ProjectionVar_IntX = ReturnKinematicVectorFromString(ProjectionVar_StrX);
  else ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  if (IsSubEventVarY) ProjectionVar_IntY = ReturnKinematicVectorFromString(ProjectionVar_StrY);
  else ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY); 

  //JM Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    const int EventSample = MCSamples[iEvent].NominalSample;
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
int SampleHandlerFD::ReturnKinematicParameterFromString(const std::string& KinematicParameterStr) const {
// ************************************************
  auto it = KinematicParameters->find(KinematicParameterStr);
  if (it != KinematicParameters->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameterStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ************************************************
std::string SampleHandlerFD::ReturnStringFromKinematicParameter(const int KinematicParameter) const {
// ************************************************
  auto it = ReversedKinematicParameters->find(KinematicParameter);
  if (it != ReversedKinematicParameters->end()) {
    return it->second;
  }

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameter);
  throw MaCh3Exception(__FILE__, __LINE__);

  return "";
}

// === JM define KinematicVector-to-string mapping functions  ===
// ************************************************
int SampleHandlerFD::ReturnKinematicVectorFromString(const std::string& KinematicVectorStr) const {
// ************************************************
  auto it = KinematicVectors->find(KinematicVectorStr);
  if (it != KinematicVectors->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Vector: {}", KinematicVectorStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ************************************************
std::string SampleHandlerFD::ReturnStringFromKinematicVector(const int KinematicVector) const {
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
std::vector<double> SampleHandlerFD::ReturnKinematicParameterBinning(const int Sample, const std::string& KinematicParameter) const {
// ************************************************
  // If x or y variable return used binning
  if(KinematicParameter == GetXBinVarName(Sample)) {
    return Binning->GetXBinEdges(Sample);
  } else if (KinematicParameter == GetYBinVarName(Sample)) {
    return Binning->GetYBinEdges(Sample);
  }

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

  // We first check if binning for a sample has been specified
  auto BinningConfig = M3OpenConfig(SampleManager->raw()["BinningFile"].as<std::string>());
  if(BinningConfig[GetSampleTitle(Sample)] && BinningConfig[GetSampleTitle(Sample)][KinematicParameter]){
    auto BinningVect = Get<std::vector<double>>(BinningConfig[GetSampleTitle(Sample)][KinematicParameter], __FILE__, __LINE__);
    return BinningVect;
  } else {
    auto BinningVect = Get<std::vector<double>>(BinningConfig[KinematicParameter], __FILE__, __LINE__);
    return BinningVect;
  }
}


bool SampleHandlerFD::IsSubEventVarString(const std::string& VarStr) {
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
// ===============================================================

TH1* SampleHandlerFD::Get1DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_Str,
    int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
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

  return Get1DVarHist(iSample, ProjectionVar_Str,SelectionVec,WeightStyle,Axis);
}

TH2* SampleHandlerFD::Get2DVarHistByModeAndChannel(const int iSample, const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY,
    int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
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

  return Get2DVarHist(iSample, ProjectionVar_StrX,ProjectionVar_StrY,SelectionVec,WeightStyle,AxisX,AxisY);
}

void SampleHandlerFD::PrintIntegral(const int iSample, const TString& OutputFileName, const int WeightStyle, const TString& OutputCSVFileName) {
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

  std::vector< std::vector< TH1* > > IntegralList;
  IntegralList.resize(Modes->GetNModes());

  std::vector<double> ChannelIntegral;
  ChannelIntegral.resize(GetNOscChannels(iSample));
  for (unsigned int i=0;i<ChannelIntegral.size();i++) {ChannelIntegral[i] = 0.;}

  for (int i=0;i<Modes->GetNModes();i++) {
    if (GetNDim(iSample) == 1) {
      IntegralList[i] = ReturnHistsBySelection1D(iSample, GetXBinVarName(iSample),1,i,WeightStyle);
    } else {
      IntegralList[i] = CastVector<TH2, TH1>(ReturnHistsBySelection2D(iSample, GetXBinVarName(iSample), GetYBinVarName(iSample),1,i,WeightStyle));
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
  // KS: Clean memory we could use smart pointers in future
  CleanContainer(IntegralList);
}

std::vector<TH1*> SampleHandlerFD::ReturnHistsBySelection1D(const int iSample, std::string KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
  std::vector<TH1*> hHistList;
  std::string legendEntry;

  if (THStackLeg != nullptr) {delete THStackLeg;}
  THStackLeg = new TLegend(0.1,0.1,0.9,0.9);

  int iMax = -1;
  if (Selection1 == FDPlotType::kModePlot) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == FDPlotType::kOscChannelPlot) {
    iMax = GetNOscChannels(iSample);
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1 == FDPlotType::kModePlot) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(iSample, KinematicProjection,i,Selection2,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(Modes->GetMaCh3ModeName(i)+Form(" : (%4.2f)",hHistList[i]->Integral())).c_str(),"f");

      hHistList[i]->SetFillColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
      hHistList[i]->SetLineColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
    }
    if (Selection1 == FDPlotType::kOscChannelPlot) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(iSample, KinematicProjection,Selection2,i,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(GetFlavourName(iSample, i)+Form(" | %4.2f",hHistList[i]->Integral())).c_str(),"f");
    }
  }

  return hHistList;
}
// ************************************************
std::vector<TH2*> SampleHandlerFD::ReturnHistsBySelection2D(const int iSample, std::string KinematicProjectionX, std::string KinematicProjectionY,
                                                            int Selection1, int Selection2, int WeightStyle,
                                                            TAxis* XAxis, TAxis* YAxis) {
// ************************************************
  std::vector<TH2*> hHistList;

  int iMax = -1;
  if (Selection1 == FDPlotType::kModePlot) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == FDPlotType::kOscChannelPlot) {
    iMax = GetNOscChannels(iSample);
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1 == FDPlotType::kModePlot) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(iSample, KinematicProjectionX,KinematicProjectionY,i,Selection2,WeightStyle,XAxis,YAxis));
    }
    if (Selection1 == FDPlotType::kOscChannelPlot) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(iSample, KinematicProjectionX,KinematicProjectionY,Selection2,i,WeightStyle,XAxis,YAxis));
    }
  }

  return hHistList;
}

// ************************************************
THStack* SampleHandlerFD::ReturnStackedHistBySelection1D(const int iSample, const std::string& KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
// ************************************************
  std::vector<TH1*> HistList = ReturnHistsBySelection1D(iSample, KinematicProjection, Selection1, Selection2, WeightStyle, XAxis);
  THStack* StackHist = new THStack((GetSampleTitle(iSample)+"_"+KinematicProjection+"_Stack").c_str(),"");
  for (unsigned int i=0;i<HistList.size();i++) {
    StackHist->Add(HistList[i]);
  }
  return StackHist;
}


// ************************************************
const double* SampleHandlerFD::GetPointerToOscChannel(const int iEvent) const {
// ************************************************
  auto& OscillationChannels = SampleDetails[MCSamples[iEvent].NominalSample].OscChannels;
  const int Channel = GetOscChannel(OscillationChannels, (*MCSamples[iEvent].nupdgUnosc), (*MCSamples[iEvent].nupdg));

  return &(OscillationChannels[Channel].ChannelIndex);
}

// ***************************************************************************
// Helper function to print rates for the samples with LLH
void SampleHandlerFD::PrintRates(const bool DataOnly) {
// ***************************************************************************
  if (SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("Data sample is empty!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Printing for {}", GetName());

  if (!DataOnly) {
    const std::string sep_full(71, '-');
    MACH3LOG_INFO("{}", sep_full);
    MACH3LOG_INFO("{:<40}{:<10}{:<10}{:<10}|", "Sample", "Data", "MC", "-LLH");
  } else {
    const std::string sep_data(51, '-');
    MACH3LOG_INFO("{}", sep_data);
    MACH3LOG_INFO("{:<40}{:<10}|", "Sample", "Data");
  }

  double sumData = 0.0;
  double sumMC = 0.0;
  double likelihood = 0.0;

  for (int iSample = 0; iSample < GetNsamples(); ++iSample) {
    std::string name = GetSampleTitle(iSample);
    double dataIntegral = GetDataHist(iSample, GetNDim(iSample))->Integral();
    sumData += dataIntegral;
    if (!DataOnly) {
      double mcIntegral = GetMCHist(iSample, GetNDim(iSample))->Integral();
      sumMC += mcIntegral;
      likelihood = GetSampleLikelihood(iSample);

      MACH3LOG_INFO("{:<40}{:<10.2f}{:<10.2f}{:<10.2f}|", name, dataIntegral, mcIntegral, likelihood);
    } else {
      MACH3LOG_INFO("{:<40}{:<10.2f}|", name, dataIntegral);
    }
  }
  if (!DataOnly) {
    likelihood = GetLikelihood();
    MACH3LOG_INFO("{:<40}{:<10.2f}{:<10.2f}{:<10.2f}|", "Total", sumData, sumMC, likelihood);
    const std::string sep_full(71, '-');
    MACH3LOG_INFO("{}", sep_full);
  } else {
    MACH3LOG_INFO("{:<40}{:<10.2f}|", "Total", sumData);
    const std::string sep_data(51, '-');
    MACH3LOG_INFO("{}", sep_data);
  }
}
