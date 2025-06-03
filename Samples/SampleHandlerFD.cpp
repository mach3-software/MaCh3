#include "SampleHandlerFD.h"
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include <cstddef>

#include <algorithm>
#include <memory>

// ************************************************
SampleHandlerFD::SampleHandlerFD(std::string ConfigFileName, ParameterHandlerGeneric* xsec_cov,
                                 ParameterHandlerOsc* osc_cov, const std::shared_ptr<OscillationHandler>& OscillatorObj_) : SampleHandlerBase() {
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
  if(!osc_cov){
    MACH3LOG_WARN("You have passed a nullptr to a covarianceOsc, this means I will not calculate oscillation weights");
  }
  OscParHandler = osc_cov;

  if (OscillatorObj_ != nullptr) {
    MACH3LOG_WARN("You have passed an Oscillator object through the constructor of a SampleHandlerFD object - this will be used for all oscillation channels");
    Oscillator = OscillatorObj_;
  }

  KinematicParameters = nullptr;
  ReversedKinematicParameters = nullptr;

  SampleHandlerFD_array = nullptr;
  SampleHandlerFD_data = nullptr;
  SampleHandlerFD_array_w2 = nullptr;
  SampleName = "";
  SampleManager = std::unique_ptr<manager>(new manager(ConfigFileName.c_str()));

  // Variables related to MC stat
  FirstTimeW2 = true;
  UpdateW2 = false;
}

SampleHandlerFD::~SampleHandlerFD()
{
  MACH3LOG_DEBUG("I'm deleting SampleHandlerFD");
  
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    if(SampleHandlerFD_array != nullptr){delete[] SampleHandlerFD_array[yBin];}
    delete[] SampleHandlerFD_array_w2[yBin];
    //ETA - there is a chance that you haven't added any data...
    if(SampleHandlerFD_data != nullptr){delete[] SampleHandlerFD_data[yBin];}
  }

  if(SampleHandlerFD_array != nullptr){delete[] SampleHandlerFD_array;}
  delete[] SampleHandlerFD_array_w2;
  //ETA - there is a chance that you haven't added any data...
  if(SampleHandlerFD_data != nullptr){delete[] SampleHandlerFD_data;}
 
  if(THStackLeg != nullptr) delete THStackLeg;
}

void SampleHandlerFD::ReadSampleConfig() 
{
  auto ModeName = Get<std::string>(SampleManager->raw()["MaCh3ModeConfig"], __FILE__ , __LINE__);
  Modes = new MaCh3Modes(ModeName);
  //SampleTitle has to be provided in the sample yaml otherwise this will throw an exception
  SampleTitle = Get<std::string>(SampleManager->raw()["SampleTitle"], __FILE__ , __LINE__);
  //SampleName has to be provided in the sample yaml otherwise this will throw an exception
  SampleName = Get<std::string>(SampleManager->raw()["SampleName"], __FILE__ , __LINE__);

  fTestStatistic = static_cast<TestStatistic>(SampleManager->GetMCStatLLH());
  if (CheckNodeExists(SampleManager->raw(), "LikelihoodOptions")) {
    UpdateW2 = GetFromManager<bool>(SampleManager->raw()["LikelihoodOptions"]["UpdateW2"], false);
  }
  //Binning
  nDimensions = 0;
  XVarStr = GetFromManager(SampleManager->raw()["Binning"]["XVarStr"], std::string(""));
  SampleXBins = GetFromManager(SampleManager->raw()["Binning"]["XVarBins"], std::vector<double>());
  if(XVarStr.length() > 0){
    nDimensions++;
  } else{
    MACH3LOG_ERROR("Please specify an X-variable string in sample config {}", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  YVarStr = GetFromManager(SampleManager->raw()["Binning"]["YVarStr"], std::string(""));
  SampleYBins = GetFromManager(SampleManager->raw()["Binning"]["YVarBins"], std::vector<double>());
  if(YVarStr.length() > 0){
    if(XVarStr.length() == 0){
      MACH3LOG_ERROR("Please specify an X-variable string in sample config {}. I won't work only with a Y-variable", SampleManager->GetFileName());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    nDimensions++;
  }
  
  if(nDimensions == 0){
    MACH3LOG_ERROR("Error setting up the sample binning");
    MACH3LOG_ERROR("Number of dimensions is {}", nDimensions);
    MACH3LOG_ERROR("Check that an XVarStr has been given in the sample config");
    throw MaCh3Exception(__FILE__, __LINE__);
  } else{
    MACH3LOG_INFO("Found {} dimensions for sample binning", nDimensions);
  }
  
  //Sanity check that some binning has been specified
  if(SampleXBins.size() == 0 && SampleYBins.size() == 0){
    MACH3LOG_ERROR("No binning specified for either X or Y of sample binning, please add some binning to the sample config {}", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  //FD file info
  if (!CheckNodeExists(SampleManager->raw(), "InputFiles", "mtupleprefix")){
    MACH3LOG_ERROR("InputFiles:mtupleprefix not given in {}, please add this", SampleManager->GetFileName());
  }
  auto mtupleprefix  = Get<std::string>(SampleManager->raw()["InputFiles"]["mtupleprefix"], __FILE__, __LINE__);
  auto mtuplesuffix  = Get<std::string>(SampleManager->raw()["InputFiles"]["mtuplesuffix"], __FILE__, __LINE__);
  auto splineprefix  = Get<std::string>(SampleManager->raw()["InputFiles"]["splineprefix"], __FILE__, __LINE__);
  auto splinesuffix  = Get<std::string>(SampleManager->raw()["InputFiles"]["splinesuffix"], __FILE__, __LINE__);
  
  int NChannels = static_cast<M3::int_t>(SampleManager->raw()["SubSamples"].size());
  OscChannels.reserve(NChannels);

  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    OscChannelInfo OscInfo;
    OscInfo.flavourName       = osc_channel["Name"].as<std::string>();
    OscInfo.flavourName_Latex = osc_channel["LatexName"].as<std::string>();
    OscInfo.InitPDG           = static_cast<NuPDG>(osc_channel["nutype"].as<int>());
    OscInfo.FinalPDG          = static_cast<NuPDG>(osc_channel["oscnutype"].as<int>());
    OscInfo.ChannelIndex      = static_cast<int>(OscChannels.size());

    OscChannels.push_back(std::move(OscInfo));

    mc_files.push_back(mtupleprefix+osc_channel["mtuplefile"].as<std::string>()+mtuplesuffix);
    spline_files.push_back(splineprefix+osc_channel["splinefile"].as<std::string>()+splinesuffix);
  }

  //Now grab the selection cuts from the manager
  for ( auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
    SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());
    auto TempBoundsVec = GetBounds(SelectionCuts["Bounds"]);
    KinematicCut CutObj;
    CutObj.LowerBound = TempBoundsVec[0];
    CutObj.UpperBound = TempBoundsVec[1];
    CutObj.ParamToCutOnIt = ReturnKinematicParameterFromString(SelectionCuts["KinematicStr"].as<std::string>());
    MACH3LOG_INFO("Adding cut on {} with bounds {} to {}", SelectionCuts["KinematicStr"].as<std::string>(), TempBoundsVec[0], TempBoundsVec[1]);
    StoredSelection.push_back(CutObj);
  }

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

void SampleHandlerFD::Initialise() {
  //First grab all the information from your sample config via your manager
  ReadSampleConfig();

  //Now initialise all the variables you will need
  Init();

  nEvents = SetupExperimentMC();

  InitialiseSingleFDMCObject();
  SetupFDMC();

  MACH3LOG_INFO("=============================================");
  MACH3LOG_INFO("Total number of events is: {}", GetNEvents());

  if (OscParHandler) {
    MACH3LOG_INFO("Setting up NuOscillator..");
    if (Oscillator != nullptr) {
      MACH3LOG_INFO("You have passed an OscillatorBase object through the constructor of a SampleHandlerFD object - this will be used for all oscillation channels");
      if(Oscillator->isEqualBinningPerOscChannel() != true) {
        MACH3LOG_ERROR("Trying to run shared NuOscillator without EqualBinningPerOscChannel, this will not work");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      auto OscParams = OscParHandler->GetOscParsFromSampleName(SampleName);
      if(OscParams.size() != Oscillator->GetOscParamsSize()){
        MACH3LOG_ERROR("Sample {} with {} has {} osc params, while shared NuOsc has {} osc params", GetTitle(), GetSampleName(),
                       OscParams.size(), Oscillator->GetOscParamsSize());
        MACH3LOG_ERROR("This indicate misconfiguration in your Osc yaml");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } else {
      InitialiseNuOscillatorObjects();
    }
    SetupNuOscillatorPointers();
  }

  MACH3LOG_INFO("Setting up Sample Binning..");
  SetupSampleBinning();
  MACH3LOG_INFO("Setting up Splines..");
  SetupSplines();
  MACH3LOG_INFO("Setting up Normalisation Pointers..");
  SetupNormParameters();
  MACH3LOG_INFO("Setting up Functional Pointers..");
  SetupFunctionalParameters();
  MACH3LOG_INFO("Setting up Weight Pointers..");
  SetupWeightPointers();
  MACH3LOG_INFO("Setting up Kinematic Map..");
  SetupKinematicMap();
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

void SampleHandlerFD::Fill1DHist()
{
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  _hPDF1D->Reset();
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    for (unsigned int xBin=0;xBin<(XBinEdges.size()-1);xBin++) {
      _hPDF1D->AddBinContent(xBin+1,SampleHandlerFD_array[yBin][xBin]);
    }
  }
}

void SampleHandlerFD::Fill2DHist()
{
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  _hPDF2D->Reset();
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    for (unsigned int xBin=0;xBin<(XBinEdges.size()-1);xBin++) {
      _hPDF2D->SetBinContent(xBin+1,yBin+1,SampleHandlerFD_array[yBin][xBin]);
    }
  }
}

// ************************************************
/// @function SampleHandlerFD::SetupSampleBinning()
/// @brief Function to setup the binning of your sample histograms and the underlying
/// arrays that get handled in fillArray() and fillArray_MP().
/// The SampleXBins are filled in the daughter class from the sample config file.
/// This "passing" can be removed. 
void SampleHandlerFD::SetupSampleBinning(){
// ************************************************
  MACH3LOG_INFO("Setting up Sample Binning");
  TString histname1d = (XVarStr).c_str();
  TString histname2d = (XVarStr+"_"+YVarStr).c_str();
  TString histtitle = "";

  //The binning here is arbitrary, now we get info from cfg so the
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D   = new TH1D("h"+histname1d+SampleTitle,histtitle, 1, 0, 1);
  dathist   = new TH1D("d"+histname1d+SampleTitle,histtitle, 1, 0, 1);
  _hPDF2D   = new TH2D("h"+histname2d+SampleTitle,histtitle, 1, 0, 1, 1, 0, 1);
  dathist2d = new TH2D("d"+histname2d+SampleTitle,histtitle, 1, 0, 1, 1, 0, 1);

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  XBinEdges.reserve(SampleXBins.size());
  YBinEdges.reserve(SampleYBins.size());

  //A string to store the binning for a nice print out
  std::string XBinEdgesStr = "";
  std::string YBinEdgesStr = "";

  for(auto XBinEdge : SampleXBins){
    XBinEdges.push_back(XBinEdge);
    XBinEdgesStr += std::to_string(XBinEdge);
    XBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("XBinning:");
  MACH3LOG_INFO("{}", XBinEdgesStr);
  
  //And now the YBin Edges
  for(auto YBinEdge : SampleYBins){
    YBinEdges.push_back(YBinEdge);
    YBinEdgesStr += std::to_string(YBinEdge);
    YBinEdgesStr += ", ";
  }
  MACH3LOG_INFO("YBinning:");
  MACH3LOG_INFO("{}", YBinEdgesStr);
  
  //Check whether you are setting up 1D or 2D binning
  if(nDimensions == 1){
    MACH3LOG_INFO("Setting up 1D binning with {}", XVarStr);
    Set1DBinning(SampleXBins);  
  }
  else if(nDimensions == 2){
    MACH3LOG_INFO("Setting up 2D binning with {} and {}", XVarStr, YVarStr);
    Set2DBinning(SampleXBins, SampleYBins);
  }
  else{
    MACH3LOG_ERROR("Number of dimensions is not 1 or 2, this is unsupported at the moment");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
bool SampleHandlerFD::IsEventSelected(const int iEvent) {
// ************************************************
  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {
    const double Val = ReturnKinematicParameter(Selection[iSelection].ParamToCutOnIt, iEvent);
    if ((Val < Selection[iSelection].LowerBound) || (Val >= Selection[iSelection].UpperBound)) {
      return false;
    }
  }
  //DB To avoid unnecessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

// === JM Define function to check if particle is selected ===
bool SampleHandlerFD::IsParticleSelected(const std::vector<KinematicCut> &ParticleCuts, const int iEvent, const unsigned int iParticle, size_t nparticles) {
  for (unsigned int iSelection=0;iSelection < ParticleCuts.size() ;iSelection++) {
    std::vector<double> Vec = ReturnKinematicVector(ParticleCuts[iSelection].ParamToCutOnIt, iEvent);
    if (nparticles != Vec.size()) {
      MACH3LOG_ERROR("Cannot apply kinematic cut on {} as it is of different size to plotting variable");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const double Val = Vec[iParticle];
    if ((Val < ParticleCuts[iSelection].LowerBound) || (Val >= ParticleCuts[iSelection].UpperBound)) {
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
/// @function SampleHandlerFD::fillArray()
/// Function which does the core reweighting. This assumes that oscillation weights have 
/// already been calculated and stored in SampleHandlerFD.osc_w[iEvent]. This
/// function takes advantage of most of the things called in setupSKMC to reduce reweighting time.
/// It also follows the ND code reweighting pretty closely. This function fills the SampleHandlerFD 
/// array array which is binned to match the sample binning, such that bin[1][1] is the 
/// equivalent of _hPDF2D->GetBinContent(2,2) {Noticing the offset}
void SampleHandlerFD::FillArray() {
//************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;
  
  PrepFunctionalParameters();

  for (unsigned int iEvent=0;iEvent<GetNEvents();iEvent++) {
    ApplyShifts(iEvent);

    if (!IsEventSelected(iEvent)) {
      continue;
    }

    M3::float_t splineweight = 1.0;
    M3::float_t normweight = 1.0;
    M3::float_t totalweight = 1.0;

    if(SplineHandler){
      splineweight = CalcWeightSpline(iEvent);
    }
    //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient. Do this on a spline-by-spline basis
    if (splineweight <= 0.){
      MCSamples.xsec_w[iEvent] = 0.;
      continue;
    }

    //Loop over stored normalisation and function pointers
    normweight = CalcWeightNorm(iEvent);

    //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
    if (normweight <= 0.){
      MCSamples.xsec_w[iEvent] = 0.;
      continue;
    }
    // Virtual by default does nothing
    CalcWeightFunc(iEvent);

    MCSamples.xsec_w[iEvent] = splineweight*normweight;

    //DB Total weight
    totalweight = GetEventWeight(iEvent);
    //DB Catch negative weights and skip any event with a negative event
    if (totalweight <= 0.){
      MCSamples.xsec_w[iEvent] = 0.;
      continue;
    }
    //DB Switch on BinningOpt to allow different binning options to be implemented
    //The alternative would be to have inheritance based on BinningOpt
    double XVar = *(MCSamples.x_var[iEvent]);

    //DB Find the relevant bin in the PDF for each event
    int XBinToFill = -1;
    int YBinToFill = MCSamples.NomYBin[iEvent];

    //DB - First, check to see if the event is still in the nominal bin
    if (XVar < MCSamples.rw_upper_xbinedge[iEvent] && XVar >= MCSamples.rw_lower_xbinedge[iEvent]) {
      XBinToFill = MCSamples.NomXBin[iEvent];
    }
    //DB - Second, check to see if the event is outside of the binning range and skip event if it is
    //ETA- note that nXBins is XBinEdges.size() - 1
    else if (XVar < XBinEdges[0] || XVar >= XBinEdges[nXBins]) {
      continue;
    }
    //DB - Thirdly, check the adjacent bins first as Eb+CC+EScale shifts aren't likely to move an Erec more than 1bin width
    //Shifted down one bin from the event bin at nominal
    else if (XVar < MCSamples.rw_lower_xbinedge[iEvent] && XVar >= MCSamples.rw_lower_lower_xbinedge[iEvent]) {
      XBinToFill = MCSamples.NomXBin[iEvent]-1;
    }
    //Shifted up one bin from the event bin at nominal
    else if (XVar < MCSamples.rw_upper_upper_xbinedge[iEvent] && XVar >= MCSamples.rw_upper_xbinedge[iEvent]) {
      XBinToFill = MCSamples.NomXBin[iEvent]+1;
    }
    //DB - If we end up in this loop, the event has been shifted outside of its nominal bin, but is still within the allowed binning range
    else {
      for (unsigned int iBin=0;iBin<(XBinEdges.size()-1);iBin++)
      {
        if (XVar >= XBinEdges[iBin] && XVar < XBinEdges[iBin+1]) {
          XBinToFill = iBin;
          break; // KS: Once found bin stop looping over remaining
        }
      }
    }

    //DB Fill relevant part of thread array
    if (XBinToFill != -1 && YBinToFill != -1) {
      SampleHandlerFD_array[YBinToFill][XBinToFill] += totalweight;
      if (FirstTimeW2) SampleHandlerFD_array_w2[YBinToFill][XBinToFill] += totalweight*totalweight;
    }
  }
}

#ifdef MULTITHREAD
// ************************************************ 
/// Multithreaded version of fillArray @see fillArray()
void SampleHandlerFD::FillArray_MP()  {
// ************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;

  PrepFunctionalParameters();

  //This is stored as [y][x] due to shifts only occurring in the x variable (Erec/Lep mom) - I believe this will help reduce cache misses
  double** SampleHandlerFD_array_private = nullptr;
  double** SampleHandlerFD_array_private_w2 = nullptr;
  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
#pragma omp parallel private(SampleHandlerFD_array_private, SampleHandlerFD_array_private_w2)
  {
    // private to each thread
    // ETA - maybe we can use parallel firstprivate to initialise these?
    SampleHandlerFD_array_private = new double*[nYBins];
    SampleHandlerFD_array_private_w2 = new double*[nYBins];
    for (size_t yBin=0;yBin<nYBins;yBin++) {
      SampleHandlerFD_array_private[yBin] = new double[nXBins];
      SampleHandlerFD_array_private_w2[yBin] = new double[nXBins];

      std::fill_n(SampleHandlerFD_array_private[yBin], nXBins, 0.0);
      std::fill_n(SampleHandlerFD_array_private_w2[yBin], nXBins, 0.0);
    }
    
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
    
    #pragma omp for
    for (unsigned int iEvent=0;iEvent<GetNEvents();iEvent++) {
      //ETA - generic functions to apply shifts to kinematic variables
      // Apply this before IsEventSelected is called.
      ApplyShifts(iEvent);

      //ETA - generic functions to apply shifts to kinematic variable
      //this is going to be slow right now due to string comps under the hood.
      //Need to implement a more efficient version of event-by-event cut checks
      if(!IsEventSelected(iEvent)){
        continue;
      }

      M3::float_t splineweight = 1.0;
      M3::float_t normweight = 1.0;
      M3::float_t totalweight = 1.0;

      //DB SKDet Syst
      //As weights were skdet::fParProp, and we use the non-shifted erec, we might as well cache the corresponding fParProp index for each event and the pointer to it

      if(SplineHandler){
        splineweight = CalcWeightSpline(iEvent);
      }
      //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
      if (splineweight <= 0.){
        MCSamples.xsec_w[iEvent] = 0.;
        continue;
      }

      normweight = CalcWeightNorm(iEvent);
      //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
      if (normweight <= 0.){
        MCSamples.xsec_w[iEvent] = 0.;
        continue;
      }

      // Virtual by default does nothing
      CalcWeightFunc(iEvent);

      MCSamples.xsec_w[iEvent] = splineweight*normweight;

      totalweight = GetEventWeight(iEvent);

      //DB Catch negative weights and skip any event with a negative event
      if (totalweight <= 0.){
        MCSamples.xsec_w[iEvent] = 0.;
        continue;
      }

      //DB Switch on BinningOpt to allow different binning options to be implemented
      //The alternative would be to have inheritance based on BinningOpt
      double XVar = (*(MCSamples.x_var[iEvent]));

      //DB Commented out by default but if we ever want to consider shifts in theta this will be needed
      //double YVar = MCSamples.rw_theta[iEvent];
      //ETA - this would actually be with (*(MCSamples.y_var[iEvent])) and done extremely
      //similarly to XVar now

      //DB Find the relevant bin in the PDF for each event
      int XBinToFill = -1;
      int YBinToFill = MCSamples.NomYBin[iEvent];

      //DB Check to see if momentum shift has moved bins
      //DB - First, check to see if the event is still in the nominal bin
      if (XVar < MCSamples.rw_upper_xbinedge[iEvent] && XVar >= MCSamples.rw_lower_xbinedge[iEvent]) {
        XBinToFill = MCSamples.NomXBin[iEvent];
      }
      //DB - Second, check to see if the event is outside of the binning range and skip event if it is
      else if (XVar < XBinEdges[0] || XVar >= XBinEdges[nXBins]) {
        continue;
      }
      //DB - Thirdly, check the adjacent bins first as Eb+CC+EScale shifts aren't likely to move an Erec more than 1bin width
      //Shifted down one bin from the event bin at nominal
      else if (XVar < MCSamples.rw_lower_xbinedge[iEvent] && XVar >= MCSamples.rw_lower_lower_xbinedge[iEvent]) {
        XBinToFill = MCSamples.NomXBin[iEvent]-1;
      }
      //Shifted up one bin from the event bin at nominal
      else if (XVar < MCSamples.rw_upper_upper_xbinedge[iEvent] && XVar >= MCSamples.rw_upper_xbinedge[iEvent]) {
        XBinToFill = MCSamples.NomXBin[iEvent]+1;
      }
      //DB - If we end up in this loop, the event has been shifted outside of its nominal bin, but is still within the allowed binning range
      else {
        for(unsigned int iBin=0;iBin<(XBinEdges.size()-1);iBin++) {
          if (XVar >= XBinEdges[iBin] && XVar < XBinEdges[iBin+1]) {
            XBinToFill = iBin;
          }
        }
      }

      //ETA - we can probably remove this final if check on the -1?
      //Maybe we can add an overflow bin to the array and assign any events to this bin?
      //Might save us an extra if call?
      //DB Fill relevant part of thread array
      if (XBinToFill != -1 && YBinToFill != -1) {
        SampleHandlerFD_array_private[YBinToFill][XBinToFill] += totalweight;
        SampleHandlerFD_array_private_w2[YBinToFill][XBinToFill] += totalweight*totalweight;
      }
    }
    //End of Calc Weights and fill Array
    //==================================================
    // DB Copy contents of 'SampleHandlerFD_array_private' into 'SampleHandlerFD_array' which can then be used in GetLikelihood
    for (size_t yBin = 0; yBin < nYBins; ++yBin) {
      for (size_t xBin = 0; xBin < nXBins; ++xBin) {
        #pragma omp atomic
        SampleHandlerFD_array[yBin][xBin] += SampleHandlerFD_array_private[yBin][xBin];
        if(FirstTimeW2) {
          #pragma omp atomic
          SampleHandlerFD_array_w2[yBin][xBin] += SampleHandlerFD_array_private_w2[yBin][xBin];
        }
      }
    }
    
    for (size_t yBin = 0; yBin < nYBins; ++yBin) {
      delete[] SampleHandlerFD_array_private[yBin];
      delete[] SampleHandlerFD_array_private_w2[yBin];
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
  for (size_t yBin = 0; yBin < nYBins; ++yBin) {
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (size_t xBin = 0; xBin < nXBins; ++xBin) {
      SampleHandlerFD_array[yBin][xBin] = 0.;
      SampleHandlerFD_array_w2[yBin][xBin] = 0.;
    }
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
  funcParsVec = ParHandler->GetFunctionalParametersFromSampleName(SampleName);
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
      bool ModeMatch = MatchCondition((*it).modes, static_cast<int>(std::round(*(MCSamples.mode[iEvent]))));
      if (!ModeMatch) {
        MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, *(MCSamples.mode[iEvent]), (*it).name);
        continue;
      }
      // Now check whether within kinematic bounds
      bool IsSelected = true;
      if ((*it).hasKinBounds) {
        const auto& kinVars = (*it).KinematicVarStr;
        const auto& selection = (*it).Selection;

        for (std::size_t iKinPar = 0; iKinPar < kinVars.size(); ++iKinPar) {
          const double kinVal = ReturnKinematicParameter(kinVars[iKinPar], static_cast<int>(iEvent));

          if (kinVal <= selection[iKinPar][0] || kinVal > selection[iKinPar][1]) {
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
  for (int fpEnum : funcParsGrid[iEvent]) {
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
M3::float_t SampleHandlerFD::CalcWeightSpline(const int iEvent) const {
// ***************************************************************************
  M3::float_t spline_weight = 1.0;
  //DB Xsec syst
  //Loop over stored spline pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iSpline = 0; iSpline < MCSamples.nxsec_spline_pointers[iEvent]; ++iSpline) {
    spline_weight *= *(MCSamples.xsec_spline_pointers[iEvent][iSpline]);
  }
  return spline_weight;
}

// ***************************************************************************
// Calculate the normalisation weight for an event
M3::float_t SampleHandlerFD::CalcWeightNorm(const int iEvent) const {
// ***************************************************************************
  M3::float_t xsecw = 1.0;
  //Loop over stored normalisation and function pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iParam = 0; iParam < MCSamples.nxsec_norm_pointers[iEvent]; ++iParam)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
    xsecw *= static_cast<M3::float_t>(*(MCSamples.xsec_norm_pointers[iEvent][iParam]));
#pragma GCC diagnostic pop
    #ifdef DEBUG
    if (TMath::IsNaN(xsecw)) MACH3LOG_WARN("iParam= {} xsecweight=nan from norms", iParam);
    #endif
  }
  return xsecw;
}


// ***************************************************************************
// Setup the norm parameters
void SampleHandlerFD::SetupNormParameters() {  
// ***************************************************************************
  norm_parameters = ParHandler->GetNormParsFromSampleName(GetSampleName());

  if(!ParHandler){
    MACH3LOG_ERROR("ParHandler is not setup!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Assign xsec norm bins in MCSamples tree
  CalcNormsBins();

  //DB Attempt at reducing impact of SystematicHandlerGeneric::calcReweight()
  int counter;
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    counter = 0;

    MCSamples.nxsec_norm_pointers[iEvent] = int(MCSamples.xsec_norms_bins[iEvent].size());
    MCSamples.xsec_norm_pointers[iEvent].resize(MCSamples.nxsec_norm_pointers[iEvent]);

    for(auto const & norm_bin: MCSamples.xsec_norms_bins[iEvent]) {
      MCSamples.xsec_norm_pointers[iEvent][counter] = ParHandler->RetPointer(norm_bin);
      counter += 1;
    }
  }
  // If not debugging let's clear memory
  #ifndef DEBUG
  MCSamples.xsec_norms_bins.clear();
  MCSamples.xsec_norms_bins.shrink_to_fit();
  #endif
}

// ************************************************
//A way to check whether a normalisation parameter applies to an event or not
void SampleHandlerFD::CalcNormsBins() {
// ************************************************
  FarDetectorCoreInfo *fdobj = &MCSamples;
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
      if (fdobj->isNC[iEvent] && (*fdobj->nupdg[iEvent] != *fdobj->nupdgUnosc[iEvent]) ) {
        MACH3LOG_TRACE("Event {}, missed NC/signal check", iEvent);
        continue;
      } //DB Abstract check on MaCh3Modes to determine which apply to neutral current
      for (std::vector<NormParameter>::iterator it = norm_parameters.begin(); it != norm_parameters.end(); ++it) {
        //Now check that the target of an interaction matches with the normalisation parameters
        bool TargetMatch = MatchCondition((*it).targets, *(fdobj->Target[iEvent]));
        if (!TargetMatch) {
          MACH3LOG_TRACE("Event {}, missed target check ({}) for dial {}", iEvent, *(fdobj->Target[iEvent]), (*it).name);
          continue;
        }

        //Now check that the neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourMatch = MatchCondition((*it).pdgs, *(fdobj->nupdg[iEvent]));
        if (!FlavourMatch) {
          MACH3LOG_TRACE("Event {}, missed PDG check ({}) for dial {}", iEvent,(*fdobj->nupdg[iEvent]), (*it).name);
          continue;
        }

        //Now check that the unoscillated neutrino flavour in an interaction matches with the normalisation parameters
        bool FlavourUnoscMatch = MatchCondition((*it).preoscpdgs, *(fdobj->nupdgUnosc[iEvent]));
        if (!FlavourUnoscMatch){
          MACH3LOG_TRACE("Event {}, missed FlavourUnosc check ({}) for dial {}", iEvent,(*fdobj->nupdgUnosc[iEvent]), (*it).name);
          continue;
        }

        //Now check that the mode of an interaction matches with the normalisation parameters
        bool ModeMatch = MatchCondition((*it).modes, static_cast<int>(std::round(*(fdobj->mode[iEvent]))));
        if (!ModeMatch) {
          MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, *(fdobj->mode[iEvent]), (*it).name);
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
            const double kinVal = ReturnKinematicParameter(kinVars[iKinPar], iEvent);

            if (kinVal <= selection[iKinPar][0] || kinVal > selection[iKinPar][1]) {
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
    fdobj->xsec_norms_bins[iEvent] = NormBins;
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
void SampleHandlerFD::SetupReweightArrays(const size_t numberXBins, const size_t numberYBins) {
// ************************************************
  //Set the number of X and Y bins now
  nXBins = numberXBins;
  nYBins = numberYBins;

  SampleHandlerFD_array = new double*[nYBins];
  SampleHandlerFD_array_w2 = new double*[nYBins];
  SampleHandlerFD_data = new double*[nYBins];
  for (size_t yBin=0;yBin<nYBins;yBin++) {
    SampleHandlerFD_array[yBin] = new double[nXBins];
    SampleHandlerFD_array_w2[yBin] = new double[nXBins];
    SampleHandlerFD_data[yBin] = new double[nXBins];
    for (size_t xBin=0;xBin<nXBins;xBin++) {
      SampleHandlerFD_array[yBin][xBin] = 0.;
      SampleHandlerFD_array_w2[yBin][xBin] = 0.;
      SampleHandlerFD_data[yBin][xBin] = 0.;
    }
  }
}

//ETA
//New versions of set binning functions is SampleHandlerBase
//so that we can set the values of the bin and lower/upper
//edges in the skmc_base. Hopefully we can use this to make
//fill1Dhist and fill2Dhist quicker
void SampleHandlerFD::Set1DBinning(size_t nbins, double* boundaries)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(static_cast<int>(nbins),boundaries);
  dathist->SetBins(static_cast<int>(nbins),boundaries);

  XBinEdges = std::vector<double>(nbins+1);
  for (size_t i=0;i<nbins+1;i++) {
    XBinEdges[i] = _hPDF1D->GetXaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  double YBinEdges_Arr[2];
  YBinEdges_Arr[0] = YBinEdges[0];
  YBinEdges_Arr[1] = YBinEdges[1];

  _hPDF2D->Reset();
  _hPDF2D->SetBins(static_cast<int>(nbins),boundaries,1,YBinEdges_Arr);
  dathist2d->SetBins(static_cast<int>(nbins),boundaries,1,YBinEdges_Arr);

  //Set the number of X and Y bins now
  SetupReweightArrays(XBinEdges.size() - 1, YBinEdges.size() - 1);

  FindNominalBinAndEdges1D();
}

void SampleHandlerFD::Set1DBinning(size_t nbins, double low, double high)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(static_cast<int>(nbins),low,high);
  dathist->SetBins(static_cast<int>(nbins),low,high);

  XBinEdges = std::vector<double>(nbins+1);
  for (size_t i=0;i<nbins+1;i++) {
    XBinEdges[i] = _hPDF1D->GetXaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  _hPDF2D->Reset();
  _hPDF2D->SetBins(static_cast<int>(nbins),low,high,1,YBinEdges[0],YBinEdges[1]);
  dathist2d->SetBins(static_cast<int>(nbins),low,high,1,YBinEdges[0],YBinEdges[1]);

  //Set the number of X and Y bins now
  SetupReweightArrays(XBinEdges.size() - 1, YBinEdges.size() - 1);

  FindNominalBinAndEdges1D();
}

void SampleHandlerFD::FindNominalBinAndEdges1D() {
  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(unsigned int event_i = 0; event_i < GetNEvents(); event_i++){
    //Set x_var and y_var values based on XVarStr and YVarStr
    MCSamples.x_var[event_i] = GetPointerToKinematicParameter(XVarStr, event_i);
    //Give y_var M3::_BAD_DOUBLE_ value for the 1D case since this won't be used
    MCSamples.y_var[event_i] = &(M3::_BAD_DOUBLE_);
    int bin = _hPDF1D->FindBin(*(MCSamples.x_var[event_i]));

    double low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
    if (bin == 0) {
      low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
    } else {
      low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin-1);
    }

    double low_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
    double upper_edge = _hPDF1D->GetXaxis()->GetBinUpEdge(bin);

    double upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
    if (bin<(_hPDF1D->GetNbinsX()-2)) {
      upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+2);
    } else {
      upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+1);
    }

    if ((bin-1) >= 0 && (bin-1) < int(XBinEdges.size()-1)) {
      MCSamples.NomXBin[event_i] = bin-1;
    } else {
      MCSamples.NomXBin[event_i] = -1;
      low_edge = M3::_DEFAULT_RETURN_VAL_;
      upper_edge = M3::_DEFAULT_RETURN_VAL_;
      low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
      upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
    }

    MCSamples.NomYBin[event_i] = 0;
    MCSamples.rw_lower_xbinedge[event_i] = low_edge;
    MCSamples.rw_upper_xbinedge[event_i] = upper_edge;
    MCSamples.rw_lower_lower_xbinedge[event_i] = low_lower_edge;
    MCSamples.rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
  }
}

void SampleHandlerFD::Set2DBinning(size_t nbins1, double* boundaries1, size_t nbins2, double* boundaries2)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(static_cast<int>(nbins1),boundaries1);
  dathist->SetBins(static_cast<int>(nbins1),boundaries1);

  _hPDF2D->Reset();
  _hPDF2D->SetBins(static_cast<int>(nbins1),boundaries1,static_cast<int>(nbins2),boundaries2);
  dathist2d->SetBins(static_cast<int>(nbins1),boundaries1,static_cast<int>(nbins2),boundaries2);

  XBinEdges = std::vector<double>(nbins1+1);
  for (size_t i=0;i<nbins1+1;i++) {
    XBinEdges[i] = _hPDF2D->GetXaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }
  YBinEdges = std::vector<double>(nbins2+1);
  for (size_t i=0;i<nbins2+1;i++) {
    YBinEdges[i] = _hPDF2D->GetYaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }
  
  //Set the number of X and Y bins now
  SetupReweightArrays(XBinEdges.size() - 1, YBinEdges.size() - 1);

  FindNominalBinAndEdges2D();
}

void SampleHandlerFD::Set2DBinning(size_t nbins1, double low1, double high1, size_t nbins2, double low2, double high2)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(static_cast<int>(nbins1),low1,high1);
  dathist->SetBins(static_cast<int>(nbins1),low1,high1);

  _hPDF2D->Reset();
  _hPDF2D->SetBins(static_cast<int>(nbins1),low1,high1,static_cast<int>(nbins2),low2,high2);
  dathist2d->SetBins(static_cast<int>(nbins1),low1,high1,static_cast<int>(nbins2),low2,high2);

  XBinEdges = std::vector<double>(nbins1+1);
  for (size_t i=0;i<nbins1+1;i++) {
    XBinEdges[i] = _hPDF2D->GetXaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }
  YBinEdges = std::vector<double>(nbins2+1);
  for (size_t i=0;i<nbins2+1;i++) {
    YBinEdges[i] = _hPDF2D->GetYaxis()->GetBinLowEdge(static_cast<int>(i+1));
  }

  //Set the number of X and Y bins now
  SetupReweightArrays(XBinEdges.size() - 1, YBinEdges.size() - 1);

  FindNominalBinAndEdges2D();
}

// ************************************************
void SampleHandlerFD::FindNominalBinAndEdges2D() {
// ************************************************
  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(unsigned int event_i = 0 ; event_i < GetNEvents(); event_i++){
    //Set x_var and y_var values based on XVarStr and YVarStr
    MCSamples.x_var[event_i] = GetPointerToKinematicParameter(XVarStr, event_i);
    MCSamples.y_var[event_i] = GetPointerToKinematicParameter(YVarStr, event_i);

    //Global bin number
    int bin = _hPDF2D->FindBin(*(MCSamples.x_var[event_i]), *(MCSamples.y_var[event_i]));

    int bin_x = M3::_BAD_INT_;
    int bin_y = M3::_BAD_INT_;
    int bin_z = M3::_BAD_INT_;
    _hPDF2D->GetBinXYZ(bin, bin_x, bin_y, bin_z);
    //erec is the x-axis so get GetXaxis then find the bin edges using the x bin number

    double low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
    if (bin==0) {
      low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
    } else {
      low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x-1);
    }

    double low_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
    double upper_edge = _hPDF2D->GetXaxis()->GetBinUpEdge(bin_x);

    double upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
    if (bin<(_hPDF2D->GetNbinsX()-2)) {
      upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+2);
    } else {
      upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+1);
    }

    if ((bin_x-1) >= 0 && (bin_x-1) < int(XBinEdges.size()-1)) {
      MCSamples.NomXBin[event_i] = bin_x-1;
    } else {
      MCSamples.NomXBin[event_i] = -1;
      low_edge = M3::_DEFAULT_RETURN_VAL_;
      upper_edge = M3::_DEFAULT_RETURN_VAL_;
      low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
      upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
    }
    MCSamples.NomYBin[event_i] = bin_y-1;
    if(MCSamples.NomYBin[event_i] < 0){
      MACH3LOG_INFO("Nominal YBin PROBLEM, y-bin is {}", MCSamples.NomYBin[event_i]);
    }
    MCSamples.rw_lower_xbinedge[event_i] = low_edge;
    MCSamples.rw_upper_xbinedge[event_i] = upper_edge;
    MCSamples.rw_lower_lower_xbinedge[event_i] = low_lower_edge;
    MCSamples.rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
  }
}

void SampleHandlerFD::AddData(std::vector<double> &data) {
  dathist2d = nullptr;
  dathist->Reset(); 
  
  if (GetNDim()!=1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  for (auto const& data_point : data){
    dathist->Fill(data_point);
  }

  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (size_t yBin=0;yBin<nYBins;yBin++) {
    for (size_t xBin=0;xBin<nXBins;xBin++) {
      SampleHandlerFD_data[yBin][xBin] = dathist->GetBinContent(static_cast<int>(xBin+1));
    }
  }
}

void SampleHandlerFD::AddData(std::vector< std::vector <double> > &data) {
  dathist = nullptr;
  dathist2d->Reset();                                                       

  if (GetNDim()!=2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //TODO: this assumes that std::vector is of length 2 and then both data.at(0) and data
  //ETA: I think this might just be wrong? We should probably just make this AddData(std::vector<double> data_x, std::vector<double> data_y)
  // or maybe something like AddData(std::vector<std::pair<double, double>> data)?
  for (int i = 0; i < int(data.size()); i++) {
    dathist2d->Fill(data.at(0)[i],data.at(1)[i]);
  }

  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (size_t yBin=0;yBin<nYBins;yBin++) {
    for (size_t xBin=0;xBin<nXBins;xBin++) {
      //Need to cast to an int (Int_t) for ROOT
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      SampleHandlerFD_data[yBin][xBin] = dathist2d->GetBinContent(static_cast<int>(xBin+1),static_cast<int>(yBin+1));
    }
  }
}

void SampleHandlerFD::AddData(TH1D* Data) {
  MACH3LOG_INFO("Adding 1D data histogram: {} with {:.2f} events", Data->GetTitle(), Data->Integral());
  dathist2d = nullptr;
  dathist = Data;
  
  if (GetNDim()!=1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram in a 2D sample - Quitting"); 
    MACH3LOG_ERROR("The number of dimensions for this sample is {}", GetNDim());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
    
  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (size_t yBin=0;yBin<nYBins;yBin++) {
    for (size_t xBin=0;xBin<nXBins;xBin++) {
      //Need to cast to an int (Int_t) for ROOT
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      SampleHandlerFD_data[yBin][xBin] = Data->GetBinContent(static_cast<int>(xBin+1));
    }
  }
}

void SampleHandlerFD::AddData(TH2D* Data) {
  MACH3LOG_INFO("Adding 2D data histogram: {} with {:.2f} events", Data->GetTitle(), Data->Integral());
  dathist2d = Data;
  dathist = nullptr;

  if (GetNDim()!=2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram in a 1D sample - Quitting"); 
    throw MaCh3Exception(__FILE__ , __LINE__ );}	
   
  if(SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("SampleHandlerFD_data haven't been initialised yet");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  for (size_t yBin=0;yBin<nYBins;yBin++) {
    for (size_t xBin=0;xBin<nXBins;xBin++) {
      //Need to cast to an int (Int_t) for ROOT
      //Need to do +1 for the bin, this is to be consistent with ROOTs binning scheme
      SampleHandlerFD_data[yBin][xBin] = dathist2d->GetBinContent(static_cast<int>(xBin+1),static_cast<int>(yBin+1));
    }
  }
}

// ************************************************
void SampleHandlerFD::InitialiseNuOscillatorObjects() {
// ************************************************
  if (!OscParHandler) {
    MACH3LOG_ERROR("Attempted to setup NuOscillator without ParameterHandler object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto NuOscillatorConfigFile = Get<std::string>(SampleManager->raw()["NuOsc"]["NuOscConfigFile"], __FILE__ , __LINE__);
  auto EqualBinningPerOscChannel = Get<bool>(SampleManager->raw()["NuOsc"]["EqualBinningPerOscChannel"], __FILE__ , __LINE__);

  // TN override the sample setting if not using binned oscillation
  if (EqualBinningPerOscChannel) {
    if (YAML::LoadFile(NuOscillatorConfigFile)["General"]["CalculationType"].as<std::string>() == "Unbinned") {
      MACH3LOG_WARN("Tried using EqualBinningPerOscChannel while using Unbinned oscillation calculation, changing EqualBinningPerOscChannel to false");
      EqualBinningPerOscChannel = false;
    }
  }
  std::vector<const double*> OscParams = OscParHandler->GetOscParsFromSampleName(SampleName);
  if (OscParams.empty()) {
    MACH3LOG_ERROR("OscParams is empty for sample '{}'.", GetTitle());
    MACH3LOG_ERROR("This likely indicates an error in your oscillation YAML configuration.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  Oscillator = std::make_shared<OscillationHandler>(NuOscillatorConfigFile, EqualBinningPerOscChannel, OscParams, static_cast<int>(OscChannels.size()));

  if (!EqualBinningPerOscChannel) {
    for(int iChannel = 0; iChannel < static_cast<int>(OscChannels.size()); iChannel++) {
      std::vector<M3::float_t> EnergyArray;
      std::vector<M3::float_t> CosineZArray;

      for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
        // KS: This is bit weird but we basically loop over all events and push to vector only these which are part of a given OscChannel
        const int Channel = GetOscChannel(OscChannels, (*MCSamples.nupdgUnosc[iEvent]), (*MCSamples.nupdg[iEvent]));
        //DB Remove NC events from the arrays which are handed to the NuOscillator objects
        if (!MCSamples.isNC[iEvent] && Channel == iChannel) {
          EnergyArray.push_back(M3::float_t(*(MCSamples.rw_etru[iEvent])));
        }
      }
      std::sort(EnergyArray.begin(),EnergyArray.end());

      //============================================================================
      //DB Atmospheric only part, can only happen if truecz has been initialised within the experiment specific code
      if (MCSamples.rw_truecz.size() > 0 && MCSamples.rw_truecz.size() == GetNEvents()) {
        for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
          // KS: This is bit weird but we basically loop over all events and push to vector only these which are part of a given OscChannel
          const int Channel = GetOscChannel(OscChannels, (*MCSamples.nupdgUnosc[iEvent]), (*MCSamples.nupdg[iEvent]));
          //DB Remove NC events from the arrays which are handed to the NuOscillator objects
          if (!MCSamples.isNC[iEvent] && Channel == iChannel) {
            CosineZArray.push_back(M3::float_t(*(MCSamples.rw_truecz[iEvent])));
          }
        }
        std::sort(CosineZArray.begin(),CosineZArray.end());
      }
      Oscillator->SetOscillatorBinning(iChannel, EnergyArray, CosineZArray);
    } // end loop over channels
  }
}

void SampleHandlerFD::SetupNuOscillatorPointers() {
  for (unsigned int iEvent=0;iEvent<GetNEvents();iEvent++) {
    MCSamples.osc_w_pointer[iEvent] = &M3::Unity;
    if (MCSamples.isNC[iEvent]) {
      if (*MCSamples.nupdg[iEvent] != *MCSamples.nupdgUnosc[iEvent]) {
        MCSamples.osc_w_pointer[iEvent] = &M3::Zero;
      } else {
        MCSamples.osc_w_pointer[iEvent] = &M3::Unity;
      }
    } else {
      int InitFlav = M3::_BAD_INT_;
      int FinalFlav = M3::_BAD_INT_;

      InitFlav =  MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples.nupdgUnosc[iEvent]));
      FinalFlav = MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples.nupdg[iEvent]));

      if (InitFlav == M3::_BAD_INT_ || FinalFlav == M3::_BAD_INT_) {
        MACH3LOG_ERROR("Something has gone wrong in the mapping between MCSamples.nutype and the enum used within NuOscillator");
        MACH3LOG_ERROR("MCSamples.nupdgUnosc: {}", (*MCSamples.nupdgUnosc[iEvent]));
        MACH3LOG_ERROR("InitFlav: {}", InitFlav);
        MACH3LOG_ERROR("MCSamples.nupdg: {}", (*MCSamples.nupdg[iEvent]));
        MACH3LOG_ERROR("FinalFlav: {}", FinalFlav);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      const int OscIndex = GetOscChannel(OscChannels, (*MCSamples.nupdgUnosc[iEvent]), (*MCSamples.nupdg[iEvent]));
      if (MCSamples.rw_truecz.size() > 0) { //Can only happen if truecz has been initialised within the experiment specific code
        //Atmospherics
        MCSamples.osc_w_pointer[iEvent] = Oscillator->GetNuOscillatorPointers(OscIndex, InitFlav, FinalFlav, FLOAT_T(*(MCSamples.rw_etru[iEvent])), FLOAT_T(*(MCSamples.rw_truecz[iEvent])));
      } else {
        //Beam
        MCSamples.osc_w_pointer[iEvent] = Oscillator->GetNuOscillatorPointers(OscIndex, InitFlav, FinalFlav, FLOAT_T(*(MCSamples.rw_etru[iEvent])));
      }
    } // end if NC
  } // end loop over events
}

std::string SampleHandlerFD::GetSampleName(int iSample) const {
  //ETA - this is just to suppress a warning for an unused variable
  (void)iSample;

  //ETA - extra safety to make sure SampleName is actually set
  // probably unnecessary due to the requirement for it to be in the yaml config
  if(SampleName.length() == 0){
    MACH3LOG_ERROR("No sample name provided");
    MACH3LOG_ERROR("Please provide a SampleName in your configuration file: {}", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return SampleName;
}

M3::float_t SampleHandlerFD::GetEventWeight(const int iEntry) const {
  M3::float_t totalweight = 1.0;
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iParam = 0; iParam < MCSamples.ntotal_weight_pointers[iEntry]; ++iParam) {
    totalweight *= *(MCSamples.total_weight_pointers[iEntry][iParam]);
  }
  
  return totalweight;
}

/// @func fillSplineBins()
/// @brief Finds the binned spline that an event should apply to and stored them in a
/// a vector for easy evaluation in the fillArray() function.
void SampleHandlerFD::FillSplineBins() {
  //Now loop over events and get the spline bin for each event
  for (unsigned int j = 0; j < GetNEvents(); ++j) {
    const int OscIndex = GetOscChannel(OscChannels, (*MCSamples.nupdgUnosc[j]), (*MCSamples.nupdg[j]));

    std::vector< std::vector<int> > EventSplines;
    switch(nDimensions){
      case 1:
        EventSplines = SplineHandler->GetEventSplines(GetSampleName(), OscIndex, int(*(MCSamples.mode[j])), *(MCSamples.rw_etru[j]), *(MCSamples.x_var[j]), 0.);
        break;
      case 2:
        EventSplines = SplineHandler->GetEventSplines(GetSampleName(), OscIndex, int(*(MCSamples.mode[j])), *(MCSamples.rw_etru[j]), *(MCSamples.x_var[j]), *(MCSamples.y_var[j]));
        break;
      default:
        MACH3LOG_ERROR("Error in assigning spline bins because nDimensions = {}", nDimensions);
        MACH3LOG_ERROR("MaCh3 only supports splines binned in Etrue + the sample binning");
        MACH3LOG_ERROR("Please check the sample binning you specified in your sample config ");
        throw MaCh3Exception(__FILE__, __LINE__);
        break;
    }
    MCSamples.nxsec_spline_pointers[j] = int(EventSplines.size());
    if(MCSamples.nxsec_spline_pointers[j] < 0){
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    MCSamples.xsec_spline_pointers[j].resize(MCSamples.nxsec_spline_pointers[j]);
    for(int spline = 0; spline < MCSamples.nxsec_spline_pointers[j]; spline++) {
      //Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
      MCSamples.xsec_spline_pointers[j][spline] = SplineHandler->retPointer(EventSplines[spline][0], EventSplines[spline][1],
                                                                                EventSplines[spline][2], EventSplines[spline][3],
                                                                                EventSplines[spline][4], EventSplines[spline][5],
                                                                                EventSplines[spline][6]);
    }
  }
}

// ************************************************
double SampleHandlerFD::GetLikelihood() {
// ************************************************
  if (SampleHandlerFD_data == nullptr) {
    MACH3LOG_ERROR("Data sample is empty! Can't calculate a likelihood!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
    
  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2) reduction(+:negLogL)
  #endif
  for (size_t yBin = 0; yBin < nYBins; ++yBin)
  {
    for (size_t xBin = 0; xBin < nXBins; ++xBin)
    {
      const double DataVal = SampleHandlerFD_data[yBin][xBin];
      const double MCPred = SampleHandlerFD_array[yBin][xBin];
      const double w2 = SampleHandlerFD_array_w2[yBin][xBin];
      
      //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
      negLogL += GetTestStatLLH(DataVal, MCPred, w2);
    }
  }
  return negLogL;
}

void SampleHandlerFD::InitialiseSingleFDMCObject() {
  FarDetectorCoreInfo *fdobj = &MCSamples;
  
  fdobj->x_var.resize(nEvents, &M3::Unity_D);
  fdobj->y_var.resize(nEvents, &M3::Unity_D);
  fdobj->rw_etru.resize(nEvents, &M3::Unity_D);
  fdobj->XBin.resize(nEvents, -1);
  fdobj->YBin.resize(nEvents, -1);
  fdobj->NomXBin.resize(nEvents, -1);
  fdobj->NomYBin.resize(nEvents, -1);
  fdobj->rw_lower_xbinedge.resize(nEvents, -1);
  fdobj->rw_lower_lower_xbinedge.resize(nEvents, -1);
  fdobj->rw_upper_xbinedge.resize(nEvents, -1);
  fdobj->rw_upper_upper_xbinedge.resize(nEvents, -1);
  fdobj->mode.resize(nEvents, &M3::Unity_D);
  fdobj->nxsec_norm_pointers.resize(nEvents);
  fdobj->xsec_norm_pointers.resize(nEvents);
  fdobj->xsec_norms_bins.resize(nEvents);
  fdobj->xsec_w.resize(nEvents, 1.0);
  fdobj->nupdg.resize(nEvents);
  fdobj->nupdgUnosc.resize(nEvents);
  fdobj->isNC = new bool[nEvents];
  fdobj->nxsec_spline_pointers.resize(nEvents);
  fdobj->xsec_spline_pointers.resize(nEvents);
  fdobj->ntotal_weight_pointers.resize(nEvents);
  fdobj->total_weight_pointers.resize(nEvents);
  fdobj->Target.resize(nEvents, 0);
  fdobj->osc_w_pointer.resize(nEvents, &M3::Unity);

  for(unsigned int iEvent = 0 ; iEvent < GetNEvents(); ++iEvent){
    fdobj->isNC[iEvent] = false;
  }
}

void SampleHandlerFD::InitialiseSplineObject() {
  std::vector<std::string> spline_filepaths;
  for(unsigned iChannel = 0 ; iChannel < OscChannels.size() ; iChannel++){
    spline_filepaths.push_back(spline_files[iChannel]);
  }

  //Keep a track of the spline variables
  SplineVarNames.push_back("TrueNeutrinoEnergy");
  if(XVarStr.length() > 0){
    SplineVarNames.push_back(XVarStr);
  }
  if(YVarStr.length() > 0){
    SplineVarNames.push_back(YVarStr);
  }
  
  SplineHandler->AddSample(SampleName, spline_filepaths, SplineVarNames);
  SplineHandler->CountNumberOfLoadedSplines(false, 1);
  SplineHandler->TransferToMonolith();

  MACH3LOG_INFO("--------------------------------");
  MACH3LOG_INFO("Setup Far Detector splines");

  FillSplineBins();

  SplineHandler->cleanUpMemory();
}

// === JM adjust GetNDVarHist functions to allow for particle-level plotting ===
TH1* SampleHandlerFD::Get1DVarHist(const std::string& ProjectionVar_Str, const std::vector< KinematicCut >& SelectionVec, int WeightStyle, TAxis* Axis) {
  bool isParticleHist = IsParticleVarString(ProjectionVar_Str);
  //DB Grab the associated enum with the argument string
  int ProjectionVar_Int;
  if (!isParticleHist) ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);
  else ProjectionVar_Int = ReturnKinematicVectorFromString(ProjectionVar_Str);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< KinematicCut > tmp_Selection = Selection;
  std::vector< KinematicCut > SelectionVecToApply;

  //JM Separate vector for particle-level kinematic cuts
  std::vector< KinematicCut > ParticleCuts;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    //JM If particle-level cut, append to a different vector (avoids checking this in IsEventSelected)
    if (SelectionVec[iSelec].isParticleVar) ParticleCuts.emplace_back(SelectionVec[iSelec]);
    else SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH1D* _h1DVar;
  if (Axis) {
    _h1DVar = new TH1D("","",Axis->GetNbins(),Axis->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges;
    if (!isParticleHist) xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_Str);
    else xBinEdges = ReturnKinematicVectorBinning(ProjectionVar_Str);
    _h1DVar = new TH1D("", "", int(xBinEdges.size())-1, xBinEdges.data());
  }

  //DB Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    if (IsEventSelected(iEvent)) {
      double Weight = GetEventWeight(iEvent);
      if (WeightStyle == 1) {
        Weight = 1.;
      }
      if (!isParticleHist) {
        double Var = ReturnKinematicParameter(ProjectionVar_Int,iEvent);
        _h1DVar->Fill(Var,Weight);
      }
      else {
        std::vector<double> Vec = ReturnKinematicVector(ProjectionVar_Int,iEvent);
        size_t nparticles = Vec.size();
        //JM Loop over all particles in event
        for (unsigned int iParticle = 0; iParticle < nparticles; iParticle++) {
          if (IsParticleSelected(ParticleCuts, iEvent, iParticle, nparticles)) {
            double Var = Vec[iParticle];
            _h1DVar->Fill(Var,Weight);
          }
        }
      }
    }
  }

  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h1DVar;
}

// ************************************************
TH2* SampleHandlerFD::Get2DVarHist(const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY,
    const std::vector< KinematicCut >& SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  // ************************************************
  bool isParticleXVar = IsParticleVarString(ProjectionVar_StrX);
  bool isParticleYVar = IsParticleVarString(ProjectionVar_StrY);
  bool isParticleHist = isParticleXVar || isParticleYVar;

  //DB Grab the associated enum with the argument string
  int ProjectionVar_IntX, ProjectionVar_IntY;
  
  if (!isParticleXVar) ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  else ProjectionVar_IntX = ReturnKinematicVectorFromString(ProjectionVar_StrX);
  
  if (!isParticleYVar) ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);
  else ProjectionVar_IntY = ReturnKinematicVectorFromString(ProjectionVar_StrY);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< KinematicCut > tmp_Selection = Selection;
  std::vector< KinematicCut > SelectionVecToApply;

  //JM Separate vector for particle-level kinematic cuts
  std::vector< KinematicCut > ParticleCuts;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    //JM If particle-level cut, append to a different vector (avoids checking this in IsEventSelected)
    if (SelectionVec[iSelec].isParticleVar) ParticleCuts.emplace_back(SelectionVec[iSelec]);
    else SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH2D* _h2DVar;
  if (AxisX && AxisY) {
    _h2DVar = new TH2D("","",AxisX->GetNbins(),AxisX->GetXbins()->GetArray(),AxisY->GetNbins(),AxisY->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges, yBinEdges;
    if (!isParticleHist) {
      xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrX);
      yBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrY);
    }
    else {
      xBinEdges = ReturnKinematicVectorBinning(ProjectionVar_StrX);
      yBinEdges = ReturnKinematicVectorBinning(ProjectionVar_StrY);
    }
    _h2DVar = new TH2D("", "", int(xBinEdges.size())-1, xBinEdges.data(), int(yBinEdges.size())-1, yBinEdges.data());
  }

  //DB Loop over all events
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); iEvent++) {
    if (IsEventSelected(iEvent)) {
      double VarX = M3::_BAD_DOUBLE_, VarY = M3::_BAD_DOUBLE_;
      double Weight = GetEventWeight(iEvent);
      if (WeightStyle == 1) {
        Weight = 1.;
      }
      if (!isParticleHist) {
        VarX = ReturnKinematicParameter(ProjectionVar_IntX, iEvent);
        VarY = ReturnKinematicParameter(ProjectionVar_IntY, iEvent);
        _h2DVar->Fill(VarX,VarY,Weight);
      }
      else {
        std::vector<double> VecX = {}, VecY = {};
        size_t nparticles = 0;
        // JM Three cases: particleX vs eventY || eventX vs particleY || particleX vs particleY
        if (isParticleXVar && !isParticleYVar) {
          VecX = ReturnKinematicVector(ProjectionVar_IntX, iEvent);
          VarY = ReturnKinematicParameter(ProjectionVar_IntY, iEvent);
          nparticles = VecX.size();
        }
        else if (!isParticleXVar && isParticleYVar) {
          VecY = ReturnKinematicVector(ProjectionVar_IntY, iEvent);
          VarX = ReturnKinematicParameter(ProjectionVar_IntX, iEvent);
          nparticles = VecY.size();
        }
        else {
          VecX = ReturnKinematicVector(ProjectionVar_IntX, iEvent);
          VecY = ReturnKinematicVector(ProjectionVar_IntY, iEvent);
          if (VecX.size() != VecY.size()) {
            MACH3LOG_ERROR("Cannot plot {} of size {} against {} of size {}", ProjectionVar_StrX, VecX.size(), ProjectionVar_StrY, VecY.size());
            throw MaCh3Exception(__FILE__, __LINE__);
          }
          nparticles = VecX.size();
        }
        //JM Loop over all particles in event
        for (unsigned int iParticle = 0; iParticle < nparticles; iParticle++) {
          if (IsParticleSelected(ParticleCuts, iEvent, iParticle, nparticles)) {
            if (isParticleXVar) VarX = VecX[iParticle];
            if (isParticleYVar) VarY = VecY[iParticle];
            _h2DVar->Fill(VarX,VarY,Weight);
          }
        }
      }
    }
  }

  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h2DVar;
}
// ================================================

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

bool SampleHandlerFD::IsParticleVarString(const std::string& VarStr) {
  if (KinematicVectors->count(VarStr)) {
    if (!KinematicParameters->count(VarStr)) return true;
    else {
      MACH3LOG_ERROR("Attempted to plot kinematic variable {}, but it appears in both KinematicVectors and KinematicParameters");
      throw MaCh3Exception(__FILE__,__LINE__);
    }
  }
  return false;
}
// ===============================================================

TH1* SampleHandlerFD::Get1DVarHistByModeAndChannel(const std::string& ProjectionVar_Str, 
    int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill != -1) {
    if (kChannelToFill > GetNOscChannels()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}", GetNOscChannels());
      MACH3LOG_ERROR("kChannelToFill given:{}", kChannelToFill);
      MACH3LOG_ERROR("Exiting.");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
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

  return Get1DVarHist(ProjectionVar_Str,SelectionVec,WeightStyle,Axis);
}

TH2* SampleHandlerFD::Get2DVarHistByModeAndChannel(const std::string& ProjectionVar_StrX, const std::string& ProjectionVar_StrY, 
    int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill > GetNOscChannels()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}", GetNOscChannels());
      MACH3LOG_ERROR("kChannelToFill given:{}", kChannelToFill);
      MACH3LOG_ERROR("Exiting.");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
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

  return Get2DVarHist(ProjectionVar_StrX,ProjectionVar_StrY,SelectionVec,WeightStyle,AxisX,AxisY);
}

void SampleHandlerFD::PrintIntegral(TString OutputFileName, int WeightStyle, TString OutputCSVFileName) {
  int space = 14;

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
  ChannelIntegral.resize(GetNOscChannels());
  for (unsigned int i=0;i<ChannelIntegral.size();i++) {ChannelIntegral[i] = 0.;}

  for (int i=0;i<Modes->GetNModes();i++) {
    if (GetNDim()==1) {
      IntegralList[i] = ReturnHistsBySelection1D(XVarStr,1,i,WeightStyle);
    } else {
      IntegralList[i] = CastVector<TH2, TH1>(ReturnHistsBySelection2D(XVarStr,YVarStr,1,i,WeightStyle));
    }
  }

  MACH3LOG_INFO("-------------------------------------------------");

  if (printToFile) {
    outfile << "\\begin{table}[ht]" << std::endl;
    outfile << "\\begin{center}" << std::endl;
    outfile << "\\caption{Integral breakdown for sample: " << GetTitle() << "}" << std::endl;
    outfile << "\\label{" << GetTitle() << "-EventRate}" << std::endl;

    TString nColumns;
    for (int i=0;i<GetNOscChannels();i++) {nColumns+="|c";}
    nColumns += "|c|";
    outfile << "\\begin{tabular}{|l" << nColumns.Data() << "}" << std::endl;
    outfile << "\\hline" << std::endl;
  }

  if(printToCSV){
    // HW Probably a better way but oh well, here I go making MaCh3 messy again
    outcsv<<"Integral Breakdown for sample :"<<GetTitle()<<"\n";
  }

  MACH3LOG_INFO("Integral breakdown for sample: {}", GetTitle());
  MACH3LOG_INFO("");

  if (printToFile) {outfile << std::setw(space) << "Mode:";}
  if(printToCSV) {outcsv<<"Mode,";}

  std::string table_headings = fmt::format("| {:<8} |", "Mode");
  std::string table_footline = "------------"; //Scalable table horizontal line
  for (size_t i = 0;i < OscChannels.size(); i++) {
    table_headings += fmt::format(" {:<17} |", OscChannels[i].flavourName);
    table_footline += "--------------------";
    if (printToFile) {outfile << "&" << std::setw(space) << OscChannels[i].flavourName_Latex << " ";}
    if (printToCSV)  {outcsv << OscChannels[i].flavourName << ",";}
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

      if (Integral<1e-100) {Integral=0;}

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

std::vector<TH1*> SampleHandlerFD::ReturnHistsBySelection1D(std::string KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
  std::vector<TH1*> hHistList;
  std::string legendEntry;

  if (THStackLeg != nullptr) {delete THStackLeg;}
  THStackLeg = new TLegend(0.1,0.1,0.9,0.9);

  int iMax = -1;
  if (Selection1 == 0) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == 1) {
    iMax = GetNOscChannels();
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1==0) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(KinematicProjection,i,Selection2,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(Modes->GetMaCh3ModeName(i)+Form(" : (%4.2f)",hHistList[i]->Integral())).c_str(),"f");

      hHistList[i]->SetFillColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
      hHistList[i]->SetLineColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
    }
    if (Selection1==1) {
      hHistList.push_back(Get1DVarHistByModeAndChannel(KinematicProjection,Selection2,i,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(OscChannels[i].flavourName+Form(" | %4.2f",hHistList[i]->Integral())).c_str(),"f");
    }
  }

  return hHistList;
}

std::vector<TH2*> SampleHandlerFD::ReturnHistsBySelection2D(std::string KinematicProjectionX, std::string KinematicProjectionY,
    int Selection1, int Selection2, int WeightStyle, 
    TAxis* XAxis, TAxis* YAxis) {
  std::vector<TH2*> hHistList;

  int iMax = -1;
  if (Selection1 == 0) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == 1) {
    iMax = GetNOscChannels();
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1==0) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(KinematicProjectionX,KinematicProjectionY,i,Selection2,WeightStyle,XAxis,YAxis));
    }
    if (Selection1==1) {
      hHistList.push_back(Get2DVarHistByModeAndChannel(KinematicProjectionX,KinematicProjectionY,Selection2,i,WeightStyle,XAxis,YAxis));
    }
  }

  return hHistList;
}

THStack* SampleHandlerFD::ReturnStackedHistBySelection1D(std::string KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
  std::vector<TH1*> HistList = ReturnHistsBySelection1D(KinematicProjection, Selection1, Selection2, WeightStyle, XAxis);
  THStack* StackHist = new THStack((GetTitle()+"_"+KinematicProjection+"_Stack").c_str(),"");
  for (unsigned int i=0;i<HistList.size();i++) {
    StackHist->Add(HistList[i]);
  }
  return StackHist;
}


// ************************************************
const double* SampleHandlerFD::GetPointerToOscChannel(int iEvent) const {
  // ************************************************
  const int Channel = GetOscChannel(OscChannels, (*MCSamples.nupdgUnosc[iEvent]), (*MCSamples.nupdg[iEvent]));

  return &(OscChannels[Channel].ChannelIndex);
}
