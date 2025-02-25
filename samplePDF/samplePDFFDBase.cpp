#include "samplePDFFDBase.h"
#include "samplePDF/Structs.h"

_MaCh3_Safe_Include_Start_ //{
#include "Oscillator/OscillatorFactory.h"
#include "Constants/OscillatorConstants.h"
_MaCh3_Safe_Include_End_ //}

#include <algorithm>
#include <memory>

samplePDFFDBase::samplePDFFDBase(std::string ConfigFileName, covarianceXsec* xsec_cov, covarianceOsc* osc_cov) : samplePDFBase()
{
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Creating SamplePDFFDBase object");

  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(!xsec_cov){
    MACH3LOG_ERROR("You've passed me a nullptr to a covarianceXsec... I need this to setup splines!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  XsecCov = xsec_cov;

  if(!osc_cov){
    MACH3LOG_WARN("You have passed a nullptr to a covarianceOsc, this means I will not calculate oscillation weights");
  }
  OscCov = osc_cov;
  
  KinematicParameters = nullptr;
  ReversedKinematicParameters = nullptr;

  samplePDFFD_array = nullptr;
  samplePDFFD_data = nullptr;
  SampleName = "";
  SampleManager = std::unique_ptr<manager>(new manager(ConfigFileName.c_str()));
}

samplePDFFDBase::~samplePDFFDBase()
{
  MACH3LOG_DEBUG("I'm deleting samplePDFFDBase");
  
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    if(samplePDFFD_array != nullptr){delete[] samplePDFFD_array[yBin];}
    delete[] samplePDFFD_array_w2[yBin];
    //ETA - there is a chance that you haven't added any data...
    if(samplePDFFD_data != nullptr){delete[] samplePDFFD_data[yBin];}
  }

  if(samplePDFFD_array != nullptr){delete[] samplePDFFD_array;}
  delete[] samplePDFFD_array_w2;
  //ETA - there is a chance that you haven't added any data...
  if(samplePDFFD_data != nullptr){delete[] samplePDFFD_data;}
 
  for (unsigned int iCalc=0;iCalc<NuOscProbCalcers.size();iCalc++) {
    delete NuOscProbCalcers[iCalc];
  }

  if(THStackLeg != nullptr) delete THStackLeg;
}

void samplePDFFDBase::ReadSampleConfig() 
{
  auto ModeName = Get<std::string>(SampleManager->raw()["MaCh3ModeConfig"], __FILE__ , __LINE__);
  Modes = new MaCh3Modes(ModeName);
  SampleTitle = Get<std::string>(SampleManager->raw()["SampleTitle"], __FILE__ , __LINE__);
  SampleName = Get<std::string>(SampleManager->raw()["SampleName"], __FILE__ , __LINE__);
  NuOscillatorConfigFile = Get<std::string>(SampleManager->raw()["NuOsc"]["NuOscConfigFile"], __FILE__ , __LINE__);
  EqualBinningPerOscChannel = Get<bool>(SampleManager->raw()["NuOsc"]["EqualBinningPerOscChannel"], __FILE__ , __LINE__);
  
  //Default TestStatistic is kPoisson
  //ETA: this can be configured with samplePDFBase::SetTestStatistic()
  if (CheckNodeExists(SampleManager->raw(), "TestStatistic")) {
    fTestStatistic = static_cast<TestStatistic>(SampleManager->raw()["TestStatistic"].as<int>());
  } else {
    MACH3LOG_WARN("Didn't find a TestStatistic specified in {}", SampleManager->GetFileName());
    MACH3LOG_WARN("Defaulting to using a poisson likelihood");
    fTestStatistic = kPoisson;
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
      MACH3LOG_ERROR("Please specify an Y-variable string in sample config {}", SampleManager->GetFileName());
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
  std::string mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  std::string mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  std::string splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  std::string splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();
  
  nSamples = static_cast<M3::int_t>(SampleManager->raw()["SubSamples"].size());
  MCSamples.resize(nSamples);
  
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    oscchan_flavnames.push_back(osc_channel["Name"].as<std::string>());
    oscchan_flavnames_Latex.push_back(osc_channel["LatexName"].as<std::string>());
    mc_files.push_back(mtupleprefix+osc_channel["mtuplefile"].as<std::string>()+mtuplesuffix);
    spline_files.push_back(splineprefix+osc_channel["splinefile"].as<std::string>()+splinesuffix);
    sample_nupdgunosc.push_back(static_cast<NuPDG>(osc_channel["nutype"].as<int>()));
    sample_nupdg.push_back(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>()));
  }

  //Now get selection cuts
  double low_bound = 0;
  double up_bound = 0;
  double KinematicParamter = 0;
  std::vector<double> SelectionVec;
  //Now grab the selection cuts from the manager
  for ( auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
    SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());
    SelectionBounds.push_back(SelectionCuts["Bounds"].as<std::vector<double>>());
    low_bound = SelectionBounds.back().at(0);
    up_bound = SelectionBounds.back().at(1);
    KinematicParamter = static_cast<double>(ReturnKinematicParameterFromString(SelectionCuts["KinematicStr"].as<std::string>()));
    MACH3LOG_INFO("Adding cut on {} with bounds {} to {}", SelectionCuts["KinematicStr"].as<std::string>(), SelectionBounds.back().at(0), SelectionBounds.back().at(1));
    SelectionVec = {KinematicParamter, low_bound, up_bound};
    StoredSelection.push_back(SelectionVec);
  }
  NSelections = int(SelectionStr.size());

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
  for( int iMode=0; iMode<Modes->GetNModes(); iMode++ ) {
    std::string modeStr = Modes->GetMaCh3ModeName(iMode);
    MACH3LOG_INFO("    - {}: {}", modeStr, _modeNomWeightMap.at(modeStr));
  }
}

void samplePDFFDBase::Initialise() {
  //First grab all the information from your sample config via your manager
  ReadSampleConfig();

  //Now initialise all the variables you will need
  Init();

  int TotalMCEvents = 0;
  for(M3::int_t iSample=0 ; iSample < nSamples ; iSample++){
    MACH3LOG_INFO("=============================================");
    MACH3LOG_INFO("Initialising sample: {}/{}", iSample, nSamples);
    MCSamples[iSample].nEvents = setupExperimentMC(iSample);
    MACH3LOG_INFO("Number of events processed: {}", MCSamples[iSample].nEvents);
    TotalMCEvents += MCSamples[iSample].nEvents;
    MACH3LOG_INFO("Initialising FDMC object..");
    InitialiseSingleFDMCObject(iSample, MCSamples[iSample].nEvents);
    setupFDMC(iSample);
    MACH3LOG_INFO("Initialised sample: {}/{}", iSample, nSamples);
  }
  MACH3LOG_INFO("=============================================");
  MACH3LOG_INFO("Total number of events is: {}", TotalMCEvents);

  MACH3LOG_INFO("Setting up NuOscillator..");
  SetupNuOscillator(); 
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
void samplePDFFDBase::SetupKinematicMap() {
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

void samplePDFFDBase::fill1DHist()
{
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  _hPDF1D->Reset();
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    for (unsigned int xBin=0;xBin<(XBinEdges.size()-1);xBin++) {
      _hPDF1D->AddBinContent(xBin+1,samplePDFFD_array[yBin][xBin]);
    }
  }
}

void samplePDFFDBase::fill2DHist()
{
  // DB Commented out by default - Code heading towards GetLikelihood using arrays instead of root objects
  // Wouldn't actually need this for GetLikelihood as TH objects wouldn't be filled
  _hPDF2D->Reset();
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    for (unsigned int xBin=0;xBin<(XBinEdges.size()-1);xBin++) {
      _hPDF2D->SetBinContent(xBin+1,yBin+1,samplePDFFD_array[yBin][xBin]);
    }
  }
}

// ************************************************
/// @function samplePDFFDBase::SetupSampleBinning()
/// @brief Function to setup the binning of your sample histograms and the underlying
/// arrays that get handled in fillArray() and fillArray_MP().
/// The SampleXBins are filled in the daughter class from the sample config file.
/// This "passing" can be removed. 
void samplePDFFDBase::SetupSampleBinning(){
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
    set1DBinning(SampleXBins);  
  }
  else if(nDimensions == 2){
    MACH3LOG_INFO("Setting up 2D binning with {} and {}", XVarStr, YVarStr);
    set2DBinning(SampleXBins, SampleYBins);
  }
  else{
    MACH3LOG_ERROR("Number of dimensions is not 1 or 2, this is unsupported at the moment");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ************************************************
bool samplePDFFDBase::IsEventSelected(const int iSample, const int iEvent) {
// ************************************************
  double Val;
  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {  
    Val = ReturnKinematicParameter(Selection[iSelection][0], iSample, iEvent);
    if ((Val<Selection[iSelection][1])||(Val>=Selection[iSelection][2])) {
      return false;
    }
  }
  //DB To avoid unnecessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

//************************************************
// Reweight function - Depending on Osc Calculator this function uses different CalcOsc functions
void samplePDFFDBase::reweight() {
//************************************************
  //KS: Reset the histograms before reweight 
  ResetHistograms();
  
  //You only need to do these things if OscCov has been initialised
  //if not then you're not considering oscillations
  if (OscCov) {
    std::vector<M3::float_t> OscVec(OscParams.size());
    for (size_t iPar = 0; iPar < OscParams.size(); ++iPar) {
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wuseless-cast"
      OscVec[iPar] = M3::float_t(*OscParams[iPar]);
      #pragma GCC diagnostic pop
    }

    if (EqualBinningPerOscChannel) {
      NuOscProbCalcers[0]->CalculateProbabilities(OscVec);
    } else {
      for (int iSample=0;iSample<int(MCSamples.size());iSample++) {
        NuOscProbCalcers[iSample]->CalculateProbabilities(OscVec);
      }
    }
    
  }
  
  fillArray();
}

//************************************************
/// @function samplePDFFDBase::fillArray()
/// Function which does the core reweighting. This assumes that oscillation weights have 
/// already been calculated and stored in samplePDFFDBase[iSample].osc_w[iEvent]. This 
/// function takes advantage of most of the things called in setupSKMC to reduce reweighting time.
/// It also follows the ND code reweighting pretty closely. This function fills the samplePDFFD 
/// array array which is binned to match the sample binning, such that bin[1][1] is the 
/// equivalent of _hPDF2D->GetBinContent(2,2) {Noticing the offset}
void samplePDFFDBase::fillArray() {
//************************************************
  //DB Reset which cuts to apply
  Selection = StoredSelection;
  
  // Call entirely different routine if we're running with openMP
#ifdef MULTITHREAD
  fillArray_MP();
#else
  //ETA we should probably store this in samplePDFFDBase
  size_t nXBins = int(XBinEdges.size()-1);
  //size_t nYBins = int(YBinEdges.size()-1);

  PrepFunctionalParameters();
  if(SplineHandler){
    SplineHandler->Evaluate();
  }

  for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
    for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
      applyShifts(iSample, iEvent);
      
      if (!IsEventSelected(iSample, iEvent)) { 
        continue;
      } 

      double splineweight = 1.0;
      double normweight = 1.0;
      double totalweight = 1.0;
      
      if(SplineHandler){
        splineweight = CalcWeightSpline(iSample, iEvent);
      }
      //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient. Do this on a spline-by-spline basis
      if (splineweight <= 0.){
        MCSamples[iSample].xsec_w[iEvent] = 0.;
        continue;
      }
      
      //Loop over stored normalisation and function pointers 
      normweight = CalcWeightNorm(iSample, iEvent);
      
      //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
      if (normweight <= 0.){
        MCSamples[iSample].xsec_w[iEvent] = 0.;
        continue;
      }
      // Virtual by default does nothing
      CalcWeightFunc(iSample,iEvent);
      
      MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight;
      
      //DB Total weight
      totalweight = GetEventWeight(iSample,iEvent);
      //DB Catch negative weights and skip any event with a negative event
      if (totalweight <= 0.){
        MCSamples[iSample].xsec_w[iEvent] = 0.;
        continue;
      }
      //DB Switch on BinningOpt to allow different binning options to be implemented
      //The alternative would be to have inheritance based on BinningOpt
      double XVar = *(MCSamples[iSample].x_var[iEvent]);
      
      //DB Find the relevant bin in the PDF for each event
      int XBinToFill = -1;
      int YBinToFill = MCSamples[iSample].NomYBin[iEvent];
      
      //DB - First, check to see if the event is still in the nominal bin
      if (XVar < MCSamples[iSample].rw_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_xbinedge[iEvent]) {
        XBinToFill = MCSamples[iSample].NomXBin[iEvent];
      }
      //DB - Second, check to see if the event is outside of the binning range and skip event if it is
      //ETA- note that nXBins is XBinEdges.size() - 1
      else if (XVar < XBinEdges[0] || XVar >= XBinEdges[nXBins]) {
        continue;
      }
      //DB - Thirdly, check the adjacent bins first as Eb+CC+EScale shifts aren't likely to move an Erec more than 1bin width
      //Shifted down one bin from the event bin at nominal
      else if (XVar < MCSamples[iSample].rw_lower_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_lower_xbinedge[iEvent]) {
        XBinToFill = MCSamples[iSample].NomXBin[iEvent]-1;
      }
      //Shifted up one bin from the event bin at nominal
      else if (XVar < MCSamples[iSample].rw_upper_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_upper_xbinedge[iEvent]) {
        XBinToFill = MCSamples[iSample].NomXBin[iEvent]+1;
      }
      //DB - If we end up in this loop, the event has been shifted outside of its nominal bin, but is still within the allowed binning range
      else {
        for (unsigned int iBin=0;iBin<(XBinEdges.size()-1);iBin++)
        {
          if (XVar >= XBinEdges[iBin] && XVar < XBinEdges[iBin+1]) {
            XBinToFill = iBin;
          }
        }
      }
      
      //DB Fill relevant part of thread array
      if (XBinToFill != -1 && YBinToFill != -1) {
        samplePDFFD_array[YBinToFill][XBinToFill] += totalweight;
        samplePDFFD_array_w2[YBinToFill][XBinToFill] += totalweight*totalweight;
      }
    }
  }
#endif // end the else in openMP
}

#ifdef MULTITHREAD
// ************************************************ 
/// Multithreaded version of fillArray @see fillArray()
void samplePDFFDBase::fillArray_MP()  {
// ************************************************
  size_t nXBins = int(XBinEdges.size()-1);
  size_t nYBins = int(YBinEdges.size()-1);

  PrepFunctionalParameters();
  //==================================================
  //Calc Weights and fill Array
  if(SplineHandler){
    SplineHandler->Evaluate();
  }

  //This is stored as [y][x] due to shifts only occurring in the x variable (Erec/Lep mom) - I believe this will help reduce cache misses
  double** samplePDFFD_array_private = nullptr;
  double** samplePDFFD_array_private_w2 = nullptr;
  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
#pragma omp parallel private(samplePDFFD_array_private, samplePDFFD_array_private_w2)
  {
    // private to each thread
    // ETA - maybe we can use parallel firstprivate to initialise these?
    samplePDFFD_array_private = new double*[nYBins];
    samplePDFFD_array_private_w2 = new double*[nYBins];
    for (size_t yBin=0;yBin<nYBins;yBin++) {
      samplePDFFD_array_private[yBin] = new double[nXBins];
      samplePDFFD_array_private_w2[yBin] = new double[nXBins];

      std::fill_n(samplePDFFD_array_private[yBin], nXBins, 0.0);
      std::fill_n(samplePDFFD_array_private_w2[yBin], nXBins, 0.0);
    }
    
    //DB - Brain dump of speedup ideas
    //
    //Those relevant to reweighting
    // 1. Don't bother storing and calculating NC signal events - Implemented and saves marginal s/step
    // 2. Loop over spline event weight calculation in the following event loop - Currently done in splineSKBase->calcWeight() where multi-threading won't be optimised - Implemented and saves 0.3s/step
    // 3. Inline getDiscVar or somehow include that calculation inside the multi-threading - Implemented and saves about 0.01s/step
    // 4. Include isCC inside SKMCStruct so don't have to have several 'if' statements determine if oscillation weight needs to be set to 1.0 for NC events - Implemented and saves marginal s/step
    // 5. Do explict check on adjacent bins when finding event XBin instead of looping over all BinEdge indices - Implemented but doesn't significantly affect s/step
    //
    //Other aspects
    // 1. Order minituples in Y-axis variable as this will *hopefully* reduce cache misses inside samplePDFFD_array_class[yBin][xBin]
    //
    // We will hit <0.1 s/step eventually! :D
    
    for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
      #pragma omp for
      for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
        //ETA - generic functions to apply shifts to kinematic variables
        // Apply this before IsEventSelected is called.
        applyShifts(iSample, iEvent);

        //ETA - generic functions to apply shifts to kinematic variable
        //this is going to be slow right now due to string comps under the hood.
        //Need to implement a more efficient version of event-by-event cut checks
        if(!IsEventSelected(iSample, iEvent)){
          continue;
        }

        M3::float_t splineweight = 1.0;
        M3::float_t normweight = 1.0;
        M3::float_t totalweight = 1.0;

        //DB SKDet Syst
        //As weights were skdet::fParProp, and we use the non-shifted erec, we might as well cache the corresponding fParProp index for each event and the pointer to it

        if(SplineHandler){
          splineweight = CalcWeightSpline(iSample, iEvent);
        }
        //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
        if (splineweight <= 0.){
          MCSamples[iSample].xsec_w[iEvent] = 0.;
          continue;
        }

        normweight = CalcWeightNorm(iSample, iEvent);
        //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
        if (normweight <= 0.){
          MCSamples[iSample].xsec_w[iEvent] = 0.;
          continue;
        }

        // Virtual by default does nothing
        CalcWeightFunc(iSample,iEvent);

        MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight;

        totalweight = GetEventWeight(iSample, iEvent);

        //DB Catch negative weights and skip any event with a negative event
        if (totalweight <= 0.){
          MCSamples[iSample].xsec_w[iEvent] = 0.;
          continue;
        }

        //DB Switch on BinningOpt to allow different binning options to be implemented
        //The alternative would be to have inheritance based on BinningOpt
        double XVar = (*(MCSamples[iSample].x_var[iEvent]));

        //DB Commented out by default but if we ever want to consider shifts in theta this will be needed
        //double YVar = MCSamples[iSample].rw_theta[iEvent];
        //ETA - this would actually be with (*(MCSamples[iSample].y_var[iEvent])) and done extremely
        //similarly to XVar now

        //DB Find the relevant bin in the PDF for each event
        int XBinToFill = -1;
        int YBinToFill = MCSamples[iSample].NomYBin[iEvent];

        //DB Check to see if momentum shift has moved bins
        //DB - First, check to see if the event is still in the nominal bin
        if (XVar < MCSamples[iSample].rw_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_xbinedge[iEvent]) {
          XBinToFill = MCSamples[iSample].NomXBin[iEvent];
        }
        //DB - Second, check to see if the event is outside of the binning range and skip event if it is
        else if (XVar < XBinEdges[0] || XVar >= XBinEdges[nXBins]) {
          continue;
        }
        //DB - Thirdly, check the adjacent bins first as Eb+CC+EScale shifts aren't likely to move an Erec more than 1bin width
        //Shifted down one bin from the event bin at nominal
        else if (XVar < MCSamples[iSample].rw_lower_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_lower_xbinedge[iEvent]) {
          XBinToFill = MCSamples[iSample].NomXBin[iEvent]-1;
        }
        //Shifted up one bin from the event bin at nominal
        else if (XVar < MCSamples[iSample].rw_upper_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_upper_xbinedge[iEvent]) {
          XBinToFill = MCSamples[iSample].NomXBin[iEvent]+1;
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
          samplePDFFD_array_private[YBinToFill][XBinToFill] += totalweight;
          samplePDFFD_array_private_w2[YBinToFill][XBinToFill] += totalweight*totalweight;
        }
      }
    }    
    //End of Calc Weights and fill Array
    //==================================================
    // DB Copy contents of 'samplePDFFD_array_private' into 'samplePDFFD_array' which can then be used in GetLikelihood
    for (size_t yBin = 0; yBin < nYBins; ++yBin) {
      for (size_t xBin = 0; xBin < nXBins; ++xBin) {
        #pragma omp atomic
        samplePDFFD_array[yBin][xBin] += samplePDFFD_array_private[yBin][xBin];
        #pragma omp atomic
        samplePDFFD_array_w2[yBin][xBin] += samplePDFFD_array_private_w2[yBin][xBin];
      }
    }
    
    for (size_t yBin = 0; yBin < nYBins; ++yBin) {
      delete[] samplePDFFD_array_private[yBin];
      delete[] samplePDFFD_array_private_w2[yBin];
    }
    delete[] samplePDFFD_array_private;
    delete[] samplePDFFD_array_private_w2;
  } //end of parallel region
}
#endif

// **************************************************
// Helper function to reset the data and MC histograms
void samplePDFFDBase::ResetHistograms() {
// **************************************************
  size_t nXBins = int(XBinEdges.size()-1);
  size_t nYBins = int(YBinEdges.size()-1);
  
  //DB Reset values stored in PDF array to 0.
  for (size_t yBin = 0; yBin < nYBins; ++yBin) {
    for (size_t xBin = 0; xBin < nXBins; ++xBin) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }
} // end function

// ***************************************************************************
// Calculate the spline weight for one event
M3::float_t samplePDFFDBase::CalcWeightSpline(const int iSample, const int iEvent) const {
// ***************************************************************************
  M3::float_t xsecw = 1.0;
  //DB Xsec syst
  //Loop over stored spline pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iSpline = 0; iSpline < MCSamples[iSample].nxsec_spline_pointers[iEvent]; ++iSpline) {
    xsecw *= *(MCSamples[iSample].xsec_spline_pointers[iEvent][iSpline]);
  }
  return xsecw;
}

// ***************************************************************************
// Calculate the normalisation weight for one event
M3::float_t samplePDFFDBase::CalcWeightNorm(const int iSample, const int iEvent) const {
// ***************************************************************************
  M3::float_t xsecw = 1.0;
  //Loop over stored normalisation and function pointers
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iParam = 0; iParam < MCSamples[iSample].nxsec_norm_pointers[iEvent]; ++iParam)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
    xsecw *= static_cast<M3::float_t>(*(MCSamples[iSample].xsec_norm_pointers[iEvent][iParam]));
#pragma GCC diagnostic pop
    #ifdef DEBUG
    if (TMath::IsNaN(xsecw)) MACH3LOG_WARN("iParam= {} xsecweight=nan from norms", iParam);
    #endif
  }
  return xsecw;
}

void samplePDFFDBase::SetupNormParameters() {  
  xsec_norms = XsecCov->GetNormParsFromSampleName(SampleName);

  if(!XsecCov){
    MACH3LOG_ERROR("XsecCov is not setup!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Assign xsec norm bins in MCSamples tree
  for (unsigned int iSample = 0; iSample < MCSamples.size(); ++iSample) {
    CalcXsecNormsBins(iSample);
  }

  //DB
  //Attempt at reducing impact of covarianceXsec::calcReweight()
  int counter;

  for (int iSample = 0; iSample < int(MCSamples.size()); ++iSample) {
    for (int iEvent = 0; iEvent < MCSamples[iSample].nEvents; ++iEvent) {
      counter = 0;

      MCSamples[iSample].nxsec_norm_pointers[iEvent] = int(MCSamples[iSample].xsec_norms_bins[iEvent].size());
      MCSamples[iSample].xsec_norm_pointers[iEvent].resize(MCSamples[iSample].nxsec_norm_pointers[iEvent]);

      for(auto const & norm_bin: MCSamples[iSample].xsec_norms_bins[iEvent]) {
        MCSamples[iSample].xsec_norm_pointers[iEvent][counter] = XsecCov->retPointer(norm_bin);
        counter += 1;
      }
    }
  }
  // If not debugging let's clear memory
  #ifndef DEBUG
  for (size_t iSample = 0; iSample < MCSamples.size(); ++iSample) {
    MCSamples[iSample].xsec_norms_bins.clear();
    MCSamples[iSample].xsec_norms_bins.shrink_to_fit();
  }
  #endif
}

// ************************************************
//A way to check whether a normalisation parameter applies to an event or not
void samplePDFFDBase::CalcXsecNormsBins(int iSample) {
// ************************************************
  FarDetectorCoreInfo *fdobj = &MCSamples[iSample];
  #ifdef DEBUG
  std::vector<int> VerboseCounter(xsec_norms.size(), 0);
  #endif
  for(int iEvent = 0; iEvent < fdobj->nEvents; ++iEvent){
    std::vector< int > XsecBins = {};
    if (XsecCov) {
      // Skip oscillated NC events
      // Not strictly needed, but these events don't get included in oscillated predictions, so
      // no need to waste our time calculating and storing information about xsec parameters
      // that will never be used.
      if (fdobj->isNC[iEvent] && (*fdobj->nupdg[iEvent] != *fdobj->nupdgUnosc[iEvent]) ) {
        MACH3LOG_TRACE("Event {}, missed NC/signal check", iEvent);
        continue;
      } //DB Abstract check on MaCh3Modes to determine which apply to neutral current
      for (std::vector<XsecNorms4>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it) {
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
        bool IsSelected = true;
        if ((*it).hasKinBounds) {
          for (unsigned int iKinematicParameter = 0 ; iKinematicParameter < (*it).KinematicVarStr.size() ; ++iKinematicParameter ) {
            if (ReturnKinematicParameter((*it).KinematicVarStr[iKinematicParameter], iSample, iEvent) <= (*it).Selection[iKinematicParameter][0]) {
              IsSelected = false;
              MACH3LOG_TRACE("Event {}, missed Kinematic var check ({}) for dial {}", iEvent, (*it).KinematicVarStr[iKinematicParameter], (*it).name);
              continue;
            }
            else if (ReturnKinematicParameter((*it).KinematicVarStr[iKinematicParameter], iSample, iEvent) > (*it).Selection[iKinematicParameter][1]) {
              MACH3LOG_TRACE("Event {}, missed Kinematic var check ({}) for dial {}", iEvent, (*it).KinematicVarStr[iKinematicParameter], (*it).name);
              IsSelected = false;
              continue;
            }
          } 
        }
        //Need to then break the event loop 
        if(!IsSelected){
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", iEvent, (*it).name);
          continue;
        }
        // Now set 'index bin' for each normalisation parameter
        // All normalisations are just 1 bin for 2015, so bin = index (where index is just the bin for that normalisation)
        int bin = (*it).index;

        XsecBins.push_back(bin);
        MACH3LOG_TRACE("Event {}, will be affected by dial {}", iEvent, (*it).name);
        #ifdef DEBUG
        VerboseCounter[std::distance(xsec_norms.begin(), it)]++;
        #endif
        //}
      } // end iteration over xsec_norms
    } // end if (xsecCov)
    fdobj->xsec_norms_bins[iEvent] = XsecBins;
  }//end loop over events
  #ifdef DEBUG
  MACH3LOG_DEBUG("Channel {}", iSample);
  MACH3LOG_DEBUG("┌──────────────────────────────────────────────────────────┐");
  for (std::size_t i = 0; i < xsec_norms.size(); ++i) {
    const auto& norm = xsec_norms[i];
    double eventRatio = static_cast<double>(VerboseCounter[i]) / static_cast<double>(fdobj->nEvents);

    MACH3LOG_DEBUG("│ Param {:<15}, affects {:<8} events ({:>6.2f}%) │",
                  XsecCov->GetParFancyName(norm.index), VerboseCounter[i], eventRatio);
  }
  MACH3LOG_DEBUG("└──────────────────────────────────────────────────────────┘");
  #endif
}

//ETA - this is all a bit (less) stupid
void samplePDFFDBase::set1DBinning(std::vector<double> &XVec){
  _hPDF1D->Reset();
  _hPDF1D->SetBins(int(XVec.size()-1), XVec.data());
  dathist->SetBins(int(XVec.size()-1), XVec.data());

  //This will overwrite XBinEdges with whatever you pass this function
  XBinEdges = XVec;
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  _hPDF2D->Reset();
  _hPDF2D  ->SetBins(int(XVec.size()-1), XVec.data(), int(YBinEdges.size()-1), YBinEdges.data());
  dathist2d->SetBins(int(XVec.size()-1), XVec.data(), int(YBinEdges.size()-1), YBinEdges.data());

  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }
  FindNominalBinAndEdges1D();
}

//ETA - this is all a bit stupid
void samplePDFFDBase::set2DBinning(std::vector<double> &XVec, std::vector<double> &YVec)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(int(XVec.size()-1), XVec.data());
  dathist->SetBins(int(XVec.size()-1), XVec.data());

  _hPDF2D->Reset();
  _hPDF2D->SetBins(int(XVec.size()-1), XVec.data(), int(YVec.size()-1), YVec.data());
  dathist2d->SetBins(int(XVec.size()-1), XVec.data(), int(YVec.size()-1), YVec.data());

  //ETA - maybe need to be careful here
  int nXBins = int(XVec.size()-1);
  int nYBins = int(YVec.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges2D();
}

//ETA
//New versions of set binning functions is samplePDFBase
//so that we can set the values of the bin and lower/upper
//edges in the skmc_base. Hopefully we can use this to make
//fill1Dhist and fill2Dhist quicker
void samplePDFFDBase::set1DBinning(int nbins, double* boundaries)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins,boundaries);
  dathist->SetBins(nbins,boundaries);

  XBinEdges = std::vector<double>(nbins+1);
  for (int i=0;i<nbins+1;i++) {
    XBinEdges[i] = _hPDF1D->GetXaxis()->GetBinLowEdge(i+1);
  }
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  double YBinEdges_Arr[2];
  YBinEdges_Arr[0] = YBinEdges[0];
  YBinEdges_Arr[1] = YBinEdges[1];

  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins,boundaries,1,YBinEdges_Arr);
  dathist2d->SetBins(nbins,boundaries,1,YBinEdges_Arr);

  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges1D();
}

void samplePDFFDBase::set1DBinning(int nbins, double low, double high)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins,low,high);
  dathist->SetBins(nbins,low,high);

  XBinEdges = std::vector<double>(nbins+1);
  for (int i=0;i<nbins+1;i++) {
    XBinEdges[i] = _hPDF1D->GetXaxis()->GetBinLowEdge(i+1);
  }
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins,low,high,1,YBinEdges[0],YBinEdges[1]);
  dathist2d->SetBins(nbins,low,high,1,YBinEdges[0],YBinEdges[1]);

  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }
  FindNominalBinAndEdges1D();
}

void samplePDFFDBase::FindNominalBinAndEdges1D() {
  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(int mc_i = 0 ; mc_i < int(MCSamples.size()) ; mc_i++){
    for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){
      
      //Set x_var and y_var values based on XVarStr and YVarStr
      MCSamples[mc_i].x_var[event_i] = GetPointerToKinematicParameter(XVarStr, mc_i, event_i);
      //Give y_var M3::_BAD_DOUBLE_ value for the 1D case since this won't be used
      MCSamples[mc_i].y_var[event_i] = &(M3::_BAD_DOUBLE_);
      int bin = _hPDF1D->FindBin(*(MCSamples[mc_i].x_var[event_i]));
      
      double low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
      if (bin==0) {
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
		  MCSamples[mc_i].NomXBin[event_i] = bin-1;
	  } else {
		  MCSamples[mc_i].NomXBin[event_i] = -1;
		  low_edge = M3::_DEFAULT_RETURN_VAL_;
		  upper_edge = M3::_DEFAULT_RETURN_VAL_;
		  low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
		  upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
	  }
      MCSamples[mc_i].NomYBin[event_i] = 0;
      
      MCSamples[mc_i].rw_lower_xbinedge[event_i] = low_edge;
      MCSamples[mc_i].rw_upper_xbinedge[event_i] = upper_edge;
      MCSamples[mc_i].rw_lower_lower_xbinedge[event_i] = low_lower_edge;
      MCSamples[mc_i].rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
    }
  }
}

void samplePDFFDBase::set2DBinning(int nbins1, double* boundaries1, int nbins2, double* boundaries2)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins1,boundaries1);
  dathist->SetBins(nbins1,boundaries1);

  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins1,boundaries1,nbins2,boundaries2);
  dathist2d->SetBins(nbins1,boundaries1,nbins2,boundaries2);

  XBinEdges = std::vector<double>(nbins1+1);
  for (int i=0;i<nbins1+1;i++) {
    XBinEdges[i] = _hPDF2D->GetXaxis()->GetBinLowEdge(i+1);
  }
  YBinEdges = std::vector<double>(nbins2+1);
  for (int i=0;i<nbins2+1;i++) {
    YBinEdges[i] = _hPDF2D->GetYaxis()->GetBinLowEdge(i+1);
  }
  
  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges2D();
}

void samplePDFFDBase::set2DBinning(int nbins1, double low1, double high1, int nbins2, double low2, double high2)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins1,low1,high1);
  dathist->SetBins(nbins1,low1,high1);

  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins1,low1,high1,nbins2,low2,high2);
  dathist2d->SetBins(nbins1,low1,high1,nbins2,low2,high2);

  XBinEdges = std::vector<double>(nbins1+1);
  for (int i=0;i<nbins1+1;i++) {
    XBinEdges[i] = _hPDF2D->GetXaxis()->GetBinLowEdge(i+1);
  }
  YBinEdges = std::vector<double>(nbins2+1);
  for (int i=0;i<nbins2+1;i++) {
    YBinEdges[i] = _hPDF2D->GetYaxis()->GetBinLowEdge(i+1);
  }

  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_array = new double*[nYBins];
  samplePDFFD_array_w2 = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    samplePDFFD_array_w2[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges2D();
}

// ************************************************
void samplePDFFDBase::FindNominalBinAndEdges2D() {
// ************************************************
  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(int mc_i = 0 ; mc_i < int(MCSamples.size()) ; mc_i++){
    for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){
      
      //Set x_var and y_var values based on XVarStr and YVarStr   
      MCSamples[mc_i].x_var[event_i] = GetPointerToKinematicParameter(XVarStr, mc_i, event_i);
      MCSamples[mc_i].y_var[event_i] = GetPointerToKinematicParameter(YVarStr, mc_i, event_i);
      
      //Global bin number
      int bin = _hPDF2D->FindBin(*(MCSamples[mc_i].x_var[event_i]), *(MCSamples[mc_i].y_var[event_i]));
      
      int bin_x = -999;
      int bin_y = -999;
      int bin_z = -999;
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
        MCSamples[mc_i].NomXBin[event_i] = bin_x-1;
      } else {
        MCSamples[mc_i].NomXBin[event_i] = -1;
        low_edge = M3::_DEFAULT_RETURN_VAL_;
        upper_edge = M3::_DEFAULT_RETURN_VAL_;
        low_lower_edge = M3::_DEFAULT_RETURN_VAL_;
        upper_upper_edge = M3::_DEFAULT_RETURN_VAL_;
      }
      MCSamples[mc_i].NomYBin[event_i] = bin_y-1; 
      if(MCSamples[mc_i].NomYBin[event_i] < 0){
        MACH3LOG_INFO("Nominal YBin PROBLEM, y-bin is {}", MCSamples[mc_i].NomYBin[event_i]);
      }
      MCSamples[mc_i].rw_lower_xbinedge[event_i] = low_edge;
      MCSamples[mc_i].rw_upper_xbinedge[event_i] = upper_edge;
      MCSamples[mc_i].rw_lower_lower_xbinedge[event_i] = low_lower_edge;
      MCSamples[mc_i].rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
    }
  }
}

void samplePDFFDBase::addData(std::vector<double> &data) {
  dathist2d = nullptr;
  dathist->Reset(); 
  
  if (GetNDim()!=1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  for (int i = 0; i < int(data.size()); i++) {
    dathist->Fill(data.at(i));
  }
  
  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);
  
  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist->GetBinContent(xBin+1);
    }
  }
}

void samplePDFFDBase::addData(std::vector< std::vector <double> > &data) {
  dathist = nullptr;
  dathist2d->Reset();                                                       

  if (GetNDim()!=2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i = 0; i < int(data.size()); i++) {
    dathist2d->Fill(data.at(0)[i],data.at(1)[i]);
  }

  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist2d->GetBinContent(xBin+1,yBin+1);
    }
  }
}

void samplePDFFDBase::addData(TH1D* Data) {
  MACH3LOG_INFO("Adding 1D data histogram: {} with {:.2f} events", Data->GetName(), Data->Integral());
  dathist2d = nullptr;
  dathist = Data;
  
  if (GetNDim()!=1) {
    MACH3LOG_ERROR("Trying to set a 1D 'data' histogram in a 2D sample - Quitting"); 
    throw MaCh3Exception(__FILE__ , __LINE__ );}
  
  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);
  
  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = Data->GetBinContent(xBin+1);
    }
  }
}

void samplePDFFDBase::addData(TH2D* Data) {
  MACH3LOG_INFO("Adding 2D data histogram: {} with {:.2f} events", Data->GetName(), Data->Integral());
  dathist2d = Data;
  dathist = nullptr;

  if (GetNDim()!=2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram in a 1D sample - Quitting"); 
    throw MaCh3Exception(__FILE__ , __LINE__ );}	
  
  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist2d->GetBinContent(xBin+1,yBin+1);
    }
  }
}

// ************************************************
void samplePDFFDBase::SetupNuOscillator() {
// ************************************************
  OscillatorFactory* OscillFactory = new OscillatorFactory();
  if (!OscCov) {
    MACH3LOG_WARN("Attempted to setup NuOscillator without covarianceOsc object");
    return;
  }

  //DB's explanation of EqualBinningPerOscChannel:
  //In the situation where we are applying binning oscillation probabilities to a samplePDF object, it maybe the case that there is identical binning per oscillation channel
  //In which case, and remembering that each NuOscillator::Oscillator object calculate the oscillation probabilities for all channels, we just have to create one Oscillator object and use the results from that
  //This means that we can get upto a factor of 12 reduction in the calculation time of the oscillation probabilities, because we don't need to repeat the operation per oscillation channel

  if (EqualBinningPerOscChannel) {
    NuOscProbCalcers = std::vector<OscillatorBase*>(1);
    LoggerPrint("NuOscillator",
		[](const std::string& message) { MACH3LOG_INFO("{}", message); },
		[this, &OscillFactory]() {
		  this->NuOscProbCalcers[0] = OscillFactory->CreateOscillator(this->NuOscillatorConfigFile);
		});

    if (!NuOscProbCalcers[0]->EvalPointsSetInConstructor()) {
      MACH3LOG_ERROR("Attempted to use equal binning per oscillation channel, but not binning has been set in the NuOscillator::Oscillator object");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    NuOscProbCalcers[0]->Setup();
  } else {
    NuOscProbCalcers = std::vector<OscillatorBase*>(int(MCSamples.size()));
    for (size_t iSample=0;iSample<MCSamples.size();iSample++) {
      MACH3LOG_INFO("Setting up NuOscillator::Oscillator object in OscillationChannel: {}/{}", iSample, MCSamples.size());
      
      LoggerPrint("NuOscillator",
		  [](const std::string& message) { MACH3LOG_INFO("{}", message); },
		  [this, iSample, &OscillFactory]() {
		    this->NuOscProbCalcers[iSample] = OscillFactory->CreateOscillator(this->NuOscillatorConfigFile);
		  });
      
      if (!NuOscProbCalcers[iSample]->EvalPointsSetInConstructor()) {
	std::vector<M3::float_t> EnergyArray;
	for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
	  //DB Remove NC events from the arrays which are handed to the NuOscillator objects
	  if (!MCSamples[iSample].isNC[iEvent]) {
	    EnergyArray.push_back(M3::float_t(*(MCSamples[iSample].rw_etru[iEvent])));
	  }
	}
	std::sort(EnergyArray.begin(),EnergyArray.end());
	NuOscProbCalcers[iSample]->SetEnergyArrayInCalcer(EnergyArray);
	
	//============================================================================
	//DB Atmospheric only part
	if (MCSamples[iSample].rw_truecz.size() > 0 && int(MCSamples[iSample].rw_truecz.size()) == MCSamples[iSample].nEvents) { //Can only happen if truecz has been initialised within the experiment specific code
	  std::vector<M3::float_t> CosineZArray;
	  for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
	    //DB Remove NC events from the arrays which are handed to the NuOscillator objects
	    if (!MCSamples[iSample].isNC[iEvent]) {
	      CosineZArray.push_back(M3::float_t(*(MCSamples[iSample].rw_truecz[iEvent])));
	    }
	  }
	  std::sort(CosineZArray.begin(),CosineZArray.end());
	  
	  NuOscProbCalcers[iSample]->SetCosineZArrayInCalcer(CosineZArray);
	}
      }
      NuOscProbCalcers[iSample]->Setup();
    }
  }
  
  for (size_t iSample=0;iSample<MCSamples.size();iSample++) {

    for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
      // KS: Sry but if we use low memory we need to point to float not double...
#ifdef _LOW_MEMORY_STRUCTS_
      MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Unity_F;
#else
      MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Unity;
#endif
      if (MCSamples[iSample].isNC[iEvent]) {
        if (*MCSamples[iSample].nupdg[iEvent] != *MCSamples[iSample].nupdgUnosc[iEvent]) {
#ifdef _LOW_MEMORY_STRUCTS_
          MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Zero_F;
#else
          MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Zero;
#endif
        } else {
#ifdef _LOW_MEMORY_STRUCTS_
          MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Unity_F;
#else
          MCSamples[iSample].osc_w_pointer[iEvent] = &M3::Unity;
#endif
        }
      } else {
        int InitFlav = M3::_BAD_INT_;
        int FinalFlav = M3::_BAD_INT_;

        InitFlav =  MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples[iSample].nupdgUnosc[iEvent]));
        FinalFlav = MaCh3Utils::PDGToNuOscillatorFlavour((*MCSamples[iSample].nupdg[iEvent]));

        if (InitFlav == M3::_BAD_INT_ || FinalFlav == M3::_BAD_INT_) {
          MACH3LOG_ERROR("Something has gone wrong in the mapping between MCSamples[iSample].nutype and the enum used within NuOscillator");
          MACH3LOG_ERROR("MCSamples[iSample].nupdgUnosc: {}", (*MCSamples[iSample].nupdgUnosc[iEvent]));
          MACH3LOG_ERROR("InitFlav: {}", InitFlav);
          MACH3LOG_ERROR("MCSamples[iSample].nupdg: {}", (*MCSamples[iSample].nupdg[iEvent]));
          MACH3LOG_ERROR("FinalFlav: {}", FinalFlav);
          throw MaCh3Exception(__FILE__, __LINE__);
        }

	int Index = 0;
	if (!EqualBinningPerOscChannel) {
	  Index = static_cast<int>(iSample);
	}
	if (MCSamples[iSample].rw_truecz.size() > 0) { //Can only happen if truecz has been initialised within the experiment specific code
	  //Atmospherics
	  MCSamples[iSample].osc_w_pointer[iEvent] = NuOscProbCalcers[Index]->ReturnWeightPointer(InitFlav,FinalFlav,FLOAT_T(*(MCSamples[iSample].rw_etru[iEvent])),FLOAT_T(*(MCSamples[iSample].rw_truecz[iEvent])));
	} else {
	  //Beam
	  MCSamples[iSample].osc_w_pointer[iEvent] = NuOscProbCalcers[Index]->ReturnWeightPointer(InitFlav,FinalFlav,FLOAT_T(*(MCSamples[iSample].rw_etru[iEvent])));
        }

	
      } // end if NC
    } // end loop over events
  }// end loop over channels
  delete OscillFactory;

  OscParams = OscCov->GetOscParsFromSampleName(SampleName);
}

M3::float_t samplePDFFDBase::GetEventWeight(const int iSample, const int iEntry) const {
  M3::float_t totalweight = 1.0;
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int iParam = 0; iParam < MCSamples[iSample].ntotal_weight_pointers[iEntry]; ++iParam) {
    totalweight *= *(MCSamples[iSample].total_weight_pointers[iEntry][iParam]);
  }
  
  return totalweight;
}

/// @func fillSplineBins()
/// @brief Finds the binned spline that an event should apply to and stored them in a
/// a vector for easy evaluation in the fillArray() function.
void samplePDFFDBase::fillSplineBins() {
  for (int i = 0; i < int(MCSamples.size()); ++i) {
    //Now loop over events and get the spline bin for each event
    for (int j = 0; j < MCSamples[i].nEvents; ++j) {
      std::vector< std::vector<int> > EventSplines;
      switch(nDimensions){
        case 1:
          EventSplines = SplineHandler->GetEventSplines(GetTitle(), i, int(*(MCSamples[i].mode[j])), *(MCSamples[i].rw_etru[j]), *(MCSamples[i].x_var[j]), 0.);
          break;
        case 2:
          EventSplines = SplineHandler->GetEventSplines(GetTitle(), i, int(*(MCSamples[i].mode[j])), *(MCSamples[i].rw_etru[j]), *(MCSamples[i].x_var[j]), *(MCSamples[i].y_var[j]));
          break;
        default:
          MACH3LOG_ERROR("Error in assigning spline bins because nDimensions = {}", nDimensions);
          MACH3LOG_ERROR("MaCh3 only supports splines binned in Etrue + the sample binning");
          MACH3LOG_ERROR("Please check the sample binning you specified in your sample config ");
          throw MaCh3Exception(__FILE__, __LINE__);
          break;
      }
      MCSamples[i].nxsec_spline_pointers[j] = int(EventSplines.size());
      if(MCSamples[i].nxsec_spline_pointers[j] < 0){
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MCSamples[i].xsec_spline_pointers[j].resize(MCSamples[i].nxsec_spline_pointers[j]);
      for(int spline=0; spline<MCSamples[i].nxsec_spline_pointers[j]; spline++){
        //Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
        MCSamples[i].xsec_spline_pointers[j][spline] = SplineHandler->retPointer(EventSplines[spline][0], EventSplines[spline][1],
                                                                                 EventSplines[spline][2], EventSplines[spline][3],
                                                                                 EventSplines[spline][4], EventSplines[spline][5],
                                                                                 EventSplines[spline][6]);
      }
    }
  }
}

// ************************************************
double samplePDFFDBase::GetLikelihood() {
// ************************************************
  if (samplePDFFD_data == nullptr) {
    MACH3LOG_ERROR("Data sample is empty! Can't calculate a likelihood!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  //This can be done only once and stored
  int nXBins = int(XBinEdges.size()-1);
  int nYBins = int(YBinEdges.size()-1);
  
  double negLogL = 0.;
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2) reduction(+:negLogL)
  #endif
  for (int xBin = 0; xBin < nXBins; ++xBin)
  {
    for (int yBin = 0; yBin < nYBins; ++yBin)
    {
      const double DataVal = samplePDFFD_data[yBin][xBin];
      const double MCPred = samplePDFFD_array[yBin][xBin];
      const double w2 = samplePDFFD_array_w2[yBin][xBin];
      
      //KS: Calculate likelihood using Barlow-Beeston Poisson or even IceCube
      negLogL += getTestStatLLH(DataVal, MCPred, w2);
    }
  }
  return negLogL;
}

void samplePDFFDBase::InitialiseSingleFDMCObject(int iSample, int nEvents_) {
  if (iSample < 0 || iSample >= nSamples) {
    MACH3LOG_ERROR("Invalid iSample index in InitialiseSingleFDMCObject");
    MACH3LOG_ERROR("Index given is {} and only {} samples found", iSample, nSamples);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  FarDetectorCoreInfo *fdobj = &MCSamples[iSample];
  
  fdobj->nEvents = nEvents_;
  fdobj->flavourName = oscchan_flavnames[iSample];
  fdobj->flavourName_Latex = oscchan_flavnames_Latex[iSample];
  fdobj->ChannelIndex = iSample;
  
  int nEvents = fdobj->nEvents;
  fdobj->x_var.resize(nEvents, &M3::Unity);
  fdobj->y_var.resize(nEvents, &M3::Unity);
  fdobj->rw_etru.resize(nEvents, &M3::Unity);
  fdobj->XBin.resize(nEvents, -1);
  fdobj->YBin.resize(nEvents, -1);
  fdobj->NomXBin.resize(nEvents, -1);
  fdobj->NomYBin.resize(nEvents, -1);
  fdobj->rw_lower_xbinedge.resize(nEvents, -1);
  fdobj->rw_lower_lower_xbinedge.resize(nEvents, -1);
  fdobj->rw_upper_xbinedge.resize(nEvents, -1);
  fdobj->rw_upper_upper_xbinedge.resize(nEvents, -1);
  fdobj->mode.resize(nEvents, &M3::Unity);
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
#ifdef _LOW_MEMORY_STRUCTS_
  fdobj->osc_w_pointer.resize(nEvents, &M3::Unity_F);
#else
  fdobj->osc_w_pointer.resize(nEvents, &M3::Unity);
#endif

  for(int iEvent = 0 ; iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->isNC[iEvent] = false;
  }
}

void samplePDFFDBase::InitialiseSplineObject() {
  std::vector<std::string> spline_filepaths;
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    spline_filepaths.push_back(spline_files[iSample]);
  }

  //Keep a track of the spline variables
  SplineVarNames.push_back("TrueNeutrinoEnergy");
  if(XVarStr.length() > 0){
    SplineVarNames.push_back(XVarStr);
  }
  if(YVarStr.length() > 0){
    SplineVarNames.push_back(YVarStr);
  }
  
  SplineHandler->AddSample(SampleTitle, SampleName, spline_filepaths, SplineVarNames);
  SplineHandler->CountNumberOfLoadedSplines(false, 1);
  SplineHandler->TransferToMonolith();

  MACH3LOG_INFO("--------------------------------");
  MACH3LOG_INFO("Setup Far Detector splines");

  fillSplineBins();

  SplineHandler->cleanUpMemory();
}

TH1* samplePDFFDBase::get1DVarHist(std::string ProjectionVar_Str, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector<double> > tmp_Selection = Selection;
  std::vector< std::vector<double> > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Check the formatting of all requested cuts, should be [cutPar,lBound,uBound]
  for (size_t iSelec=0;iSelec<SelectionVecToApply.size();iSelec++) {
    if (SelectionVecToApply[iSelec].size()!=3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect size == 3, given: {}",iSelec,SelectionVecToApply[iSelec].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH1D* _h1DVar;
  if (Axis) {
    _h1DVar = new TH1D("","",Axis->GetNbins(),Axis->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_Str);
    _h1DVar = new TH1D("", "", int(xBinEdges.size())-1, xBinEdges.data());
  }

  //DB Loop over all events
  for (int iSample=0;iSample<getNMCSamples();iSample++) {
    for (int iEvent=0;iEvent<getNEventsInSample(iSample);iEvent++) {
      if (IsEventSelected(iSample,iEvent)) {
        double Weight = GetEventWeight(iSample,iEvent);
        if (WeightStyle==1) {
          Weight = 1.;
        }
        double Var = ReturnKinematicParameter(ProjectionVar_Int,iSample,iEvent);
        _h1DVar->Fill(Var,Weight);
      }
    }
  }
  
  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h1DVar;
}

TH2* samplePDFFDBase::get2DVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  int ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector<double> > tmp_Selection = Selection;
  std::vector< std::vector<double> > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Check the formatting of all requested cuts, should be [cutPar,lBound,uBound]
  for (size_t iSelec=0;iSelec<SelectionVecToApply.size();iSelec++) {
    if (SelectionVecToApply[iSelec].size()!=3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect size == 3, given: {}",iSelec,SelectionVecToApply[iSelec].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH2D* _h2DVar;
  if (AxisX && AxisY) {
    _h2DVar = new TH2D("","",AxisX->GetNbins(),AxisX->GetXbins()->GetArray(),AxisY->GetNbins(),AxisY->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrX);
    std::vector<double> yBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrY);
    _h2DVar = new TH2D("", "", int(xBinEdges.size())-1, xBinEdges.data(), int(yBinEdges.size())-1, yBinEdges.data());
  }

  //DB Loop over all events
  for (int iSample=0;iSample<getNMCSamples();iSample++) {
    for (int iEvent=0;iEvent<getNEventsInSample(iSample);iEvent++) {
      if (IsEventSelected(iSample,iEvent)) {
        double Weight = GetEventWeight(iSample,iEvent);
        if (WeightStyle==1) {
          Weight = 1.;
        }
        double VarX = ReturnKinematicParameter(ProjectionVar_IntX,iSample,iEvent);
	double VarY = ReturnKinematicParameter(ProjectionVar_IntY,iSample,iEvent);
        _h2DVar->Fill(VarX,VarY,Weight);
      }
    }
  }
  
  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h2DVar;
}

// ************************************************
int samplePDFFDBase::ReturnKinematicParameterFromString(const std::string& KinematicParameterStr) const {
// ************************************************
  auto it = KinematicParameters->find(KinematicParameterStr);
  if (it != KinematicParameters->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameterStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return -999;
}

// ************************************************
std::string samplePDFFDBase::ReturnStringFromKinematicParameter(const int KinematicParameter) const {
// ************************************************
  auto it = ReversedKinematicParameters->find(KinematicParameter);
  if (it != ReversedKinematicParameters->end()) {
    return it->second;
  }

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameter);
  throw MaCh3Exception(__FILE__, __LINE__);

  return "";
}

TH1* samplePDFFDBase::get1DVarHistByModeAndChannel(std::string ProjectionVar_Str, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>getNMCSamples()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}", getNMCSamples());
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

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = ReturnKinematicParameterFromString("Mode");
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = ReturnKinematicParameterFromString("OscillationChannel");
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(ProjectionVar_Str,SelectionVec,WeightStyle,Axis);
}

TH2* samplePDFFDBase::get2DVarHistByModeAndChannel(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>getNMCSamples()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}", getNMCSamples());
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

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = ReturnKinematicParameterFromString("Mode");
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = ReturnKinematicParameterFromString("OscillationChannel");
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get2DVarHist(ProjectionVar_StrX,ProjectionVar_StrY,SelectionVec,WeightStyle,AxisX,AxisY);
}

void samplePDFFDBase::PrintIntegral(TString OutputFileName, int WeightStyle, TString OutputCSVFileName) {
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
  ChannelIntegral.resize(getNMCSamples());
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
    for (int i=0;i<getNMCSamples();i++) {nColumns+="|c";}
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
  for (int i=0;i<getNMCSamples();i++) {
    table_headings += fmt::format(" {:<17} |", MCSamples[i].flavourName);
    table_footline += "--------------------";
    if (printToFile) {outfile << "&" << std::setw(space) << MCSamples[i].flavourName_Latex << " ";}
    if (printToCSV)  {outcsv << MCSamples[i].flavourName << ",";}
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

std::vector<TH1*> samplePDFFDBase::ReturnHistsBySelection1D(std::string KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
  std::vector<TH1*> hHistList;
  std::string legendEntry;

  if (THStackLeg != nullptr) {delete THStackLeg;}
  THStackLeg = new TLegend(0.1,0.1,0.9,0.9);

  int iMax = -1;
  if (Selection1 == 0) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == 1) {
    iMax = getNMCSamples();
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1==0) {
      hHistList.push_back(get1DVarHistByModeAndChannel(KinematicProjection,i,Selection2,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(Modes->GetMaCh3ModeName(i)+Form(" : (%4.2f)",hHistList[i]->Integral())).c_str(),"f");

      hHistList[i]->SetFillColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
      hHistList[i]->SetLineColor(static_cast<Color_t>(Modes->GetMaCh3ModePlotColor(i)));
    }
    if (Selection1==1) {
      hHistList.push_back(get1DVarHistByModeAndChannel(KinematicProjection,Selection2,i,WeightStyle,XAxis));
      THStackLeg->AddEntry(hHistList[i],(MCSamples[i].flavourName+Form(" | %4.2f",hHistList[i]->Integral())).c_str(),"f");
    }
  }

  return hHistList;
}

std::vector<TH2*> samplePDFFDBase::ReturnHistsBySelection2D(std::string KinematicProjectionX, std::string KinematicProjectionY, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis, TAxis* YAxis) {
  std::vector<TH2*> hHistList;

  int iMax = -1;
  if (Selection1 == 0) {
    iMax = Modes->GetNModes();
  }
  if (Selection1 == 1) {
    iMax = getNMCSamples();
  }
  if (iMax == -1) {
    MACH3LOG_ERROR("You've passed me a Selection1 which was not implemented in ReturnHistsBySelection1D. Selection1 and Selection2 are counters for different indexable quantities");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i=0;i<iMax;i++) {
    if (Selection1==0) {
      hHistList.push_back(get2DVarHistByModeAndChannel(KinematicProjectionX,KinematicProjectionY,i,Selection2,WeightStyle,XAxis,YAxis));
    }
    if (Selection1==1) {
      hHistList.push_back(get2DVarHistByModeAndChannel(KinematicProjectionX,KinematicProjectionY,Selection2,i,WeightStyle,XAxis,YAxis));
    }
  }

  return hHistList;
}

THStack* samplePDFFDBase::ReturnStackedHistBySelection1D(std::string KinematicProjection, int Selection1, int Selection2, int WeightStyle, TAxis* XAxis) {
  std::vector<TH1*> HistList = ReturnHistsBySelection1D(KinematicProjection, Selection1, Selection2, WeightStyle, XAxis);
  THStack* StackHist = new THStack((GetTitle()+"_"+KinematicProjection+"_Stack").c_str(),"");
  for (unsigned int i=0;i<HistList.size();i++) {
    StackHist->Add(HistList[i]);
  }
  return StackHist;
}
