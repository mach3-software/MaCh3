#include "samplePDFFDBase.h"

#include "Oscillator/OscillatorFactory.h"
#include "Constants/OscillatorConstants.h"

#include<algorithm>


samplePDFFDBase::samplePDFFDBase(double pot, std::string mc_version_, covarianceXsec* xsec_cov)
  : samplePDFBase()
//DB Throughout constructor and init, pot is livetime for atmospheric samples
{
  (void) pot;


  std::cout << "-------------------------------------------------------------------" <<std::endl;
  MACH3LOG_INFO("Ceating SamplePDFFDBase object");
    
  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){
	MACH3LOG_ERROR("You've passed me a NULL xsec covariance matrix... I need this to setup splines!");
   	throw MaCh3Exception(__FILE__, __LINE__);
  }
  SetXsecCov(xsec_cov);
  
  samplePDFFD_array = NULL;
  samplePDFFD_data = NULL;
 
  char* sample_char = (char*)mc_version_.c_str();
  SampleManager = new manager(sample_char);
}

samplePDFFDBase::~samplePDFFDBase()
{
  MACH3LOG_INFO("I'm deleting samplePDFFDBase");

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

  if(oscpars != nullptr) delete[] oscpars;
}

void samplePDFFDBase::ReadSampleConfig() 
{
   
  if (CheckNodeExists(SampleManager->raw(), "SampleName")) {
    samplename = SampleManager->raw()["SampleName"].as<std::string>();
  } else{
    MACH3LOG_ERROR("SampleName not defined in {}, please add this!", SampleManager->GetFileName());
  }
  
  if (CheckNodeExists(SampleManager->raw(), "NSubSamples")) {
    nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  } else{
    MACH3LOG_ERROR("NSubSamples not defined in {}, please add this!", SampleManager->GetFileName());
  }
  
  if (CheckNodeExists(SampleManager->raw(), "DetID")) {
    SampleDetID = SampleManager->raw()["DetID"].as<int>();
  } else{
    MACH3LOG_ERROR("ID not defined in {}, please add this!", SampleManager->GetFileName());
  }
  
  if (CheckNodeExists(SampleManager->raw(), "NuOsc", "NuOscConfigFile")) {
    NuOscillatorConfigFile = SampleManager->raw()["NuOsc"]["NuOscConfigFile"].as<std::string>();
  } else {
    MACH3LOG_ERROR("NuOsc::NuOscConfigFile is not defined in {}, please add this!", SampleManager->GetFileName());
  }
  
  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

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
	  MACH3LOG_ERROR("Please specify an X-variable string in sample config {}", SampleManager->GetFileName());
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

  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    mtuple_files.push_back(mtupleprefix+osc_channel["mtuplefile"].as<std::string>()+mtuplesuffix);
    spline_files.push_back(splineprefix+osc_channel["splinefile"].as<std::string>()+splinesuffix);
    sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
    sample_nutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
    sample_oscnutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
    sample_signal.push_back(osc_channel["signal"].as<bool>());
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
  NSelections = SelectionStr.size();

  return;
}

void samplePDFFDBase::Initialise() {

  //First grab all the information from your sample config via your manager
  ReadSampleConfig();

  //Now initialise all the variables you will need
  Init();

  int TotalMCEvents = 0;
  for(unsigned iSample=0 ; iSample < nSamples ; iSample++){
    std::cout << "=============================================" << std::endl;
    std::cout << "Initialising sample: " << iSample << "/" << nSamples << std::endl;
    MCSamples[iSample].nEvents = setupExperimentMC(iSample);
    std::cout << "Number of events processed: " << MCSamples[iSample].nEvents << std::endl;
    TotalMCEvents += MCSamples[iSample].nEvents;
    std::cout << "Initialising FDMC object.." << std::endl;
    InitialiseSingleFDMCObject(iSample, MCSamples[iSample].nEvents);
    setupFDMC(iSample);
    std::cout << "Initialised sample: " << iSample << "/" << nSamples << std::endl;
  }
  std::cout << "=============================================" << std::endl;
  std::cout << "Total number of events is: " << TotalMCEvents << std::endl;

  std::cout << "Setting up NuOscillator.." << std::endl;
  SetupNuOscillator();
  std::cout << "Setting up Sample Binning.." << std::endl;
  SetupSampleBinning();
  std::cout << "Setting up Splines.." << std::endl;
  SetupSplines();
  std::cout << "Setting up Normalisation Pointers.." << std::endl;
  SetupNormParameters();
  std::cout << "Setting up Weight Pointers.." << std::endl;
  SetupWeightPointers();

  std::cout << "=======================================================" << std::endl;
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
  return;
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
  return;
}

// ************************************************
/// @function samplePDFFDBase::SetupSampleBinning()
/// @brief Function to setup the binning of your sample histograms and the underlying 
/// arrays that get handled in fillArray() and fillArray_MP().
/// The SampleXBins are filled in the daughter class from the sample config file.
/// This "passing" can be removed. 
void samplePDFFDBase::SetupSampleBinning(){
// ************************************************
  std::cout << "Setting up Sample Binning " << std::endl;
  TString histname1d = (XVarStr).c_str();
  TString histname2d = (XVarStr+"_"+YVarStr).c_str();
  TString histtitle = "";

  //The binning here is arbitrary, now we get info from cfg so the
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D   = new TH1D("h"+histname1d+samplename,histtitle, 1, 0, 1);
  dathist   = new TH1D("d"+histname1d+samplename,histtitle, 1, 0, 1);
  _hPDF2D   = new TH2D("h"+histname2d+samplename,histtitle, 1, 0, 1, 1, 0, 1);
  dathist2d = new TH2D("d"+histname2d+samplename,histtitle, 1, 0, 1, 1, 0, 1);

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

  return; 
}

// ************************************************
bool samplePDFFDBase::IsEventSelected(const int iSample, const int iEvent) {
// ************************************************

  double Val;
  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {  
    Val = ReturnKinematicParameter(Selection[iSelection][0], iSample, iEvent);
    //DB If multiple return values, it will consider each value seperately
    //DB Already checked that Selection vector is correctly sized  
	//std::cout << Val << " " << Selection[iSelection][1] << " " << Selection[iSelection][2] << std::endl;	
    //DB In the case where Selection[0].size()==3, Only Events with Val >= Selection[iSelection][1] and Val < Selection[iSelection][2] are considered Passed
    if ((Val<Selection[iSelection][1])||(Val>=Selection[iSelection][2])) {
      return false;
    }
  }

  //DB To avoid unneccessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

// ************************************************
bool samplePDFFDBase::IsEventSelected(const std::vector< std::string >& ParameterStr,
                                      const int iSample, const int iEvent) {
// ************************************************

  double Val;

  for (unsigned int iSelection=0;iSelection<ParameterStr.size();iSelection++) {
 
	Val = ReturnKinematicParameter(ParameterStr[iSelection], iSample, iEvent);
	//ETA - still need to support other method of you specifying the cut you want in ReturnKinematicParameter
	//like in Dan's version below from T2K
    //Val = ReturnKinematicParameter(static_cast<KinematicTypes>(Selection[iSelection][0]),iSample,iEvent);
    //DB If multiple return values, it will consider each value seperately
    //DB Already checked that Selection vector is correctly sized
	//DB In the case where Selection[0].size()==3, Only Events with Val >= Selection[iSelection][1] and Val < Selection[iSelection][2] are considered Passed

    if ((Val<SelectionBounds[iSelection][0])||(Val>=SelectionBounds[iSelection][1])) {
	  return false;
    }
  }

  //DB To avoid unneccessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

// ************************************************
//Same as the function above but just acts on the vector and the event
bool samplePDFFDBase::IsEventSelected(const std::vector< std::string >& ParameterStr,
                                      const std::vector< std::vector<double> > &SelectionCuts,
                                      const int iSample, const int iEvent) {
// ************************************************

  double Val;

  for (unsigned int iSelection=0;iSelection<ParameterStr.size();iSelection++) {
    
    Val = ReturnKinematicParameter(ParameterStr[iSelection], iSample, iEvent);
    //DB If multiple return values, it will consider each value seperately
    //DB Already checked that SelectionCuts vector is correctly sized
    
    //DB In the case where SelectionCuts[0].size()==3, Only Events with Val >= SelectionCuts[iSelection][1] and Val < SelectionCuts[iSelection][2] are considered Passed
    //ETA - also check whether we're actually applying a lower or upper cut by checking they aren't -999
    if(Val >= SelectionCuts[iSelection][1] && SelectionCuts[iSelection][0] != -999){
      return false;
    }
    else if(Val < SelectionCuts[iSelection][0] && SelectionCuts[iSelection][1] != -999){
      return false;
    }
  }
  
  //DB To avoid unneccessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

void samplePDFFDBase::reweight() // Reweight function - Depending on Osc Calculator this function uses different CalcOsc functions
{
  //KS: Reset the histograms before reweight 
  ResetHistograms();

  std::vector<_float_> OscVec(OscCov->GetNumParams());
  for (int iPar=0;iPar<OscCov->GetNumParams();iPar++) {
    OscVec[iPar] = OscCov->getParProp(iPar);
  } 
  for (int iSample=0;iSample<(int)MCSamples.size();iSample++) {
    NuOscProbCalcers[iSample]->CalculateProbabilities(OscVec);
  }

  fillArray();

  return;
}

//************************************************
/// @function samplePDFFDBase::fillArray()
/// Function which does the core reweighting. This assumes that oscillation weights have 
/// already been calculated and stored in samplePDFFDBase[iSample].osc_w[iEvent]. This 
/// function takes advantage of most of the things called in setupSKMC to reduce reweighting time.
/// It also follows the ND code reweighting pretty closely. This function fills the samplePDFFD 
/// array array which is binned to match the sample binning, such that bin[1][1] is the 
/// equivalent of _hPDF2D->GetBinContent(2,2) {Noticing the offset}
//************************************************
void samplePDFFDBase::fillArray() {
  //DB Reset which cuts to apply
  Selection = StoredSelection;

  // Call entirely different routine if we're running with openMP
#ifdef MULTITHREAD
  fillArray_MP();
#else

  //ETA we should probably store this in samplePDFFDBase
  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  //DB Reset values stored in PDF array to 0.
  for (int yBin=0;yBin<nYBins;yBin++) {
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }

  PrepFunctionalParameters();
  if(splineFile){
	splineFile->Evaluate();
  }

  for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
    for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
      applyShifts(iSample, iEvent);
      
      if (!IsEventSelected(iSample, iEvent)) { 
	continue;
      } 

      double splineweight = 1.0;
      double normweight = 1.0;
      double funcweight = 1.0;
      double totalweight = 1.0;
      
      if(splineFile){
	splineweight *= CalcXsecWeightSpline(iSample, iEvent);
      }
      //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient. Do this on a spline-by-spline basis
      if (splineweight <= 0.){
	MCSamples[iSample].xsec_w[iEvent] = 0.;
	continue;
      }
      
      //Loop over stored normalisation and function pointers 
      normweight *= CalcXsecWeightNorm(iSample, iEvent);
      
      //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
      if (normweight <= 0.){
	MCSamples[iSample].xsec_w[iEvent] = 0.;
	continue;
      }
      
      funcweight = CalcXsecWeightFunc(iSample,iEvent);
      //DB Catch negative func weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
      if (funcweight <= 0.){          
	MCSamples[iSample].xsec_w[iEvent] = 0.;
	continue;
      }
      
      MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight*funcweight;
      
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
	for (unsigned int iBin=0;iBin<(XBinEdges.size()-1);iBin++) {
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
  return;
}

#ifdef MULTITHREAD
// ************************************************ 
/// Multithreaded version of fillArray @see fillArray()
// ************************************************ 
void samplePDFFDBase::fillArray_MP() 
{
  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;
  
  //DB Reset values stored in PDF array to 0.
  for (int yBin=0;yBin<nYBins;yBin++) {
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
      samplePDFFD_array_w2[yBin][xBin] = 0.;
    }
  }
  
  //reconfigureFuncPars();
  
  //This is stored as [y][x] due to shifts only occuring in the x variable (Erec/Lep mom) - I believe this will help reduce cache misses 
  double** samplePDFFD_array_private = NULL;
  double** samplePDFFD_array_private_w2 = NULL;
  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
#pragma omp parallel private(samplePDFFD_array_private, samplePDFFD_array_private_w2)
  {
    // private to each thread
    // ETA - maybe we can use parallel firstprivate to initialise these?
    samplePDFFD_array_private = new double*[nYBins];
    samplePDFFD_array_private_w2 = new double*[nYBins];
    for (int yBin=0;yBin<nYBins;yBin++) {
      samplePDFFD_array_private[yBin] = new double[nXBins];
      samplePDFFD_array_private_w2[yBin] = new double[nXBins];
      for (int xBin=0;xBin<nXBins;xBin++) {
	samplePDFFD_array_private[yBin][xBin] = 0.;
        samplePDFFD_array_private_w2[yBin][xBin] = 0.;
      }
    }
    
    //DB - Brain dump of speedup ideas
    //
    //Those relevant to reweighting
    // 1. Don't bother storing and calculating NC signal events - Implemented and saves marginal s/step
    // 2. Loop over spline event weight calculation in the following event loop - Currently done in splineSKBase->calcWeight() where multi-threading won't be optmised - Implemented and saves 0.3s/step
    // 3. Inline getDiscVar or somehow include that calculation inside the multi-threading - Implemented and saves about 0.01s/step
    // 4. Include isCC inside SKMCStruct so don't have to have several 'if' statements determine if oscillation weight needs to be set to 1.0 for NC events - Implemented and saves marginal s/step
    // 5. Do explict check on adjacent bins when finding event XBin instead of looping over all BinEdge indicies - Implemented but doesn't significantly affect s/step
    //
    //Other aspects
    // 1. Order minituples in Y-axis variable as this will *hopefully* reduce cache misses inside samplePDFFD_array_class[yBin][xBin]
    //
    // We will hit <0.1 s/step eventually! :D
    
    //ETA - does these three calls need to be inside the omp parrallel region? 
    //I don't think this will affect anything but maybe should check.
    PrepFunctionalParameters();
    //==================================================
    //Calc Weights and fill Array
    if(splineFile){
      splineFile->Evaluate();
    }
    
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
	
	double splineweight = 1.0;
	double normweight = 1.0;
	double funcweight = 1.0;
	double totalweight = 1.0;
	
	//DB SKDet Syst
	//As weights were skdet::fParProp, and we use the non-shifted erec, we might as well cache the corresponding fParProp index for each event and the pointer to it
	
	if(splineFile){
	  splineweight *= CalcXsecWeightSpline(iSample, iEvent);
	}
	//DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
	if (splineweight <= 0.){
	  MCSamples[iSample].xsec_w[iEvent] = 0.;
	  continue;
	}
	
        normweight *= CalcXsecWeightNorm(iSample, iEvent);
	//DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
	if (normweight <= 0.){
	  MCSamples[iSample].xsec_w[iEvent] = 0.;
	  continue;
	}
	
	funcweight = CalcXsecWeightFunc(iSample,iEvent);
	//DB Catch negative func weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
	if (funcweight <= 0.){
	  MCSamples[iSample].xsec_w[iEvent] = 0.;
	  continue;
	}
	
	MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight*funcweight;
	
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
    for (int yBin=0;yBin<nYBins;yBin++) {
      for (int xBin=0;xBin<nXBins;xBin++) {
#pragma omp atomic
	samplePDFFD_array[yBin][xBin] += samplePDFFD_array_private[yBin][xBin];
#pragma omp atomic    
	samplePDFFD_array_w2[yBin][xBin] += samplePDFFD_array_private_w2[yBin][xBin];
      }
    }
    
    for (int yBin=0;yBin<nYBins;yBin++) {
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
  
  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;
  
  //DB Reset values stored in PDF array to 0.
  for (int yBin = 0; yBin < nYBins; yBin++) {
    for (int xBin = 0; xBin < nXBins; xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
    }
  }
} // end function

// ***************************************************************************
// Calculate the spline weight for one event
double samplePDFFDBase::CalcXsecWeightSpline(const int iSample, const int iEvent) {
// ***************************************************************************

  double xsecw = 1.0;
  //DB Xsec syst
  //Loop over stored spline pointers
  for (int iSpline=0;iSpline<MCSamples[iSample].nxsec_spline_pointers[iEvent];iSpline++) {
    xsecw *= *(MCSamples[iSample].xsec_spline_pointers[iEvent][iSpline]);
  }
  return xsecw;
}

// ***************************************************************************
// Calculate the normalisation weight for one event
double samplePDFFDBase::CalcXsecWeightNorm(const int iSample, const int iEvent) {
// ***************************************************************************

  double xsecw = 1.0;
  //Loop over stored normalisation and function pointers
  for (int iParam = 0;iParam < MCSamples[iSample].nxsec_norm_pointers[iEvent]; iParam++)
  {
    xsecw *= *(MCSamples[iSample].xsec_norm_pointers[iEvent][iParam]);
    #ifdef DEBUG
    if (TMath::IsNaN(xsecw)) std::cout << "iParam=" << iParam << "xsecweight=nan from norms" << std::endl;
    #endif
  }
  return xsecw;
}

void samplePDFFDBase::SetXsecCov(covarianceXsec *xsec){

  MACH3LOG_INFO("SETTING UP XSEC COV!!");
  XsecCov = xsec;

  // Get the map between the normalisation parameters index, their name, what mode they should apply to, and what target
  //This information is used later in CalcXsecNormsBins to decide if a parameter applies to an event

  //DB Now get this information using the DetID from the config
  xsec_norms = XsecCov->GetNormParsFromDetID(SampleDetID);
  nFuncParams = XsecCov->GetNumFuncParamsFromDetID(SampleDetID);
  funcParsNames = XsecCov->GetFuncParsNamesFromDetID(SampleDetID);
  funcParsIndex = XsecCov->GetFuncParsIndexFromDetID(SampleDetID);

  MACH3LOG_INFO("Found {} normalisation parameters", xsec_norms.size());
  MACH3LOG_INFO("Found {} functional parameters", funcParsNames.size());
	
  return;
}

void samplePDFFDBase::SetOscCov(covarianceOsc* osc_cov){
  OscCov = osc_cov;
  int nOscPars = OscCov->GetNumParams();
  oscpars = new const double*[nOscPars];
  for(auto osc_par_i = 0; osc_par_i < nOscPars ; ++osc_par_i){
    oscpars[osc_par_i] = OscCov->retPointer(osc_par_i);
  } 
  return;
}

void samplePDFFDBase::SetupNormParameters(){

  if(!XsecCov){
	MACH3LOG_ERROR("XsecCov is not setup!");
	throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Assign xsec norm bins in MCSamples tree
  for (int iSample = 0; iSample < (int)MCSamples.size(); ++iSample) {
	CalcXsecNormsBins(iSample);
  }

  //DB
  //Attempt at reducing impact of covarianceXsec::calcReweight()
  int counter;

  for (int iSample = 0; iSample < (int)MCSamples.size(); ++iSample) {
	for (int iEvent = 0; iEvent < MCSamples[iSample].nEvents; ++iEvent) {
	  counter = 0;

	  MCSamples[iSample].nxsec_norm_pointers[iEvent] = MCSamples[iSample].xsec_norms_bins[iEvent].size();
	  MCSamples[iSample].xsec_norm_pointers[iEvent] = new const double*[MCSamples[iSample].nxsec_norm_pointers[iEvent]];

	  for(std::list< int >::iterator lit = MCSamples[iSample].xsec_norms_bins[iEvent].begin();lit!=MCSamples[iSample].xsec_norms_bins[iEvent].end();lit++) {
		MCSamples[iSample].xsec_norm_pointers[iEvent][counter] = XsecCov->retPointer(*lit);
		counter += 1;
	  }

	}
  }

  return;
}

//A way to check whether a normalisation parameter applies to an event or not
void samplePDFFDBase::CalcXsecNormsBins(int iSample){

  fdmc_base *fdobj = &MCSamples[iSample];

  for(int iEvent=0; iEvent < fdobj->nEvents; ++iEvent){
    std::list< int > XsecBins = {};
	if (XsecCov) {
	  for (std::vector<XsecNorms4>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it) {
		// Skip oscillated NC events
		// Not strictly needed, but these events don't get included in oscillated predictions, so
		// no need to waste our time calculating and storing information about xsec parameters
		// that will never be used.
		if (fdobj->isNC[iEvent] && fdobj->signal) {continue;} //DB Abstract check on MaCh3Modes to determine which apply to neutral current

		//Now check that the target of an interaction matches with the normalisation parameters
		bool TargetMatch=false;
		//If no target specified then apply to all modes
		if ((*it).targets.size()==0) {
		  TargetMatch=true;
		} else {
		  for (unsigned iTarget=0;iTarget<(*it).targets.size();iTarget++) {
			if ((*it).targets.at(iTarget)== *(fdobj->Target[iEvent])) {
			  TargetMatch=true;
			}
		  }
		}
		if (!TargetMatch) {continue;}

		//Now check that the neutrino flavour in an interaction matches with the normalisation parameters
		bool FlavourMatch=false;
		//If no mode specified then apply to all modes
		if ((*it).pdgs.size()==0) {
		  FlavourMatch=true;
		} else {
		  for (unsigned iPDG=0;iPDG<(*it).pdgs.size();iPDG++) {
			if ((*it).pdgs.at(iPDG)== fdobj->nupdg) {
			  FlavourMatch=true;
			}
		  }
		}
		if (!FlavourMatch){continue;}

		//Now check that the unoscillated neutrino flavour in an interaction matches with the normalisation parameters
		bool FlavourUnoscMatch=false;
		//If no mode specified then apply to all modes
		if ((*it).preoscpdgs.size()==0) {
		  FlavourUnoscMatch=true;
		} else {
		  for (unsigned iPDG=0;iPDG<(*it).preoscpdgs.size();iPDG++) {
			if ((*it).preoscpdgs.at(iPDG) == fdobj->nupdgUnosc) {
			  FlavourUnoscMatch=true;
			}
		  }
		}
		if (!FlavourUnoscMatch){continue;}
		
		//Now check that the mode of an interaction matches with the normalisation parameters
		bool ModeMatch=false;
		//If no mode specified then apply to all modes
		if ((*it).modes.size()==0) {
		  ModeMatch=true;
		} else {
		  for (unsigned imode=0;imode<(*it).modes.size();imode++) {
			if ((*it).modes.at(imode)== *(fdobj->mode[iEvent])) {
			  ModeMatch=true;
			}
		  }
		}
		if (!ModeMatch) {continue;}

		//Now check whether the norm has kinematic bounds
		//i.e. does it only apply to events in a particular kinematic region?
		bool IsSelected = true;
		if ((*it).hasKinBounds) {
		  for (unsigned int iKinematicParameter = 0 ; iKinematicParameter < (*it).KinematicVarStr.size() ; ++iKinematicParameter ) {
			if (ReturnKinematicParameter((*it).KinematicVarStr[iKinematicParameter], iSample, iEvent) <= (*it).Selection[iKinematicParameter][0]) { 
			  IsSelected = false;
			  continue;
			}
			else if (ReturnKinematicParameter((*it).KinematicVarStr[iKinematicParameter], iSample, iEvent) > (*it).Selection[iKinematicParameter][1]) {
			  IsSelected = false;
			  continue;
			}
		  } 
		}
		//Need to then break the event loop 
		if(!IsSelected){
		  continue;
		}
	
		// Now set 'index bin' for each normalisation parameter
		// All normalisations are just 1 bin for 2015, so bin = index (where index is just the bin for that normalisation)
		int bin = (*it).index;

		//If syst on applies to a particular detector
		if ((XsecCov->GetParDetID(bin) & SampleDetID)==SampleDetID) {
		  XsecBins.push_back(bin);
		}
	  } // end iteration over xsec_norms
	} // end if (xsecCov)
	fdobj->xsec_norms_bins[iEvent]=XsecBins;
  }//end loop over events
  return;
}

//ETA - this is all a bit (less) stupid
void samplePDFFDBase::set1DBinning(std::vector<double> &XVec){

  _hPDF1D->Reset();
  _hPDF1D->SetBins(XVec.size()-1, XVec.data());
  dathist->SetBins(XVec.size()-1, XVec.data());

  //This will overwrite XBinEdges with whatever you pass this function
  XBinEdges = XVec;
  YBinEdges = std::vector<double>(2);
  YBinEdges[0] = -1e8;
  YBinEdges[1] = 1e8;

  _hPDF2D->Reset();
  _hPDF2D  ->SetBins(XVec.size()-1, XVec.data(), YBinEdges.size()-1, YBinEdges.data());
  dathist2d->SetBins(XVec.size()-1, XVec.data(), YBinEdges.size()-1, YBinEdges.data());

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

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
  _hPDF1D->SetBins(XVec.size()-1, XVec.data());
  dathist->SetBins(XVec.size()-1, XVec.data());

  _hPDF2D->Reset();
  _hPDF2D->SetBins(XVec.size()-1, XVec.data(), YVec.size()-1, YVec.data());
  dathist2d->SetBins(XVec.size()-1, XVec.data(), YVec.size()-1, YVec.data());

  //XBinEdges = XVec;
  //YBinEdges = YVec;

  //ETA - maybe need to be careful here
  int nXBins = XVec.size()-1;
  int nYBins = YVec.size()-1;

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
//New versions of set binning funcitons is samplePDFBase
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

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

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

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

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
  for(int mc_i = 0 ; mc_i < (int)MCSamples.size() ; mc_i++){
	for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){

	  //Set x_var and y_var values based on XVarStr and YVarStr
      MCSamples[mc_i].x_var[event_i] = ReturnKinematicParameterByReference(XVarStr, mc_i, event_i);
	  //Give y)_var a dummy value
      MCSamples[mc_i].y_var[event_i] = &(MCSamples[mc_i].dummy_value);
	  int bin = _hPDF1D->FindBin(*(MCSamples[mc_i].x_var[event_i]));

	  double low_lower_edge = _DEFAULT_RETURN_VAL_;
	  if (bin==0) {
		low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
	  } else {
		low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin-1);
	  }

	  double low_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
	  double upper_edge = _hPDF1D->GetXaxis()->GetBinUpEdge(bin);

	  double upper_upper_edge = _DEFAULT_RETURN_VAL_;
	  if (bin<(_hPDF1D->GetNbinsX()-2)) {
		upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+2);
	  } else {
		upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+1);
	  }

	  if ((bin-1) >= 0 && (bin-1) < int(XBinEdges.size()-1)) {
		MCSamples[mc_i].NomXBin[event_i] = bin-1;
	  } else {
		MCSamples[mc_i].NomXBin[event_i] = -1;
		low_edge = _DEFAULT_RETURN_VAL_;
		upper_edge = _DEFAULT_RETURN_VAL_;
		low_lower_edge = _DEFAULT_RETURN_VAL_;
		upper_upper_edge = _DEFAULT_RETURN_VAL_;
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

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

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

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

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

void samplePDFFDBase::FindNominalBinAndEdges2D() {

  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(int mc_i = 0 ; mc_i < (int)MCSamples.size() ; mc_i++){
	for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){

   	  //Set x_var and y_var values based on XVarStr and YVarStr   
      MCSamples[mc_i].x_var[event_i] = ReturnKinematicParameterByReference(XVarStr, mc_i, event_i);
      MCSamples[mc_i].y_var[event_i] = ReturnKinematicParameterByReference(YVarStr, mc_i, event_i);

	  //Global bin number
	  int bin = _hPDF2D->FindBin(*(MCSamples[mc_i].x_var[event_i]), *(MCSamples[mc_i].y_var[event_i]));

	  int bin_x = -999;
	  int bin_y = -999;
	  int bin_z = -999;
	  _hPDF2D->GetBinXYZ(bin, bin_x, bin_y, bin_z);
	  //erec is the x-axis so get GetXaxis then find the bin edges using the x bin number

	  double low_lower_edge = _DEFAULT_RETURN_VAL_;
	  if (bin==0) {
		low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
	  } else {
		low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x-1);
	  }

	  double low_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
	  double upper_edge = _hPDF2D->GetXaxis()->GetBinUpEdge(bin_x);

	  double upper_upper_edge = _DEFAULT_RETURN_VAL_;
	  if (bin<(_hPDF2D->GetNbinsX()-2)) {
		upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+2);
	  } else {
		upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+1);
	  }

	  if ((bin_x-1) >= 0 && (bin_x-1) < int(XBinEdges.size()-1)) {
		MCSamples[mc_i].NomXBin[event_i] = bin_x-1;
	  } else {
		MCSamples[mc_i].NomXBin[event_i] = -1;
		low_edge = _DEFAULT_RETURN_VAL_;
		upper_edge = _DEFAULT_RETURN_VAL_;
		low_lower_edge = _DEFAULT_RETURN_VAL_;
		upper_upper_edge = _DEFAULT_RETURN_VAL_;
	  }
	  MCSamples[mc_i].NomYBin[event_i] = bin_y-1; 
	  if(MCSamples[mc_i].NomYBin[event_i] < 0){ 
		std::cout << "Nominal YBin PROBLEM, y-bin is " << MCSamples[mc_i].NomYBin[event_i] << std::endl;
	  }
	  MCSamples[mc_i].rw_lower_xbinedge[event_i] = low_edge;
	  MCSamples[mc_i].rw_upper_xbinedge[event_i] = upper_edge;
	  MCSamples[mc_i].rw_lower_lower_xbinedge[event_i] = low_lower_edge;
	  MCSamples[mc_i].rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
	}
  }
  return;
}

void samplePDFFDBase::addData(std::vector<double> &data) {
  dataSample = new std::vector<double>(data);
  dataSample2D = NULL;
  dathist2d = NULL;
  dathist->Reset(); 

  if (GetNDim()!=1) {
	MACH3LOG_ERROR("Trying to set a 1D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
	MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
	throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i = 0; i < int(dataSample->size()); i++) {
    dathist->Fill(dataSample->at(i));
  }

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist->GetBinContent(xBin+1);
    }
  }

  return;
}

void samplePDFFDBase::addData(std::vector< std::vector <double> > &data) {
  dataSample2D = new std::vector< std::vector <double> >(data);
  dataSample = NULL;
  dathist = NULL;
  dathist2d->Reset();                                                       

  if (GetNDim()!=2) {
    MACH3LOG_ERROR("Trying to set a 2D 'data' histogram when the number of dimensions for this sample is {}", GetNDim());
    MACH3LOG_ERROR("This won't work, please specify the correct dimensions in your sample config with the X and Y variables");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (int i = 0; i < int(dataSample2D->size()); i++) {
    dathist2d->Fill(dataSample2D->at(0)[i],dataSample2D->at(1)[i]);
  }

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist2d->GetBinContent(xBin+1,yBin+1);
    }
  }

  return;
}

void samplePDFFDBase::addData(TH1D* Data) {
  MACH3LOG_INFO("Adding 1D data histogram : {} with {} events", Data->GetName(), Data->Integral());
  dathist2d = NULL;
  dathist = Data;
  dataSample = NULL;
  dataSample2D = NULL;

  if (GetNDim()!=1) {
	  MACH3LOG_ERROR("Trying to set a 1D 'data' histogram in a 2D sample - Quitting"); 
	  throw MaCh3Exception(__FILE__ , __LINE__ );}

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = Data->GetBinContent(xBin+1);
    }
  }

  return;
}

void samplePDFFDBase::addData(TH2D* Data) {
  MACH3LOG_INFO("Adding 2D data histogram : {} with {} events", Data->GetName(), Data->Integral());
  dathist2d = Data;
  dathist = NULL;
  dataSample = NULL;
  dataSample2D = NULL;

  if (GetNDim()!=2) {
	  MACH3LOG_ERROR("Trying to set a 2D 'data' histogram in a 1D sample - Quitting"); 
	  throw MaCh3Exception(__FILE__ , __LINE__ );}	

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = dathist2d->GetBinContent(xBin+1,yBin+1);
    }
  }
}

void samplePDFFDBase::SetupNuOscillator() {
  OscillatorFactory* OscillFactory = new OscillatorFactory();
  
  NuOscProbCalcers = std::vector<OscillatorBase*>((int)MCSamples.size());
  for (int iSample=0;iSample<(int)MCSamples.size();iSample++) {
    std::cout << "Setting up NuOscillator::Oscillator object in OscillationChannel: " << iSample << "/" << MCSamples.size() << " ====================" << std::endl;
    NuOscProbCalcers[iSample] = OscillFactory->CreateOscillator(NuOscillatorConfigFile);

    std::vector<_float_> EnergyArray((int)MCSamples[iSample].nEvents);
    for (int iEvent=0;iEvent<(int)MCSamples[iSample].nEvents;iEvent++) {
      EnergyArray[iEvent] = *(MCSamples[iSample].rw_etru[iEvent]);
    }
    std::sort(EnergyArray.begin(),EnergyArray.end());

    if (!NuOscProbCalcers[iSample]->EvalPointsSetInConstructor()) {
	  std::cout << "Made it through EvalPointsSetInConstructor check" << std::endl;
      NuOscProbCalcers[iSample]->SetEnergyArrayInCalcer(EnergyArray);

      //============================================================================
      //DB Atmospheric only part
      if (MCSamples[iSample].rw_truecz != NULL) { //Can only happen if truecz has been initialised within the experiment specific code
	std::vector<_float_> CosineZArray((int)MCSamples[iSample].nEvents);
	for (int iEvent=0;iEvent<(int)MCSamples[iSample].nEvents;iEvent++) {
	  CosineZArray[iEvent] = *(MCSamples[iSample].rw_truecz[iEvent]);
	}
	std::sort(CosineZArray.begin(),CosineZArray.end());
	
	NuOscProbCalcers[iSample]->SetCosineZArrayInCalcer(CosineZArray);
      }
      //============================================================================
      
    }
    
    NuOscProbCalcers[iSample]->Setup();

    for (int iEvent=0;iEvent<(int)MCSamples[iSample].nEvents;iEvent++) {
      MCSamples[iSample].osc_w_pointer[iEvent] = &Unity;
      if (MCSamples[iSample].isNC[iEvent]) {
	if (MCSamples[iSample].signal) {
	  MCSamples[iSample].osc_w_pointer[iEvent] = &Zero;
	} else {
	  MCSamples[iSample].osc_w_pointer[iEvent] = &Unity;
	}
      } else {
	
	int InitFlav = NULL;
	int FinalFlav = NULL;
	
	if (std::abs(MCSamples[iSample].nutype) == 1) {
	  InitFlav = NuOscillator::kElectron;
	} else if (std::abs(MCSamples[iSample].nutype) == 2) {
	  InitFlav = NuOscillator::kMuon;
	} else if (std::abs(MCSamples[iSample].nutype) == 3) {
	  InitFlav = NuOscillator::kTau;
	}

	if (std::abs(MCSamples[iSample].oscnutype) == 1) {
	  FinalFlav = NuOscillator::kElectron;
	} else if (std::abs(MCSamples[iSample].oscnutype) == 2) {
	  FinalFlav = NuOscillator::kMuon;
	} else if (std::abs(MCSamples[iSample].oscnutype) == 3) {
	  FinalFlav = NuOscillator::kTau;
	}

	if (InitFlav == NULL || FinalFlav == NULL) {
	  std::cerr << "Something has gone wrong in the mapping between MCSamples[iSample].nutype and the enum used within NuOscillator" << std::endl;
	  std::cerr << "MCSamples[iSample].nutype:" << MCSamples[iSample].nutype << std::endl;
	  std::cerr << "InitFlav:" << InitFlav << std::endl;
	  std::cerr << "MCSamples[iSample].oscnutype:" << MCSamples[iSample].oscnutype << std::endl;
	  std::cerr << "FinalFlav:" << FinalFlav << std::endl;
	  throw;
	}

	//Assuming that if the generated neutrino is antineutrino, the detected flavour will also be antineutrino
	if (MCSamples[iSample].nutype<0) {
	  InitFlav *= -1;
	  FinalFlav *= -1;
	}
	
	if (MCSamples[iSample].rw_truecz != NULL) { //Can only happen if truecz has been initialised within the experiment specific code
	  //Atmospherics
	  MCSamples[iSample].osc_w_pointer[iEvent] = NuOscProbCalcers[iSample]->ReturnWeightPointer(InitFlav,FinalFlav,*(MCSamples[iSample].rw_etru[iEvent]),*(MCSamples[iSample].rw_truecz[iEvent]));
	} else {
	  //Beam
	  MCSamples[iSample].osc_w_pointer[iEvent] = NuOscProbCalcers[iSample]->ReturnWeightPointer(InitFlav,FinalFlav,*(MCSamples[iSample].rw_etru[iEvent]));
	}
	
      }
    }
    
  }
}

double samplePDFFDBase::GetEventWeight(int iSample, int iEntry) {
  double totalweight = 1.0;
  for (int iParam=0;iParam<MCSamples[iSample].ntotal_weight_pointers[iEntry];iParam++) {
    totalweight *= *(MCSamples[iSample].total_weight_pointers[iEntry][iParam]);
  }
  return totalweight;
}

/// @func fillSplineBins()
/// @brief Finds the binned spline that an event should apply to and stored them in a
/// a vector for easy evaluation in the fillArray() function.
void samplePDFFDBase::fillSplineBins() {
  for (int i = 0; i < (int)MCSamples.size(); ++i) {
    //Now loop over events and get the spline bin for each event
    for (int j = 0; j < MCSamples[i].nEvents; ++j) {
      std::vector< std::vector<int> > EventSplines;
      switch(nDimensions){
      case 1:
	EventSplines = splineFile->GetEventSplines(GetName(), i, *(MCSamples[i].mode[j]), *(MCSamples[i].rw_etru[j]), *(MCSamples[i].x_var[j]), 0.);
	break;
      case 2:
	EventSplines = splineFile->GetEventSplines(GetName(), i, *(MCSamples[i].mode[j]), *(MCSamples[i].rw_etru[j]), *(MCSamples[i].x_var[j]), *(MCSamples[i].y_var[j]));
	break;
      default:
	MACH3LOG_ERROR("Error in assigning spline bins because nDimensions = {}", nDimensions);
	MACH3LOG_ERROR("MaCh3 only supports splines binned in Etrue + the sample binning");
	MACH3LOG_ERROR("Please check the sample binning you specified in your sample config ");
	break;
      }
      MCSamples[i].nxsec_spline_pointers[j] = EventSplines.size();
	  if(MCSamples[i].nxsec_spline_pointers[j] < 0){
		throw MaCh3Exception(__FILE__, __LINE__);
	  }
      MCSamples[i].xsec_spline_pointers[j] = new const double*[MCSamples[i].nxsec_spline_pointers[j]];
      for(int spline=0; spline<MCSamples[i].nxsec_spline_pointers[j]; spline++){          
	//Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
	MCSamples[i].xsec_spline_pointers[j][spline] = splineFile->retPointer(EventSplines[spline][0], EventSplines[spline][1], EventSplines[spline][2], 
									      EventSplines[spline][3], EventSplines[spline][4], EventSplines[spline][5], EventSplines[spline][6]);
      }
      
    }
  }

  return;
}

double samplePDFFDBase::GetLikelihood()
{
  if (samplePDFFD_data == NULL) {

	MACH3LOG_ERROR("Data sample is empty! Can't calculate a likelihood!");
	throw MaCh3Exception(__FILE__, __LINE__);
  }

  //This can be done only once and stored
  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  int xBin;
  int yBin;

  double negLogL = 0.;

#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:negLogL) private(xBin, yBin)
#endif
  for (xBin = 0; xBin < nXBins; xBin++) 
  {
    for (yBin = 0; yBin < nYBins; yBin++) 
    {
        double DataVal = samplePDFFD_data[yBin][xBin];
        double MCPred = samplePDFFD_array[yBin][xBin];
        double w2 = samplePDFFD_array_w2[yBin][xBin];

        //KS: Calcaualte likelihood using Barlow-Beestion Poisson or even IceCube
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
  
  MCSamples[iSample] = fdmc_base();
  fdmc_base *fdobj = &MCSamples[iSample];
  
  fdobj->nEvents = nEvents_;
  fdobj->nutype = -9;
  fdobj->oscnutype = -9;
  fdobj->signal = false;
  fdobj->Unity = 1.;
  fdobj->Unity_Int = 1.;
  
  fdobj->x_var = new const double*[fdobj->nEvents];
  fdobj->y_var = new const double*[fdobj->nEvents];
  fdobj->rw_etru = new const double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
  fdobj->XBin = new int [fdobj->nEvents];
  fdobj->YBin = new int [fdobj->nEvents];;	 
  fdobj->rw_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_lower_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->mode = new int*[fdobj->nEvents];
  fdobj->nxsec_norm_pointers = new int[fdobj->nEvents];
  fdobj->xsec_norm_pointers = new const double**[fdobj->nEvents];
  fdobj->xsec_norms_bins = new std::list< int >[fdobj->nEvents];
  fdobj->xsec_w = new double[fdobj->nEvents];
  fdobj->isNC = new bool[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new const double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];
  fdobj->osc_w_pointer = new const double*[fdobj->nEvents];
  //fdobj->rw_truecz = new const double*[fdobj->nEvents];
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &fdobj->Unity;
    //fdobj->rw_truecz[iEvent] = &fdobj->Unity;
    fdobj->mode[iEvent] = &fdobj->Unity_Int;
    fdobj->Target[iEvent] = 0;
    fdobj->NomXBin[iEvent] = -1;
    fdobj->NomYBin[iEvent] = -1;
    fdobj->XBin[iEvent] = -1;
    fdobj->YBin[iEvent] = -1;	 
    fdobj->rw_lower_xbinedge[iEvent] = -1;
    fdobj->rw_lower_lower_xbinedge[iEvent] = -1;
    fdobj->rw_upper_xbinedge[iEvent] = -1;
    fdobj->rw_upper_upper_xbinedge[iEvent] = -1;
    fdobj->xsec_w[iEvent] = 1.0;
    fdobj->isNC[iEvent] = false;
    fdobj->SampleDetID = -1;
    fdobj->osc_w_pointer[iEvent] = &(fdobj->Unity);
    fdobj->x_var[iEvent] = &fdobj->Unity;
    fdobj->y_var[iEvent] = &fdobj->Unity;
  }

}

void samplePDFFDBase::InitialiseSplineObject() {
  std::vector<std::string> spline_filepaths;
  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    spline_filepaths.push_back(splineprefix+spline_files[iSample]+splinesuffix);
  }

  //Keep a track of the spline variables
  SplineVarNames.push_back("TrueNeutrinoEnergy");
  if(XVarStr.length() > 0){
	SplineVarNames.push_back(XVarStr);
  }
  if(YVarStr.length() > 0){
	SplineVarNames.push_back(YVarStr);
  }
  
  splineFile->AddSample(samplename, SampleDetID, spline_filepaths, SplineVarNames);
  splineFile->PrintArrayDimension();
  splineFile->CountNumberOfLoadedSplines(false, 1);
  splineFile->TransferToMonolith();

  std::cout << "--------------------------------" <<std::endl;
  MACH3LOG_INFO("Setup Far Detector splines");

  fillSplineBins();
}
