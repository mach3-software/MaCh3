#include "samplePDFFDBase.h"

// Constructors for erec-binned errors

samplePDFFDBase::samplePDFFDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov)
  : samplePDFBase()
//DB Throughout constructor and init, pot is livetime for atmospheric samples
{
  (void) pot;
  (void) mc_version;
  MACH3LOG_INFO("-------------------------------------------------------------------");	
  MACH3LOG_INFO("Creating samplePDFFDBase object..");
	
  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){
	  MACH3LOG_ERROR("[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!");
	  throw MaCh3Exception(__FILE__ , __LINE__ );}

  samplePDFFD_array = nullptr;
  samplePDFFD_data = nullptr;
  oscpars = nullptr;

  //KS: For now FD support only one sample
  nSamples = 1;
  SampleName.push_back("FDsample");

  //Default TestStatistic is kPoisson
  //ETA: this can be configured with samplePDFBase::SetTestStatistic()
  fTestStatistic = kPoisson;

  //Default values for oscillation-related things
  doubled_angle = true;
  osc_binned = false;
  Osc = NULL;
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
  TString histname1d = (XVarStr).c_str();
  TString histname2d = (XVarStr+YVarStr).c_str();
  TString histtitle = "";

  //The binning here is arbitrary, now we get info from cfg so the
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D   = new TH1D("h"+histname1d+samplename,histtitle, 1, 0, 1);
  dathist   = new TH1D("d"+histname1d+samplename,histtitle, 1, 0, 1);
  _hPDF2D   = new TH2D("h"+histname2d+samplename,histtitle, 1, 0, 1, 1, 0, 1);
  dathist2d = new TH2D("d"+histname2d+samplename,histtitle, 1, 0, 1, 1, 0, 1);

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  XBinEdges.reserve(SampleNXBins);
  YBinEdges.reserve(SampleNYBins);
  std::cout << "XBinning: " << std::endl;
  for(unsigned XBin_i = 0 ; XBin_i < SampleNXBins ; XBin_i++){
	XBinEdges.push_back(SampleXBins[XBin_i]);
	std::cout << SampleXBins[XBin_i] << ", ";
  }
  std::cout << "\n" << std::endl;

  //And now the YBin Edges
  std::cout << "YBinning: " << std::endl;
  for(unsigned YBin_i = 0 ; YBin_i < SampleNYBins ; YBin_i++){
	YBinEdges.push_back(SampleYBins[YBin_i]);
	std::cout << SampleYBins[YBin_i] << ", ";
  }
  std::cout << "\n" << std::endl;
 
  if(XVarStr.length() > 0 && YVarStr.length() == 0){
	set1DBinning(SampleXBins);  
  }
  else if(XVarStr.length() > 0 && YVarStr.length() > 0){
	MACH3LOG_INFO("Setting Up 2D binning");
	MACH3LOG_INFO("{} : {}", XVarStr, YVarStr);  
	set2DBinning(SampleXBins, SampleYBins);
  }
}

void samplePDFFDBase::UseBinnedOscReweighting(bool ans) 
{
  osc_binned = ans;
  if(ans == true) {
    std::cout << "WARNING: you are using binned oscillation weight without specifying the binning. It will use E_reco binning to set the energy (recommended to use set1Dbinning first !). If you want another binning use : \n useBinnedOscReweighting(true, int, double*) \nwhere int is the number of bins and double* is an array with the bins boundaries." << std::endl ;

    /// Get the binning from the MC histogram
    const int nb_bins = _hPDF1D -> GetXaxis() -> GetNbins() ;      
    const double* osc_bins = _hPDF1D -> GetXaxis() -> GetXbins() -> GetArray();
    osc_binned_axis = new TAxis(nb_bins, osc_bins) ;
  }

  return;
}

void samplePDFFDBase::UseBinnedOscReweighting(bool ans, int nbins, double *osc_bins) 
{
  osc_binned = ans;
  if(ans == true) {
    osc_binned_axis = new TAxis(nbins, osc_bins) ;
  }

  return;
}

// ************************************************
bool samplePDFFDBase::IsEventSelected(const int iSample, const int iEvent) {
// ************************************************

  double Val;

  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {
    
    Val = ReturnKinematicParameter(Selection[iSelection][0], iSample, iEvent);
    //std::cout << "Val returned for selection " << iSelection << " is " << Val << std::endl;
    //DB If multiple return values, it will consider each value seperately
    //DB Already checked that Selection vector is correctly sized
    
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
	  //std::cout << "Cutting event as " << Val << " is greater than " << SelectionCuts[iSelection][1]
	  return false;
	}
	else if(Val < SelectionCuts[iSelection][0] && SelectionCuts[iSelection][1] != -999){
	  return false;
	}
  }

  //DB To avoid unneccessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

//CalcOsc for Prob3++ CPU
#if defined (USE_PROB3) && defined (CPU_ONLY)
double samplePDFFDBase::calcOscWeights(int sample, int nutype, int oscnutype, double en)
{
  MCSamples[sample].Oscillator->SetMNS(*oscpars[0], *oscpars[2], *oscpars[1], *oscpars[3], *oscpars[4], *oscpars[5], en, doubled_angle, nutype);
  MCSamples[sample].Oscillator->propagateLinear(nutype , *oscpars[6], *oscpars[7]);

  return MCSamples[sample].Oscillator->GetProb(nutype, oscnutype);
}
#endif

//CalcOsc for Prob3++ GPU (ProbGpu)
#if defined (USE_PROB3) && not defined (CPU_ONLY)
extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw); 

void samplePDFFDBase::calcOscWeights(int nutype, int oscnutype, double *en, double *w, int num)
{
  setMNS(*oscpars[0], *oscpars[2], *oscpars[1], *oscpars[3], *oscpars[4], *oscpars[5], doubled_angle);
  GetProb(nutype, oscnutype, *oscpar[6], *oscpar[7], en, num, w);

  if (std::isnan(w[10]))
  {
    MACH3LOG_ERROR("WARNING: ProbGPU oscillation weight returned NaN! {}", w[10]);
  }
}
#endif

//CalcOsc for CUDAProb3 CPU/GPU
#if not defined (USE_PROB3)
void samplePDFFDBase::calcOscWeights(int sample, int nutype, double *w)
{
  MCSamples[sample].Oscillator->setMNSMatrix(asin(sqrt(*oscpars[0])),asin(sqrt(*oscpars[2])), asin(sqrt(*oscpars[1])), (*oscpars[5]), nutype);
  MCSamples[sample].Oscillator->setNeutrinoMasses(*oscpars[3], *oscpars[4]);
  MCSamples[sample].Oscillator->calculateProbabilities(MCSamples[sample].NeutrinoType);
  MCSamples[sample].Oscillator->getProbabilityArr(w, MCSamples[sample].ProbType);
}
#endif 

void samplePDFFDBase::reweight() // Reweight function - Depending on Osc Calculator this function uses different CalcOsc functions
{

  if (Osc!=NULL) {
	std::cout << "Osc is not NULL!! i.e. doing atm oscillations " << std::endl;
    //DB Currently hardcoded to assume rho_electrons = rho_matter/2, 25km production height
    Osc->FillOscillogram(oscpars,25.0,0.5);
    for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
      for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
		MCSamples[iSample].osc_w[iEvent] = *(MCSamples[iSample].osc_w_pointer[iEvent]);
      }
    }
  } else {
    for (int i=0; i< (int)MCSamples.size(); ++i) {
      
#if defined (USE_PROB3) && defined (CPU_ONLY)
      //Prob3 CPU needs to loop through events too
      for(int j = 0; j < MCSamples[i].nEvents; ++j) {
		MCSamples[i].osc_w[j] = calcOscWeights(i, MCSamples[i].nutype, MCSamples[i].oscnutype, *(MCSamples[i].rw_etru[j]));
      } //event loop
#endif
      
#if defined (USE_PROB3) && not defined (CPU_ONLY)
      calcOscWeights(MCSamples[i].nutype, MCSamples[i].oscnutype, *(MCSamples[i].rw_etru), MCSamples[i].osc_w, MCSamples[i].nEvents);
#endif
      
#if not defined (USE_PROB3)
      calcOscWeights(i, MCSamples[i].nutype, MCSamples[i].osc_w);
#endif
    }// Sample loop
  }

  //KS: Reset the histograms before reweight 
  ResetHistograms();

  fillArray();

  return;
}

// ************************************************
//DB Function which does the core reweighting. This assumes that oscillation weights have already been calculated and stored in samplePDFFDBase[iSample].osc_w[iEvent]
//This function takes advantage of most of the things called in setupSKMC to reduce reweighting time
//It also follows the ND code reweighting pretty closely
//This function fills the samplePDFFD_array array which is binned to match the sample binning, such that bin[1][1] is the equivalent of _hPDF2D->GetBinContent(2,2) {Noticing the offset}
void samplePDFFDBase::fillArray() {
// ************************************************
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
      
      //DB Set oscillation weights for NC events to 1.0
      //DB Another speedup - Why bother storing NC signal events and calculating the oscillation weights when we just throw them out anyway? Therefore they are skipped in setupMC
	  //
	  //LW Checking if NC event is signal (oscillated or not), if yes: osc_w = 0 || if no: osc_w = 1.0
      if (MCSamples[iSample].isNC[iEvent] && MCSamples[iSample].signal) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
	  	MCSamples[iSample].osc_w[iEvent] = 0.0;
	  	continue;
      }

	  //ETA - I need to check that this doesn't cause a problem for atmospherics and or tau samples
      if (MCSamples[iSample].isNC[iEvent] && !MCSamples[iSample].signal) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
	  	MCSamples[iSample].osc_w[iEvent] = 1.0;
      }
	  //DB Set oscillation weights for NC events to 1.0
	  //DB Another speedup - Why bother storing NC signal events and calculating the oscillation weights when we just throw them out anyway? Therefore they are skipped in setupSKMC
	  if (MCSamples[iSample].isNC[iEvent]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
		MCSamples[iSample].osc_w[iEvent] = 1.0;	
	  }
      
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
		//std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		samplePDFFD_array[YBinToFill][XBinToFill] += totalweight;
		samplePDFFD_array_w2[YBinToFill][XBinToFill] += totalweight*totalweight;
      }
    }
  }

#endif // end the else in openMP
  return;
}

#ifdef MULTITHREAD
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
		//std::cout << "norm weight is " << normweight << std::endl;
		if (normweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		funcweight = CalcXsecWeightFunc(iSample,iEvent);
		//DB Catch negative func weights and skip any event with a negative event. Previously we would set weight to zero and continue but that is inefficient
		//std::cout << "Func weight is " << funcweight << std::endl;
		if (funcweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight*funcweight;

		//DB Set oscillation weights for NC events to 1.0
		//DB Another speedup - Why bother storing NC signal events and calculating the oscillation weights when we just throw them out anyway? Therefore they are skipped in setupSKMC
		//LW Checking if NC event is signal (oscillated or not), if yes: osc_w = 0 || if no: osc_w = 1.0
		if (MCSamples[iSample].isNC[iEvent] && MCSamples[iSample].signal) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
		  MCSamples[iSample].osc_w[iEvent] = 0.0;
		  continue;
		}
		if (MCSamples[iSample].isNC[iEvent] && !MCSamples[iSample].signal) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
		  MCSamples[iSample].osc_w[iEvent] = 1.0;
		}

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

		//std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		//std::cout << "XVar is " << XVar << " and rw_upper_bin_edge is " << MCSamples[iSample].rw_upper_xbinedge[iEvent] << " and rw_lower_xbinedge is " << MCSamples[iSample].rw_lower_xbinedge[iEvent] << std::endl;
		//DB Check to see if momentum shift has moved bins
		//DB - First, check to see if the event is still in the nominal bin	
		if (XVar < MCSamples[iSample].rw_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_xbinedge[iEvent]) {
		  XBinToFill = MCSamples[iSample].NomXBin[iEvent];
		  //std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
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
		else{
		  //std::cout << "Not filled samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
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

//DB Adding in Oscillator class support for smeared oscillation probabilities
void samplePDFFDBase::SetOscillator(Oscillator* Osc_) {
#if defined (USE_PROB3)
  MACH3LOG_ERROR("Atmospheric Oscillator only defined using CUDAProb3 - USE_PROB3 is defined and indicates that Prob3++/probGPU is being used");	
  throw MaCh3Exception(__FILE__ , __LINE__ );
#endif

  Osc = Osc_;
  MACH3LOG_INFO("Set Oscillator");
	
  FindEventOscBin();
}

void samplePDFFDBase::FindEventOscBin() {
  for(int i = 0; i < getNMCSamples(); i++) {
    for (int j = 0;j < getNEventsInSample(i); j++) {
      MCSamples[i].osc_w_pointer[j] = Osc->retPointer(MCSamples[i].nutype,MCSamples[i].oscnutype,*(MCSamples[i].rw_etru[j]),MCSamples[i].rw_truecz[j]);
    }
  }
  MACH3LOG_INFO("Set all oscillation pointers to Oscillator");
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

  MACH3LOG_INFO("FD Object has {} events", fdobj->nEvents);

  for(int iEvent=0; iEvent < fdobj->nEvents; ++iEvent){
    std::list< int > XsecBins = {};
	if (XsecCov) {
	  //std::cout << "Xsec norms is of size " << xsec_norms.size() << std::endl;
	  for (std::vector<XsecNorms4>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it) {
		//std::cout << "Looping over systematic " << (*it).name << std::endl;
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
			  //std::cout << "DID MATCH " << fdobj->nupdgUnosc << " with " << (*it).preoscpdgs.at(iPDG) << std::endl;
			}
			//else{std::cout << "Didn't match " << fdobj->nupdgUnosc << " with " << (*it).preoscpdgs.at(iPDG) << std::endl;}
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
			  //if((*it).name.find("b_") != std::string::npos){
			  //  std::cout << "Failed because " << ReturnKinematicParameter((*it).KinematicVarStr[iKinematicParameter], iSample, iEvent) << " is less than " << (*it).Selection[iKinematicParameter][0] << std::endl; 
			  //}
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

//LW 
//Setup chosen oscillation calculator for each subsample
//Default Baseline Implementation
//Add your own implementation in experiment specific SamplePDF Class if necessary!!
// ETA - pass the yaml config used in the executable. This will include the cov osc
// and other information. Need to double check that this is sensible
void samplePDFFDBase::SetupOscCalc(double PathLength, double Density)
{

  for (int iSample=0; iSample < (int)MCSamples.size(); iSample++) {

#if defined (USE_PROB3) && defined (CPU_ONLY)
// if we're using Prob3++ CPU then initialise BargerPropagator object
// if we're using Prob3++ GPU then we don't need to do this since event information gets passed straight to ProbGpu.cu in CalcOscWeights
    MCSamples[iSample].Oscillator = new BargerPropagator();
    MCSamples[iSample].Oscillator->UseMassEigenstates(false);
    MCSamples[iSample].Oscillator->SetOneMassScaleMode(false);
    MCSamples[iSample].Oscillator->SetWarningSuppression(true);
#endif

#if not defined (USE_PROB3)
//if we're using CUDAProb3 then make vector of energies and convert to CUDAProb3 structs
    std::vector<double> etruVector(*(MCSamples[iSample].rw_etru), *(MCSamples[iSample].rw_etru) + MCSamples[iSample].nEvents);
    MCSamples[iSample].ProbType = SwitchToCUDAProbType(GetCUDAProbFlavour(MCSamples[iSample].nutype, MCSamples[iSample].oscnutype));
	// CUDAProb3 takes probType and antineutrino/neutrino separately
    if (MCSamples[iSample].nutype < 0) {MCSamples[iSample].NeutrinoType = cudaprob3::NeutrinoType::Antineutrino;}
    else {MCSamples[iSample].NeutrinoType = cudaprob3::NeutrinoType::Neutrino;}
#if defined (CPU_ONLY) || defined (USE_FPGA)
//if we just want to use CUDAProb3 CPU then setup BeamCpuPropagator object
#if defined (MULTITHREAD)
//if we want to multithread then get number of threads from OMP_NUM_THREADS env variable
    MCSamples[iSample].Oscillator = new cudaprob3::BeamCpuPropagator<double>(MCSamples[iSample].nEvents, omp_get_max_threads());
  MCSamples[iSample].Oscillator->setPathLength(PathLength);
  MCSamples[iSample].Oscillator->setDensity(Density);
#else
//if we're not mulithreading then just set it to 1
    MCSamples[iSample].Oscillator = new cudaprob3::BeamCpuPropagator<double>(MCSamples[iSample].nEvents, 1);
  MCSamples[iSample].Oscillator->setPathLength(PathLength);
  MCSamples[iSample].Oscillator->setDensity(Density);
#endif //MULTITHREAD
#else
//if we want to use CUDAProb3 GPU then setup BeamCudaPropagator object
    MCSamples[iSample].Oscillator = new cudaprob3::BeamCudaPropagatorSingle(0, MCSamples[iSample].nEvents);
    MCSamples[iSample].Oscillator->setPathLength(PathLength);
    MCSamples[iSample].Oscillator->setDensity(Density);
#endif // CPU_ONLY
    MCSamples[iSample].Oscillator->setEnergyList(etruVector);
#endif // USE_PROB3
  }
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
      MCSamples[mc_i].y_var[event_i] = &(MCSamples[mc_i].dummy_value);//-9999.;// = dummy;//ReturnKinematicParameterByReference(XVarStr, mc_i, event_i);

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

//ETA - this can be changed quite easily to check the number of XBins and YBins.
//We can slowly but surely remove any trace of BinningOpt
/*int samplePDFFDBase::GetNDim() {
  switch(BinningOpt) {
  case 0: 
  case 1:
    return 1;
  case 2:
  case 3:
  case 4: 
    return 2;
  default:
	std::cerr << "Error, unrecognsied BinningOpt!!" << std::endl;
	throw;
    return 0;
  }  
}
*/

void samplePDFFDBase::addData(std::vector<double> &data) {
  dataSample = new std::vector<double>(data);
  dataSample2D = NULL;
  dathist2d = NULL;
  dathist->Reset(); 

  if (GetNDim()!=1) {std::cerr << "Trying to set a 1D 'data' histogram in a 2D sample - Quitting" << std::endl; throw;}

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
	  MACH3LOG_ERROR("Trying to set a 2D 'data' histogram in a 1D sample - Quitting");
	  throw MaCh3Exception(__FILE__ , __LINE__ );}

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

double samplePDFFDBase::GetEventWeight(int iSample, int iEntry) {
  //HW : DON'T EDIT THIS!!!! (Pls make a weights pointer instead ^_^)
  double totalweight = 1.0;
  for (int iParam=0;iParam<MCSamples[iSample].ntotal_weight_pointers[iEntry];iParam++) {
	//std::cout << "Weight " << iParam << " is " <<  *(MCSamples[iSample].total_weight_pointers[iEntry][iParam]) << std::endl;
    totalweight *= *(MCSamples[iSample].total_weight_pointers[iEntry][iParam]);
  }
  return totalweight;
}

/// @func fillSplineBins()
/// @brief Finds the binned spline that an event should apply to and stored them in a
/// a vector for easy evaluation in the fillArray() function.
void samplePDFFDBase::fillSplineBins() {

  std::cout << "Now in fillSplineBins" << std::endl;
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
      MCSamples[i].xsec_spline_pointers[j] = new const double*[MCSamples[i].nxsec_spline_pointers[j]];
 
	  for(int spline=0; spline<MCSamples[i].nxsec_spline_pointers[j]; spline++){          
		//Event Splines indexed as: sample name, oscillation channel, syst, mode, etrue, var1, var2 (var2 is a dummy 0 for 1D splines)
		MCSamples[i].xsec_spline_pointers[j][spline] = splineFile->retPointer(EventSplines[spline][0], EventSplines[spline][1], EventSplines[spline][2], 
			EventSplines[spline][3], EventSplines[spline][4], EventSplines[spline][5], EventSplines[spline][6]);
	  }

	}
  }
  std::cout << "Filled spline bins" << std::endl;
  return;
}

double samplePDFFDBase::GetLikelihood()
{
  if (samplePDFFD_data == NULL) {
      MACH3LOG_ERROR("data sample is empty!");
      return -1;
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
 
#ifndef USE_PROB3 
// ************************************************
// Switch from MaCh3 CUDAProb flavour to CUDAProb Probtype
inline cudaprob3::ProbType samplePDFFDBase::SwitchToCUDAProbType(CUDAProb_nu CUDAProb_nu) {
//*************************************************  
  switch(CUDAProb_nu)
  {
    case CUDAProb_nu::e_e : return cudaprob3::ProbType::e_e;
    case CUDAProb_nu::e_m : return cudaprob3::ProbType::e_m;
    case CUDAProb_nu::e_t : return cudaprob3::ProbType::e_m;
    case CUDAProb_nu::m_e : return cudaprob3::ProbType::m_e;
    case CUDAProb_nu::m_m : return cudaprob3::ProbType::m_m;
    case CUDAProb_nu::m_t : return cudaprob3::ProbType::m_t;
    case CUDAProb_nu::t_e : return cudaprob3::ProbType::t_e;
    case CUDAProb_nu::t_m : return cudaprob3::ProbType::t_m;
    case CUDAProb_nu::t_t : return cudaprob3::ProbType::t_t;
    default:
      MACH3LOG_ERROR("Unknown CUDAProbType!");
      throw;
  }
}
#endif

