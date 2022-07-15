#include <TROOT.h>
#include "manager/manager.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "samplePDFFDBase.h"

// Constructors for erec-binned errors

samplePDFFDBase::samplePDFFDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov)
  : samplePDFBase(pot)
//DB Throughout constructor and init, pot is livetime for atmospheric samples
{
  std::cout << "-------------------------------------------------------------------" <<std::endl;
  std::cout << "Creating samplePDFFDBase object.." << "\n" << std::endl;
  std::cout << "- Using SK sample config in this file " << mc_version << std::endl;

  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!" << std::endl; throw;}

  init(pot, mc_version, xsec_cov);          

  bNu = new BargerPropagator();
  bNu->UseMassEigenstates(false);
  bNu->SetOneMassScaleMode(false);
  bNu->SetWarningSuppression(true);
  
  Beta=1;
  useBeta=false;
  applyBetaNue=false;
  applyBetaDiag=false;

  //ETA - leave this out for now, need to fix and make things nice and configurable
  //EnergyScale *energy_first = new EnergyScale();
  //energy_first->SetUncertainty(1.2);
  //ShiftFunctors.push_back(energy_first);
}


samplePDFFDBase::~samplePDFFDBase()
{

}

void samplePDFFDBase::fill1DHist()
{
  // DB Commented out by default - Code heading towards getLikelihood using arrays instead of root objects 
  // Wouldn't actually need this for getLikelihood as TH objects wouldn't be filled   
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
  // DB Commented out by default - Code heading towards getLikelihood using arrays instead of root objects 
  // Wouldn't actually need this for getLikelihood as TH objects wouldn't be filled   
  _hPDF2D->Reset();
  for (unsigned int yBin=0;yBin<(YBinEdges.size()-1);yBin++) {
    for (unsigned int xBin=0;xBin<(XBinEdges.size()-1);xBin++) {
	  //std::cout << "Filling _hPDF2D with " << samplePDFFD_array[yBin][xBin] << std::endl;
      _hPDF2D->SetBinContent(xBin+1,yBin+1,samplePDFFD_array[yBin][xBin]);
    }
  }

  return;
}


void samplePDFFDBase::useBinnedOscReweighting(bool ans) 
{
  osc_binned = ans;
  if(ans == true) {
    std::cout << "WARNING: you are using binned oscillation weight without specifying the binning. It will use E_reco binning to set the energy (recommended to use set1Dbinning first !). If you want another binning use : \n useBinnedOscReweighting(true, int, double*) \nwhere int is the number of bins and double* is an array with the bins boundaries." << std::endl ;

    /// Get the binning from the MC histogram

    const int nb_bins = _hPDF1D -> GetXaxis() -> GetNbins() ;      
    const double* osc_bins = _hPDF1D -> GetXaxis() -> GetXbins() -> GetArray();
    osc_binned_axis = new TAxis(nb_bins, osc_bins) ;
  }
}

void samplePDFFDBase::useBinnedOscReweighting(bool ans, int nbins, double *osc_bins) 
{
  osc_binned = ans;
  if(ans == true) {
    osc_binned_axis = new TAxis(nbins, osc_bins) ;
  }
}

bool samplePDFFDBase::IsEventSelected(int iSample, int iEvent) {
  
  double Val;

  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
 
    Val = ReturnKinematicParameter(static_cast<KinematicTypes>(Selection[iSelection][0]),iSample,iEvent);
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

//Same as the function above but just acts on the vector and the event
bool samplePDFFDBase::IsEventSelected(std::vector< std::vector<double> > &Selection, int iSample, int iEvent) {
  
  double Val;

  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {

    
    Val = ReturnKinematicParameter(static_cast<KinematicTypes>(Selection[iSelection][0]),iSample,iEvent);
    //DB If multiple return values, it will consider each value seperately
    //DB Already checked that Selection vector is correctly sized
    
    //DB In the case where Selection[0].size()==3, Only Events with Val >= Selection[iSelection][1] and Val < Selection[iSelection][2] are considered Passed
	//ETA - also check whether we're actually applying a lower or upper cut by checking they aren't -999
	if(Val >= Selection[iSelection][2] && Selection[iSelection][2] != -999){
	  //std::cout << "Cutting event as " << Val << " is greater than " << Selection[iSelection][1]
	  return false;
	}
	else if(Val < Selection[iSelection][1] && Selection[iSelection][1] != -999){
	  return false;
	}
  }

  //DB To avoid unneccessary checks, now return false rather than setting bool to true and continuing to check
  return true;
}

void samplePDFFDBase::reweight(double *oscpar) // Reweight function (this should be different depending on whether you use one or 2 sets of oscpars)
{

  if (MCSamples.size()==6) {
	reweight(oscpar, oscpar);
	return;
  }
  else {

	for (int i=0; i< (int)MCSamples.size(); ++i) {
	  for(int j = 0; j < MCSamples[i].nEvents; ++j) {
		MCSamples[i].osc_w[j] = calcOscWeights(MCSamples[i].nutype, MCSamples[i].oscnutype, *(MCSamples[i].rw_etru[j]), oscpar);
	  }
	}

  }
    
  fillArray();

  return;
}

void samplePDFFDBase::reweight(double *oscpar_nub, double *oscpar_nu) // Reweight function (this should be different for one vs 2 sets of oscpars)
{
  for(int i = 0; i < (int)MCSamples.size(); ++i) {
    for(int j = 0; j < MCSamples[i].nEvents; ++j) {
      MCSamples[i].osc_w[j] = calcOscWeights(MCSamples[i].nutype, MCSamples[i].oscnutype, *(MCSamples[i].rw_etru[j]), oscpar_nub, oscpar_nu);
    }
  }

  fillArray();
}

double samplePDFFDBase::calcOscWeights(int nutype, int oscnutype, double en, double *oscpar)
{
  bNu->SetMNS(oscpar[0], oscpar[2], oscpar[1], oscpar[3], oscpar[4], oscpar[5], en, doubled_angle);
  bNu->propagateLinear(nutype , oscpar[7], oscpar[8]); 

  return bNu->GetProb(nutype, oscnutype);
}

double samplePDFFDBase::calcOscWeights(int nutype, int oscnutype, double en, double *oscpar_nub, double *oscpar_nu)
{
  if (nutype < 0) // is antinu
    {
      bNu->SetMNS(oscpar_nub[0], oscpar_nub[2], oscpar_nub[1], oscpar_nub[3], oscpar_nub[4], oscpar_nub[5], en, doubled_angle, nutype);
      bNu->propagateLinear(nutype , oscpar_nub[7], oscpar_nub[8]);
      return bNu->GetProb(nutype, oscnutype);
    }
  else // is nu
    {
      bNu->SetMNS(oscpar_nu[0], oscpar_nu[2], oscpar_nu[1], oscpar_nu[3], oscpar_nu[4], oscpar_nu[5], en, doubled_angle); 
      bNu->propagateLinear(nutype , oscpar_nu[7], oscpar_nu[8]);
      return bNu->GetProb(nutype, oscnutype);
    }
}

//DB Function which does the core reweighting. This assumes that oscillation weights have already been calculated and stored in samplePDFFDBase[iSample].osc_w[iEvent]
//This function takes advantage of most of the things called in setupSKMC to reduce reweighting time
//It also follows the ND code reweighting pretty closely
//This function fills the samplePDFFD_array array which is binned to match the sample binning, such that bin[1][1] is the equivalent of _hPDF2D->GetBinContent(2,2) {Noticing the offset}
void samplePDFFDBase::fillArray() {

  //DB Reset which cuts to apply
  Selection = StoredSelection;

  // Call entirely different routine if we're running with openMP
#ifdef MULTITHREAD
  fillArray_MP();
#else

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  //DB Reset values stored in PDF array to 0.
  for (int yBin=0;yBin<nYBins;yBin++) {
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
    }
  }

  reconfigureFuncPars();

  for (int iSample=0;iSample<(int)MCSamples.size();iSample++) {
    MCSamples[iSample].splineFile->FindSplineSegment();
    MCSamples[iSample].splineFile->calcWeights();
  }

  for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
    for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {
      
     applyShifts(iSample, iEvent);
     if(!IsEventSelected(iSample, iEvent)){
	   continue;
	 }

#if USEBETA == 1
      MCSamples[iSample].osc_w[iEvent] = ApplyBetaWeights(MCSamples[iSample].osc_w[iEvent],iSample);
#endif

      double splineweight = 1.0;
      double normweight = 1.0;
      double funcweight = 1.0;
      double totalweight = 1.0;
      
      //DB SKDet Syst
      //As weights were skdet::fParProp, and we use the non-shifted erec, we might as well cache the corresponding fParProp index for each event and the pointer to it
      //MCSamples[iSample].skdet_w[iEvent] = *(MCSamples[iSample].skdet_pointer[iEvent]);

      //DB Xsec syst
      //Loop over stored spline pointers
      for (int iSpline=0;iSpline<MCSamples[iSample].nxsec_spline_pointers[iEvent];iSpline++) {
		splineweight *= *(MCSamples[iSample].xsec_spline_pointers[iEvent][iSpline]);
      }
      //DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient. Do this on a spline-by-spline basis
      if (splineweight <= 0.){
		MCSamples[iSample].xsec_w[iEvent] = 0.;
	   	continue;
	  }

      //Loop over stored normalisation and function pointers 
      for (int iParam=0;iParam<MCSamples[iSample].nxsec_norm_pointers[iEvent];iParam++) {
	normweight *= *(MCSamples[iSample].xsec_norm_pointers[iEvent][iParam]);
      }
      //DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
      if (normweight <= 0.){
		MCSamples[iSample].xsec_w[iEvent] = 0.;
	   	continue;
	  }

      funcweight = calcFuncSystWeight(iSample,iEvent);
      //DB Catch negative func weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
      if (funcweight <= 0.){
		MCSamples[iSample].xsec_w[iEvent] = 0.;
	   	continue;
	  }

      MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight*funcweight;
      
      //DB Set oscillation weights for NC events to 1.0
      //DB Another speedup - Why bother storing NC signal events and calculating the oscillation weights when we just throw them out anyway? Therefore they are skipped in setupMC
      if (MCSamples[iSample].isNC[iEvent]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
		MCSamples[iSample].osc_w[iEvent] = 1.0;
		continue;
      }

      //DB Total weight
      for (int iParam=0;iParam<MCSamples[iSample].ntotal_weight_pointers[iEvent];iParam++) {
	totalweight *= *(MCSamples[iSample].total_weight_pointers[iEvent][iParam]);
      }
      //DB Catch negative weights and skip any event with a negative event
      if (totalweight <= 0.){
		MCSamples[iSample].xsec_w[iEvent] = 0.;
	   	continue;
	  }
      
      //DB Switch on BinningOpt to allow different binning options to be implemented
      //The alternative would be to have inheritance based on BinningOpt
      double XVar = *(MCSamples[iSample].x_var[iEvent]);

      //DB Commented out by default but if we ever want to consider shifts in theta this will be needed
      //double YVar = MCSamples[iSample].rw_theta[iEvent];

      //DB Find the relevant bin in the PDF for each event
      int XBinToFill = -1;
      int YBinToFill = MCSamples[iSample].NomYBin[iEvent];

      //DB Check to see if momentum shift has moved bins
      //
      //DB - First, check to see if the event is still in the nominal bin
      if (XVar < MCSamples[iSample].rw_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_xbinedge[iEvent]) {
	XBinToFill = MCSamples[iSample].NomXBin[iEvent];
      }
      //DB - Second, check to see if the event is outside of the binning range and skip event if it is
      else if (XVar < XBinEdges[0] || XVar >= up_bnd) {
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

      //DB Fill relevant part of thread array
      if (XBinToFill != -1 && YBinToFill != -1) {
		std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
	samplePDFFD_array[YBinToFill][XBinToFill] += totalweight;
      }
	  else{
		std::cout << "Not filled samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
	  }

    }
  }

#endif // end the else in openMP
  return;
}

#ifdef MULTITHREAD
void samplePDFFDBase::fillArray_MP() {

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  //DB Reset values stored in PDF array to 0.
  for (int yBin=0;yBin<nYBins;yBin++) {
	for (int xBin=0;xBin<nXBins;xBin++) {
	  samplePDFFD_array[yBin][xBin] = 0.;
	}
  }

  reconfigureFuncPars();

  //DB Unfortunately need to have this because OMP doesn't allow member variables to be used in the shared environment
  //At end of function, copy contents of 'samplePDFFD_array_class' into 'samplePDFFD_array' which can then be used in getLikelihood
  //
  //This is stored as [y][x] due to shifts only occuring in the x variable (Erec/Lep mom) - I believe this will help reduce cache misses 
  double** samplePDFFD_array_class = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
	samplePDFFD_array_class[yBin] = new double[nXBins];
	for (int xBin=0;xBin<nXBins;xBin++) {
	  samplePDFFD_array_class[yBin][xBin] = 0.;
	}
  }

  double** samplePDFFD_array_private = NULL;

  // Declare the omp parallel region
  // The parallel region needs to stretch beyond the for loop!
#pragma omp parallel private(samplePDFFD_array_private)
  {
	// private to each thread
	samplePDFFD_array_private = new double*[nYBins];
	for (int yBin=0;yBin<nYBins;yBin++) {
	  samplePDFFD_array_private[yBin] = new double[nXBins];
	  for (int xBin=0;xBin<nXBins;xBin++) {
		samplePDFFD_array_private[yBin][xBin] = 0.;
	  }
	}

	//DB From Clarence's suggestion, moved spline weight calculation into this OMP parallel region but this did not reduce s/step
	//Maybe more efficient to OMP inside splineFile->FindSplineSegment and splineFile->calcWeights
#pragma omp for
	for (int iSample=0;iSample<(int)MCSamples.size();iSample++) {
	  MCSamples[iSample].splineFile->FindSplineSegment();
	  MCSamples[iSample].splineFile->calcWeights();
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
	// 1. Multi-thread getLikelihood() - Implemented and doesn't significantly affect s/step
	// 2. Don't use TH objects in getLikelihood() and instead just use arrays - Implemented and doesn't significantly affect s/step
	// 3. Order minituples in Y-axis variable as this will *hopefully* reduce cache misses inside samplePDFFD_array_class[yBin][xBin]
	//
	// We will hit <0.1 s/step eventually! :D

	//==================================================
	//Calc Weights and fill Array

	for (unsigned int iSample=0;iSample<MCSamples.size();iSample++) {
#pragma omp for
	  for (int iEvent=0;iEvent<MCSamples[iSample].nEvents;iEvent++) {

        //ETA - generic functions to apply shifts to kinematic variables
		applyShifts(iSample, iEvent);

        //ETA - generic functions to apply shifts to kinematic variable
		if(!IsEventSelected(iSample, iEvent)){
		  continue;
		}

#if USEBETA == 1
		MCSamples[iSample].osc_w[iEvent] = ApplyBetaWeights(MCSamples[iSample].osc_w[iEvent],iSample);
#endif

		double splineweight = 1.0;
		double normweight = 1.0;
		double funcweight = 1.0;
		double totalweight = 1.0;

		//DB SKDet Syst
		//As weights were skdet::fParProp, and we use the non-shifted erec, we might as well cache the corresponding fParProp index for each event and the pointer to it
		//ETA - removing this,  
		//MCSamples[iSample].skdet_w[iEvent] = *(MCSamples[iSample].skdet_pointer[iEvent]);

		//DB Xsec syst
		//Loop over stored spline pointers
		for (int iSpline=0;iSpline<MCSamples[iSample].nxsec_spline_pointers[iEvent];iSpline++) {
		  splineweight *= *(MCSamples[iSample].xsec_spline_pointers[iEvent][iSpline]);
		}
		//std::cout << "Spline weight is " << splineweight << std::endl;
		//DB Catch negative spline weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
		if (splineweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		//Loop over stored normalisation and function pointers
		for (int iParam=0;iParam<MCSamples[iSample].nxsec_norm_pointers[iEvent];iParam++) {
		  normweight *= *(MCSamples[iSample].xsec_norm_pointers[iEvent][iParam]);
		}
		//DB Catch negative norm weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
		//std::cout << "norm weight is " << normweight << std::endl;
		if (normweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		funcweight = calcFuncSystWeight(iSample,iEvent);
		//DB Catch negative func weights and skip any event with a negative event. Previously we would set weight to zere and continue but that is inefficient
		//std::cout << "Func weight is " << funcweight << std::endl;
		if (funcweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		MCSamples[iSample].xsec_w[iEvent] = splineweight*normweight*funcweight;

		//DB Set oscillation weights for NC events to 1.0
		//DB Another speedup - Why bother storing NC signal events and calculating the oscillation weights when we just throw them out anyway? Therefore they are skipped in setupSKMC
		if (MCSamples[iSample].isNC[iEvent]) { //DB Abstract check on MaCh3Modes to determine which apply to neutral current
		  MCSamples[iSample].osc_w[iEvent] = 1.0;	
		}

		//DB Total weight
		for (int iParam=0;iParam<MCSamples[iSample].ntotal_weight_pointers[iEvent];iParam++) {
		  totalweight *= *(MCSamples[iSample].total_weight_pointers[iEvent][iParam]);
		}

		//std::cout << "Oscillation weight is " << MCSamples[iSample].osc_w[iEvent] << std::endl;

		//std::cout << "Total weight is " << totalweight << std::endl;
		//DB Catch negative weights and skip any event with a negative event
		if (totalweight <= 0.){
		  MCSamples[iSample].xsec_w[iEvent] = 0.;
		  continue;
		}

		//DB Switch on BinningOpt to allow different binning options to be implemented
		//The alternative would be to have inheritance based on BinningOpt
		double XVar = (*(MCSamples[iSample].x_var[iEvent]));

		//DB Commented out by default but if we ever want to consider shifts in theta this will be needed
		//double YVar = skmcSamples[iSample].rw_theta[iEvent];

		//DB Find the relevant bin in the PDF for each event
		int XBinToFill = -1;
		int YBinToFill = MCSamples[iSample].NomYBin[iEvent];

		//std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		//std::cout << "XVar is " << XVar << " and rw_upper_bin_edge is " << MCSamples[iSample].rw_upper_xbinedge[iEvent] << " and rw_lower_xbinedge is " << MCSamples[iSample].rw_lower_xbinedge[iEvent] << std::endl;
		//DB Check to see if momentum shift has moved bins
		//
		//DB - First, check to see if the event is still in the nominal bin	
		if (XVar < MCSamples[iSample].rw_upper_xbinedge[iEvent] && XVar >= MCSamples[iSample].rw_lower_xbinedge[iEvent]) {
		  XBinToFill = MCSamples[iSample].NomXBin[iEvent];
		  //std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		}
		//DB - Second, check to see if the event is outside of the binning range and skip event if it is
		else if (XVar < XBinEdges[0] || XVar >= up_bnd) {
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

		//DB Fill relevant part of thread array
		if (XBinToFill != -1 && YBinToFill != -1) {
		  //std::cout << "Filling samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		  samplePDFFD_array_private[YBinToFill][XBinToFill] += totalweight;
		}
		else{
		  //std::cout << "Not filled samplePDFFD_array at YBin: " << YBinToFill << " and XBin: " << XBinToFill << std::endl;
		}

	  }
	}    

	//End of Calc Weights and fill Array
	//==================================================

#pragma omp parallel shared(samplePDFFD_array_class)
	{
	  for (int yBin=0;yBin<nYBins;yBin++) {
		for (int xBin=0;xBin<nXBins;xBin++) {
#pragma omp atomic
		  samplePDFFD_array_class[yBin][xBin] += samplePDFFD_array_private[yBin][xBin];
		}
	  }
	}

	for (int yBin=0;yBin<nYBins;yBin++) {
	  delete[] samplePDFFD_array_private[yBin];
	}
	delete[] samplePDFFD_array_private;
  }

  // DB Copy contents of 'samplePDFFD_array_class' into 'samplePDFFD_array' which can then be used in getLikelihood
  for (int yBin=0;yBin<nYBins;yBin++) {
	for (int xBin=0;xBin<nXBins;xBin++) {
	  samplePDFFD_array[yBin][xBin] = samplePDFFD_array_class[yBin][xBin];
	}
  }

  for (int yBin=0;yBin<nYBins;yBin++) {
	delete[] samplePDFFD_array_class[yBin];
  }
  delete[] samplePDFFD_array_class;
}
#endif

//ETA
void samplePDFFDBase::setXsecCov(covarianceXsec *xsec){

  xsecCov = xsec;

  // Get the map between the normalisation parameters index, their name, what mode they should apply to, and what target
  //This information is used later in calcXsecNormsBins to decide if a parameter applies to an event

  //DB Now get this information using the DetID from the config
  xsec_norms = xsecCov->GetNormParsFromDetID(SampleDetID);
  nFuncParams = xsecCov->GetNumFuncParamsFromDetID(SampleDetID);
  funcParsNames = xsecCov->GetFuncParsNamesFromDetID(SampleDetID);
  funcParsIndex = xsecCov->GetFuncParsIndexFromDetID(SampleDetID);

  return;
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
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
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
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges1D();
}

void samplePDFFDBase::FindNominalBinAndEdges1D() {

  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(int mc_i = 0 ; mc_i < (int)MCSamples.size() ; mc_i++){
    for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){
      int bin = _hPDF1D->FindBin(*(MCSamples[mc_i].x_var[event_i]));//ETA - TODO get rid of 0.001... this is only for SK
	  //std::cout << "FOUND XNOMBIN AT " << bin << " from " << *(MCSamples[mc_i].x_var[event_i]) << std::endl;
      
      double low_lower_edge = __DEFAULT_RETURN_VAL__;
      if (bin==0) {
	low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
      } else {
	low_lower_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin-1);
      }
      
      double low_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin);
      double upper_edge = _hPDF1D->GetXaxis()->GetBinUpEdge(bin);
	  //std::cout << "FINDING EDGES" << std::endl; 
	  //std::cout << "Low edge is " << low_edge << std::endl;
	  //std::cout << "Upper edge is " << upper_edge << std::endl;
      
      double upper_upper_edge = __DEFAULT_RETURN_VAL__;
      if (bin<(_hPDF1D->GetNbinsX()-2)) {
	upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+2);
      } else {
	upper_upper_edge = _hPDF1D->GetXaxis()->GetBinLowEdge(bin+1);
      }

      if ((bin-1) > 0 && (bin-1) < int(XBinEdges.size()-1)) {

	MCSamples[mc_i].NomXBin[event_i] = bin-1;
      } else {
	MCSamples[mc_i].NomXBin[event_i] = -1;
	low_edge = __DEFAULT_RETURN_VAL__;
	upper_edge = __DEFAULT_RETURN_VAL__;
	low_lower_edge = __DEFAULT_RETURN_VAL__;
	upper_upper_edge = __DEFAULT_RETURN_VAL__;
      }
      MCSamples[mc_i].NomYBin[event_i] = 0;
      
	  /*std::cout << "FOUND EDGES TO BE " << std::endl;
	  std::cout << "Low edges " << low_edge << std::endl;
	  std::cout << "Upper edges " << upper_edge << std::endl;
	  std::cout << "Lower lower edge " << low_lower_edge << std::endl;
	  std::cout << "Upper uper edge " << upper_upper_edge << std::endl;
	  */
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
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
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
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_array[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_array[yBin][xBin] = 0.;
    }
  }

  FindNominalBinAndEdges2D();
}


void samplePDFFDBase::FindNominalBinAndEdges2D() {

  //Set rw_pdf_bin and rw_upper_xbinedge and rw_lower_xbinedge for each skmc_base
  for(int mc_i = 0 ; mc_i < (int)MCSamples.size() ; mc_i++){
    for(int event_i = 0 ; event_i < MCSamples[mc_i].nEvents ; event_i++){
      //Global bin number

      int bin = _hPDF2D->FindBin(*(MCSamples[mc_i].x_var[event_i]), *(MCSamples[mc_i].y_var[event_i]));//ETA - TODO get rid of 0.001... this is only for SK
      
      int bin_x = -999;
      int bin_y = -999;
      int bin_z = -999;
      _hPDF2D->GetBinXYZ(bin, bin_x, bin_y, bin_z);
      //erec is the x-axis so get GetXaxis then find the bin edges using the x bin number
      
      double low_lower_edge = __DEFAULT_RETURN_VAL__;
      if (bin==0) {
	low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
      } else {
	low_lower_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x-1);
      }
      
      double low_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x);
      double upper_edge = _hPDF2D->GetXaxis()->GetBinUpEdge(bin_x);
      
      double upper_upper_edge = __DEFAULT_RETURN_VAL__;
      if (bin<(_hPDF2D->GetNbinsX()-2)) {
	upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+2);
      } else {
	upper_upper_edge = _hPDF2D->GetXaxis()->GetBinLowEdge(bin_x+1);
      }
      
      if ((bin_x-1) > 0 && (bin_x-1) < int(XBinEdges.size()-1)) {
	    MCSamples[mc_i].NomXBin[event_i] = bin_x-1;
      } else {
	MCSamples[mc_i].NomXBin[event_i] = -1;
	low_edge = __DEFAULT_RETURN_VAL__;
	upper_edge = __DEFAULT_RETURN_VAL__;
	low_lower_edge = __DEFAULT_RETURN_VAL__;
	upper_upper_edge = __DEFAULT_RETURN_VAL__;
      }
      MCSamples[mc_i].NomYBin[event_i] = bin_y-1; 
      MCSamples[mc_i].rw_lower_xbinedge[event_i] = low_edge;
      MCSamples[mc_i].rw_upper_xbinedge[event_i] = upper_edge;
      MCSamples[mc_i].rw_lower_lower_xbinedge[event_i] = low_lower_edge;
      MCSamples[mc_i].rw_upper_upper_xbinedge[event_i] = upper_upper_edge;
    }
  }

}

int samplePDFFDBase::getNDim() {
  switch(BinningOpt) {
  case 0: 
  case 1:
    return 1;
  case 2:
  case 3:
  case 4: 
    return 2;
  default:
    return 0;
  }
  
}

void samplePDFFDBase::addData(std::vector<double> &data) {
  dataSample = new std::vector<double>(data);
  dataSample2D = NULL;
  dathist2d = NULL;
  dathist->Reset(); 

  if (getNDim()!=1) {std::cerr << "Trying to set a 1D 'data' histogram in a 2D sample - Quitting" << std::endl; throw;}

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
}

void samplePDFFDBase::addData(std::vector< vector <double> > &data) {
  dataSample2D = new std::vector< vector <double> >(data);
  dataSample = NULL;
  dathist = NULL;
  dathist2d->Reset();                                                       

  if (getNDim()!=2) {std::cerr << "Trying to set a 2D 'data' histogram in a 1D sample - Quitting" << std::endl; throw;}

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
}

void samplePDFFDBase::addData(TH1D* Data) {
  std::cout << "adding 1D data histogram : " << Data->GetName() << " with " << Data->Integral() << " events" << std::endl;
  dathist2d = NULL;
  dathist = Data;
  dataSample = NULL;
  dataSample2D = NULL;

  if (getNDim()!=1) {std::cerr << "Trying to set a 1D 'data' histogram in a 2D sample - Quitting" << std::endl; throw;}

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  samplePDFFD_data = new double*[nYBins];
  for (int yBin=0;yBin<nYBins;yBin++) {
    samplePDFFD_data[yBin] = new double[nXBins];
    for (int xBin=0;xBin<nXBins;xBin++) {
      samplePDFFD_data[yBin][xBin] = Data->GetBinContent(xBin+1);
    }
  }
}

void samplePDFFDBase::addData(TH2D* Data) {
  std::cout << "adding 2D data histogram : " << Data->GetName() << " with " << Data->Integral() << " events" << std::endl;
  dathist2d = Data;
  dathist = NULL;
  dataSample = NULL;
  dataSample2D = NULL;

  if (getNDim()!=2) {std::cerr << "Trying to set a 2D 'data' histogram in a 1D sample - Quitting" << std::endl; throw;}

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

void samplePDFFDBase::setupSplines(fdmc_base *fdobj, const char *splineFile, int nutype, int signal) {
  // Only works with 2017 splines
  std::string str_splineFile(splineFile);
  int nevents = fdobj->nEvents;

  std::cout << "##################" << std::endl;
  std::cout << "I AM LOOKING FOR " << splineFile << std::endl;
  std::cout << "FDOBJ has " << nevents << std::endl;
  std::cout << "##################" << std::endl;

  std::cout << "Using 2020 splines" << std::endl;
  if (BinningOpt == 0){ // Splines binned in erec
	if (signal){
	  fdobj->splineFile = new splineFDBase((char*)splineFile, nutype, nevents, SampleDetID, xsecCov);
	}
	else{
	  fdobj->splineFile = new splineFDBase((char*)splineFile, nutype, nevents, SampleDetID, xsecCov);
	  if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)){
		std::cerr << "problem setting up splines in erec" << std::endl;
	  }
	}
  }
  else{ // Splines binned in 2D
	std::cout << "Creating splineFDBase" << std::endl;
	fdobj->splineFile = new splineFDBase((char*)splineFile, nutype, nevents, (double)BinningOpt, SampleDetID, xsecCov);
	if (!(nutype==1 || nutype==-1 || nutype==2 || nutype==-2)) {
	  std::cerr << "problem setting up splines in erec" << std::endl;
	} 
  }

  return;
}

double samplePDFFDBase::getEventWeight(int iSample, int iEntry)
{
  double totalweight = 1.0;
  for (int iParam=0;iParam<MCSamples[iSample].ntotal_weight_pointers[iEntry];iParam++) {
    totalweight *= *(MCSamples[iSample].total_weight_pointers[iEntry][iParam]);
  }

  return totalweight;
}

void samplePDFFDBase::fillSplineBins()
{
  for (int i = 0; i < (int)MCSamples.size(); ++i) {
    for (int j = 0; j < MCSamples[i].nEvents; ++j) {

      std::vector< std::vector<int> > EventSplines;
	
	  double erec = *(MCSamples[i].x_var[j]) * 1000;
	  switch (BinningOpt) {
		//ETA - factor of 1000 needed here because for some reason Erec is in units of MeV for Etrue-Erec splines and GeV for Etrue-Erec-theta splines?!?
		case 0: // splines binned in erec
		  MCSamples[i].splineFile->GetSplineBins(MCSamples[i].nutype, MCSamples[i].signal, *(MCSamples[i].rw_etru[j]), erec, MCSamples[i].enu_s_bin[j], MCSamples[i].xvar_s_bin[j]);
		  EventSplines = MCSamples[i].splineFile->getEventSplines(j, *(MCSamples[i].mode[j]), MCSamples[i].enu_s_bin[j], MCSamples[i].xvar_s_bin[j]); 
		  break;
		case 2:
		  MCSamples[i].splineFile->GetSplineBins(MCSamples[i].nutype, MCSamples[i].signal, *(MCSamples[i].rw_etru[j]), *(MCSamples[i].x_var[j]), *(MCSamples[i].y_var[j]), MCSamples[i].enu_s_bin[j], MCSamples[i].xvar_s_bin[j], MCSamples[i].yvar_s_bin[j]);
		  EventSplines = MCSamples[i].splineFile->getEventSplines(j, *(MCSamples[i].mode[j]), MCSamples[i].enu_s_bin[j], MCSamples[i].xvar_s_bin[j], MCSamples[i].yvar_s_bin[j]); 
		  break;
		case 3: // splines binned in erec-Q2 (not implemented yet!)
		  std::cerr << "ERROR : splines not available for erec-Q2 ! Code it in samplePDFSKBase.cpp if you have them " << std::endl ;
		default:
		  std::cout << "Error in assigning spline bins because BinningOpt = " << BinningOpt << std::endl;
		  break;
	  }

      MCSamples[i].nxsec_spline_pointers[j] = EventSplines.size();
      MCSamples[i].xsec_spline_pointers[j] = new const double*[MCSamples[i].nxsec_spline_pointers[j]];

      //Event Splines indexed as: param, mode, etrue, var1, [var2]
	  for (int spline=0;spline<MCSamples[i].nxsec_spline_pointers[j];++spline) {
		if (EventSplines[spline].size()==4) {
		  MCSamples[i].xsec_spline_pointers[j][spline] = MCSamples[i].splineFile->retPointer(EventSplines[spline][0],EventSplines[spline][1],EventSplines[spline][2],EventSplines[spline][3]);
		} else {
		  MCSamples[i].xsec_spline_pointers[j][spline] = MCSamples[i].splineFile->retPointer(EventSplines[spline][0],EventSplines[spline][1],EventSplines[spline][2],EventSplines[spline][3],EventSplines[spline][4]);
		}
	  }

	}
  }

  std::cout << "Filled spline bins" << std::endl;

  return;
}

double samplePDFFDBase::getLikelihood()
{
  if (samplePDFFD_data==NULL) {
      std::cerr << "data sample is empty!" << std::endl;
      return -1;
  }

  int nXBins = XBinEdges.size()-1;
  int nYBins = YBinEdges.size()-1;

  int xBin;
  int yBin;

  double negLogL = 0;

#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:negLogL) private(xBin, yBin)
#endif

  for (xBin=0;xBin<nXBins;xBin++) {
    double negLogLsample = 0.;

    for (yBin=0;yBin<nYBins;yBin++) {

      double MCPred = samplePDFFD_array[yBin][xBin];
      double DataVal = samplePDFFD_data[yBin][xBin];
      
      negLogLsample += getTestStatLLH(DataVal, MCPred);
    }

    negLogL += negLogLsample;
  }

  return negLogL;
}
