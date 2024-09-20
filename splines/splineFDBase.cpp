#include "splineFDBase.h"

//****************************************
splineFDBase::splineFDBase(covarianceXsec *xsec_)
              : SplineBase() {
//****************************************
  if (xsec_ == NULL) {
    MACH3LOG_ERROR("Trying to create splineSKBase with NULL covariance object");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  xsec = xsec_;

  // Keep these in class scope, important for using 1 monolith/sample!
  MonolithIndex = 0; //Keeps track of the monolith index we're on when filling arrays (declared here so we can have multiple FillSampleArray calls)
  CoeffIndex = 0; //Keeps track of our indexing the coefficient arrays [x, ybcd]
}
//****************************************
splineFDBase::~splineFDBase(){
//****************************************
  if(manycoeff_arr != nullptr) delete[] manycoeff_arr;
  if(xcoeff_arr != nullptr) delete[] xcoeff_arr;

}
//****************************************
void splineFDBase::cleanUpMemory() {
//****************************************

  //Call once everything's been allocated in samplePDFSKBase, cleans up junk from memory!
  //Not a huge saving but it's better than leaving everything up to the compiler
  MACH3LOG_INFO("Cleaning up spline memory");

  indexvec.clear();
  indexvec.shrink_to_fit();
  SplineFileParPrefixNames.clear();
  SplineFileParPrefixNames.shrink_to_fit();
  SplineBinning.clear();
  SplineBinning.shrink_to_fit();
  GlobalSystIndex.clear();
  GlobalSystIndex.shrink_to_fit();
  UniqueSystNames.clear();
  UniqueSystNames.shrink_to_fit();
  splinevec_Monolith.clear();
  splinevec_Monolith.shrink_to_fit();
  if(isflatarray != nullptr) delete isflatarray;
}

//****************************************
bool splineFDBase::AddSample(std::string SampleName, int NSplineDimensions, int DetID, std::vector<std::string> OscChanFileNames, std::vector<std::string> SplineVarNames)
//Adds samples to the large array
//****************************************
{
  SampleNames.push_back(SampleName);
  Dimensions.push_back(NSplineDimensions);
  DimensionLabels.push_back(SplineVarNames);
  DetIDs.push_back(DetID);

  int nSplineParam = xsec->GetNumSplineParamsFromDetID(DetID);
  nSplineParams.push_back(nSplineParam);

  //This holds the global index of the spline i.e. 0 -> _fNumPar
  std::vector<int> GlobalSystIndex_Sample = xsec->GetGlobalSystIndexFromDetID(DetID, kSpline);
  //Keep track of this for all the samples
  GlobalSystIndex.push_back(GlobalSystIndex_Sample);

  //std::vector<int> SplineParsIndex_Sample_temp = xsec->GetSplineParsIndexFromDetID(DetID);

  std::vector<std::string> SplineFileParPrefixNames_Sample = xsec->GetSplineParsNamesFromDetID(DetID);
  SplineFileParPrefixNames.push_back(SplineFileParPrefixNames_Sample);

  std::vector<std::vector<int>> SplineModeVecs_Sample = StripDuplicatedModes(xsec->GetSplineModeVecFromDetID(DetID));
  SplineModeVecs.push_back(SplineModeVecs_Sample);

  int nOscChan = OscChanFileNames.size();
  nOscChans.push_back(nOscChan);

  PrintSampleDetails(SampleName);

  std::vector<std::vector<TAxis *>> SampleBinning(nOscChan);
  for (int iOscChan = 0; iOscChan < nOscChan; iOscChan++)
  {
    std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;
    SampleBinning[iOscChan] = FindSplineBinning(OscChanFileNames[iOscChan], SampleName);
  }
  std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;
  SplineBinning.push_back(SampleBinning);

  BuildSampleIndexingArray(SampleName);
  PrintArrayDetails(SampleName);
  std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;

  FillSampleArray(SampleName, OscChanFileNames);
  std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;

  return true;
}

//****************************************
void splineFDBase::TransferToMonolith()
//****************************************
{
  PrepForReweight(); 
  MonolithSize = CountNumberOfLoadedSplines(false, 1);

  if(MonolithSize!=MonolithIndex){
    MACH3LOG_ERROR("Something's gone wrong when we tried to get the size of your monolith");
    MACH3LOG_ERROR("MonolithSize is {}", MonolithSize);
    MACH3LOG_ERROR("MonolithIndex is {}", MonolithIndex);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  MACH3LOG_INFO("Now transferring splines to a monolith if size {}", MonolithSize);

  uniquesplinevec_Monolith.resize(MonolithSize);
  weightvec_Monolith.resize(MonolithSize);
  isflatarray = new bool[MonolithSize];
  
  xcoeff_arr = new _float_[CoeffIndex];
  manycoeff_arr = new _float_[CoeffIndex*4];

  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  { // Loop over sample
    for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
    { // Loop over oscillation channels
      for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++)
      { // Loop over systematics
        for (unsigned int iMode = 0; iMode < indexvec[iSample][iOscChan][iSyst].size(); iMode++)
        { // Loop over modes
          for (unsigned int iVar1 = 0; iVar1 < indexvec[iSample][iOscChan][iSyst][iMode].size(); iVar1++)
          { // Loop over first dimension
            for (unsigned int iVar2 = 0; iVar2 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(); iVar2++)
            { // Loop over second dimension
              for (unsigned int iVar3 = 0; iVar3 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(); iVar3++)
              {
                int splineindex=indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3];
                weightvec_Monolith[splineindex] = 1.0;

                bool foundUniqueSpline = false;
                for (int iUniqueSyst = 0; iUniqueSyst < nUniqueSysts; iUniqueSyst++)
                {
                  if (SplineFileParPrefixNames[iSample][iSyst] == UniqueSystNames[iUniqueSyst])
                  {
                    uniquesplinevec_Monolith[splineindex] = iUniqueSyst;
                    foundUniqueSpline = true;
                  }
                }//unique syst loop end

                if (!foundUniqueSpline)
                {
                  MACH3LOG_ERROR("Unique spline index not found");
				  MACH3LOG_ERROR("For Spline {}", SplineFileParPrefixNames[iSample][iSyst]);
				  MACH3LOG_ERROR("Couldn't match {} with any of the following {} systs:", SplineFileParPrefixNames[iSample][iSyst], nUniqueSysts);
				  for (int iUniqueSyst = 0; iUniqueSyst < nUniqueSysts; iUniqueSyst++)
				  {
					MACH3LOG_ERROR("{},", UniqueSystNames.at(iUniqueSyst));
				  }//unique syst loop end				
                  throw MaCh3Exception(__FILE__ , __LINE__ );
                }

                int splineKnots;
                if(splinevec_Monolith[splineindex]!=NULL){
                  isflatarray[splineindex]=false;
                  splineKnots=splinevec_Monolith[splineindex]->GetNp();

                  //Now to fill up our coefficient arrayss
                  _float_* tmpXCoeffArr = new _float_[splineKnots];
                  _float_* tmpManyCoeffArr = new _float_[splineKnots*4];

                  int iCoeff=coeffindexvec[splineindex];
                  getSplineCoeff_SepMany(splineindex, tmpXCoeffArr, tmpManyCoeffArr);

                  #ifdef MULTITHREAD
                  #pragma omp parallel for
                  #endif
                  for(int i=0; i<splineKnots; i++){

                    if(tmpXCoeffArr[i]==-999){
                      std::cerr<<"ERROR : looks like we've got a bad X, index = "<<i<<std::endl;
                      throw MaCh3Exception(__FILE__ , __LINE__ );
                    }
                    xcoeff_arr[iCoeff+i]=tmpXCoeffArr[i];

                    for(int j=0; j<4; j++){
                      if(tmpManyCoeffArr[i*4+j]==-999){
                        std::cerr<<"Bad ybcd, index : "<<i<<", "<<j<<std::endl;
                        std::cerr<<"Param Values : "<<tmpManyCoeffArr[i*4]<<", "<<tmpManyCoeffArr[i*4+1]<<", "<<tmpManyCoeffArr[i*4+2]<<", "<<tmpManyCoeffArr[i*4+3]<<std::endl;
                        throw MaCh3Exception(__FILE__ , __LINE__ );
                      }
                    manycoeff_arr[(iCoeff+i)*4+j]=tmpManyCoeffArr[i*4+j];
                    }
                  }
                  delete tmpXCoeffArr;
                  delete tmpManyCoeffArr;
                }
                else {
                    isflatarray[splineindex]=true;
                }
              }//3d loop end
            }//2d loop end
          }//1d loop end
        }//mode loop end
      }//syst2 loop end
    }//osc loop end
  }//syst1 loop end

  return;
}

// *****************************************
void splineFDBase::Evaluate() {
// *****************************************

  // There's a parameter mapping that goes from spline parameter to a global parameter index
  // Find the spline segments
  FindSplineSegment();

  //KS: Huge MP loop over all valid splines
  CalcSplineWeights();

  //KS: Huge MP loop over all events calculating total weight
  ModifyWeights();

  return;
}

//****************************************
// ETA - find the spline segment that the current parameter
// value is in. This is now extremely similar to the
// function in SplineMonolith.cpp
void splineFDBase::FindSplineSegment()
//****************************************
{
  //HW okay let's try this, we delete+refill a new array which we'll fill with x-s for our segment
  #ifdef MULTITHREAD
  #pragma omp parallel //for schedule(dynamic)
  #endif
  for (int iSyst = 0; iSyst < nUniqueSysts; iSyst++)
  {
    int nPoints = UniqueSystNKnots[iSyst];
    std::vector<_float_> xArray = UniqueSystXPts[iSyst];

    // Get the variation for this reconfigure for the ith parameter
    int GlobalIndex = UniqueSystIndices[iSyst];

    _float_ xvar=_float_(xsec->getParProp(GlobalIndex));

    xVarArray[iSyst]=xvar;
    
    _int_ segment = 0;
	_int_ kHigh = nPoints - 1;

    //KS: We expect new segment is very close to previous
    const _int_ PreviousSegment = UniqueSystCurrSegment[iSyst];
    //KS: It is quite probable the new segment is same as in previous step so try to avoid binary search
    if( xArray[PreviousSegment+1] > xvar && xvar >= xArray[PreviousSegment] ){segment = PreviousSegment;}
    // If the variation is below the lowest saved spline point
	else if (xvar <= xArray[0]) {
	  segment = 0;
	  // If the variation is above the highest saved spline point
	} else if (xvar >= xArray[nPoints-1]) {
	  //CW: Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
	  segment = kHigh;
	  //KS: It is quite probable the new segment is same as in previous step so try to avoid binary search
	} else {
      // The top point we've got
      _int_ kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point

    if (segment >= nPoints-1 && nPoints > 1){segment = nPoints-2;}
    UniqueSystCurrSegment[iSyst] = segment; 
      
//#ifdef DEBUG
//    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
//      std::cerr << "Found a segment which is _ABOVE_ the variation!" << std::endl;
//      std::cerr << "IT SHOULD ALWAYS BE BELOW! (except when segment 0)" << std::endl;
//      std::cerr << "Spline: "<< i << std::endl;
//
//      std::cerr << "Found segment   = " << segment << std::endl;
//      std::cerr << "Doing variation = " << xvar << std::endl;
//      std::cerr << "x in spline     = " << SplineInfoArray[i].xPts[segment] << std::endl;
//      for (__int__ j = 0; j < SplineInfoArray[j].nPts; ++j) {
//        std::cerr << "    " << j << " = " << SplineInfoArray[i].xPts[j] << std::endl;
//      }
//      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
//      throw;
//    }
//#endif
  } //end loop over params
}

//****************************************
void splineFDBase::CalcSplineWeights()
//****************************************
{
  #ifdef MULTITHREAD
  #pragma omp parallel for simd
  #endif
  for (unsigned int iCoeff = 0; iCoeff < uniquecoeffindices.size(); iCoeff++)
  {

    int iSpline = uniquecoeffindices[iCoeff];
    short int uniqueIndex=uniquesplinevec_Monolith[iSpline];
    short int currentsegment=UniqueSystCurrSegment[uniqueIndex];

    int segCoeff = coeffindexvec[iSpline]+currentsegment;

    // These are what we can extract from the TSpline3
    _float_ x = xcoeff_arr[segCoeff];
    _float_ y = manycoeff_arr[(segCoeff)*4+kCoeffY];
    _float_ b = manycoeff_arr[(segCoeff)*4+kCoeffB];
    _float_ c = manycoeff_arr[(segCoeff)*4+kCoeffC];
    _float_ d = manycoeff_arr[(segCoeff)*4+kCoeffD];

    // Get the variation for this reconfigure for the ith parameter
    _float_ xvar = xVarArray[uniqueIndex];
    // The Delta(x)
    _float_ dx = xvar - x;

    //Speedy 1% time boost https://en.cppreference.com/w/c/numeric/math/fma (see ND code!)
    _float_ weight = fmaf(dx, fmaf(dx, fmaf(dx, d, c), b), y);
    //This is the speedy version of writing dx^3+b*dx^2+c*dx+d


    //ETA - do we need this? We check later for negative weights and I wonder if this is even
    //possible with the fmaf line above?
    if(weight<0){weight=0;}  //Stops is getting negative weights

    weightvec_Monolith[iSpline]=double(weight);
  }
}

//****************************************
void splineFDBase::BuildSampleIndexingArray(std::string SampleName)
//Creates an array to be filled with monolith indexes for each sample (allows for indexing between 7D binning and 1D Vector)
//Only need 1 indexing array everything else interfaces with this to get binning properties
//****************************************
{  
  int iSample = getSampleIndex(SampleName);
  int nSplineSysts = nSplineParams[iSample];
  int nOscChannels = nOscChans[iSample];
  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>> indexvec_OscChan;
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++)
  { // Loop over oscillation channels
    std::vector<std::vector<std::vector<std::vector<std::vector<int >>>>> indexvec_Syst;
    for (int iSyst = 0; iSyst < nSplineSysts; iSyst++)
    { // Loop over systematics
      std::vector<std::vector<std::vector<std::vector<int >>>> indexvec_Mode;
      int nModesInSyst = SplineModeVecs[iSample][iSyst].size();
      for (int iMode = 0; iMode < nModesInSyst; iMode++)
      { // Loop over modes
        std::vector<std::vector<std::vector<int >>> indexvec_Var1;
        for (int iVar1 = 0; iVar1 < (SplineBinning[iSample][iOscChan][0])->GetNbins(); iVar1++)
        { // Loop over first dimension
          std::vector<std::vector<int >> indexvec_Var2;
          for (int iVar2 = 0; iVar2 < (SplineBinning[iSample][iOscChan][1])->GetNbins(); iVar2++)
          { // Loop over second dimension
            std::vector<int> indexvec_Var3;
            for (int iVar3 = 0; iVar3 < (SplineBinning[iSample][iOscChan][2])->GetNbins(); iVar3++)
            { // Loop over third dimension
              indexvec_Var3.push_back(0); //Don't start counting yet!
            } // end iVar3 loop
            indexvec_Var2.push_back(indexvec_Var3);
          } // end iVar2 loop
          indexvec_Var1.push_back(indexvec_Var2);
        } // end iVar1 loop
        indexvec_Mode.push_back(indexvec_Var1);
      } // end of iMode loop
      indexvec_Syst.push_back(indexvec_Mode);
    } // end of iOscChan loop
    indexvec_OscChan.push_back(indexvec_Syst);
  } // end of iSyst loop
  indexvec.push_back(indexvec_OscChan);
}

//****************************************
std::vector<TAxis *> splineFDBase::FindSplineBinning(std::string FileName, std::string SampleName)
//****************************************
{
  std::vector<TAxis *> ReturnVec;
  int iSample=getSampleIndex(SampleName);

  TFile *File = new TFile(FileName.c_str());
  if (!File || File->IsZombie())
  {
    std::cerr << "File " << FileName << " not found" << std::endl;
    std::cerr << "This is caused by something here! "<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  std::cout << "Finding binning for:" << std::endl;
  std::cout << FileName << std::endl;

  bool isHist2D = false;
  bool isHist3D = false;

  TH2F *Hist2D = NULL;
  TH3F *Hist3D = NULL;

  TObject *Obj = File->Get("dev_tmp_0_0");
  if (!Obj)
  {
    Obj = File->Get("dev_tmp.0.0");
    if (!Obj)
    {
      std::cerr << "Error: could not find dev_tmp_0_0 in spline file. Spline binning will not be set!" << std::endl;
      std::cerr << "FileName: " << FileName << std::endl;
      std::cerr << "0_0, I'm here! "<<__FILE__<<" : "<<__LINE__<<std::endl;
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  if (Obj->IsA() == TH2F::Class())
  {
    isHist2D = true;
  }

  //For T2K annoyingly all objects are TH3Fs
  if (Obj->IsA() == TH3F::Class())
  {
    isHist3D = true;
  }

  if (!isHist2D && !isHist3D)
  {
    std::cerr << "Object doesn't inherit from either TH2D and TH3D - Odd A" << std::endl;
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (isHist2D)
  {
    if (Dimensions[iSample] != 2)
    {
      std::cerr << "Trying to load a 2D spline template when nDim=" << Dimensions[iSample] << std::endl;
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    Hist2D = (TH2F *)File->Get("dev_tmp_0_0");
  }

  if (isHist3D)
  {

    if (Dimensions[iSample] != 3 && Hist3D->GetZaxis()->GetNbins() != 1)
    {
      std::cerr << "Trying to load a 3D spline template when nDim=" << Dimensions[iSample] << std::endl;
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    Hist3D = (TH3F *)Obj->Clone();
  }

  int nDummyBins = 1;
  double *DummyEdges = new double[2];
  DummyEdges[0] = -1e15;
  DummyEdges[1] = 1e15;
  TAxis *DummyAxis = new TAxis(nDummyBins, DummyEdges);

  if (Dimensions[iSample] == 2)
  {
	if(isHist2D){
	  ReturnVec.push_back((TAxis *)(Hist2D->GetXaxis())->Clone());
	  ReturnVec.push_back((TAxis *)(Hist2D->GetYaxis())->Clone());
	  ReturnVec.push_back((TAxis *)(DummyAxis)->Clone());
	}
	else if(isHist3D){
	  ReturnVec.push_back((TAxis *)(Hist3D->GetXaxis())->Clone());
	  ReturnVec.push_back((TAxis *)(Hist3D->GetYaxis())->Clone());
	  ReturnVec.push_back((TAxis *)(DummyAxis)->Clone());
	}
  }
  else if (Dimensions[iSample] == 3)
  {
    ReturnVec.push_back((TAxis *)(Hist3D->GetXaxis())->Clone());
    ReturnVec.push_back((TAxis *)(Hist3D->GetYaxis())->Clone());
    ReturnVec.push_back((TAxis *)(Hist3D->GetZaxis())->Clone());
  }
  else
  {
    std::cerr << "Number of dimensions not valid! Given:" << Dimensions[iSample] << std::endl;
    std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  for (unsigned int iAxis = 0; iAxis < ReturnVec.size(); ++iAxis)
  {
    PrintBinning(ReturnVec[iAxis]);
  }

  MACH3LOG_INFO("Left PrintBinning now tidying up");
  //This could be NULL if 2D
  if(isHist2D){
	delete Hist2D;
  } else {
    delete Hist3D;
  }

  File->Close();
  delete File;
  delete DummyAxis;
  return ReturnVec;
}

//****************************************
int splineFDBase::CountNumberOfLoadedSplines(bool NonFlat, int Verbosity)
//****************************************
{
  int SampleCounter_NonFlat = 0;
  int SampleCounter_All = 0;
  int FullCounter_NonFlat = 0;
  int FullCounter_All = 0;

  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  { // Loop over systematics
    SampleCounter_NonFlat = 0;
    SampleCounter_All = 0;
    std::string SampleName = SampleNames[iSample];
    for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
    { // Loop over oscillation channels
      for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++)
      { // Loop over systematics
   		for (unsigned int iMode = 0; iMode < indexvec[iSample][iOscChan][iSyst].size(); iMode++)
        { // Loop over modes
          for (unsigned int iVar1 = 0; iVar1 < indexvec[iSample][iOscChan][iSyst][iMode].size(); iVar1++)
          { // Loop over first dimension
            for (unsigned int iVar2 = 0; iVar2 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(); iVar2++)
            { // Loop over second dimension
              for (unsigned int iVar3 = 0; iVar3 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(); iVar3++)
              { // Loop over third dimension
                if (isValidSplineIndex(SampleName, iOscChan, iSyst, iMode, iVar1, iVar2, iVar3))
                {
                  int splineindex = indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3];
                  if (splinevec_Monolith[splineindex] != NULL)
                  {
                    SampleCounter_NonFlat += 1;
                  }
                  SampleCounter_All += 1;
                }
              }
            }
          }
        }
      }
    }
    if (Verbosity > 0)
    {
      std::cout << std::setw(10) << SampleNames[iSample] << " has " << std::setw(10) << SampleCounter_All << " splines, of which " << std::setw(10) << SampleCounter_NonFlat << " are not flat" << std::endl;
    }

    FullCounter_NonFlat += SampleCounter_NonFlat;
    FullCounter_All += SampleCounter_All;
  }

  if (Verbosity > 0)
  {
    MACH3LOG_INFO("Total number of splines loaded: {}", FullCounter_All);
    MACH3LOG_INFO("Total number of non-flat splines loaded: {}", FullCounter_NonFlat);
  }

  if (NonFlat)
  {
    return FullCounter_NonFlat;
  }
  else
  {
    return FullCounter_All;
  }
}

//****************************************
void splineFDBase::PrepForReweight()
//****************************************
{

  std::vector<TSpline3_red*> UniqueSystSplines;

  // DB Find all the Unique systs across each sample and oscillation channel
  //    This assumes that each occurence of the same systematic spline has the same knot spacing
  //    Which is a reasonable assumption for the current implementation of spline evaluations
  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  { // Loop over systematics
    for (unsigned int iSyst = 0; iSyst < indexvec[iSample][0].size(); iSyst++)
    { // Loop over systematics
      std::string SystName = SplineFileParPrefixNames[iSample][iSyst];
      bool FoundSyst = false;

	  //ETA - this always seems to be empty to begin with??
	  //so this loop never gets used?
      for (unsigned int iFoundSyst = 0; iFoundSyst < UniqueSystNames.size(); iFoundSyst++)
      {
        if (SystName == UniqueSystNames[iFoundSyst])
        {
          FoundSyst = true;
        }
      }
      if (!FoundSyst)
      {
        bool FoundNonFlatSpline = false;
        //This is all basically to check that you have some resposne from a spline
        //systematic somewhere
        for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
        { // Loop over oscillation channels
          for (unsigned int iMode = 0; iMode < indexvec[iSample][iOscChan][iSyst].size(); iMode++)
          { // Loop over modes
            for (unsigned int iVar1 = 0; iVar1 < indexvec[iSample][iOscChan][iSyst][iMode].size(); iVar1++)
            { // Loop over first dimension
              for (unsigned int iVar2 = 0; iVar2 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(); iVar2++)
              { // Loop over second dimension
                for (unsigned int iVar3 = 0; iVar3 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(); iVar3++)
                { // Loop over third dimension
                  int splineindex=indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3];
                  if (splinevec_Monolith[splineindex] != NULL)
                  {
                    UniqueSystSplines.push_back(splinevec_Monolith[splineindex]);
                    UniqueSystIndices.push_back(GlobalSystIndex[iSample][iSyst]);
                    FoundNonFlatSpline = true;
                  }
                  if (FoundNonFlatSpline)
                  {
                    break;
                  }
                }//3D loop end
                if (FoundNonFlatSpline)
                {
                  break;
                }
              }//2D loop end
              if (FoundNonFlatSpline)
              {
                break;
              }
            }//1D loop end
            if (FoundNonFlatSpline)
            {
              break;
            }
          }//mode loop end
          if (FoundNonFlatSpline)
          {
            break;
          }
        }//osc loop end
		 //ETA - only push back unique name if a non-flat response has been found
		if(FoundNonFlatSpline){
		  UniqueSystNames.push_back(SystName);
		}

        if (!FoundNonFlatSpline)
        {
		  MACH3LOG_INFO("{} syst has no response in sample {}", SystName, iSample);
		  MACH3LOG_INFO("Whilst this isn't neccessarily a problem, it seems odd");
          continue;
        }
      }
    }//Syst loop end
  }
  
  nUniqueSysts = UniqueSystSplines.size();

  // DB Find the number of splines knots which assumes each instance of the syst has the same number of knots
  UniqueSystNKnots.resize(nUniqueSysts);
  UniqueSystCurrSegment.resize(nUniqueSysts);
  UniqueSystXPts.resize(nUniqueSysts);
  xVarArray=new _float_[nUniqueSysts];

  for (int iSpline = 0; iSpline < nUniqueSysts; iSpline++)
  {
    UniqueSystNKnots[iSpline] = UniqueSystSplines[iSpline]->GetNp();
    UniqueSystXPts[iSpline].resize(UniqueSystNKnots[iSpline]);
    for (int iKnot = 0; iKnot < UniqueSystNKnots[iSpline]; iKnot++)
    {
      _float_ xPoint = -999;
      _float_ yPoint = -999;
      UniqueSystSplines[iSpline]->GetKnot(iKnot, xPoint, yPoint);
      if (xPoint == -999 || yPoint == -999)
      {
        std::cerr << "Something has gone wrong in the knot finding" << std::endl;
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      UniqueSystXPts[iSpline][iKnot] = xPoint;
    }
	//ETA - let this just be set as the first segment by default
    UniqueSystCurrSegment[iSpline] = 0;
    xVarArray[iSpline]=0;
  }
  

  std::cout << "nUniqueSysts:" << nUniqueSysts << " -----------------" << std::endl;
  std::cout << std::endl;

  std::cout << std::setw(15) << "Spline Index"
            << " | " << std::setw(20) << "Syst Name"
            << " | " << std::setw(6) << "nKnots" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  for (int iUniqueSyst = 0; iUniqueSyst < nUniqueSysts; iUniqueSyst++)
  {
    std::cout << std::setw(15) << iUniqueSyst << " | " << std::setw(20) << UniqueSystNames[iUniqueSyst] << " | " << std::setw(6) << UniqueSystNKnots[iUniqueSyst] << std::endl;
  }
  std::cout << std::endl;

  //ETA
  //Isn't this just doing what CountNumberOfLoadedSplines() does?
  int nCombinations_FlatSplines = 0;
  int nCombinations_All = 0;
  // DB Now actually loop over splines to determine which are all NULL
  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  { // Loop over systematics
	for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
	{ // Loop over oscillation channels
	  for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++)
	  { // Loop over systematics
		for (unsigned int iMode = 0; iMode < indexvec[iSample][iOscChan][iSyst].size(); iMode++)
		{ // Loop over modes
		  //bool isFlat = false;
		  for (unsigned int iVar1 = 0; iVar1 < indexvec[iSample][iOscChan][iSyst][iMode].size(); iVar1++)
		  { // Loop over first dimension
			for (unsigned int iVar2 = 0; iVar2 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(); iVar2++)
			{ // Loop over second dimension
			  for (unsigned int iVar3 = 0; iVar3 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(); iVar3++)
			  { // Loop over third dimension
				int splineindex=indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3];
				if (splinevec_Monolith[splineindex] != NULL)
				{
				  nCombinations_All += 1;
				} else{
				  nCombinations_FlatSplines += 1;
				  nCombinations_All += 1;
				}
			  }
			}
		  }
		}
	  }
	}
  }
  // We need to grab the maximum number of knots
  MACH3LOG_INFO("Number of combinations of Sample, OscChan, Syst and Mode which have entirely flat response: {} / {}", nCombinations_FlatSplines, nCombinations_All);
  
}

//****************************************
// Rather work with spline coefficients in the splines, let's copy ND and use coefficient arrays
void splineFDBase::getSplineCoeff_SepMany(int splineindex, _float_* &xArray, _float_* &manyArray){
//****************************************
  // Initialise all arrays to 1.0
  int nPoints;
  //No point evaluating a flat spline
  nPoints = splinevec_Monolith[splineindex]->GetNp();

  for (int i = 0; i < nPoints; i++) {
    xArray[i] = 1.0;
    for (int j = 0; j < 4; j++) {
      manyArray[i*4+j] = 1.0;
    }
  }

  for(int i=0; i<nPoints; i++){
    // Spline coefficients to be
    _float_ x = -999.99;
    _float_ y = -999.99;
    _float_ b = -999.99;
    _float_ c = -999.99;
    _float_ d = -999.99;
    splinevec_Monolith[splineindex]->GetCoeff(i, x, y, b, c, d);
    //Let's save some memory and store them as floats! (It's a surprise tool that will help with GPU later)
    xArray[i]=_float_(x);

    //Might as well copy ND here and 
    xArray[i] = _float_(x);
    manyArray[i*4] = _float_(y); // 4 because manyArray stores y,b,c,d
    manyArray[i*4+1] = _float_(b);
    manyArray[i*4+2] = _float_(c);
    manyArray[i*4+3] = _float_(d);
    
    if((xArray[i] == -999) | (manyArray[i*4] == -999) | (manyArray[i*4+1] == -999) | (manyArray[i*4+2] == -999) | (manyArray[i*4+3] == -999)){
      MACH3LOG_ERROR("*********** Bad params in getSplineCoeff_SepMany() ************");
      MACH3LOG_ERROR("pre cast to _float_ (x, y, b, c, d) = {}, {}, {}, {}, {}",x, y, b, c, d);
      MACH3LOG_ERROR("post cast to float (x, y, b, c, d) = {}, {}, {}, {}, {}",xArray[i], manyArray[i*4], manyArray[i*4+1], manyArray[i*4+2], manyArray[i*4+3]);	    
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }
  //We now clean up the splines!
  delete splinevec_Monolith[splineindex];
  splinevec_Monolith[splineindex] = NULL;
}

//****************************************
//ETA - this may need to be virtual and then we can define this in the experiment.
//Equally though could just use KinematicVariable to map back
std::string splineFDBase::getDimLabel(int iSample, unsigned int Axis)
//****************************************
{
  if(Axis > DimensionLabels[iSample].size()){
    MACH3LOG_ERROR("The spline Axis you are trying to get the label of is larger than the number of dimensions");
    MACH3LOG_ERROR("You are trying to get axis {} but have only got {}", Axis, Dimensions[iSample]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return DimensionLabels.at(iSample).at(Axis);
}

//Returns sample index in 
int splineFDBase::getSampleIndex(std::string SampleName){
  int SampleIndex = -1;
  for (unsigned int iSample = 0; iSample < SampleNames.size(); iSample++)
  {
    if (SampleName == SampleNames[iSample])
    {
      SampleIndex = iSample;
    }
  }
  if (SampleIndex == -1)
  {
    MACH3LOG_ERROR("Sample name not found : {}", SampleName);	  
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return SampleIndex;
}

//****************************************
void splineFDBase::PrintSampleDetails(std::string SampleName)
//****************************************
{
  int iSample = getSampleIndex(SampleName);

  MACH3LOG_INFO("Details about sample: {:<20}", SampleNames[iSample]);
  MACH3LOG_INFO("\t Dimension: {:<35}", Dimensions[iSample]);
  MACH3LOG_INFO("\t DetID: {:<35}", DetIDs[iSample]);
  MACH3LOG_INFO("\t nSplineParam: {:<35}", nSplineParams[iSample]);
  MACH3LOG_INFO("\t nOscChan: {:<35}", nOscChans[iSample]);

  return;
}

//****************************************
void splineFDBase::PrintArrayDetails(std::string SampleName)
//****************************************
{
  int iSample = getSampleIndex(SampleName);
  int nOscChannels = indexvec[iSample].size();
  MACH3LOG_INFO("Sample {} has {} oscillation channels", iSample, nOscChannels);	
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++)
  {
    int nSysts = indexvec[iSample][iOscChan].size();
    MACH3LOG_INFO("Oscillation channel {} has {} systematics", iOscChan, nSysts);	  
  }
}

//****************************************
void splineFDBase::PrintArrayDimension()
//****************************************
{
  std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;
  std::cout << "Array dimensions.." << std::endl;
  std::cout << std::endl;

  std::cout << std::setw(20) << "nSamples:" << indexvec.size() << std::endl;
  std::cout << std::endl;

  std::cout << std::setw(20) << "nOscChans:";
  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  {
    std::cout << indexvec[iSample].size() << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << std::setw(20) << "nSysts:" << std::endl;
  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  {
    std::cout << "\t"
              << "Sample:" << iSample << "\t";
    for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
    {
      std::cout << indexvec[iSample][iOscChan].size() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << std::setw(20) << "nModes:" << std::endl;
  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++)
  {
    std::cout << "\t"
              << "Sample:" << iSample << "\t"
              << "--------------------------" << std::endl;
    for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++)
    {
      std::cout << "\t\t"
                << "OscChan:" << iOscChan << "\t";
      for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++)
      {
        std::cout << indexvec[iSample][iOscChan][iSyst].size() << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;
}

//****************************************
bool splineFDBase::isValidSplineIndex(std::string SampleName, int iOscChan, int iSyst, int iMode, int iVar1, int iVar2, int iVar3)
//****************************************
{

  int iSample=getSampleIndex(SampleName);
  bool isValid = true;

  if (iSample < 0 || iSample >= (int)indexvec.size())
  {
    MACH3LOG_ERROR("Sample index is invalid! 0 <= Index < {} ", indexvec.size());
    isValid = false;
  }

  if (iOscChan < 0 || iOscChan >= (int)indexvec[iSample].size())
  {
    MACH3LOG_ERROR("OscChan index is invalid! 0 <= Index < {} ", indexvec[iSample].size());
    isValid = false;
  }

  if (iSyst < 0 || iSyst >= (int)indexvec[iSample][iOscChan].size())
  {
    MACH3LOG_ERROR("Syst index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan].size());
    isValid = false;
  }

  if (iMode < 0 || iMode >= (int)indexvec[iSample][iOscChan][iSyst].size())
  {
    MACH3LOG_ERROR("Mode index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst].size());
    isValid = false;
  }

  if (iVar1 < 0 || iVar1 >= (int)indexvec[iSample][iOscChan][iSyst][iMode].size())
  {
    MACH3LOG_ERROR("Var1 index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst][iMode].size());	  
    isValid = false;
  }

  if (iVar2 < 0 || iVar2 >= (int)indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size())
  {
    MACH3LOG_ERROR("Var2 index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size());
    isValid = false;
  }

  if (iVar3 < 0 || iVar3 >= (int)indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size())
  {
    MACH3LOG_ERROR("Var3 index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size());
    isValid = false;
  }

  if (!isValid)
  { 
    MACH3LOG_ERROR("Given iSample: {}", iSample);
    MACH3LOG_ERROR("Given iOscChan: {}", iOscChan);
    MACH3LOG_ERROR("Given iSyst: {}", iSyst);
    MACH3LOG_ERROR("Given iMode: {}", iMode);
    MACH3LOG_ERROR("Given iVar1: {}", iVar1);
    MACH3LOG_ERROR("Given iVar2: {}", iVar2);
    MACH3LOG_ERROR("Given iVar3: {}", iVar3);	
    MACH3LOG_ERROR("Come visit me at : {} : {}", __FILE__, __LINE__);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  return true;
}

//****************************************
void splineFDBase::PrintBinning(TAxis *Axis)
//****************************************
{
  const int NBins = Axis->GetNbins();
  const double *BinEdges = Axis->GetXbins()->GetArray();
  
  std::cout << "\t";
  for (int iBin = 0; iBin < (NBins + 1); iBin++)
  {
    std::cout << BinEdges[iBin] << " ";
  }
  std::cout << std::endl;
  return;
}
