#include "splineFDBase.h"

//****************************************
splineFDBase::splineFDBase(covarianceXsec *xsec_)
              : SplineBase() {
//****************************************
  if (xsec_ == NULL) {
    std::cerr << "Trying to create splineSKBase with NULL covariance object" << std::endl;
    throw;
  }
  xsec = xsec_;

  // Keep these in class scope, important for using 1 monolith/sample!
  MonolithIndex=0; //Keeps track of the monolith index we're on when filling arrays (declared here so we can have multiple FillSampleArray calls)
  CoeffIndex=0; //Keeps track of our indexing the coefficient arrays [x, ybcd]
}

//****************************************
bool splineFDBase::AddSample(std::string SampleName, int BinningOpt, int DetID, std::vector<std::string> OscChanFileNames)
//Adds samples to the large array
//****************************************
{
  SampleNames.push_back(SampleName);
  BinningOpts.push_back(BinningOpt);
  Dimensions.push_back(getNDim(BinningOpt));
  DetIDs.push_back(DetID);

  int nSplineParam = xsec->GetNumSplineParamsFromDetID(DetID);
  std::cout << "FOund " << nSplineParam << " spline parameters" << std::endl;
  nSplineParams.push_back(nSplineParam);

  std::cout << "Filling with GetSplineParsIndexFromDetID" << std::endl;
  std::vector<int> SplineParsIndex_Sample = xsec->GetSplineParsIndexFromDetID(DetID);
  SplineParsIndex.push_back(SplineParsIndex_Sample);

  std::cout << "Filling with GetFDSplineFileParsNamesFromDetID" << std::endl;
  std::vector<std::string> SplineFileParPrefixNames_Sample = xsec->GetFDSplineFileParsNamesFromDetID(DetID);
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
  std::cout << "TRANSFERING TO MONOLITH!!" << std::endl;
  //FindUniqueModes();
  PrepForReweight(); 
  MonolithSize = CountNumberOfLoadedSplines();

  if(MonolithSize!=MonolithIndex){
    std::cerr<<"Something's gone wrong when we tried to get the size of your monolith"<<std::endl;
	std::cout << "MonolishSize is " << MonolithSize << std::endl;
	std::cout << "MonolithIndex is " << MonolithIndex << std::endl;
    std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }
  uniquesplinevec_Monolith.reserve(MonolithSize);
  weightvec_Monolith.reserve(MonolithSize);
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
                  std::cerr << "Unique spline index not found" << std::endl;
                  std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
                  throw;
                }
                int splineKnots;
                if(splinevec_Monolith[splineindex]!=NULL){
				  //std::cout << "Looking at splineindex " << splineindex << "and found a non-flat spline (wahoo!)" << std::endl;
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
                      std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
                      throw;
                    }
                    xcoeff_arr[iCoeff+i]=tmpXCoeffArr[i];

                    for(int j=0; j<4; j++){
                      if(tmpManyCoeffArr[i*4+j]==-999){
                        std::cerr<<"Bad ybcd, index : "<<i<<", "<<j<<std::endl;
                        std::cerr<<"Param Values : "<<tmpManyCoeffArr[i*4]<<", "<<tmpManyCoeffArr[i*4+1]<<", "<<tmpManyCoeffArr[i*4+2]<<", "<<tmpManyCoeffArr[i*4+3]<<std::endl;
                        std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
                        throw;
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

    _float_ xvar=_float_(xsec->calcReWeight(GlobalIndex));
    xVarArray[iSyst]=xvar;
    
    //Rather than starting from 0 everytime, let's use what we know!    
    int segment = UniqueSystCurrSegment[iSyst];
    int kHigh = nPoints - 1;
    // If the variation is below the lowest saved spline point
    if(segment<nPoints-1 && segment>=0 && xvar>=xArray[segment] && xvar<xArray[segment+1]) //We're still in the same place!
      {
	      continue;
      }
    //Now check to see if it's the segment below it
    else if(segment >= 1  && segment<=nPoints && xvar>=xArray[segment-1] && xvar<xArray[segment])
      {
	      segment--;
      }
    //Okay what if we're above it
    else if(segment <= nPoints-2 && segment>=0 && xvar>=xArray[segment+1] && xvar<xArray[segment+2])
      {
	      segment++;
      }
    else if (xvar <= xArray[0])
      {
	      segment=0;
	// If the variation is above the highest saved spline point
      }
    else if (xvar >= xArray[nPoints - 1])
    {
      // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search
    }
    //Now we resort to binary search!
    else
    {
      segment=0;
      // The top point we've got
      int kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1)
      {
        // Increment the half-step
        kHalf = (segment + kHigh) / 2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf])
        {
          segment = kHalf;
          // Else move kHigh down
        }
        else
        {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    }   // End the else: we've now found our point
    if (segment >= nPoints - 1 && nPoints > 1){
      segment = nPoints - 2;}
    // Save the segment for the ith parameter
    UniqueSystCurrSegment[iSyst] = segment;
  }
  // std::cout << "#----------------------------------------------------------------------------------------------------------------------------------#" << std::endl;
}

//****************************************
void splineFDBase::CalcSplineWeights()
//****************************************
{
  #ifdef MULTITHREAD
  #pragma omp parallel for// schedule(dynamic)
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
    std::cerr << "Error: could not find dev_tmp_0_0 in spline file. Spline binning will not be set!" << std::endl;
    std::cerr << "FileName: " << FileName << std::endl;
    std::cerr << "0_0, I'm here! "<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  if (Obj->IsA() == TH2F::Class())
  {
    isHist2D = true;
  }
  if (Obj->IsA() == TH3F::Class())
  {
    isHist3D = true;
  }

  if (isHist2D && isHist3D)
  {
    std::cerr << "Object inherits from both TH2D and TH3D - Odd" << std::endl;
    std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }
  if (!isHist2D && !isHist3D)
  {
    std::cerr << "Object doesn't inherit from either TH2D and TH3D - Odd A" << std::endl;
    std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  if (isHist2D)
  {
    if (Dimensions[iSample] != 2)
    {
      std::cerr << "Trying to load a 2D spline template when nDim=" << Dimensions[iSample] << std::endl;
      std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
      throw;
    }
    Hist2D = (TH2F *)File->Get("dev_tmp_0_0");
  }
  if (isHist3D)
  {
    if (Dimensions[iSample] != 3)
    {
      std::cerr << "Trying to load a 3D spline template when nDim=" << Dimensions[iSample] << std::endl;
      std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
      throw;
    }
    Hist3D = (TH3F *)File->Get("dev_tmp_0_0");
  }

  int nDummyBins = 1;
  double *DummyEdges = new double[2];
  DummyEdges[0] = -1e15;
  DummyEdges[1] = 1e15;
  TAxis *DummyAxis = new TAxis(nDummyBins, DummyEdges);

  if (Dimensions[iSample] == 2)
  {
    ReturnVec.push_back((TAxis *)(Hist2D->GetXaxis())->Clone());
    ReturnVec.push_back((TAxis *)(Hist2D->GetYaxis())->Clone());
    ReturnVec.push_back((TAxis *)(DummyAxis)->Clone());
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

  for (unsigned int iAxis = 0; iAxis < ReturnVec.size(); iAxis++)
  {
    std::cout << "Stored Var " << iAxis << " (" << getDimLabel(BinningOpts[iSample], iAxis) << ") Spline Binning for sample " << SampleNames[iSample] << ":" << std::endl;
    PrintBinning(ReturnVec[iAxis]);
  }

  //This could be NULL if 2D
  if(isHist2D){
	delete Hist2D;
  } else {
    delete Hist3D;
  }

  File->Close();
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
    std::cout << "Total number of splines loaded:" << FullCounter_All << std::endl;
    std::cout << "Total number of non-flat splines loaded:" << FullCounter_NonFlat << std::endl;
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

  std::vector<TSpline3_red *> UniqueSystSplines;

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
        UniqueSystNames.push_back(SystName);

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
                    UniqueSystIndices.push_back(SplineParsIndex[iSample][iSyst]);
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

        if (!FoundNonFlatSpline)
        {
          std::cerr << SystName << " syst has no response in sample " << iSample << std::endl;
          std::cerr << "Whilst this isn't neccessarily a problem, it seems odd" << std::endl;
          std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
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
      double xPoint = -999;
      double yPoint = -999;
      UniqueSystSplines[iSpline]->GetKnot(iKnot, xPoint, yPoint);
      if (xPoint == -999 || yPoint == -999)
      {
        std::cerr << "Something has gone wrong in the knot finding" << std::endl;
        std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
        throw;
      }
      UniqueSystXPts[iSpline][iKnot] = xPoint;
    }
    UniqueSystCurrSegment[iSpline] = -999;
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

  std::cout << "Number of combinations of Sample, OscChan, Syst and Mode which have entirely flat response:" << nCombinations_FlatSplines << " / " << nCombinations_All << std::endl;
}

// Rather work with spline coefficients in the splines, let's copy ND and use coefficient arrays
void splineFDBase::getSplineCoeff_SepMany(int splineindex, _float_* &xArray, _float_* &manyArray){

  // Initialise all arrays to 1.0
  int nPoints;
  //No point evalutating a flat spline

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
      std::cerr << "*********** Bad params in getSplineCoeff_SepMany() ************"<<std::endl;
      std::cerr << "pre cast to _float_ (x, y, b, c, d) = "<<x<<", "<<y<<", "<<b<<", "<<c<<", "<<d<<std::endl;
      std::cerr << "post cast to float (x, y, b, c, d) = "<<xArray[i]<<", "<<manyArray[i*4]<<", "<<manyArray[i*4+1]<<", "<<manyArray[i*4+2]<<", "<<manyArray[i*4+3]<<std::endl;
      std::cerr<<__FILE__<<"::"<<__LINE__<<std::endl;
      std::cerr << "***************************************************************"<<std::endl;
     throw;
      }
    }
  //We now clean up the splines!
  delete splinevec_Monolith[splineindex];
  splinevec_Monolith[splineindex] = NULL;
}

//****************************************
int splineFDBase::getNDim(int BinningOpt)
//****************************************
{
  int ReturnVal = -1;

  switch (BinningOpt)
  {
  case 0:
    ReturnVal = 2;
    break;
  case 1:
  case 2:
  case 3:
  case 4:
    ReturnVal = 3;
    break;
  default:
    std::cout << "Unrecognised BinningOpt = " << BinningOpt << std::endl;
    throw;
  }

  return ReturnVal;
}

//ETA - this may need to be virtual and then we can define this in the experiment.
//Equally though could just use KinematicVariable to map back
//****************************************
TString splineFDBase::getDimLabel(int BinningOpt, int Axis)
//****************************************
{
  if (Axis < 0 || Axis >= 3)
  {
    std::cerr << "Invalid axis:" << Axis << std::endl;
    throw;
  }

  std::string ReturnVal;
  switch (BinningOpt)
  {
  case 0:
    if (Axis == 0)
      ReturnVal = "ETrue";
    if (Axis == 1)
      ReturnVal = "Erec";
    if (Axis == 2)
      ReturnVal = "Dummy";
    break;
  case 1:
    if (Axis == 0)
      ReturnVal = "ETrue";
    if (Axis == 1)
      ReturnVal = "Momentum";
    if (Axis == 2)
      ReturnVal = "Theta";
    break;
  case 2:
    if (Axis == 0)
      ReturnVal = "ETrue";
    if (Axis == 1)
      ReturnVal = "Erec";
    if (Axis == 2)
      ReturnVal = "Theta";
    break;
  case 3:
    if (Axis == 0)
      ReturnVal = "ETrue";
    if (Axis == 1)
      ReturnVal = "Erec";
    if (Axis == 2)
      ReturnVal = "Q2";
    break;
  case 4:
    if (Axis == 0)
      ReturnVal = "ETrue";
    if (Axis == 1)
      ReturnVal = "Momentum";
    if (Axis == 2)
      ReturnVal = "CosineZenith";
    break;
  default:
    std::cout << "Unrecognised BinningOpt = " << BinningOpt << std::endl;
    throw;
  }

  return ReturnVal;
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
    std::cerr << "Sample name not found : "<<SampleName << std::endl;
    std::cerr << __FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }
  return SampleIndex;
}

//****************************************
void splineFDBase::PrintSampleDetails(std::string SampleName)
//****************************************
{

  int iSample = getSampleIndex(SampleName);

  std::cout << "Details about sample: " << std::setw(20) << SampleNames[iSample] << std::endl;
  std::cout << "\t" << std::setw(35) << "Binning Option"
            << ":" << BinningOpts[iSample] << std::endl;
  std::cout << "\t" << std::setw(35) << "Dimension"
            << ":" << Dimensions[iSample] << std::endl;
  std::cout << "\t" << std::setw(35) << "DetID"
            << ":" << DetIDs[iSample] << std::endl;
  std::cout << "\t" << std::setw(35) << "Number of Spline Params"
            << ":" << nSplineParams[iSample] << std::endl;
  std::cout << "\t" << std::setw(35) << "Number of Oscillation Channels"
            << ":" << nOscChans[iSample] << std::endl;
}

//****************************************
void splineFDBase::PrintArrayDetails(std::string SampleName)
//****************************************
{
  int iSample = getSampleIndex(SampleName);
  int nOscChans = indexvec[iSample].size();
  std::cout << "Sample " << iSample << " has " << nOscChans << " oscillation channels" << std::endl;

  for (int iOscChan = 0; iOscChan < nOscChans; iOscChan++)
  {
    int nSysts = indexvec[iSample][iOscChan].size();
    std::cout << "Oscillation channel " << iOscChan << " has " << nSysts << " systematics" << std::endl;
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
    std::cerr << "Sample index is invalid! 0 <= Index < " << indexvec.size() << std::endl;
    isValid = false;
  }

  if (iOscChan < 0 || iOscChan >= (int)indexvec[iSample].size())
  {
    std::cerr << "OscChan index is invalid! 0 <= Index < " << indexvec[iSample].size() << std::endl;
    isValid = false;
  }

  if (iSyst < 0 || iSyst >= (int)indexvec[iSample][iOscChan].size())
  {
    std::cerr << "Syst index is invalid! 0 <= Index < " << indexvec[iSample][iOscChan].size() << std::endl;
    isValid = false;
  }

  if (iMode < 0 || iMode >= (int)indexvec[iSample][iOscChan][iSyst].size())
  {
    std::cerr << "Mode index is invalid! 0 <= Index < " << indexvec[iSample][iOscChan][iSyst].size() << std::endl;
    isValid = false;
  }

  if (iVar1 < 0 || iVar1 >= (int)indexvec[iSample][iOscChan][iSyst][iMode].size())
  {
    std::cerr << "Var1 index is invalid! 0 <= Index < " << indexvec[iSample][iOscChan][iSyst][iMode].size() << std::endl;
    isValid = false;
  }

  if (iVar2 < 0 || iVar2 >= (int)indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size())
  {
    std::cerr << "Var2 index is invalid! 0 <= Index < " << indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size() << std::endl;
    isValid = false;
  }

  if (iVar3 < 0 || iVar3 >= (int)indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size())
  {
    std::cerr << "Var3 index is invalid! 0 <= Index < " << indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size() << std::endl;
    isValid = false;
  }

  if (!isValid)
  {
    std::cerr << "Given iSample:" << iSample << std::endl;
    std::cerr << "Given iOscChan:" << iOscChan << std::endl;
    std::cerr << "Given iSyst:" << iSyst << std::endl;
    std::cerr << "Given iMode:" << iMode << std::endl;
    std::cerr << "Given iVar1:" << iVar1 << std::endl;
    std::cerr << "Given iVar2:" << iVar2 << std::endl;
    std::cerr << "Given iVar3:" << iVar3 << std::endl;
    std::cerr << "Come visit me at : "<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
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
}
