#include "splineFDBase.h"
#include <memory>

#include "TROOT.h"

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"

//****************************************
splineFDBase::splineFDBase(covarianceXsec *xsec_, MaCh3Modes *Modes_)
              : SplineBase() {
//****************************************
  if (!xsec_) {
    MACH3LOG_ERROR("Trying to create splineFDBase with uninitialised covariance object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  xsec = xsec_;

  if (!Modes_) {
    MACH3LOG_ERROR("Trying to create splineFDBase with uninitialised MaCh3Modes object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  Modes = Modes_;

  // Keep these in class scope, important for using 1 monolith/sample!
  MonolithIndex = 0; //Keeps track of the monolith index we're on when filling arrays (declared here so we can have multiple FillSampleArray calls)
  CoeffIndex = 0; //Keeps track of our indexing the coefficient arrays [x, ybcd]
}
//****************************************
splineFDBase::~splineFDBase(){
//****************************************
  if(manycoeff_arr != nullptr) delete[] manycoeff_arr;
  if(xcoeff_arr != nullptr) delete[] xcoeff_arr;
  if(SplineSegments != nullptr) delete[] SplineSegments;
  if(ParamValues != nullptr) delete[] ParamValues;
}
//****************************************
void splineFDBase::cleanUpMemory() {
//****************************************
  //Call once everything's been allocated in samplePDFSKBase, cleans up junk from memory!
  //Not a huge saving but it's better than leaving everything up to the compiler
  MACH3LOG_INFO("Cleaning up spline memory");
  CleanVector(indexvec);
  CleanVector(SplineFileParPrefixNames);
  CleanVector(GlobalSystIndex);
  CleanVector(SplineModeVecs);
  CleanVector(UniqueSystNames);
  CleanVector(SplineInterpolationTypes);

  for (auto& Binning2D : SplineBinning) {
    for (auto& Binning1D : Binning2D) {
      for (TAxis* axis : Binning1D) {
        delete axis;
      }
      Binning1D.clear();
      Binning1D.shrink_to_fit();
    }
    Binning2D.clear();
    Binning2D.shrink_to_fit();
  }
  SplineBinning.clear();
  SplineBinning.shrink_to_fit();

  //Really make sure all the memory is cleared
  for(auto Spline : splinevec_Monolith){
    if(Spline){delete Spline;}
  }
  splinevec_Monolith.clear();
  splinevec_Monolith.shrink_to_fit();
  if(isflatarray) delete [] isflatarray;
}

//****************************************
bool splineFDBase::AddSample(const std::string& SampleName,
                             const std::string& DetID,
                             const std::vector<std::string>& OscChanFileNames,
                             const std::vector<std::string>& SplineVarNames)
//Adds samples to the large array
//****************************************
{
  SampleNames.push_back(SampleName);
  Dimensions.push_back(int(SplineVarNames.size()));
  DimensionLabels.push_back(SplineVarNames);
  DetIDs.push_back(DetID);

  int nSplineParam = xsec->GetNumParamsFromDetID(DetID, SystType::kSpline);
  nSplineParams.push_back(nSplineParam);

  //This holds the global index of the spline i.e. 0 -> _fNumPar
  std::vector<int> GlobalSystIndex_Sample = xsec->GetGlobalSystIndexFromDetID(DetID, SystType::kSpline);
  //Keep track of this for all the samples
  GlobalSystIndex.push_back(GlobalSystIndex_Sample);

  std::vector<SplineInterpolation> SplineInterpolation_Sample = xsec->GetSplineInterpolationFromDetID(DetID);
  // Keep track of this for all samples
  SplineInterpolationTypes.push_back(SplineInterpolation_Sample);

  std::vector<std::string> SplineFileParPrefixNames_Sample = xsec->GetSplineParsNamesFromDetID(DetID);
  SplineFileParPrefixNames.push_back(SplineFileParPrefixNames_Sample);

  MACH3LOG_INFO("Create SplineModeVecs_Sample");
  std::vector<std::vector<int>> SplineModeVecs_Sample = StripDuplicatedModes(xsec->GetSplineModeVecFromDetID(DetID));
  MACH3LOG_INFO("SplineModeVecs_Sample is of size {}", SplineModeVecs_Sample.size());
  SplineModeVecs.push_back(SplineModeVecs_Sample);

  MACH3LOG_INFO("SplineModeVecs is of size {}", SplineModeVecs.size());

  int nOscChan = int(OscChanFileNames.size());
  nOscChans.push_back(nOscChan);

  PrintSampleDetails(SampleName);

  std::vector<std::vector<TAxis *>> SampleBinning(nOscChan);
  for (int iOscChan = 0; iOscChan < nOscChan; iOscChan++)
  {
    SampleBinning[iOscChan] = FindSplineBinning(OscChanFileNames[iOscChan], SampleName);
  }
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
  SplineBinning.push_back(SampleBinning);

  BuildSampleIndexingArray(SampleName);
  PrintArrayDetails(SampleName);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");

  FillSampleArray(SampleName, OscChanFileNames);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");

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
  
  xcoeff_arr = new M3::float_t[CoeffIndex];
  manycoeff_arr = new M3::float_t[CoeffIndex*4];

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
                for (int iUniqueSyst = 0; iUniqueSyst < nParams; iUniqueSyst++)
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
                  MACH3LOG_ERROR("Couldn't match {} with any of the following {} systs:", SplineFileParPrefixNames[iSample][iSyst], nParams);
                  for (int iUniqueSyst = 0; iUniqueSyst < nParams; iUniqueSyst++)
                  {
                    MACH3LOG_ERROR("{},", UniqueSystNames.at(iUniqueSyst));
                  }//unique syst loop end
                  throw MaCh3Exception(__FILE__ , __LINE__ );
                }

                int splineKnots;
                if(splinevec_Monolith[splineindex]){
                  isflatarray[splineindex]=false;
                  splineKnots=splinevec_Monolith[splineindex]->GetNp();

                  //Now to fill up our coefficient arrayss
                  M3::float_t* tmpXCoeffArr = new M3::float_t[splineKnots];
                  M3::float_t* tmpManyCoeffArr = new M3::float_t[splineKnots*4];

                  int iCoeff=coeffindexvec[splineindex];
                  getSplineCoeff_SepMany(splineindex, tmpXCoeffArr, tmpManyCoeffArr);

                  for(int i = 0; i < splineKnots; i++){
                    xcoeff_arr[iCoeff+i]=tmpXCoeffArr[i];

                    for(int j=0; j<4; j++){
                      manycoeff_arr[(iCoeff+i)*4+j]=tmpManyCoeffArr[i*4+j];
                    }
                  }
                  delete[] tmpXCoeffArr;
                  delete[] tmpManyCoeffArr;
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
}

//****************************************
void splineFDBase::CalcSplineWeights()
//****************************************
{
  #ifdef MULTITHREAD
  #pragma omp parallel for simd
  #endif
  for (size_t iCoeff = 0; iCoeff < uniquecoeffindices.size(); ++iCoeff)
  {
    const int iSpline = uniquecoeffindices[iCoeff];
    const short int uniqueIndex = short(uniquesplinevec_Monolith[iSpline]);
    const short int currentsegment = short(SplineSegments[uniqueIndex]);

    const int segCoeff = coeffindexvec[iSpline]+currentsegment;
    const int coeffOffset = segCoeff * 4;
    // These are what we can extract from the TSpline3
    const M3::float_t y = manycoeff_arr[coeffOffset+kCoeffY];
    const M3::float_t b = manycoeff_arr[coeffOffset+kCoeffB];
    const M3::float_t c = manycoeff_arr[coeffOffset+kCoeffC];
    const M3::float_t d = manycoeff_arr[coeffOffset+kCoeffD];

    // Get the variation for this reconfigure for the ith parameter
    /// @todo KS: Once could use "ParamValues" but this will result in tiny bit different results due to floating point precision
    const M3::float_t xvar = (*SplineInfoArray[uniqueIndex].splineParsPointer);
    // The Delta(x) = xvar - x
    const M3::float_t dx = xvar - xcoeff_arr[segCoeff];

    //Speedy 1% time boost https://en.cppreference.com/w/c/numeric/math/fma (see ND code!)
    M3::float_t weight = M3::fmaf_t(dx, M3::fmaf_t(dx, M3::fmaf_t(dx, d, c), b), y);
    //This is the speedy version of writing dx^3+b*dx^2+c*dx+d

    //ETA - do we need this? We check later for negative weights and I wonder if this is even
    //possible with the fmaf line above?
    if(weight < 0){weight = 0.;}  //Stops is getting negative weights

    weightvec_Monolith[iSpline] = weight;
  }
}

//****************************************
void splineFDBase::BuildSampleIndexingArray(const std::string& SampleName)
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
      int nModesInSyst = int(SplineModeVecs[iSample][iSyst].size());
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

  //Try declaring these outside of TFile so they aren't owned by File
  int nDummyBins = 1;
  double DummyEdges[2];
  DummyEdges[0] = -1e15;
  DummyEdges[1] = 1e15;
  TAxis* DummyAxis = new TAxis(nDummyBins, DummyEdges);
  TH2F* Hist2D = nullptr;
  TH3F* Hist3D = nullptr;

  auto File = std::unique_ptr<TFile>(TFile::Open(FileName.c_str(), "READ"));
  if (!File || File->IsZombie())
  {
    MACH3LOG_ERROR("File {} not found", FileName);
    MACH3LOG_ERROR("This is caused by something here! {} : {}", __FILE__, __LINE__);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  MACH3LOG_INFO("Finding binning for:");
  MACH3LOG_INFO("{}", FileName);

  bool isHist2D = false;
  bool isHist3D = false;

  std::string TemplateName = "dev_tmp_0_0";
  TObject *Obj = File->Get(TemplateName.c_str());
  //If you can't find dev_tmp_0_0 then this will cause a problem
  if (!Obj)
  {
    TemplateName = "dev_tmp.0.0";
    Obj = File->Get(TemplateName.c_str());
    if (!Obj)
    {
      MACH3LOG_ERROR("Error: could not find dev_tmp_0_0 in spline file. Spline binning cannot be set!");
      MACH3LOG_ERROR("FileName: {}", FileName);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  //Now check if dev_tmp_0_0 is a TH2 i.e. specifying the dimensions of the splines is 2D
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
    MACH3LOG_ERROR("Object doesn't inherit from either TH2D and TH3D - Odd A");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (isHist2D)
  {
    if (Dimensions[iSample] != 2)
    {
      MACH3LOG_ERROR("Trying to load a 2D spline template when nDim={}", Dimensions[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    //Hist2D = std::unique_ptr<TH2F>(File->Get<TH2F>("dev_tmp_0_0"));
    Hist2D = File->Get<TH2F>(TemplateName.c_str());
  }

  if (isHist3D)
  {
    Hist3D = File->Get<TH3F>((TemplateName.c_str()));

    if (Dimensions[iSample] != 3 && Hist3D->GetZaxis()->GetNbins() != 1)
    {
      MACH3LOG_ERROR("Trying to load a 3D spline template when nDim={}", Dimensions[iSample]);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    Hist3D = File->Get<TH3F>(TemplateName.c_str());
  }

  if (Dimensions[iSample] == 2) {
    if(isHist2D){
      ReturnVec.push_back(static_cast<TAxis*>(Hist2D->GetXaxis()->Clone()));
      ReturnVec.push_back(static_cast<TAxis*>(Hist2D->GetYaxis()->Clone()));
      ReturnVec.push_back(static_cast<TAxis*>(DummyAxis->Clone()));
    } else if (isHist3D) {
      ReturnVec.push_back(static_cast<TAxis*>(Hist3D->GetXaxis()->Clone()));
      ReturnVec.push_back(static_cast<TAxis*>(Hist3D->GetYaxis()->Clone()));
      ReturnVec.push_back(static_cast<TAxis*>(DummyAxis->Clone()));
    }
  } else if (Dimensions[iSample] == 3) {
    ReturnVec.push_back(static_cast<TAxis*>(Hist3D->GetXaxis()->Clone()));
    ReturnVec.push_back(static_cast<TAxis*>(Hist3D->GetYaxis()->Clone()));
    ReturnVec.push_back(static_cast<TAxis*>(Hist3D->GetZaxis()->Clone()));
  } else {
    MACH3LOG_ERROR("Number of dimensions not valid! Given: {}", Dimensions[iSample]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  for (unsigned int iAxis = 0; iAxis < ReturnVec.size(); ++iAxis) {
    PrintBinning(ReturnVec[iAxis]);
  }

  MACH3LOG_INFO("Left PrintBinning now tidying up");
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
                  if (splinevec_Monolith[splineindex])
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
    MACH3LOG_DEBUG("{:<10} has {:<10} splines, of which {:<10} are not flat", SampleNames[iSample], SampleCounter_All, SampleCounter_NonFlat);

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
void splineFDBase::PrepForReweight() {
//****************************************
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
                  if (splinevec_Monolith[splineindex])
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
  
  nParams = static_cast<short int>(UniqueSystSplines.size());

  // DB Find the number of splines knots which assumes each instance of the syst has the same number of knots
  SplineSegments = new short int[nParams]();
  ParamValues = new float[nParams]();
  SplineInfoArray.resize(nParams);
  for (int iSpline = 0; iSpline < nParams; iSpline++)
  {
    SplineInfoArray[iSpline].nPts = static_cast<M3::int_t>(UniqueSystSplines[iSpline]->GetNp());
    SplineInfoArray[iSpline].xPts.resize(SplineInfoArray[iSpline].nPts);
    SplineInfoArray[iSpline].splineParsPointer = xsec->retPointer(UniqueSystIndices[iSpline]);
    for (int iKnot = 0; iKnot < SplineInfoArray[iSpline].nPts; iKnot++)
    {
      M3::float_t xPoint;
      M3::float_t yPoint;
      UniqueSystSplines[iSpline]->GetKnot(iKnot, xPoint, yPoint);
      SplineInfoArray[iSpline].xPts[iKnot] = xPoint;
    }
    //ETA - let this just be set as the first segment by default
    SplineSegments[iSpline] = 0;
    ParamValues[iSpline] = 0.;
  }
  
  MACH3LOG_INFO("nUniqueSysts: {}", nParams);
  MACH3LOG_INFO("{:<15} | {:<20} | {:<6}", "Spline Index", "Syst Name", "nKnots");
  for (int iUniqueSyst = 0; iUniqueSyst < nParams; iUniqueSyst++)
  {
    MACH3LOG_INFO("{:<15} | {:<20} | {:<6}", iUniqueSyst, UniqueSystNames[iUniqueSyst], SplineInfoArray[iUniqueSyst].nPts);
  }

  //ETA
  //Isn't this just doing what CountNumberOfLoadedSplines() does?
  int nCombinations_FlatSplines = 0;
  int nCombinations_All = 0;
  // DB Now actually loop over splines to determine which are all null i.e. flat
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
				if (splinevec_Monolith[splineindex])
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
void splineFDBase::getSplineCoeff_SepMany(int splineindex, M3::float_t* &xArray, M3::float_t* &manyArray){
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

  for(int i=0; i<nPoints; i++) {
    // Spline coefficients to be
    // M3::float_t type is defined by the LOW_MEMORY_STRUCTS compiler flag
    // so M3::float_t can be double or float depending on this
    M3::float_t x = M3::float_t(-999.99);
    M3::float_t y = M3::float_t(-999.99);
    M3::float_t b = M3::float_t(-999.99);
    M3::float_t c = M3::float_t(-999.99);
    M3::float_t d = M3::float_t(-999.99);
    splinevec_Monolith[splineindex]->GetCoeff(i, x, y, b, c, d);

    // Store the coefficients for each knot contiguously in memory
    // 4 because manyArray stores y,b,c,d
    xArray[i] = x;
    manyArray[i*4] = y; 
    manyArray[i*4+1] = b;
    manyArray[i*4+2] = c;
    manyArray[i*4+3] = d;    
  }

  //We now clean up the splines!
  delete splinevec_Monolith[splineindex];
  splinevec_Monolith[splineindex] = nullptr;
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
int splineFDBase::getSampleIndex(const std::string& SampleName){
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
void splineFDBase::PrintSampleDetails(const std::string& SampleName)
//****************************************
{
  int iSample = getSampleIndex(SampleName);

  MACH3LOG_INFO("Details about sample: {:<20}", SampleNames[iSample]);
  MACH3LOG_INFO("\t Dimension: {:<35}", Dimensions[iSample]);
  MACH3LOG_INFO("\t DetID: {:<35}", DetIDs[iSample]);
  MACH3LOG_INFO("\t nSplineParam: {:<35}", nSplineParams[iSample]);
  MACH3LOG_INFO("\t nOscChan: {:<35}", nOscChans[iSample]);
}

//****************************************
void splineFDBase::PrintArrayDetails(const std::string& SampleName)
//****************************************
{
  int iSample = getSampleIndex(SampleName);
  int nOscChannels = int(indexvec[iSample].size());
  MACH3LOG_INFO("Sample {} has {} oscillation channels", SampleName, nOscChannels);	
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++)
  {
    int nSysts = int(indexvec[iSample][iOscChan].size());
    MACH3LOG_INFO("Oscillation channel {} has {} systematics", iOscChan, nSysts);	  
  }

  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
  MACH3LOG_INFO("Printing no. of modes affected by each systematic for each oscillation channel");
  for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++) {
    std::string modes = fmt::format("OscChan: {}\t", iOscChan);
    for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++) {
      modes += fmt::format("{} ", indexvec[iSample][iOscChan][iSyst].size());
    }
    MACH3LOG_INFO("{}", modes);
  }
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
}

//****************************************
bool splineFDBase::isValidSplineIndex(const std::string& SampleName, int iOscChan, int iSyst, int iMode, int iVar1, int iVar2, int iVar3)
//****************************************
{
  int iSample=getSampleIndex(SampleName);
  bool isValid = true;

  if (iSample < 0 || iSample >= int(indexvec.size()))
  {
    MACH3LOG_ERROR("Sample index is invalid! 0 <= Index < {} ", indexvec.size());
    isValid = false;
  }

  if (iOscChan < 0 || iOscChan >= int(indexvec[iSample].size()))
  {
    MACH3LOG_ERROR("OscChan index is invalid! 0 <= Index < {} ", indexvec[iSample].size());
    isValid = false;
  }

  if (iSyst < 0 || iSyst >= int(indexvec[iSample][iOscChan].size()))
  {
    MACH3LOG_ERROR("Syst index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan].size());
    isValid = false;
  }

  if (iMode < 0 || iMode >= int(indexvec[iSample][iOscChan][iSyst].size()))
  {
    MACH3LOG_ERROR("Mode index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst].size());
    isValid = false;
  }

  if (iVar1 < 0 || iVar1 >= int(indexvec[iSample][iOscChan][iSyst][iMode].size()))
  {
    MACH3LOG_ERROR("Var1 index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst][iMode].size());	  
    isValid = false;
  }

  if (iVar2 < 0 || iVar2 >= int(indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size()))
  {
    MACH3LOG_ERROR("Var2 index is invalid! 0 <= Index < {} ", indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size());
    isValid = false;
  }

  if (iVar3 < 0 || iVar3 >= int(indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size()))
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
  std::string text = "";
  for (int iBin = 0; iBin <= NBins; iBin++) {
    text += fmt::format("{} ", Axis->GetXbins()->GetAt(iBin));
  }
  MACH3LOG_INFO("{}", text);
}


// checks if there are multiple modes with the same SplineSuffix
// (for example if CCRES and CCCoherent are treated as one spline mode)
std::vector< std::vector<int> > splineFDBase::StripDuplicatedModes(std::vector< std::vector<int> > InputVector) {

  //ETA - this is of size nPars from the xsec model                                                                                                                                                                   
  size_t InputVectorSize = InputVector.size();
  std::vector< std::vector<int> > ReturnVec(InputVectorSize);

  //ETA - loop over all systematics                                                                                                                                                                                   
  for (size_t iSyst=0;iSyst<InputVectorSize;iSyst++) {
    std::vector<int> TmpVec;
    std::vector<std::string> TestVec;

    //Loop over the modes that we've listed in xsec cov                                                                                                                                                               
    for (unsigned int iMode = 0 ; iMode < InputVector[iSyst].size() ; iMode++) {
      int Mode = InputVector[iSyst][iMode];
      std::string ModeName = Modes->GetSplineSuffixFromMaCh3Mode(Mode);

      bool IncludeMode = true;
      for (auto TestString : TestVec) {
        if (ModeName == TestString) {
          IncludeMode = false;
          break;
	}
      }

      if (IncludeMode) {
        TmpVec.push_back(Mode);
        TestVec.push_back(ModeName);
      }
    }

    ReturnVec[iSyst] = TmpVec;
  }

  return ReturnVec;
}

//****************************************
std::vector< std::vector<int> > splineFDBase::GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val)
//****************************************
{
  std::vector<std::vector<int>> ReturnVec;
  int SampleIndex = -1;
  for (unsigned int iSample = 0; iSample < SampleNames.size(); iSample++) {
    if (SampleName == SampleNames[iSample]) {
      SampleIndex = iSample;
    }
  }

  if (SampleIndex == -1)
  {
          MACH3LOG_ERROR("Sample not found: {}", SampleName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  size_t nSplineSysts = indexvec[SampleIndex][iOscChan].size();
  int Mode = EventMode;

  int Var1Bin = SplineBinning[SampleIndex][iOscChan][0]->FindBin(Var1Val)-1;
  if (Var1Bin < 0 || Var1Bin >= SplineBinning[SampleIndex][iOscChan][0]->GetNbins()){
    return ReturnVec;
  }

  int Var2Bin = SplineBinning[SampleIndex][iOscChan][1]->FindBin(Var2Val)-1;
  if (Var2Bin < 0 || Var2Bin >= SplineBinning[SampleIndex][iOscChan][1]->GetNbins()){
    return ReturnVec;
  }

  int Var3Bin = SplineBinning[SampleIndex][iOscChan][2]->FindBin(Var3Val)-1;

  if (Var3Bin < 0 || Var3Bin >= SplineBinning[SampleIndex][iOscChan][2]->GetNbins()){
    return ReturnVec;
  }

  for(size_t iSyst=0; iSyst<nSplineSysts; iSyst++){
    std::vector<int> spline_modes = SplineModeVecs[SampleIndex][iSyst];
    size_t nSampleModes = spline_modes.size();

    //ETA - look here at the length of spline_modes and what you're actually comparing against                                                                                                                        
    for(size_t iMode = 0; iMode<nSampleModes ; iMode++){
      //Only consider if the event mode (Mode) matches ones of the spline modes                                                                                                                                       
      if (Mode == spline_modes[iMode]) {
        std::vector<int> event_vec(7);
        event_vec[0]=SampleIndex;
        event_vec[1]=iOscChan;
        event_vec[2]=static_cast<int>(iSyst);
        event_vec[3]=static_cast<int>(iMode);
        event_vec[4]=Var1Bin;
        event_vec[5]=Var2Bin;
        event_vec[6]=Var3Bin;
        int splineID=indexvec[SampleIndex][iOscChan][iSyst][iMode][Var1Bin][Var2Bin][Var3Bin];
        //Also check that the spline isn't flat                                                                                                                                                                       
        if(!isflatarray[splineID]){
          ReturnVec.push_back(event_vec);
        }
      }
    }
  }
  return ReturnVec;
}

void splineFDBase::FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames)
{
  int iSample = getSampleIndex(SampleName);
  
  int nOscChannels = nOscChans[iSample];
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
    MACH3LOG_INFO("Processing: {}", OscChanFileNames[iOscChan]);
    
    TSpline3* mySpline = nullptr;
    TSpline3_red* Spline = nullptr;
    TString Syst, Mode;
    int nKnots, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin = M3::_BAD_INT_;
    double x,y, Eval = M3::_BAD_DOUBLE_;
    bool isFlat = true;

    auto File = std::unique_ptr<TFile>(TFile::Open(OscChanFileNames[iOscChan].c_str()));

    if (!File || File->IsZombie()) {
      MACH3LOG_ERROR("File {} not found", OscChanFileNames[iOscChan]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (auto k : *File->GetListOfKeys()) {
      auto Key = static_cast<TKey*>(k);
      TClass *Class = gROOT->GetClass(Key->GetClassName(), false);
      if(!Class->InheritsFrom("TSpline3")) {
        continue;
      }

      TString FullSplineName = TString(Key->GetName());
      std::vector<std::string> Tokens = GetTokensFromSplineName(std::string(FullSplineName));

      if (Tokens.size() != kNTokens) {
	std::cerr << "Invalid tokens from spline name - Expected " << kNTokens << " tokens. Check implementation in GetTokensFromSplineName()" << std::endl;
	throw;
      }

      Syst = Tokens[kSystToken];
      Mode = Tokens[kModeToken];
      Var1Bin = std::stoi(Tokens[kVar1BinToken]);
      Var2Bin = std::stoi(Tokens[kVar2BinToken]);
      Var3Bin = std::stoi(Tokens[kVar3BinToken]);

      SystNum = -1;
      for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++) {
        if (strcmp(Syst, SplineFileParPrefixNames[iSample][iSyst].c_str()) == 0) {
          SystNum = iSyst;
          break;
        }
      }

      // If the syst doesn't match any of the spline names then skip it
      if (SystNum == -1){
        MACH3LOG_DEBUG("Couldn't Match any systematic name in xsec yaml with spline name: {}" , FullSplineName.Data());
        continue;
      }

      ModeNum = -1;
      for (unsigned int iMode = 0; iMode < SplineModeVecs[iSample][SystNum].size(); iMode++) {
        if (strcmp(Mode, Modes->GetSplineSuffixFromMaCh3Mode(SplineModeVecs[iSample][SystNum][iMode]).c_str()) == 0) {
          ModeNum = iMode;
          break;
        }
      }

      if (ModeNum == -1) {
        MACH3LOG_ERROR("Couldn't find mode for {} in {}. Problem Spline is : {} ", Mode, Syst, FullSplineName);
        throw;
      }

      mySpline = Key->ReadObject<TSpline3>();

      if (isValidSplineIndex(SampleName, iOscChan, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin)) {
	// loop over all the spline knots and check their value
        // if the value is 1 then set the flat bool to false
        nKnots = mySpline->GetNp();
        isFlat = true;
	for (int iKnot = 0; iKnot < nKnots; iKnot++)
          {
            mySpline->GetKnot(iKnot, x, y);

            Eval = mySpline->Eval(x);
            if (Eval < 0.99999 || Eval > 1.00001) {
              isFlat = false;
              break;
            }
          }

	//Rather than keeping a mega vector of splines then converting, this should just keep everything nice in memory
        indexvec[iSample][iOscChan][SystNum][ModeNum][Var1Bin][Var2Bin][Var3Bin]=MonolithIndex;
        coeffindexvec.push_back(CoeffIndex);
        // Should save memory rather saving [x_i_0 ,... x_i_maxknots] for every spline!
        if (isFlat) {
          splinevec_Monolith.push_back(nullptr);
          delete mySpline;
        } else {
          Spline = new TSpline3_red(mySpline, SplineInterpolationTypes[iSample][SystNum]);
          delete mySpline;
          splinevec_Monolith.push_back(Spline);
          uniquecoeffindices.push_back(MonolithIndex); //So we can get the unique coefficients and skip flat splines later on
          CoeffIndex+=nKnots;
        }
      	//Incrementing MonolithIndex to keep track of number of valid spline indices
        MonolithIndex+=1;
      } else {
        //Potentially you are not a valid spline index
	delete mySpline;
      }
    }//End of loop over all TKeys in file
    //A bit of clean up                                                                                                                                                                                               
    File->Delete("*");
    File->Close();
  } //End of oscillation channel loop
}
