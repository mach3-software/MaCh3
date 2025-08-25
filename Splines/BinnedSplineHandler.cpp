#include "BinnedSplineHandler.h"
#include <memory>

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"

_MaCh3_Safe_Include_Start_ //{
#include "TROOT.h"
#include "TKey.h"
_MaCh3_Safe_Include_End_ //}

//****************************************
BinnedSplineHandler::BinnedSplineHandler(ParameterHandlerGeneric *xsec_, MaCh3Modes *Modes_) : SplineBase() {
//****************************************
  if (!xsec_) {
    MACH3LOG_ERROR("Trying to create BinnedSplineHandler with uninitialized covariance object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  xsec = xsec_;

  if (!Modes_) {
    MACH3LOG_ERROR("Trying to create BinnedSplineHandler with uninitialized MaCh3Modes object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  Modes = Modes_;

  // Keep these in class scope, important for using 1 monolith/sample!
  MonolithIndex = 0; //Keeps track of the monolith index we're on when filling arrays (declared here so we can have multiple FillSampleArray calls)
  CoeffIndex = 0; //Keeps track of our indexing the coefficient arrays [x, ybcd]
  isflatarray = nullptr;
}

//****************************************
BinnedSplineHandler::~BinnedSplineHandler(){
//****************************************
  if(manycoeff_arr != nullptr) delete[] manycoeff_arr;
  if(xcoeff_arr != nullptr) delete[] xcoeff_arr;
  if(SplineSegments != nullptr) delete[] SplineSegments;
  if(ParamValues != nullptr) delete[] ParamValues;
}
//****************************************
void BinnedSplineHandler::cleanUpMemory() {
//****************************************
  //Call once everything's been allocated in SampleHandlerFDBase, cleans up junk from memory!
  //Not a huge saving but it's better than leaving everything up to the compiler
  MACH3LOG_INFO("Cleaning up spline memory");
  CleanVector(indexvec);
  CleanVector(SplineFileParPrefixNames);
  CleanVector(GlobalSystIndex);
  CleanVector(SplineModeVecs);
  CleanVector(UniqueSystNames);
  CleanVector(SplineInterpolationTypes);
  CleanVector(nOscChans);
  CleanVector(nSplineParams);
  CleanVector(DimensionLabels);
  CleanVector(SampleNames);
  CleanVector(SampleTittles);
  CleanVector(Dimensions);
  CleanContainer(splinevec_Monolith);
  CleanContainer(SplineBinning);
  CleanVector(UniqueSystIndices);
  if(isflatarray) delete [] isflatarray;
}

//****************************************
void BinnedSplineHandler::AddSample(const std::string& SampleName,
                                    const std::string& SampleTittle,
                                    const std::vector<std::string>& OscChanFileNames,
                                    const std::vector<std::string>& SplineVarNames)
//Adds samples to the large array
//****************************************
{
  SampleNames.push_back(SampleName);
  SampleTittles.push_back(SampleTittle);
  Dimensions.push_back(static_cast<int>(SplineVarNames.size()));
  DimensionLabels.push_back(SplineVarNames);

  int nSplineParam = xsec->GetNumParamsFromSampleName(SampleName, SystType::kSpline);
  nSplineParams.push_back(nSplineParam);

  //This holds the global index of the spline i.e. 0 -> _fNumPar
  std::vector<int> GlobalSystIndex_Sample = xsec->GetGlobalSystIndexFromSampleName(SampleName, SystType::kSpline);
  //Keep track of this for all the samples
  GlobalSystIndex.push_back(GlobalSystIndex_Sample);

  std::vector<SplineInterpolation> SplineInterpolation_Sample = xsec->GetSplineInterpolationFromSampleName(SampleName);
  // Keep track of this for all samples
  SplineInterpolationTypes.push_back(SplineInterpolation_Sample);

  std::vector<std::string> SplineFileParPrefixNames_Sample = xsec->GetSplineParsNamesFromSampleName(SampleName);
  SplineFileParPrefixNames.push_back(SplineFileParPrefixNames_Sample);

  MACH3LOG_INFO("Create SplineModeVecs_Sample");
  std::vector<std::vector<int>> SplineModeVecs_Sample = StripDuplicatedModes(xsec->GetSplineModeVecFromSampleName(SampleName));
  MACH3LOG_INFO("SplineModeVecs_Sample is of size {}", SplineModeVecs_Sample.size());
  SplineModeVecs.push_back(SplineModeVecs_Sample);

  MACH3LOG_INFO("SplineModeVecs is of size {}", SplineModeVecs.size());

  int nOscChan = int(OscChanFileNames.size());
  nOscChans.push_back(nOscChan);

  PrintSampleDetails(SampleTittle);

  std::vector<std::vector<TAxis *>> SampleBinning(nOscChan);
  for (int iOscChan = 0; iOscChan < nOscChan; iOscChan++)
  {
    SampleBinning[iOscChan] = FindSplineBinning(OscChanFileNames[iOscChan], SampleTittle);
  }
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
  SplineBinning.push_back(SampleBinning);

  BuildSampleIndexingArray(SampleTittle);
  PrintArrayDetails(SampleTittle);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");

  FillSampleArray(SampleTittle, OscChanFileNames);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
}

//****************************************
void BinnedSplineHandler::InvestigateMissingSplines() const {
//****************************************
  // Map: iSample → iSyst → modeSuffix → {totalSplines, zeroCount}
  std::map<unsigned int, std::map<unsigned int, std::map<std::string, std::pair<unsigned int, unsigned int>>>> systZeroCounts;

  for (unsigned int iSample = 0; iSample < indexvec.size(); iSample++) {
    std::string SampleName = SampleNames[iSample];

    // Get list of systematic names for this sample
    std::vector<std::string> SplineFileParPrefixNames_Sample =
    xsec->GetParsNamesFromSampleName(SampleName, kSpline);

    for (unsigned int iOscChan = 0; iOscChan < indexvec[iSample].size(); iOscChan++) {
      for (unsigned int iSyst = 0; iSyst < indexvec[iSample][iOscChan].size(); iSyst++) {
        unsigned int zeroCount = 0;
        unsigned int totalSplines = 0;

        for (unsigned int iMode = 0; iMode < indexvec[iSample][iOscChan][iSyst].size(); iMode++) {
          // Get the mode suffix string
          std::string modeSuffix =
          Modes->GetSplineSuffixFromMaCh3Mode(SplineModeVecs[iSample][iSyst][iMode]);

          for (unsigned int iVar1 = 0; iVar1 < indexvec[iSample][iOscChan][iSyst][iMode].size(); iVar1++) {
            for (unsigned int iVar2 = 0; iVar2 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(); iVar2++) {
              for (unsigned int iVar3 = 0; iVar3 < indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(); iVar3++) {
                totalSplines++;
                if (indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3] == 0) {
                  zeroCount++;
                  if(zeroCount > 1){
                    systZeroCounts[iSample][iSyst][modeSuffix] = {totalSplines, zeroCount};
                  }
                  MACH3LOG_DEBUG(
                    "Sample '{}' | OscChan {} | Syst '{}' | Mode '{}' | Var1 {} | Var2 {} | Var3 {} => Value: {}",
                    SampleTittles[iSample],
                    iOscChan,
                    SplineFileParPrefixNames_Sample[iSyst],
                    modeSuffix,
                    iVar1,
                    iVar2,
                    iVar3,
                    indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2][iVar3]
                  );
                }
              }
            }
          }
        }
      }
    }
  }

  // KS: Let's print this atrocious mess...
  for (const auto& samplePair : systZeroCounts) {
    unsigned int iSample = samplePair.first;
    std::vector<std::string> SplineFileParPrefixNames_Sample = xsec->GetParsNamesFromSampleName(SampleNames[iSample], kSpline);
    for (const auto& systPair : samplePair.second) {
      unsigned int iSyst = systPair.first;
      const auto& systName = SplineFileParPrefixNames_Sample[iSyst];
      for (const auto& modePair : systPair.second) {
        const auto& modeSuffix = modePair.first;
        const auto& counts = modePair.second;
        MACH3LOG_CRITICAL(
          "Sample '{}': Systematic '{}' has missing splines in mode '{}'. Expected Splines: {}, Missing Splines: {}",
          SampleTittles[iSample],
          systName,
          modeSuffix,
          counts.first,
          counts.second
        );
      }
    }
  }
}

//****************************************
void BinnedSplineHandler::TransferToMonolith()
//****************************************
{
  PrepForReweight(); 
  MonolithSize = CountNumberOfLoadedSplines(false, 1);

  if(MonolithSize!=MonolithIndex){
    InvestigateMissingSplines();
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
  manycoeff_arr = new M3::float_t[CoeffIndex*_nCoeff_];

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
                  M3::float_t* tmpManyCoeffArr = new M3::float_t[splineKnots*_nCoeff_];

                  int iCoeff=coeffindexvec[splineindex];
                  getSplineCoeff_SepMany(splineindex, tmpXCoeffArr, tmpManyCoeffArr);

                  for(int i = 0; i < splineKnots; i++){
                    xcoeff_arr[iCoeff+i]=tmpXCoeffArr[i];

                    for(int j = 0; j < _nCoeff_; j++){
                      manycoeff_arr[(iCoeff+i)*_nCoeff_+j]=tmpManyCoeffArr[i*_nCoeff_+j];
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
void BinnedSplineHandler::Evaluate() {
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
void BinnedSplineHandler::CalcSplineWeights()
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
    const int coeffOffset = segCoeff * _nCoeff_;
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
//Creates an array to be filled with monolith indexes for each sample (allows for indexing between 7D binning and 1D Vector)
//Only need 1 indexing array everything else interfaces with this to get binning properties
void BinnedSplineHandler::BuildSampleIndexingArray(const std::string& SampleTittle)
//****************************************
{  
  int iSample = getSampleIndex(SampleTittle);
  int nSplineSysts = nSplineParams[iSample];
  int nOscChannels = nOscChans[iSample];

  // Resize the main indexing structure
  indexvec.emplace_back(nOscChannels);

  for (int iOscChan = 0; iOscChan < nOscChannels; ++iOscChan)
  {
    indexvec.back()[iOscChan].resize(nSplineSysts);
    for (int iSyst = 0; iSyst < nSplineSysts; ++iSyst)
    {
      int nModesInSyst = static_cast<int>(SplineModeVecs[iSample][iSyst].size());
      indexvec.back()[iOscChan][iSyst].resize(nModesInSyst);

      for (int iMode = 0; iMode < nModesInSyst; ++iMode)
      {
        const int nBins1 = SplineBinning[iSample][iOscChan][0]->GetNbins();
        const int nBins2 = SplineBinning[iSample][iOscChan][1]->GetNbins();
        const int nBins3 = SplineBinning[iSample][iOscChan][2]->GetNbins();

        indexvec.back()[iOscChan][iSyst][iMode]
                .resize(nBins1,std::vector<std::vector<int>>(nBins2, std::vector<int>(nBins3, 0)));
      }
    } // end of iSyst loop
  } // end of iOscChan loop
}

//****************************************
std::vector<TAxis *> BinnedSplineHandler::FindSplineBinning(const std::string& FileName, const std::string& SampleTittle)
//****************************************
{
  int iSample=getSampleIndex(SampleTittle);

  //Try declaring these outside of TFile so they aren't owned by File
  constexpr int nDummyBins = 1;
  constexpr double DummyEdges[2] = {-1e15, 1e15};
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

  std::string TemplateName = "dev_tmp_0_0";
  TObject *Obj = File->Get(TemplateName.c_str());
  //If you can't find dev_tmp_0_0 then this will cause a problem
  if (!Obj)
  {
    TemplateName = "dev_tmp.0.0";
    Obj = File->Get(TemplateName.c_str());
    if (!Obj)
    {
      MACH3LOG_ERROR("Could not find dev_tmp_0_0 in spline file. Spline binning cannot be set!");
      MACH3LOG_ERROR("FileName: {}", FileName);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  //Now check if dev_tmp_0_0 is a TH2 i.e. specifying the dimensions of the splines is 2D
  bool isHist2D = Obj->IsA() == TH2F::Class();
  //For T2K annoyingly all objects are TH3Fs
  bool isHist3D = Obj->IsA() == TH3F::Class();
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

  std::vector<TAxis*> ReturnVec;
  // KS: Resize to reduce impact of push back and memory fragmentation
  ReturnVec.resize(3);
  if (Dimensions[iSample] == 2) {
    if (isHist2D) {
      ReturnVec[0] = static_cast<TAxis*>(Hist2D->GetXaxis()->Clone());
      ReturnVec[1] = static_cast<TAxis*>(Hist2D->GetYaxis()->Clone());
      ReturnVec[2] = static_cast<TAxis*>(DummyAxis->Clone());
    } else if (isHist3D) {
      ReturnVec[0] = static_cast<TAxis*>(Hist3D->GetXaxis()->Clone());
      ReturnVec[1] = static_cast<TAxis*>(Hist3D->GetYaxis()->Clone());
      ReturnVec[2] = static_cast<TAxis*>(DummyAxis->Clone());
    }
  } else if (Dimensions[iSample] == 3) {
    ReturnVec[0] = static_cast<TAxis*>(Hist3D->GetXaxis()->Clone());
    ReturnVec[1] = static_cast<TAxis*>(Hist3D->GetYaxis()->Clone());
    ReturnVec[2] = static_cast<TAxis*>(Hist3D->GetZaxis()->Clone());
  } else {
    MACH3LOG_ERROR("Number of dimensions not valid! Given: {}", Dimensions[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (unsigned int iAxis = 0; iAxis < ReturnVec.size(); ++iAxis) {
    PrintBinning(ReturnVec[iAxis]);
  }

  MACH3LOG_INFO("Left PrintBinning now tidying up");
  delete DummyAxis;

  return ReturnVec;
}

//****************************************
int BinnedSplineHandler::CountNumberOfLoadedSplines(bool NonFlat, int Verbosity)
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
    std::string SampleTittle = SampleTittles[iSample];
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
                if (isValidSplineIndex(SampleTittle, iOscChan, iSyst, iMode, iVar1, iVar2, iVar3))
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
    MACH3LOG_DEBUG("{:<10} has {:<10} splines, of which {:<10} are not flat", SampleTittles[iSample], SampleCounter_All, SampleCounter_NonFlat);

    FullCounter_NonFlat += SampleCounter_NonFlat;
    FullCounter_All += SampleCounter_All;
  }

  if (Verbosity > 0)
  {
    MACH3LOG_INFO("Total number of splines loaded: {}", FullCounter_All);
    MACH3LOG_INFO("Total number of non-flat splines loaded: {}", FullCounter_NonFlat);
  }

  if (NonFlat) {
    return FullCounter_NonFlat;
  } else {
    return FullCounter_All;
  }
}

//****************************************
void BinnedSplineHandler::PrepForReweight() {
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
                  if (FoundNonFlatSpline) { break;}
                }//3D loop end
                if (FoundNonFlatSpline) { break;}
              }//2D loop end
              if (FoundNonFlatSpline){ break; }
            }//1D loop end
            if (FoundNonFlatSpline){ break; }
          }//mode loop end
          if (FoundNonFlatSpline) { break; }
        }//osc loop end
        //ETA - only push back unique name if a non-flat response has been found
        if(FoundNonFlatSpline){
          UniqueSystNames.push_back(SystName);
        }

        if (!FoundNonFlatSpline)
        {
          MACH3LOG_INFO("{} syst has no response in sample {}", SystName, iSample);
          MACH3LOG_INFO("Whilst this isn't necessarily a problem, it seems odd");
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
    SplineInfoArray[iSpline].splineParsPointer = xsec->RetPointer(UniqueSystIndices[iSpline]);
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
void BinnedSplineHandler::getSplineCoeff_SepMany(int splineindex, M3::float_t* &xArray, M3::float_t* &manyArray){
//****************************************
  //No point evaluating a flat spline
  int nPoints = splinevec_Monolith[splineindex]->GetNp();

  for (int i = 0; i < nPoints; i++) {
    xArray[i] = 1.0;
    for (int j = 0; j < _nCoeff_; j++) {
      manyArray[i*_nCoeff_+j] = 1.0;
    }
  }

  for(int i=0; i<nPoints; i++) {
    M3::float_t x = M3::float_t(-999.99);
    M3::float_t y = M3::float_t(-999.99);
    M3::float_t b = M3::float_t(-999.99);
    M3::float_t c = M3::float_t(-999.99);
    M3::float_t d = M3::float_t(-999.99);
    splinevec_Monolith[splineindex]->GetCoeff(i, x, y, b, c, d);

    // Store the coefficients for each knot contiguously in memory
    // 4 because manyArray stores y,b,c,d
    xArray[i] = x;
    manyArray[i * _nCoeff_ + kCoeffY] = y;
    manyArray[i * _nCoeff_ + kCoeffB] = b;
    manyArray[i * _nCoeff_ + kCoeffC] = c;
    manyArray[i * _nCoeff_ + kCoeffD] = d;
  }

  //We now clean up the splines!
  delete splinevec_Monolith[splineindex];
  splinevec_Monolith[splineindex] = nullptr;
}

//****************************************
//Equally though could just use KinematicVariable to map back
std::string BinnedSplineHandler::getDimLabel(const int iSample, const unsigned int Axis) const
//****************************************
{
  if(Axis > DimensionLabels[iSample].size()){
    MACH3LOG_ERROR("The spline Axis you are trying to get the label of is larger than the number of dimensions");
    MACH3LOG_ERROR("You are trying to get axis {} but have only got {}", Axis, Dimensions[iSample]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return DimensionLabels.at(iSample).at(Axis);
}

//****************************************
//Returns sample index in
int BinnedSplineHandler::getSampleIndex(const std::string& SampleTittle) const{
//****************************************
  for (size_t iSample = 0; iSample < SampleTittles.size(); ++iSample) {
    if (SampleTittle == SampleTittles[iSample]) {
      return static_cast<int>(iSample);
    }
  }
  MACH3LOG_ERROR("Sample name not found: {}", SampleTittle);
  throw MaCh3Exception(__FILE__, __LINE__);
}

//****************************************
void BinnedSplineHandler::PrintSampleDetails(const std::string& SampleTittle) const
//****************************************
{
  const int iSample = getSampleIndex(SampleTittle);

  MACH3LOG_INFO("Details about sample: {:<20}", SampleTittles[iSample]);
  MACH3LOG_INFO("\t Dimension: {:<35}", Dimensions[iSample]);
  MACH3LOG_INFO("\t nSplineParam: {:<35}", nSplineParams[iSample]);
  MACH3LOG_INFO("\t nOscChan: {:<35}", nOscChans[iSample]);
}

//****************************************
void BinnedSplineHandler::PrintArrayDetails(const std::string& SampleTittle) const
//****************************************
{
  int iSample = getSampleIndex(SampleTittle);
  int nOscChannels = int(indexvec[iSample].size());
  MACH3LOG_INFO("Sample {} has {} oscillation channels", SampleTittle, nOscChannels);
  
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
bool BinnedSplineHandler::isValidSplineIndex(const std::string& SampleTittle, int iOscChan, int iSyst, int iMode, int iVar1, int iVar2, int iVar3)
//****************************************
{
  int iSample = getSampleIndex(SampleTittle);
  bool isValid = true;

  // Lambda to check if an index is valid for a specific dimension
  auto checkIndex = [&isValid](int index, size_t size, const std::string& name) {
    if (index < 0 || index >= int(size)) {
      MACH3LOG_ERROR("{} index is invalid! 0 <= Index < {} ", name, size);
      isValid = false;
    }
  };

  checkIndex(iSample, indexvec.size(), "Sample");
  if (isValid) checkIndex(iOscChan, indexvec[iSample].size(), "OscChan");
  if (isValid) checkIndex(iSyst, indexvec[iSample][iOscChan].size(), "Syst");
  if (isValid) checkIndex(iMode, indexvec[iSample][iOscChan][iSyst].size(), "Mode");
  if (isValid) checkIndex(iVar1, indexvec[iSample][iOscChan][iSyst][iMode].size(), "Var1");
  if (isValid) checkIndex(iVar2, indexvec[iSample][iOscChan][iSyst][iMode][iVar1].size(), "Var2");
  if (isValid) checkIndex(iVar3, indexvec[iSample][iOscChan][iSyst][iMode][iVar1][iVar2].size(), "Var3");

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
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return true;
}

//****************************************
void BinnedSplineHandler::PrintBinning(TAxis *Axis) const
//****************************************
{
  const int NBins = Axis->GetNbins();
  std::string text = "";
  for (int iBin = 0; iBin <= NBins; iBin++) {
    text += fmt::format("{} ", Axis->GetXbins()->GetAt(iBin));
  }
  MACH3LOG_INFO("{}", text);
}

//****************************************
std::vector< std::vector<int> > BinnedSplineHandler::GetEventSplines(const std::string& SampleTittle, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val)
//****************************************
{
  std::vector<std::vector<int>> ReturnVec;
  int SampleIndex = -1;
  for (unsigned int iSample = 0; iSample < SampleTittles.size(); iSample++) {
    if (SampleTittle == SampleTittles[iSample]) {
      SampleIndex = iSample;
    }
  }

  if (SampleIndex == -1) {
    MACH3LOG_ERROR("Sample not found: {}", SampleTittle);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  int nSplineSysts = static_cast<int>(indexvec[SampleIndex][iOscChan].size());


  int Mode = -1;
  std::string SuffixForEventMode = Modes->GetSplineSuffixFromMaCh3Mode(EventMode);
  for (int iMode=0;iMode<Modes->GetNModes();iMode++) {
    if (SuffixForEventMode == Modes->GetSplineSuffixFromMaCh3Mode(iMode)) {
      Mode = iMode;
      break;
    }
  }
  if (Mode == -1) {
    return ReturnVec;
  }

  int Var1Bin = SplineBinning[SampleIndex][iOscChan][0]->FindBin(Var1Val)-1;
  if (Var1Bin < 0 || Var1Bin >= SplineBinning[SampleIndex][iOscChan][0]->GetNbins()) {
    return ReturnVec;
  }

  int Var2Bin = SplineBinning[SampleIndex][iOscChan][1]->FindBin(Var2Val)-1;
  if (Var2Bin < 0 || Var2Bin >= SplineBinning[SampleIndex][iOscChan][1]->GetNbins()) {
    return ReturnVec;
  }

  int Var3Bin = SplineBinning[SampleIndex][iOscChan][2]->FindBin(Var3Val)-1;
  if (Var3Bin < 0 || Var3Bin >= SplineBinning[SampleIndex][iOscChan][2]->GetNbins()){
    return ReturnVec;
  }

  for(int iSyst=0; iSyst<nSplineSysts; iSyst++){
    std::vector<int> spline_modes = SplineModeVecs[SampleIndex][iSyst];
    int nSampleModes = static_cast<int>(spline_modes.size());

    //ETA - look here at the length of spline_modes and what you're actually comparing against
    for(int iMode = 0; iMode<nSampleModes ; iMode++){
      //Only consider if the event mode (Mode) matches ones of the spline modes
      if (Mode == spline_modes[iMode]) {
        int splineID=indexvec[SampleIndex][iOscChan][iSyst][iMode][Var1Bin][Var2Bin][Var3Bin];
        //Also check that the spline isn't flat
        if(!isflatarray[splineID]){
          ReturnVec.push_back({SampleIndex, iOscChan, iSyst, iMode, Var1Bin, Var2Bin, Var3Bin});
        }
      }
    }
  }
  
  return ReturnVec;
}

// checks if there are multiple modes with the same SplineSuffix
// (for example if CCRES and CCCoherent are treated as one spline mode)
std::vector< std::vector<int> > BinnedSplineHandler::StripDuplicatedModes(const std::vector< std::vector<int> >& InputVector) {

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

void BinnedSplineHandler::FillSampleArray(std::string SampleTittle, std::vector<std::string> OscChanFileNames)
{
  int iSample = getSampleIndex(SampleTittle);
  int nOscChannels = nOscChans[iSample];
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
    MACH3LOG_INFO("Processing: {}", OscChanFileNames[iOscChan]);
    TSpline3* mySpline = nullptr;
    TSpline3_red* Spline = nullptr;
    TString Syst, Mode;
    int nKnots, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin = M3::_BAD_INT_;
    double x,y = M3::_BAD_DOUBLE_;
    bool isFlat = true;

    std::set<std::string> SplineFileNames;

    auto File = std::unique_ptr<TFile>(TFile::Open(OscChanFileNames[iOscChan].c_str()));

    if (!File || File->IsZombie()) {
      MACH3LOG_ERROR("File {} not found", OscChanFileNames[iOscChan]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    //This is the MC specific part of the code
    //i.e. we always assume that the splines are just store in  single TDirectory and they're all in there as single objects   
    for (auto k : *File->GetListOfKeys()) {
      auto Key = static_cast<TKey*>(k);
      TClass *Class = gROOT->GetClass(Key->GetClassName(), false);
      if(!Class->InheritsFrom("TSpline3")) {
        continue;
      }

      std::string FullSplineName = std::string(Key->GetName());

      if (SplineFileNames.count(FullSplineName) > 0) {
        MACH3LOG_CRITICAL("Skipping spline - Found a spline whose name has already been encountered before: {}", FullSplineName);
        continue;
      }
      SplineFileNames.insert(FullSplineName);

      std::vector<std::string> Tokens = GetTokensFromSplineName(FullSplineName);

      if (Tokens.size() != kNTokens) {
        MACH3LOG_ERROR("Invalid tokens from spline name - Expected {} tokens. Check implementation in GetTokensFromSplineName()", static_cast<int>(kNTokens));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      
      Syst = Tokens[kSystToken];
      Mode = Tokens[kModeToken];
      Var1Bin = std::stoi(Tokens[kVar1BinToken]);
      Var2Bin = std::stoi(Tokens[kVar2BinToken]);
      Var3Bin = std::stoi(Tokens[kVar3BinToken]);

      SystNum = -1;
      for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++) {
        if (Syst == SplineFileParPrefixNames[iSample][iSyst]) {
          SystNum = iSyst;
          break;
        }
      }

      // If the syst doesn't match any of the spline names then skip it
      if (SystNum == -1){
        MACH3LOG_DEBUG("Couldn't match!!");
        MACH3LOG_DEBUG("Couldn't Match any systematic name in ParameterHandler with spline name: {}" , FullSplineName);
        continue;
      }

      ModeNum = -1;
      for (unsigned int iMode = 0; iMode < SplineModeVecs[iSample][SystNum].size(); iMode++) {
        if (Mode == Modes->GetSplineSuffixFromMaCh3Mode(SplineModeVecs[iSample][SystNum][iMode])) {
          ModeNum = iMode;
          break;
        }
      }

      if (ModeNum == -1) {
      //DB - If you have splines in the root file that you don't want to use (e.g. removing a mode from a syst), this will cause a throw
      //     Therefore include as debug warning and continue instead
        MACH3LOG_DEBUG("Couldn't find mode for {} in {}. Problem Spline is : {} ", Mode, Syst, FullSplineName);
        continue;
      }

      mySpline = Key->ReadObject<TSpline3>();

      if (isValidSplineIndex(SampleTittle, iOscChan, SystNum, ModeNum, Var1Bin, Var2Bin, Var3Bin)) { // loop over all the spline knots and check their value
        MACH3LOG_DEBUG("Pushed back monolith for spline {}", FullSplineName);
        // if the value is 1 then set the flat bool to false
        nKnots = mySpline->GetNp();
        isFlat = true;
        for (int iKnot = 0; iKnot < nKnots; iKnot++) {
          mySpline->GetKnot(iKnot, x, y);
          if (y < 0.99999 || y > 1.00001)
          {
            isFlat = false;
            break;
          }
        }

        //Rather than keeping a mega vector of splines then converting, this should just keep everything nice in memory!
        indexvec[iSample][iOscChan][SystNum][ModeNum][Var1Bin][Var2Bin][Var3Bin]=MonolithIndex;
        coeffindexvec.push_back(CoeffIndex);
        // Should save memory rather saving [x_i_0 ,... x_i_maxknots] for every spline!
        if (isFlat) {
          splinevec_Monolith.push_back(nullptr);
          delete mySpline;
        } else {
          ApplyKnotWeightCapTSpline3(mySpline, SystNum, xsec);
          Spline = new TSpline3_red(mySpline, SplineInterpolationTypes[iSample][SystNum]);
          if(mySpline) delete mySpline;

          splinevec_Monolith.push_back(Spline);
          uniquecoeffindices.push_back(MonolithIndex); //So we can get the unique coefficients and skip flat splines later on!
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

// *****************************************
// Load SplineMonolith from ROOT file
void BinnedSplineHandler::LoadSplineFile(std::string FileName) {
// *****************************************
  if (std::getenv("MACH3") != nullptr) {
    FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
  }
  // Check for spaces in the filename
  size_t pos = FileName.find(' ');
  if (pos != std::string::npos) {
    MACH3LOG_WARN("Filename ({}) contains spaces. Replacing spaces with underscores.", FileName);
    while ((pos = FileName.find(' ')) != std::string::npos) {
      FileName[pos] = '_';
    }
  }
  auto SplineFile = std::make_unique<TFile>(FileName.c_str(), "OPEN");
  LoadSettingsDir(SplineFile);
  LoadMonolithDir(SplineFile);
  LoadIndexDir(SplineFile);
  LoadFastSplineInfoDir(SplineFile);

  for (int iSpline = 0; iSpline < nParams; iSpline++) {
    SplineInfoArray[iSpline].splineParsPointer = xsec->RetPointer(UniqueSystIndices[iSpline]);
  }
  SplineFile->Close();
}

// *****************************************
// KS: Prepare Fast Spline Info within SplineFile
void BinnedSplineHandler::LoadSettingsDir(std::unique_ptr<TFile>& SplineFile) {
// *****************************************
  TTree *Settings = SplineFile->Get<TTree>("Settings");
  int CoeffIndex_temp, MonolithSize_temp;
  short int nParams_temp;
  Settings->SetBranchAddress("CoeffIndex", &CoeffIndex_temp);
  Settings->SetBranchAddress("MonolithSize", &MonolithSize_temp);
  Settings->SetBranchAddress("nParams", &nParams_temp);
  int indexvec_sizes[7];
  for (int i = 0; i < 7; ++i) {
    Settings->SetBranchAddress(("indexvec_size" + std::to_string(i+1)).c_str(), &indexvec_sizes[i]);
  }

  int SplineBinning_size1, SplineBinning_size2, SplineBinning_size3;
  Settings->SetBranchAddress("SplineBinning_size1", &SplineBinning_size1);
  Settings->SetBranchAddress("SplineBinning_size2", &SplineBinning_size2);
  Settings->SetBranchAddress("SplineBinning_size3", &SplineBinning_size3);
  int SplineModeVecs_size1, SplineModeVecs_size2, SplineModeVecs_size3;
  Settings->SetBranchAddress("SplineModeVecs_size1", &SplineModeVecs_size1);
  Settings->SetBranchAddress("SplineModeVecs_size2", &SplineModeVecs_size2);
  Settings->SetBranchAddress("SplineModeVecs_size3", &SplineModeVecs_size3);
  std::vector<std::string>* SampleNames_temp = nullptr;
  Settings->SetBranchAddress("SampleNames", &SampleNames_temp);
  std::vector<std::string>* SampleTittles_temp = nullptr;
  Settings->SetBranchAddress("SampleTittles", &SampleTittles_temp);
  Settings->GetEntry(0);

  CoeffIndex = CoeffIndex_temp;
  MonolithSize = MonolithSize_temp;
  SampleNames = *SampleNames_temp;
  SampleTittles = *SampleTittles_temp;

  nParams = nParams_temp;

  SplineSegments = new short int[nParams]();
  ParamValues = new float[nParams]();

  // Resize indexvec according to saved dimensions
  indexvec.resize(indexvec_sizes[0]);
  for (int i = 0; i < indexvec_sizes[0]; ++i) {
    indexvec[i].resize(indexvec_sizes[1]);
    for (int j = 0; j < indexvec_sizes[1]; ++j) {
      indexvec[i][j].resize(indexvec_sizes[2]);
      for (int k = 0; k < indexvec_sizes[2]; ++k) {
        indexvec[i][j][k].resize(indexvec_sizes[3]);
        for (int l = 0; l < indexvec_sizes[3]; ++l) {
          indexvec[i][j][k][l].resize(indexvec_sizes[4]);
          for (int m = 0; m < indexvec_sizes[4]; ++m) {
            indexvec[i][j][k][l][m].resize(indexvec_sizes[5]);
            for (int n = 0; n < indexvec_sizes[5]; ++n) {
              indexvec[i][j][k][l][m][n].resize(indexvec_sizes[6]);
            }
          }
        }
      }
    }
  }

  auto Resize3D = [](auto& vec, int d1, int d2, int d3) {
    vec.resize(d1);
    for (int i = 0; i < d1; ++i) {
      vec[i].resize(d2);
      for (int j = 0; j < d2; ++j) {
        vec[i][j].resize(d3);
      }
    }
  };

  Resize3D(SplineBinning, SplineBinning_size1, SplineBinning_size2, SplineBinning_size3);
  Resize3D(SplineModeVecs, SplineModeVecs_size1, SplineModeVecs_size2, SplineModeVecs_size3);
}

// *****************************************
// KS: Prepare Fast Spline Info within SplineFile
void BinnedSplineHandler::LoadMonolithDir(std::unique_ptr<TFile>& SplineFile) {
// *****************************************
  TTree *MonolithTree = SplineFile->Get<TTree>("MonolithTree");

  manycoeff_arr = new M3::float_t[CoeffIndex * _nCoeff_];
  MonolithTree->SetBranchAddress("manycoeff", manycoeff_arr);
  isflatarray = new bool[MonolithSize];
  weightvec_Monolith.resize(MonolithSize);
  MonolithTree->SetBranchAddress("isflatarray", isflatarray);

  // Load vectors
  std::vector<int>* coeffindexvec_temp = nullptr;
  MonolithTree->SetBranchAddress("coeffindexvec", &coeffindexvec_temp);
  std::vector<int>* uniquecoeffindices_temp = nullptr;
  MonolithTree->SetBranchAddress("uniquecoeffindices", &uniquecoeffindices_temp);
  std::vector<int>* uniquesplinevec_Monolith_temp = nullptr;
  MonolithTree->SetBranchAddress("uniquesplinevec_Monolith", &uniquesplinevec_Monolith_temp);
  std::vector<int>* UniqueSystIndices_temp = nullptr;
  MonolithTree->SetBranchAddress("UniqueSystIndices", &UniqueSystIndices_temp);

  // Allocate and load xcoeff_arr
  xcoeff_arr = new M3::float_t[CoeffIndex];
  MonolithTree->SetBranchAddress("xcoeff", xcoeff_arr);

  MonolithTree->GetEntry(0);

  coeffindexvec       = *coeffindexvec_temp;
  uniquecoeffindices  = *uniquecoeffindices_temp;
  uniquesplinevec_Monolith = *uniquesplinevec_Monolith_temp;
  UniqueSystIndices = *UniqueSystIndices_temp;
}

// *****************************************
// KS: Prepare Fast Spline Info within SplineFile
void BinnedSplineHandler::LoadIndexDir(std::unique_ptr<TFile>& SplineFile) {
// *****************************************
  TTree *IndexTree = SplineFile->Get<TTree>("IndexVec");

  std::vector<int> Dim(7);
  int value;
  for (int d = 0; d < 7; ++d) {
    IndexTree->SetBranchAddress(Form("dim%d", d+1), &Dim[d]);
  }
  IndexTree->SetBranchAddress("value", &value);

  // Fill indexvec with data from IndexTree
  for (Long64_t i = 0; i < IndexTree->GetEntries(); ++i) {
    IndexTree->GetEntry(i);
    indexvec[Dim[0]][Dim[1]][Dim[2]][Dim[3]][Dim[4]][Dim[5]][Dim[6]] = value;
  }

  // Load SplineBinning data
  TTree *SplineBinningTree = SplineFile->Get<TTree>("SplineBinningTree");
  std::vector<int> indices(3);
  SplineBinningTree->SetBranchAddress("i", &indices[0]);
  SplineBinningTree->SetBranchAddress("j", &indices[1]);
  SplineBinningTree->SetBranchAddress("k", &indices[2]);
  TAxis* axis = nullptr;
  SplineBinningTree->SetBranchAddress("axis", &axis);

  // Reconstruct TAxis objects
  for (Long64_t entry = 0; entry < SplineBinningTree->GetEntries(); ++entry) {
    SplineBinningTree->GetEntry(entry);
    int i = indices[0];
    int j = indices[1];
    int k = indices[2];
    SplineBinning[i][j][k] = static_cast<TAxis*>(axis->Clone());
  }

  std::vector<int> indices_mode(3);
  int mode_value;
  TTree *SplineModeTree = SplineFile->Get<TTree>("SplineModeTree");
  SplineModeTree->SetBranchAddress("i", &indices_mode[0]);
  SplineModeTree->SetBranchAddress("j", &indices_mode[1]);
  SplineModeTree->SetBranchAddress("k", &indices_mode[2]);
  SplineModeTree->SetBranchAddress("value", &mode_value);

  // Fill SplineModeVecs with values from the tree
  for (Long64_t entry = 0; entry < SplineModeTree->GetEntries(); ++entry) {
    SplineModeTree->GetEntry(entry);
    int i = indices_mode[0];
    int j = indices_mode[1];
    int k = indices_mode[2];
    SplineModeVecs[i][j][k] = mode_value;
  }
}

// *****************************************
// Save SplineMonolith into ROOT file
void BinnedSplineHandler::PrepareSplineFile(std::string FileName) {
// *****************************************
  if (std::getenv("MACH3") != nullptr) {
    FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
  }
  // Check for spaces in the filename
  size_t pos = FileName.find(' ');
  if (pos != std::string::npos) {
    MACH3LOG_WARN("Filename ({}) contains spaces. Replacing spaces with underscores.", FileName);
    while ((pos = FileName.find(' ')) != std::string::npos) {
      FileName[pos] = '_';
    }
  }

  auto SplineFile = std::make_unique<TFile>(FileName.c_str(), "recreate");
  PrepareSettingsDir(SplineFile);
  PrepareMonolithDir(SplineFile);
  PrepareIndexDir(SplineFile);
  PrepareOtherInfoDir(SplineFile);
  PrepareFastSplineInfoDir(SplineFile);

  SplineFile->Close();
}

// *****************************************
void BinnedSplineHandler::PrepareSettingsDir(std::unique_ptr<TFile>& SplineFile) const {
// *****************************************
  TTree *Settings = new TTree("Settings", "Settings");
  int CoeffIndex_temp = CoeffIndex;
  int MonolithSize_temp = MonolithSize;
  short int nParams_temp = nParams;

  Settings->Branch("CoeffIndex", &CoeffIndex_temp, "CoeffIndex/I");
  Settings->Branch("MonolithSize", &MonolithSize_temp, "MonolithSize/I");
  Settings->Branch("nParams", &nParams_temp, "nParams/S");

  int indexvec_sizes[7];
  indexvec_sizes[0] = static_cast<int>(indexvec.size());
  indexvec_sizes[1] = (indexvec_sizes[0] > 0) ? static_cast<int>(indexvec[0].size()) : 0;
  indexvec_sizes[2] = (indexvec_sizes[1] > 0) ? static_cast<int>(indexvec[0][0].size()) : 0;
  indexvec_sizes[3] = (indexvec_sizes[2] > 0) ? static_cast<int>(indexvec[0][0][0].size()) : 0;
  indexvec_sizes[4] = (indexvec_sizes[3] > 0) ? static_cast<int>(indexvec[0][0][0][0].size()) : 0;
  indexvec_sizes[5] = (indexvec_sizes[4] > 0) ? static_cast<int>(indexvec[0][0][0][0][0].size()) : 0;
  indexvec_sizes[6] = (indexvec_sizes[5] > 0) ? static_cast<int>(indexvec[0][0][0][0][0][0].size()) : 0;

  for (int i = 0; i < 7; ++i) {
    Settings->Branch(("indexvec_size" + std::to_string(i+1)).c_str(),
                     &indexvec_sizes[i], ("indexvec_size" + std::to_string(i+1) + "/I").c_str());
  }

  int SplineBinning_size1 = static_cast<int>(SplineBinning.size());
  int SplineBinning_size2 = (SplineBinning_size1 > 0) ? static_cast<int>(SplineBinning[0].size()) : 0;
  int SplineBinning_size3 = (SplineBinning_size2 > 0) ? static_cast<int>(SplineBinning[0][0].size()) : 0;

  Settings->Branch("SplineBinning_size1", &SplineBinning_size1, "SplineBinning_size1/I");
  Settings->Branch("SplineBinning_size2", &SplineBinning_size2, "SplineBinning_size2/I");
  Settings->Branch("SplineBinning_size3", &SplineBinning_size3, "SplineBinning_size3/I");

  int SplineModeVecs_size1 = static_cast<int>(SplineModeVecs.size());
  int SplineModeVecs_size2 = (SplineModeVecs_size1 > 0) ? static_cast<int>(SplineModeVecs[0].size()) : 0;
  int SplineModeVecs_size3 = (SplineModeVecs_size2 > 0) ? static_cast<int>(SplineModeVecs[0][0].size()) : 0;

  Settings->Branch("SplineModeVecs_size1", &SplineModeVecs_size1, "SplineModeVecs_size1/I");
  Settings->Branch("SplineModeVecs_size2", &SplineModeVecs_size2, "SplineModeVecs_size2/I");
  Settings->Branch("SplineModeVecs_size3", &SplineModeVecs_size3, "SplineModeVecs_size3/I");

  std::vector<std::string> SampleNames_temp = SampleNames;
  Settings->Branch("SampleNames", &SampleNames_temp);
  std::vector<std::string> SampleTittles_temp = SampleTittles;
  Settings->Branch("SampleTittles", &SampleNames_temp);

  Settings->Fill();
  SplineFile->cd();
  Settings->Write();
  delete Settings;
}

// *****************************************
void BinnedSplineHandler::PrepareMonolithDir(std::unique_ptr<TFile>& SplineFile) const {
// *****************************************
  TTree *MonolithTree = new TTree("MonolithTree", "MonolithTree");
  MonolithTree->Branch("manycoeff", manycoeff_arr, Form("manycoeff[%d]/%s", CoeffIndex * _nCoeff_, M3::float_t_str));
  MonolithTree->Branch("isflatarray", isflatarray, Form("isflatarray[%d]/O", MonolithSize));

  std::vector<int> coeffindexvec_temp = coeffindexvec;
  MonolithTree->Branch("coeffindexvec", &coeffindexvec_temp);
  std::vector<int> uniquecoeffindices_temp = uniquecoeffindices;
  MonolithTree->Branch("uniquecoeffindices", &uniquecoeffindices_temp);
  std::vector<int> uniquesplinevec_Monolith_temp = uniquesplinevec_Monolith;
  MonolithTree->Branch("uniquesplinevec_Monolith", &uniquesplinevec_Monolith_temp);
  std::vector<int> UniqueSystIndices_temp = UniqueSystIndices;
  MonolithTree->Branch("UniqueSystIndices", &UniqueSystIndices_temp);
  MonolithTree->Branch("xcoeff", xcoeff_arr, Form("xcoeff[%d]/%s", CoeffIndex, M3::float_t_str));

  MonolithTree->Fill();
  SplineFile->cd();
  MonolithTree->Write();
  delete MonolithTree;
}

// *****************************************
void BinnedSplineHandler::PrepareIndexDir(std::unique_ptr<TFile>& SplineFile) const {
// *****************************************
  // Create a TTree to store the data
  TTree *IndexTree = new TTree("IndexVec", "IndexVec");

  // Vector holding the 7 dims
  std::vector<int> Dim(7);
  int value;

  // Create branches for each dimension
  for (int d = 0; d < 7; ++d) {
    IndexTree->Branch(Form("dim%d", d+1), &Dim[d], Form("dim%d/I", d+1));
  }
  IndexTree->Branch("value", &value, "value/I");

  // Fill the tree
  for (size_t i = 0; i < indexvec.size(); ++i) {
    for (size_t j = 0; j < indexvec[i].size(); ++j) {
      for (size_t k = 0; k < indexvec[i][j].size(); ++k) {
        for (size_t l = 0; l < indexvec[i][j][k].size(); ++l) {
          for (size_t m = 0; m < indexvec[i][j][k][l].size(); ++m) {
            for (size_t n = 0; n < indexvec[i][j][k][l][m].size(); ++n) {
              for (size_t p = 0; p < indexvec[i][j][k][l][m][n].size(); ++p) {
                Dim[0] = static_cast<int>(i);
                Dim[1] = static_cast<int>(j);
                Dim[2] = static_cast<int>(k);
                Dim[3] = static_cast<int>(l);
                Dim[4] = static_cast<int>(m);
                Dim[5] = static_cast<int>(n);
                Dim[6] = static_cast<int>(p);
                value = static_cast<int>(indexvec[i][j][k][l][m][n][p]);
                IndexTree->Fill();
              }
            }
          }
        }
      }
    }
  }

  SplineFile->cd();
  // Write the tree to the file
  IndexTree->Write();
  delete IndexTree;
}

// *****************************************
void BinnedSplineHandler::PrepareOtherInfoDir(std::unique_ptr<TFile>& SplineFile) const {
// *****************************************
  // Create a new tree for SplineBinning data
  TTree *SplineBinningTree = new TTree("SplineBinningTree", "SplineBinningTree");
  std::vector<int> indices(3); // To store the 3D indices
  TAxis* axis = nullptr;
  SplineBinningTree->Branch("i", &indices[0], "i/I");
  SplineBinningTree->Branch("j", &indices[1], "j/I");
  SplineBinningTree->Branch("k", &indices[2], "k/I");
  SplineBinningTree->Branch("axis", "TAxis", &axis);

  // Fill the SplineBinningTree
  for (size_t i = 0; i < SplineBinning.size(); ++i) {
    for (size_t j = 0; j < SplineBinning[i].size(); ++j) {
      for (size_t k = 0; k < SplineBinning[i][j].size(); ++k) {
        axis = SplineBinning[i][j][k];
        indices[0] = static_cast<int>(i);
        indices[1] = static_cast<int>(j);
        indices[2] = static_cast<int>(k);
        SplineBinningTree->Fill();
      }
    }
  }
  SplineFile->cd();
  SplineBinningTree->Write();
  delete SplineBinningTree;

  std::vector<int> indices_mode(3); // to store 3D indices
  int mode_value;

  TTree *SplineModeTree = new TTree("SplineModeTree", "SplineModeTree");
  // Create branches for indices and value
  SplineModeTree->Branch("i", &indices_mode[0], "i/I");
  SplineModeTree->Branch("j", &indices_mode[1], "j/I");
  SplineModeTree->Branch("k", &indices_mode[2], "k/I");
  SplineModeTree->Branch("value", &mode_value, "value/I");

  // Fill the tree
  for (size_t i = 0; i < SplineModeVecs.size(); ++i) {
    for (size_t j = 0; j < SplineModeVecs[i].size(); ++j) {
      for (size_t k = 0; k < SplineModeVecs[i][j].size(); ++k) {
        indices_mode[0] = static_cast<int>(i);
        indices_mode[1] = static_cast<int>(j);
        indices_mode[2] = static_cast<int>(k);
        mode_value = SplineModeVecs[i][j][k];
        SplineModeTree->Fill();
      }
    }
  }
  // Write the tree to the file
  SplineFile->cd();
  SplineModeTree->Write();
  delete SplineModeTree;
}
