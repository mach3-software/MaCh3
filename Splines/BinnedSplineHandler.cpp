#include "BinnedSplineHandler.h"
#include <memory>

#pragma GCC diagnostic ignored "-Wuseless-cast"

_MaCh3_Safe_Include_Start_ //{
#include "TROOT.h"
#include "TKey.h"
#include "TH3F.h"
_MaCh3_Safe_Include_End_ //}

//****************************************
BinnedSplineHandler::BinnedSplineHandler(ParameterHandlerGeneric *ParHandler_, MaCh3Modes *Modes_) : SplineBase() {
//****************************************
  if (!ParHandler_) {
    MACH3LOG_ERROR("Trying to create BinnedSplineHandler with uninitialized covariance object");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  ParHandler = ParHandler_;

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
void BinnedSplineHandler::CleanUpMemory() {
//****************************************
  //Call once everything's been allocated in SampleHandlerBase, cleans up junk from memory!
  //Not a huge saving but it's better than leaving everything up to the compiler
  MACH3LOG_INFO("Cleaning up spline memory");
  CleanVector(IndexVect);
  IndexVectMap.clear();
  CleanVector(SplineFileParPrefixNames);
  CleanVector(GlobalSystIndex);
  CleanVector(SplineModeVecs);
  CleanVector(UniqueSystNames);
  CleanVector(SplineInterpolationTypes);
  CleanVector(nOscChans);
  CleanVector(nSplineParams);
  CleanVector(DimensionLabels);
  CleanVector(SampleNames);
  CleanVector(SampleTitles);
  CleanVector(Dimensions);
  CleanContainer(splinevec_Monolith);
  CleanContainer(SplineBinning);
  CleanVector(UniqueSystIndices);
  if(isflatarray) delete [] isflatarray;
}

//****************************************
//Adds samples to the large array
void BinnedSplineHandler::AddSample(const std::string& SampleName,
                                    const std::string& SampleTitle,
                                    const std::vector<std::string>& OscChanFileNames,
                                    const std::vector<std::string>& SplineVarNames) {
//****************************************
  SampleNames.push_back(SampleName);
  SampleTitles.push_back(SampleTitle);
  Dimensions.push_back(static_cast<int>(SplineVarNames.size()));
  DimensionLabels.push_back(SplineVarNames);

  int nSplineParam = ParHandler->GetNumParamsFromSampleName(SampleName, SystType::kSpline);
  nSplineParams.push_back(nSplineParam);

  //This holds the global index of the spline i.e. 0 -> _fNumPar
  std::vector<int> GlobalSystIndex_Sample = ParHandler->GetGlobalSystIndexFromSampleName(SampleName, SystType::kSpline);
  //Keep track of this for all the samples
  GlobalSystIndex.push_back(GlobalSystIndex_Sample);

  std::vector<SplineInterpolation> SplineInterpolation_Sample = ParHandler->GetSplineInterpolationFromSampleName(SampleName);
  // Keep track of this for all samples
  SplineInterpolationTypes.push_back(SplineInterpolation_Sample);

  std::vector<std::string> SplineFileParPrefixNames_Sample = ParHandler->GetSplineParsNamesFromSampleName(SampleName);
  SplineFileParPrefixNames.push_back(SplineFileParPrefixNames_Sample);

  MACH3LOG_INFO("Create SplineModeVecs_Sample");
  std::vector<std::vector<int>> SplineModeVecs_Sample = StripDuplicatedModes(ParHandler->GetSplineModeVecFromSampleName(SampleName));
  MACH3LOG_INFO("SplineModeVecs_Sample is of size {}", SplineModeVecs_Sample.size());
  SplineModeVecs.push_back(SplineModeVecs_Sample);

  MACH3LOG_INFO("SplineModeVecs is of size {}", SplineModeVecs.size());

  int nOscChan = int(OscChanFileNames.size());
  nOscChans.push_back(nOscChan);

  PrintSampleDetails(SampleTitle);

  std::vector<std::vector<TAxis *>> SampleBinning(nOscChan);
  for (int iOscChan = 0; iOscChan < nOscChan; iOscChan++)
  {
    SampleBinning[iOscChan] = FindSplineBinning(OscChanFileNames[iOscChan], SampleTitle);
  }
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
  SplineBinning.push_back(SampleBinning);

  BuildSampleIndexingArray(SampleTitle);
  PrintArrayDetails(SampleTitle);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");

  FillSampleArray(SampleTitle, OscChanFileNames);
  MACH3LOG_INFO("#----------------------------------------------------------------------------------------------------------------------------------#");
}

//****************************************
void BinnedSplineHandler::InvestigateMissingSplines() const {
//****************************************
  // Map: iSample → iSyst → modeSuffix → {totalSplines, zeroCount}
  std::map<unsigned int, std::map<unsigned int, std::map<std::string, std::pair<unsigned int, unsigned int>>>> systZeroCounts;
  for (const auto& entry : IndexVect) {
    unsigned int iSample = entry.iSample;
    unsigned int iSyst   = entry.iSyst;

    std::string SampleName = SampleNames[iSample];
    // Get the mode suffix string
    std::string modeSuffix = Modes->GetSplineSuffixFromMaCh3Mode(SplineModeVecs[iSample][entry.iSyst][entry.iMode]);
    auto& counts = systZeroCounts[iSample][iSyst][modeSuffix];

    counts.first++; // totalSplines
    if (entry.value == -1)
    {
      counts.second++; // zeroCount
      if (counts.second > 1)
      {
        MACH3LOG_DEBUG(
          "Sample '{}' | OscChan {} | Syst '{}' | Mode '{}' | Var {} => Value: {}",
          SampleTitles[iSample],
          entry.iOscChan,
          SplineFileParPrefixNames[iSample][iSyst],
          modeSuffix,
          fmt::join(entry.iVar, " "),
          entry.value
        );
      }
    }
  }

  // KS: Let's print this atrocious mess...
  for (const auto& samplePair : systZeroCounts) {
    unsigned int iSample = samplePair.first;
    std::vector<std::string> SplineFileParPrefixNames_Sample = ParHandler->GetParsNamesFromSampleName(SampleNames[iSample], kSpline);
    for (const auto& systPair : samplePair.second) {
      unsigned int iSyst = systPair.first;
      const auto& systName = SplineFileParPrefixNames_Sample[iSyst];
      for (const auto& modePair : systPair.second) {
        const auto& modeSuffix = modePair.first;
        const auto& counts = modePair.second;
        MACH3LOG_CRITICAL(
          "Sample '{}': Systematic '{}' has missing splines in mode '{}'. Expected Splines: {}, Missing Splines: {}",
          SampleTitles[iSample],
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
void BinnedSplineHandler::TransferToMonolith() {
//****************************************
  PrepForReweight(); 
  MonolithSize = CountNumberOfLoadedSplines(false, 1);

  if(MonolithSize != MonolithIndex){
    InvestigateMissingSplines();
    MACH3LOG_ERROR("Something's gone wrong when we tried to get the size of your monolith");
    MACH3LOG_ERROR("MonolithSize is {}", MonolithSize);
    MACH3LOG_ERROR("MonolithIndex is {}", MonolithIndex);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  MACH3LOG_INFO("Now transferring splines to a monolith if size {}", MonolithSize);
  // Maps single spline object with single parameter
  uniquesplinevec_Monolith.resize(MonolithSize);
  weightvec_Monolith.resize(MonolithSize);
  isflatarray = new bool[MonolithSize];
  
  xcoeff_arr = new M3::float_t[CoeffIndex];
  manycoeff_arr = new M3::float_t[CoeffIndex*_nCoeff_];

  for (const auto& entry : IndexVect) {
    int splineindex = entry.value;
    weightvec_Monolith[splineindex] = 1.0;

    bool foundUniqueSpline = false;
    // We are trying to match Spline Object with single parameter (like MAQE)
    for (int iUniqueSyst = 0; iUniqueSyst < nParams; iUniqueSyst++)
    {
      if (SplineFileParPrefixNames[entry.iSample][entry.iSyst] == UniqueSystNames[iUniqueSyst])
      {
        uniquesplinevec_Monolith[splineindex] = iUniqueSyst;
        foundUniqueSpline = true;
        break;
      }
    } //unique syst loop end

    // If current spline object hasn't been matched with actual parameter this means misconfiguration
    if (!foundUniqueSpline)
    {
      MACH3LOG_ERROR("Unique spline index not found");
      MACH3LOG_ERROR("For Spline {}", SplineFileParPrefixNames[entry.iSample][entry.iSyst]);
      MACH3LOG_ERROR("Couldn't match {} with any of the following {} systs:", SplineFileParPrefixNames[entry.iSample][entry.iSyst], nParams);
      for (int iUniqueSyst = 0; iUniqueSyst < nParams; iUniqueSyst++)
      {
        MACH3LOG_ERROR("{},", UniqueSystNames.at(iUniqueSyst));
      }//unique syst loop end
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }

    int splineKnots;
    if(splinevec_Monolith[splineindex]){
      isflatarray[splineindex] = false;
      splineKnots=splinevec_Monolith[splineindex]->GetNp();

      //Now to fill up our coefficient arrayss
      M3::float_t* tmpXCoeffArr = new M3::float_t[splineKnots];
      M3::float_t* tmpManyCoeffArr = new M3::float_t[splineKnots*_nCoeff_];

      int iCoeff=coeffindexvec[splineindex];
      GetSplineCoeff_SepMany(splineindex, tmpXCoeffArr, tmpManyCoeffArr);

      for(int i = 0; i < splineKnots; i++){
        xcoeff_arr[iCoeff+i] = tmpXCoeffArr[i];

        for(int j = 0; j < _nCoeff_; j++){
          manycoeff_arr[(iCoeff+i)*_nCoeff_+j]=tmpManyCoeffArr[i*_nCoeff_+j];
        }
      }
      delete[] tmpXCoeffArr;
      delete[] tmpManyCoeffArr;
    } else {
      isflatarray[splineindex]=true;
    }
  }
}

// *****************************************
void BinnedSplineHandler::Evaluate() {
// *****************************************
  // There's a parameter mapping that goes from spline parameter to a global parameter index
  // Find the spline segments
  FindSplineSegment();

  //KS: Huge MP loop over all valid splines
  CalcSplineWeights();
}

//****************************************
void BinnedSplineHandler::CalcSplineWeights() {
//****************************************
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

    // Get the variation for this reconfigure for the i-th parameter
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
void BinnedSplineHandler::BuildSampleIndexingArray(const std::string& SampleTitle) {
//****************************************
  int iSample = GetSampleIndex(SampleTitle);
  for (int iOscChan = 0; iOscChan < nOscChans[iSample]; ++iOscChan)
  {
    for (int iSyst = 0; iSyst < nSplineParams[iSample]; ++iSyst)
    {
      int nModesInSyst = static_cast<int>(SplineModeVecs[iSample][iSyst].size());
      for (int iMode = 0; iMode < nModesInSyst; ++iMode)
      {
        const int nBins1 = SplineBinning[iSample][iOscChan][0]->GetNbins();
        const int nBins2 = SplineBinning[iSample][iOscChan][1]->GetNbins();
        const int nBins3 = SplineBinning[iSample][iOscChan][2]->GetNbins();
        for (int iVar1 = 0; iVar1 < nBins1; ++iVar1) {
          for (int iVar2 = 0; iVar2 < nBins2; ++iVar2) {
            for (int iVar3 = 0; iVar3 < nBins3; ++iVar3) {
              SplineIndex entry;
              entry.value    = -1;
              entry.iSample  = iSample;
              entry.iOscChan = iOscChan;
              entry.iSyst    = iSyst;
              entry.iMode    = iMode;
              entry.iVar = {iVar1, iVar2, iVar3};
              IndexVect.push_back(entry);
              IndexVectMap[std::make_tuple(iSample, iOscChan, iSyst, iMode, std::vector<int>{iVar1, iVar2, iVar3})] = static_cast<int>(IndexVect.size() - 1);
            }
          }
        }
      }
    }
  }
}

//****************************************
std::vector<TAxis *> BinnedSplineHandler::FindSplineBinning(const std::string& FileName, const std::string& SampleTitle) {
//****************************************
  int iSample = GetSampleIndex(SampleTitle);

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
const M3::float_t* BinnedSplineHandler::RetPointer(const SplineIndex& Variables) const {
//****************************************
  int Index = IndexVectMap.at(std::make_tuple(Variables.iSample, Variables.iOscChan, Variables.iSyst,
                                              Variables.iMode, Variables.iVar));
  return &weightvec_Monolith[IndexVect[Index].value];
}

//****************************************
int BinnedSplineHandler::CountNumberOfLoadedSplines(bool NonFlat, int Verbosity) const {
//****************************************
  std::vector<int> SampleAll(SampleTitles.size(), 0);
  std::vector<int> SampleNonFlat(SampleTitles.size(), 0);

  int FullCounter_All = 0;
  int FullCounter_NonFlat = 0;

  for (unsigned int index = 0; index < IndexVect.size(); index++)
  {
    const auto& entry = IndexVect[index];

    int iSample = entry.iSample;

    std::string SampleTitle = SampleTitles[iSample];

    if (!isValidSplineIndex(SampleTitle, entry.iOscChan, entry.iSyst,
         entry.iMode, entry.iVar)) {
      continue;
    }
    SampleAll[iSample]++;

    if (splinevec_Monolith[entry.value]) {
      SampleNonFlat[iSample]++;
    }
  }

  // Print per-sample summary
  for (size_t iSample = 0; iSample < SampleTitles.size(); iSample++)
  {
    MACH3LOG_DEBUG("{:<10} has {:<10} splines, of which {:<10} are not flat",
                   SampleTitles[iSample],
                   SampleAll[iSample],
                   SampleNonFlat[iSample]);

    FullCounter_All += SampleAll[iSample];
    FullCounter_NonFlat += SampleNonFlat[iSample];
  }

  if (Verbosity > 0) {
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
  bool FoundNonFlatSpline = false;
  int SampleCounter = M3::_BAD_INT_;
  int SystCounter = M3::_BAD_INT_;
  // DB Find all the Unique systs across each sample and oscillation channel
  //    This assumes that each occurrence of the same systematic spline has the same knot spacing
  //    Which is a reasonable assumption for the current implementation of spline evaluations
  for (const auto& entry : IndexVect)
  {
    int splineindex = entry.value;
    // KS: reset if we moved to another sample or syst
    if(SampleCounter != entry.iSample || SystCounter != entry.iSyst) {
      SampleCounter = entry.iSample;
      SystCounter = entry.iSyst;
      FoundNonFlatSpline = false;
    }
    std::string SystName = SplineFileParPrefixNames[entry.iSample][entry.iSyst];
    bool FoundSyst = false;
    for (unsigned int iFoundSyst = 0; iFoundSyst < UniqueSystNames.size(); iFoundSyst++) {
      if (SystName == UniqueSystNames[iFoundSyst]) {
        FoundSyst = true;
      }
    }
    if (FoundSyst) continue;

    if (splinevec_Monolith[splineindex])
    {
      UniqueSystSplines.push_back(splinevec_Monolith[splineindex]);
      UniqueSystIndices.push_back(GlobalSystIndex[entry.iSample][entry.iSyst]);

      FoundNonFlatSpline = true;
    }
    if (FoundNonFlatSpline) {
      UniqueSystNames.push_back(SystName);
    } else {
      MACH3LOG_DEBUG("{} syst has no response in sample {}", SystName, entry.iSample);
      MACH3LOG_DEBUG("Whilst this isn't necessarily a problem, it seems odd");
    }
  } // end loop over indices
  nParams = static_cast<short int>(UniqueSystSplines.size());

  // DB Find the number of splines knots which assumes each instance of the syst has the same number of knots
  SplineSegments = new short int[nParams]();
  ParamValues = new float[nParams]();
  SplineInfoArray.resize(nParams);
  for (int iSpline = 0; iSpline < nParams; iSpline++)
  {
    SplineInfoArray[iSpline].nPts = static_cast<M3::int_t>(UniqueSystSplines[iSpline]->GetNp());
    SplineInfoArray[iSpline].xPts.resize(SplineInfoArray[iSpline].nPts);
    SplineInfoArray[iSpline].splineParsPointer = ParHandler->RetPointer(UniqueSystIndices[iSpline]);
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

  int nCombinations_FlatSplines = 0;
  int nCombinations_All = 0;
  // DB Now actually loop over splines to determine which are all null i.e. flat
  for (unsigned int index = 0; index < IndexVect.size(); index++)
  {
    int splineindex = IndexVect[index].value;;
    nCombinations_All++;
    if (!splinevec_Monolith[splineindex]) {
      nCombinations_FlatSplines++;
    }
  }

  // We need to grab the maximum number of knots
  MACH3LOG_INFO("Number of combinations of Sample, OscChan, Syst and Mode which have entirely flat response: {} / {}", nCombinations_FlatSplines, nCombinations_All);
}

//****************************************
// Rather work with spline coefficients in the splines, let's copy ND and use coefficient arrays
void BinnedSplineHandler::GetSplineCoeff_SepMany(int splineindex, M3::float_t* &xArray, M3::float_t* &manyArray) {
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
//Returns sample index in
int BinnedSplineHandler::GetSampleIndex(const std::string& SampleTitle) const {
//****************************************
  for (size_t iSample = 0; iSample < SampleTitles.size(); ++iSample) {
    if (SampleTitle == SampleTitles[iSample]) {
      return static_cast<int>(iSample);
    }
  }
  MACH3LOG_ERROR("Sample name not found: {}", SampleTitle);
  throw MaCh3Exception(__FILE__, __LINE__);
}

//****************************************
void BinnedSplineHandler::PrintSampleDetails(const std::string& SampleTitle) const {
//****************************************
  const int iSample = GetSampleIndex(SampleTitle);

  MACH3LOG_INFO("Details about sample: {:<20}", SampleTitles[iSample]);
  MACH3LOG_INFO("\t Dimension: {:<35}", Dimensions[iSample]);
  MACH3LOG_INFO("\t nSplineParam: {:<35}", nSplineParams[iSample]);
  MACH3LOG_INFO("\t nOscChan: {:<35}", nOscChans[iSample]);
}

//****************************************
void BinnedSplineHandler::PrintArrayDetails(const std::string& SampleTitle) const {
//****************************************
  int iSample = GetSampleIndex(SampleTitle);
  // count oscillation channels
  std::map<int, std::set<int>> OscToSysts;

  for (const auto& entry : IndexVect) {
    if (entry.iSample != iSample) continue;
    OscToSysts[entry.iOscChan].insert(entry.iSyst);
  }
  MACH3LOG_INFO("Sample {} has {} oscillation channels", SampleTitle, OscToSysts.size());

  for (const auto& OscPair : OscToSysts) {
    int osc = OscPair.first;
    const std::set<int>& SystSet = OscPair.second;
    MACH3LOG_INFO("Oscillation channel {} has {} systematics", osc, SystSet.size());
  }
}

//****************************************
bool BinnedSplineHandler::isValidSplineIndex(const std::string& SampleTitle, int iOscChan,
                                             int iSyst, int iMode, const std::vector<int>& iVar) const {
//****************************************
  int iSample = GetSampleIndex(SampleTitle);

  bool found = IndexVectMap.find(std::make_tuple(iSample, iOscChan, iSyst, iMode, iVar)) != IndexVectMap.end();

  if (!found)
  {
    MACH3LOG_ERROR("Given iSample: {}", iSample);
    MACH3LOG_ERROR("Given iOscChan: {}", iOscChan);
    MACH3LOG_ERROR("Given iSyst: {}", iSyst);
    MACH3LOG_ERROR("Given iMode: {}", iMode);
    for (size_t i = 0; i < iVar.size(); ++i) {
      MACH3LOG_ERROR("Given iVar{}: {}", i, iVar[i]);
    }
    MACH3LOG_ERROR("Come visit me at : {} : {}", __FILE__, __LINE__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return true;
}

//****************************************
void BinnedSplineHandler::PrintBinning(TAxis *Axis) const {
//****************************************
  const int NBins = Axis->GetNbins();
  std::string text = "";
  for (int iBin = 0; iBin <= NBins; iBin++) {
    text += fmt::format("{} ", Axis->GetXbins()->GetAt(iBin));
  }
  MACH3LOG_INFO("{}", text);
}

//****************************************
std::vector<SplineIndex> BinnedSplineHandler::GetEventSplines(const std::string& SampleTitle,
                                                              int iOscChan, int EventMode, const std::vector<double>& VarVals) {
//****************************************
  std::vector<SplineIndex> ReturnVec;
  int SampleIndex = GetSampleIndex(SampleTitle);

  int Mode = -1;
  std::string SuffixForEventMode = Modes->GetSplineSuffixFromMaCh3Mode(EventMode);
  for (int iMode = 0; iMode< Modes->GetNModes(); iMode++) {
    if (SuffixForEventMode == Modes->GetSplineSuffixFromMaCh3Mode(iMode)) {
      Mode = iMode;
      break;
    }
  }
  if (Mode == -1) {
    return ReturnVec;
  }

  std::vector<int> var_bins;
  for (size_t i = 0; i < VarVals.size(); ++i) {
    int bin = SplineBinning[SampleIndex][iOscChan][i]->FindBin(VarVals[i]) - 1;
    if (bin < 0 || bin >= SplineBinning[SampleIndex][iOscChan][i]->GetNbins()) {
      return ReturnVec;
    }
    var_bins.push_back(bin);
  }

  for(int iSyst=0; iSyst < nSplineParams[SampleIndex]; iSyst++){
    std::vector<int> spline_modes = SplineModeVecs[SampleIndex][iSyst];
    int nSampleModes = static_cast<int>(spline_modes.size());

    //ETA - look here at the length of spline_modes and what you're actually comparing against
    for(int iMode = 0; iMode<nSampleModes; iMode++) {
      //Only consider if the event mode (Mode) matches ones of the spline modes
      if (Mode == spline_modes[iMode]) {
        int index = IndexVectMap.at(std::make_tuple(SampleIndex, iOscChan, iSyst, iMode,
                                                    var_bins));
        int splineID = IndexVect[index].value;
        //Also check that the spline isn't flat
        if(!isflatarray[splineID]) {
          SplineIndex idx;
          idx.iSample  = SampleIndex;
          idx.iOscChan = iOscChan;
          idx.iSyst    = iSyst;
          idx.iMode    = iMode;
          idx.iVar     = var_bins;
          ReturnVec.push_back(idx);
        }
      }
    }
  }
  return ReturnVec;
}

//****************************************
// checks if there are multiple modes with the same SplineSuffix
// (for example if CCRES and CCCoherent are treated as one spline mode)
std::vector< std::vector<int> > BinnedSplineHandler::StripDuplicatedModes(const std::vector< std::vector<int> >& InputVector) const {
//****************************************
  //ETA - this is of size nPars from the syst model
  size_t InputVectorSize = InputVector.size();
  std::vector< std::vector<int> > ReturnVec(InputVectorSize);

  //ETA - loop over all systematics
  for (size_t iSyst=0;iSyst<InputVectorSize;iSyst++) {
    std::vector<int> TmpVec;
    std::vector<std::string> TestVec;

    //Loop over the modes that we've listed in ParHandler
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

void BinnedSplineHandler::FillSampleArray(const std::string& SampleTitle, const std::vector<std::string>& OscChanFileNames)
{
  int iSample = GetSampleIndex(SampleTitle);
  int nOscChannels = nOscChans[iSample];
  
  for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
    MACH3LOG_INFO("Processing: {}", OscChanFileNames[iOscChan]);
    TSpline3* mySpline = nullptr;
    TSpline3_red* Spline = nullptr;

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
      
      TString Syst = Tokens[kSystToken];
      TString Mode = Tokens[kModeToken];
      std::vector<int> VarBins;
      VarBins.reserve(Tokens.size() - kVarBinToken);
      for (std::size_t i = kVarBinToken; i < Tokens.size(); ++i) {
        VarBins.push_back(std::stoi(Tokens.at(i)));
      }

      int SystNum = -1;
      for (unsigned iSyst = 0; iSyst < SplineFileParPrefixNames[iSample].size(); iSyst++) {
        if (Syst == SplineFileParPrefixNames[iSample][iSyst]) {
          SystNum = iSyst;
          break;
        }
      }

      // If the syst doesn't match any of the spline names then skip it
      if (SystNum == -1){
        MACH3LOG_DEBUG("Couldn't Match any systematic name in ParameterHandler with spline name: {}" , FullSplineName);
        continue;
      }

      int ModeNum = -1;
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
      // loop over all the spline knots and check their value
      if (isValidSplineIndex(SampleTitle, iOscChan, SystNum, ModeNum, VarBins)) {
        MACH3LOG_TRACE("Pushed back monolith for spline {}", FullSplineName);
        // if the value is 1 then set the flat bool to false
        int nKnots = mySpline->GetNp();
        bool isFlat = true;
        for (int iKnot = 0; iKnot < nKnots; iKnot++) {
          double x, y = M3::_BAD_DOUBLE_;
          mySpline->GetKnot(iKnot, x, y);
          if (y < 0.99999 || y > 1.00001)
          {
            isFlat = false;
            break;
          }
        }

        //Rather than keeping a mega vector of splines then converting, this should just keep everything nice in memory!
        int index = IndexVectMap.at(std::make_tuple(iSample, iOscChan, SystNum, ModeNum, VarBins));
        IndexVect[index].value = MonolithIndex;
        coeffindexvec.push_back(CoeffIndex);
        // Should save memory rather saving [x_i_0 ,... x_i_maxknots] for every spline!
        if (isFlat) {
          splinevec_Monolith.push_back(nullptr);
          delete mySpline;
        } else {
          ApplyKnotWeightCapTSpline3(mySpline, SystNum, ParHandler);
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
  M3::AddPath(FileName);

  // Check for spaces in the filename
  size_t pos = FileName.find(' ');
  if (pos != std::string::npos) {
    MACH3LOG_WARN("Filename ({}) contains spaces. Replacing spaces with underscores.", FileName);
    while ((pos = FileName.find(' ')) != std::string::npos) {
      FileName[pos] = '_';
    }
  }
  auto SplineFile = std::make_unique<TFile>(FileName.c_str(), "OPEN");

  TMacro *ConfigCov = SplineFile->Get<TMacro>("ParameterHandler");
  // Config which was in MCMC from which we are starting
  YAML::Node CovSettings = TMacroToYAML(*ConfigCov);
  // Config from currently used cov object
  YAML::Node ConfigCurrent = ParHandler->GetConfig();

  if (!compareYAMLNodes(CovSettings, ConfigCurrent))
  {
    MACH3LOG_ERROR("Loading precomputed spline file, however encountered different YAML config, please regenerate input");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  LoadSettingsDir(SplineFile);
  LoadMonolithDir(SplineFile);
  LoadIndexDir(SplineFile);
  LoadFastSplineInfoDir(SplineFile);

  for (int iSpline = 0; iSpline < nParams; iSpline++) {
    SplineInfoArray[iSpline].splineParsPointer = ParHandler->RetPointer(UniqueSystIndices[iSpline]);
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
  std::vector<std::string>* SampleTitles_temp = nullptr;
  Settings->SetBranchAddress("SampleTitles", &SampleTitles_temp);
  std::vector<int>* nSplineParams_temp = nullptr;
  Settings->SetBranchAddress("nSplineParams", &nSplineParams_temp);
  Settings->GetEntry(0);

  CoeffIndex = CoeffIndex_temp;
  MonolithSize = MonolithSize_temp;
  SampleNames = *SampleNames_temp;
  SampleTitles = *SampleTitles_temp;
  nSplineParams = *nSplineParams_temp;

  nParams = nParams_temp;

  SplineSegments = new short int[nParams]();
  ParamValues = new float[nParams]();

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

  SplineIndex* IndexTemp = nullptr;
  IndexTree->SetBranchAddress("SplineIndex", &IndexTemp);
  IndexVect.resize(IndexTree->GetEntries());
  // Fill indexvec with data from IndexTree
  for (Long64_t iEntry = 0; iEntry < IndexTree->GetEntries(); ++iEntry) {
    IndexTree->GetEntry(iEntry);
    IndexVect[iEntry] = *IndexTemp;

    auto key = std::make_tuple(IndexTemp->iSample, IndexTemp->iOscChan, IndexTemp->iSyst,
                               IndexTemp->iMode, IndexTemp->iVar);
    IndexVectMap[key] = static_cast<int>(iEntry);
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
  M3::AddPath(FileName);
  // Check for spaces in the filename
  size_t pos = FileName.find(' ');
  if (pos != std::string::npos) {
    MACH3LOG_WARN("Filename ({}) contains spaces. Replacing spaces with underscores.", FileName);
    while ((pos = FileName.find(' ')) != std::string::npos) {
      FileName[pos] = '_';
    }
  }
  // Save ROOT File
  auto SplineFile = std::make_unique<TFile>(FileName.c_str(), "recreate");
  YAML::Node ConfigCurrent = ParHandler->GetConfig();
  TMacro ConfigSave = YAMLtoTMacro(ConfigCurrent, "ParameterHandler");
  ConfigSave.Write();

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
  std::vector<std::string> SampleTitles_temp = SampleTitles;
  Settings->Branch("SampleTitles", &SampleTitles_temp);
  std::vector<int> nSplineParams_temp = nSplineParams;
  Settings->Branch("nSplineParams", &nSplineParams_temp);

  Settings->Fill();
  SplineFile->cd();
  Settings->Write();
  delete Settings;
}

// *****************************************
void BinnedSplineHandler::PrepareMonolithDir(std::unique_ptr<TFile>& SplineFile) const {
// *****************************************
  TTree *MonolithTree = new TTree("MonolithTree", "MonolithTree");
  MonolithTree->Branch("manycoeff", manycoeff_arr, Form("manycoeff[%d]/%s", CoeffIndex * _nCoeff_, M3::float_t_str_repr));
  MonolithTree->Branch("isflatarray", isflatarray, Form("isflatarray[%d]/O", MonolithSize));

  std::vector<int> coeffindexvec_temp = coeffindexvec;
  MonolithTree->Branch("coeffindexvec", &coeffindexvec_temp);
  std::vector<int> uniquecoeffindices_temp = uniquecoeffindices;
  MonolithTree->Branch("uniquecoeffindices", &uniquecoeffindices_temp);
  std::vector<int> uniquesplinevec_Monolith_temp = uniquesplinevec_Monolith;
  MonolithTree->Branch("uniquesplinevec_Monolith", &uniquesplinevec_Monolith_temp);
  std::vector<int> UniqueSystIndices_temp = UniqueSystIndices;
  MonolithTree->Branch("UniqueSystIndices", &UniqueSystIndices_temp);
  MonolithTree->Branch("xcoeff", xcoeff_arr, Form("xcoeff[%d]/%s", CoeffIndex, M3::float_t_str_repr));

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
  SplineIndex entry;
  IndexTree->Branch("SplineIndex", &entry);
  for (const auto& e : IndexVect) {
    entry = e;
    IndexTree->Fill();
  }

  SplineFile->cd();
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
