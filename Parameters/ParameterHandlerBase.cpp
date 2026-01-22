#include "Parameters/ParameterHandlerBase.h"

#include "Samples/HistogramUtils.h"

#include <regex>

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string name, std::string file, double threshold, int FirstPCA, int LastPCA)
                     : inputFile(file), pca(true) {
// ********************************************
  MACH3LOG_DEBUG("Constructing instance of ParameterHandler");
  doSpecialStepProposal = false;
  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("NOTE: {} {}", name, file);
    MACH3LOG_INFO("Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }
  Init(name, file);

  // Call the innocent helper function
  if (pca) ConstructPCA(threshold, FirstPCA, LastPCA);
}
// ********************************************
ParameterHandlerBase::ParameterHandlerBase(const std::vector<std::string>& YAMLFile, std::string name, double threshold, int FirstPCA, int LastPCA)
                     : inputFile(YAMLFile[0].c_str()), matrixName(name), pca(true) {
// ********************************************
  MACH3LOG_INFO("Constructing instance of ParameterHandler using");
  doSpecialStepProposal = false;
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    MACH3LOG_INFO("{}", YAMLFile[i]);
  }
  MACH3LOG_INFO("as an input");

  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }

  Init(YAMLFile);
  // Call the innocent helper function
  if (pca) ConstructPCA(threshold, FirstPCA, LastPCA);
}

// ********************************************
//Destructor
ParameterHandlerBase::~ParameterHandlerBase(){
// ********************************************
  delete[] randParams;
  delete[] corr_throw;

  if (covMatrix != nullptr) delete covMatrix;
  if (invCovMatrix != nullptr) delete invCovMatrix;
  if (throwMatrix != nullptr) delete throwMatrix;
  for(int i = 0; i < _fNumPar; i++) {
    delete[] throwMatrixCholDecomp[i];
  }
  delete[] throwMatrixCholDecomp;
}

// ********************************************
void ParameterHandlerBase::ConstructPCA(const double eigen_threshold, int FirstPCAdpar, int LastPCAdpar) {
// ********************************************
  if(AdaptiveHandler) {
    MACH3LOG_ERROR("Adaption has been enabled and now trying to enable PCA. Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  PCAObj = std::make_unique<PCAHandler>();
  //Check whether first and last pcadpar are set and if not just PCA everything
  if(FirstPCAdpar == -999 || LastPCAdpar == -999){
    if(FirstPCAdpar == -999 && LastPCAdpar == -999){
      FirstPCAdpar = 0;
      LastPCAdpar = covMatrix->GetNrows()-1;
    }
    else{
      MACH3LOG_ERROR("You must either leave FirstPCAdpar and LastPCAdpar at -999 or set them both to something");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }

  PCAObj->ConstructPCA(covMatrix, FirstPCAdpar, LastPCAdpar, eigen_threshold, _fNumPar);
  PCAObj->SetupPointers(&_fCurrVal, &_fPropVal);
  // Make a note that we have now done PCA
  pca = true;
}

// ********************************************
void ParameterHandlerBase::Init(const std::string& name, const std::string& file) {
// ********************************************
  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  TFile *infile = new TFile(file.c_str(), "READ");
  if (infile->IsZombie()) {
    MACH3LOG_ERROR("Could not open input covariance ROOT file {} !!!", file);
    MACH3LOG_ERROR("Was about to retrieve matrix with name {}", name);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  TMatrixDSym *CovMat = static_cast<TMatrixDSym*>(infile->Get(name.c_str()));

  if (!CovMat) {
    MACH3LOG_ERROR("Could not find covariance matrix name {} in file {}", name, file);
    MACH3LOG_ERROR("Are you really sure {} exists in the file?", name);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  PrintLength = 35;

  const int nThreads = M3::GetNThreads();
  //KS: set Random numbers for each thread so each thread has different seed
  //or for one thread if without MULTITHREAD
  random_number.reserve(nThreads);
  for (int iThread = 0; iThread < nThreads; iThread++) {
    random_number.emplace_back(std::make_unique<TRandom3>(0));
  }
  // Not using adaptive by default
  use_adaptive = false;
  // Set the covariance matrix
  _fNumPar = CovMat->GetNrows();

  InvertCovMatrix.resize(_fNumPar, std::vector<double>(_fNumPar, 0.0));
  throwMatrixCholDecomp = new double*[_fNumPar]();
  // Set the defaults to true
  for(int i = 0; i < _fNumPar; i++) {
    throwMatrixCholDecomp[i] = new double[_fNumPar]();
    for (int j = 0; j < _fNumPar; j++) {
      throwMatrixCholDecomp[i][j] = 0.;
    }
  }
  SetName(name);
  MakePosDef(CovMat);
  SetCovMatrix(CovMat);
  if (_fNumPar <= 0) {
    MACH3LOG_CRITICAL("Covariance matrix {} has {} entries!", GetName(), _fNumPar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  ReserveMemory(_fNumPar);

  infile->Close();

  MACH3LOG_INFO("Created covariance matrix named: {}", GetName());
  MACH3LOG_INFO("from file: {}", file);
  delete infile;
}

// ********************************************
// ETA An init function for the YAML constructor
// All you really need from the YAML file is the number of Systematics
void ParameterHandlerBase::Init(const std::vector<std::string>& YAMLFile) {
// ********************************************

  std::map<std::pair<int, int>, std::unique_ptr<TMatrixDSym>> ThrowSubMatrixOverrides;
  int running_num_file_pars = 0;

  _fYAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    YAML::Node YAMLDocTemp = M3OpenConfig(YAMLFile[i]);

    if (YAMLDocTemp["ThrowMatrixOverride"]) { // LP: this allows us to put in
                                              // proposal matrix overrides per
                                              // parameter-containing file, add
                                              // the block diagonal proposal
                                              // matrix to a list and overwrite
                                              // the throw matrix after set up.
      auto filename =
          YAMLDocTemp["ThrowMatrixOverride"]["file"].as<std::string>();
      TFile *submatrix_file = M3::Open(filename, "OPEN", __FILE__, __LINE__);

      auto matrixname =
          YAMLDocTemp["ThrowMatrixOverride"]["matrix"].as<std::string>();
      std::unique_ptr<TMatrixDSym> submatrix{
          submatrix_file->Get<TMatrixDSym>(matrixname.c_str())};
      if (!submatrix) {
        MACH3LOG_CRITICAL("Covariance matrix {} doesn't exist in file: {}",
                          matrixname, filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      auto numrows = submatrix->GetNrows();
      // LP: the -1 here is because we specify the last index for consistency
      // with PCAHandler, not the first index after the end as is more common
      // throughout computer science...
      ThrowSubMatrixOverrides[{running_num_file_pars,
                               running_num_file_pars + (numrows - 1)}] =
          std::move(submatrix);

      // LP: check names by default, but have option to disable check if you
      // know what you're doing
      if (!bool(YAMLDocTemp["ThrowMatrixOverride"]["check_names"]) ||
          YAMLDocTemp["ThrowMatrixOverride"]["check_names"].as<bool>()) {
        auto nametree = submatrix_file->Get<TTree>("param_names");
        if (!nametree) {
          MACH3LOG_CRITICAL("TTree param_names doesn't exist in file: {}. Set "
                            "ThrowMatrixOverride: {{ check_names: False }} to "
                            "disable this check.",
                            filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        std::string *param_name = nullptr;
        nametree->SetBranchAddress("name", &param_name);

        if (nametree->GetEntries() != int(YAMLDocTemp["Systematics"].size())) {
          MACH3LOG_CRITICAL("TTree param_names in file: {} has {} entries, but "
                            "the corresponding yaml file only declares {} "
                            "parameters. Set ThrowMatrixOverride: {{ "
                            "check_names: False }} to disable this check.",
                            filename, nametree->GetEntries(),
                            YAMLDocTemp["Systematics"].size());
          throw MaCh3Exception(__FILE__, __LINE__);
        }

        int pit = 0;
        for (const auto &param : YAMLDocTemp["Systematics"]) {
          nametree->GetEntry(pit++);
          auto yaml_pname = Get<std::string>(
              param["Systematic"]["Names"]["FancyName"], __FILE__, __LINE__);
          if ((*param_name) != yaml_pname) {
            MACH3LOG_CRITICAL(
                "TTree param_names in file: {} at entry {} has parameter {}, "
                "but "
                "the corresponding yaml parameter is named {}. Set "
                "ThrowMatrixOverride: {{ "
                "check_names: False }} to disable this check.",
                filename, pit, (*param_name), yaml_pname);
            throw MaCh3Exception(__FILE__, __LINE__);
          }
        }
      }
      submatrix_file->Close();
    }

    for (const auto& item : YAMLDocTemp["Systematics"]) {
      _fYAMLDoc["Systematics"].push_back(item);
      running_num_file_pars++;
    }
  }

  const int nThreads = M3::GetNThreads();
  //KS: set Random numbers for each thread so each thread has different seed
  //or for one thread if without MULTITHREAD
  random_number.reserve(nThreads);
  for (int iThread = 0; iThread < nThreads; iThread++) {
    random_number.emplace_back(std::make_unique<TRandom3>(0));
  }
  PrintLength = 35;

  // Set the covariance matrix
  _fNumPar = int(_fYAMLDoc["Systematics"].size());

  use_adaptive = false;

  InvertCovMatrix.resize(_fNumPar, std::vector<double>(_fNumPar, 0.0));
  throwMatrixCholDecomp = new double*[_fNumPar]();
  for(int i = 0; i < _fNumPar; i++) {
    throwMatrixCholDecomp[i] = new double[_fNumPar]();
    for (int j = 0; j < _fNumPar; j++) {
      throwMatrixCholDecomp[i][j] = 0.;
    }
  }
  ReserveMemory(_fNumPar);

  TMatrixDSym* _fCovMatrix = new TMatrixDSym(_fNumPar);
  int i = 0;
  std::vector<std::map<std::string,double>> Correlations(_fNumPar);
  std::map<std::string, int> CorrNamesMap;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPar for Names,
  //PreFitValues etc etc.

  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fFancyNames[i] = Get<std::string>(param["Systematic"]["Names"]["FancyName"], __FILE__ , __LINE__);
    _fPreFitValue[i] = Get<double>(param["Systematic"]["ParameterValues"]["PreFitValue"], __FILE__ , __LINE__);
    _fIndivStepScale[i] = Get<double>(param["Systematic"]["StepScale"]["MCMC"], __FILE__ , __LINE__);
    _fError[i] = Get<double>(param["Systematic"]["Error"], __FILE__ , __LINE__);
    _fSampleNames[i] = GetFromManager<std::vector<std::string>>(param["Systematic"]["SampleNames"], {}, __FILE__, __LINE__);
    if(_fError[i] <= 0) {
      MACH3LOG_ERROR("Error for param {}({}) is negative and equal to {}", _fFancyNames[i], i, _fError[i]);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    //ETA - a bit of a fudge but works
    auto TempBoundsVec = GetBounds(param["Systematic"]["ParameterBounds"]);
    _fLowBound[i] = TempBoundsVec[0];
    _fUpBound[i] = TempBoundsVec[1];

    //ETA - now for parameters which are optional and have default values
    _fFlatPrior[i] = GetFromManager<bool>(param["Systematic"]["FlatPrior"], false, __FILE__ , __LINE__);

    // Allow to fix param, this setting should be used only for params which are permanently fixed like baseline, please use global config for fixing param more flexibly
    if(GetFromManager<bool>(param["Systematic"]["FixParam"], false, __FILE__ , __LINE__)) {
      ToggleFixParameter(_fFancyNames[i]);
    }

    if(param["Systematic"]["SpecialProposal"]) {
      EnableSpecialProposal(param["Systematic"]["SpecialProposal"], i);
    }

    //Fill the map to get the correlations later as well
    CorrNamesMap[param["Systematic"]["Names"]["FancyName"].as<std::string>()]=i;

    //Also loop through the correlations
    if(param["Systematic"]["Correlations"]) {
      for(unsigned int Corr_i = 0; Corr_i < param["Systematic"]["Correlations"].size(); ++Corr_i){
        for (YAML::const_iterator it = param["Systematic"]["Correlations"][Corr_i].begin(); it!=param["Systematic"]["Correlations"][Corr_i].end();++it) {
          Correlations[i][it->first.as<std::string>()] = it->second.as<double>();
        }
      }
    }
    i++;
  } // end loop over para
  if(i != _fNumPar) {
    MACH3LOG_CRITICAL("Inconsistent number of params in Yaml  {} vs {}, this indicate wrong syntax", i, i, _fNumPar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  // ETA Now that we've been through all systematic let's fill the covmatrix
  //This makes the root TCov from YAML
  for(int j = 0; j < _fNumPar; j++) {
    (*_fCovMatrix)(j, j) = _fError[j]*_fError[j];
    //Get the map of parameter name to correlation from the Correlations object
    for (auto const& pair : Correlations[j]) {
      auto const& key = pair.first;
      auto const& val = pair.second;
      int index = -1;
      //If you found the parameter name then get the index
      if (CorrNamesMap.find(key) != CorrNamesMap.end()) {
        index = CorrNamesMap[key];
      } else {
        MACH3LOG_ERROR("Parameter {} not in list! Check your spelling?", key);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      double Corr1 = val;
      double Corr2 = 0;
      if(Correlations[index].find(_fFancyNames[j]) != Correlations[index].end()) {
        Corr2 = Correlations[index][_fFancyNames[j]];
        //Do they agree to better than float precision?
        if(std::abs(Corr2 - Corr1) > FLT_EPSILON) {
          MACH3LOG_ERROR("Correlations are not equal between {} and {}", _fFancyNames[j], key);
          MACH3LOG_ERROR("Got : {} and {}", Corr2, Corr1);
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
      } else {
        MACH3LOG_ERROR("Correlation does not appear reciprocally between {} and {}", _fFancyNames[j], key);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      (*_fCovMatrix)(j, index)= (*_fCovMatrix)(index, j) = Corr1*_fError[j]*_fError[index];
    }
  }

  //Now make positive definite
  MakePosDef(_fCovMatrix);
  SetCovMatrix(_fCovMatrix);

  if (_fNumPar <= 0) {
    MACH3LOG_ERROR("ParameterHandler object has {} systematics!", _fNumPar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  for(auto const & matovr : ThrowSubMatrixOverrides){
    SetSubThrowMatrix(matovr.first.first, matovr.first.second, *matovr.second);
  }

  Tunes = std::make_unique<ParameterTunes>(_fYAMLDoc["Systematics"]);

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for(const auto &file : YAMLFile){
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total", _fNumPar);
  MACH3LOG_INFO("----------------");
}

// ********************************************
void ParameterHandlerBase::EnableSpecialProposal(const YAML::Node& param, const int Index){
// ********************************************
  doSpecialStepProposal = true;

  bool CircEnabled = false;
  bool FlipEnabled = false;

  if (param["CircularBounds"]) {
    CircEnabled = true;
  }

  if (param["FlipParameter"]) {
    FlipEnabled = true;
  }

  if (!CircEnabled && !FlipEnabled) {
    MACH3LOG_ERROR("None of Special Proposal were enabled even though param {}, has SpecialProposal entry in Yaml", GetParFancyName(Index));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (CircEnabled) {
    CircularBoundsIndex.push_back(Index);
    CircularBoundsValues.push_back(Get<std::pair<double, double>>(param["CircularBounds"], __FILE__, __LINE__));
    MACH3LOG_INFO("Enabling CircularBounds for parameter {} with range [{}, {}]",
                  GetParFancyName(Index),
                  CircularBoundsValues.back().first,
                  CircularBoundsValues.back().second);
    // KS: Make sure circular bounds are within physical bounds. If we are outside of physics bound MCMC will never explore such phase space region
    if (CircularBoundsValues.back().first < _fLowBound.at(Index) || CircularBoundsValues.back().second > _fUpBound.at(Index)) {
      MACH3LOG_ERROR("Circular bounds [{}, {}] for parameter {} exceed physical bounds [{}, {}]",
                     CircularBoundsValues.back().first, CircularBoundsValues.back().second,
                     GetParFancyName(Index),
                     _fLowBound.at(Index), _fUpBound.at(Index));
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  if (FlipEnabled) {
    FlipParameterIndex.push_back(Index);
    FlipParameterPoint.push_back(Get<double>(param["FlipParameter"], __FILE__, __LINE__));
    MACH3LOG_INFO("Enabling Flipping for parameter {} with value {}",
                  GetParFancyName(Index),
                  FlipParameterPoint.back());
  }

  if (CircEnabled && FlipEnabled) {
    if (FlipParameterPoint.back() < CircularBoundsValues.back().first || FlipParameterPoint.back() > CircularBoundsValues.back().second) {
      MACH3LOG_ERROR("FlipParameter value {} for parameter {} is outside the CircularBounds [{}, {}]",
                     FlipParameterPoint.back(), GetParFancyName(Index), CircularBoundsValues.back().first, CircularBoundsValues.back().second);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    const double low = CircularBoundsValues.back().first;
    const double high = CircularBoundsValues.back().second;

    // Sanity check: ensure flipping any x in [low, high] keeps the result in [low, high]
    const double flipped_low = 2 * FlipParameterPoint.back() - low;
    const double flipped_high = 2 * FlipParameterPoint.back() - high;
    const double min_flip = std::min(flipped_low, flipped_high);
    const double max_flip = std::max(flipped_low, flipped_high);

    if (min_flip < low || max_flip > high) {
      MACH3LOG_ERROR("Flipping about point {} for parameter {} would leave circular bounds [{}, {}]",
                     FlipParameterPoint.back(), GetParFancyName(Index), low, high);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// ********************************************
// Set the covariance matrix for this class
void ParameterHandlerBase::SetCovMatrix(TMatrixDSym *cov) {
// ********************************************
  if (cov == nullptr) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to {}", __func__ );
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  covMatrix = cov;

  invCovMatrix = static_cast<TMatrixDSym *>(cov->Clone());
  invCovMatrix->Invert();
  //KS: ROOT has bad memory management, using standard double means we can decrease most operation by factor 2 simply due to cache hits
  for (int i = 0; i < _fNumPar; i++)
  {
    for (int j = 0; j < _fNumPar; ++j)
    {
      InvertCovMatrix[i][j] = (*invCovMatrix)(i,j);
    }
  }

  SetThrowMatrix(cov);
}
// ********************************************
void ParameterHandlerBase::ReserveMemory(const int SizeVec) {
// ********************************************
  _fNames = std::vector<std::string>(SizeVec);
  _fFancyNames = std::vector<std::string>(SizeVec);
  _fPreFitValue = std::vector<double>(SizeVec);
  _fError = std::vector<double>(SizeVec);
  _fCurrVal = std::vector<double>(SizeVec);
  _fPropVal = std::vector<double>(SizeVec);
  _fLowBound = std::vector<double>(SizeVec);
  _fUpBound = std::vector<double>(SizeVec);
  _fFlatPrior = std::vector<bool>(SizeVec);
  _fIndivStepScale = std::vector<double>(SizeVec);
  _fSampleNames = std::vector<std::vector<std::string>>(_fNumPar);

  corr_throw = new double[SizeVec];
  // set random parameter vector (for correlated steps)
  randParams = new double[SizeVec];

  // Set the defaults to true
  for(int i = 0; i < SizeVec; i++) {
    _fPreFitValue.at(i) = 1.;
    _fError.at(i) = 1.;
    _fCurrVal.at(i) = 0.;
    _fPropVal.at(i) = 0.;
    _fLowBound.at(i) = -999.99;
    _fUpBound.at(i) = 999.99;
    _fFlatPrior.at(i) = false;
    _fIndivStepScale.at(i) = 1.;
    corr_throw[i] = 0.0;
    randParams[i] = 0.0;
  }

  _fGlobalStepScale = 1.0;
}

// ********************************************
// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void ParameterHandlerBase::SetPar(const int i , const double val) {
// ********************************************
  MACH3LOG_DEBUG("Over-riding {}: _fPropVal ({}), _fCurrVal ({}), _fPreFitValue ({}) to ({})",
                 GetParFancyName(i), _fPropVal[i], _fCurrVal[i], _fPreFitValue[i], val);

  _fPropVal[i] = val;
  _fCurrVal[i] = val;
  _fPreFitValue[i] = val;

  // Transfer the parameter values to the PCA basis
  if (pca) PCAObj->TransferToPCA();
}

// ********************************************
std::vector<double> ParameterHandlerBase::GetProposed() const {
// ********************************************
  std::vector<double> props(_fNumPar);
  for (int i = 0; i < _fNumPar; ++i) props[i] = _fPropVal[i];
  return props;
}

// *************************************
// Throw the parameters according to the covariance matrix
// This shouldn't be used in MCMC code ase it can break Detailed Balance;
void ParameterHandlerBase::ThrowParameters() {
// *************************************
  // First draw new randParams
  Randomize();

  M3::MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

  // KS: We use PCA very rarely on top PCA functionality isn't implemented for this function.
  // Use __builtin_expect to give compiler a hint which option is more likely, which should help
  // with better optimisation. This isn't critical but more to have example
  if (__builtin_expect(!pca, 1)) {
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < _fNumPar; ++i) {
      // Check if parameter is fixed first: if so don't randomly throw
      if (IsParameterFixed(i)) continue;

      _fPropVal[i] = _fPreFitValue[i] + corr_throw[i];

      int throws = 0;
      // Try again if we the initial parameter proposal falls outside of the range of the parameter
      while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
        randParams[i] = random_number[M3::GetThreadIndex()]->Gaus(0, 1);
        const double corr_throw_single = M3::MatrixVectorMultiSingle(throwMatrixCholDecomp, randParams, _fNumPar, i);
        _fPropVal[i] = _fPreFitValue[i] + corr_throw_single;
        if (throws > 10000)
        {
          //KS: Since we are multithreading there is danger that those messages
          //will be all over the place, small price to pay for faster code
          MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws, i);
          MACH3LOG_WARN("Matrix: {}", matrixName);
          MACH3LOG_WARN("Param: {}", _fNames[i]);
          MACH3LOG_WARN("Setting _fPropVal:  {} to {}", _fPropVal[i], _fPreFitValue[i]);
          MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
          _fPropVal[i] = _fPreFitValue[i];
          //throw MaCh3Exception(__FILE__ , __LINE__ );
        }
        throws++;
      }
      _fCurrVal[i] = _fPropVal[i];
    }
  }
  else
  {
    PCAObj->ThrowParameters(random_number, throwMatrixCholDecomp,
                            randParams, corr_throw,
                            _fPreFitValue, _fLowBound, _fUpBound, _fNumPar);
  } // end if pca

  // KS: At the end once we are happy with proposal do special proposal
  SpecialStepProposal();
}

// *************************************
// Throw each parameter within their 1 sigma range
// Used to start the chain in different states
void ParameterHandlerBase::RandomConfiguration() {
// *************************************
  // Have the 1 sigma for each parameter in each covariance class, sweet!
  // Don't want to change the prior array because that's what determines our likelihood
  // Want to change the _fPropVal, _fCurrVal, _fPreFitValue
  // _fPreFitValue and the others will already be set
  for (int i = 0; i < _fNumPar; ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (IsParameterFixed(i)) continue;
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    const double paramrange = _fUpBound[i] - _fLowBound[i];
    const double sigma = sqrt((*covMatrix)(i,i));
    double throwrange = sigma;
    if (paramrange < sigma) throwrange = paramrange;

    _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
    // Try again if we the initial parameter proposal falls outside of the range of the parameter
    int throws = 0;
    while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
      if (throws > 1000) {
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws, i);
        MACH3LOG_WARN("Matrix: {}", matrixName);
        MACH3LOG_WARN("Param: {}", _fNames[i]);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
      throws++;
    }
    MACH3LOG_INFO("Setting current step in {} param {} = {} from {}", matrixName, i, _fPropVal[i], _fCurrVal[i]);
    _fCurrVal[i] = _fPropVal[i];
  }
  if (pca) PCAObj->TransferToPCA();

  // KS: At the end once we are happy with proposal do special proposal
  SpecialStepProposal();
}

// *************************************
// Set a single parameter
void ParameterHandlerBase::SetSingleParameter(const int parNo, const double parVal) {
// *************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  MACH3LOG_DEBUG("Setting {} (parameter {}) to {})", GetParFancyName(parNo),  parNo, parVal);
  if (pca) PCAObj->TransferToPCA();
}

// ********************************************
void ParameterHandlerBase::SetParCurrProp(const int parNo, const double parVal) {
// ********************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  MACH3LOG_DEBUG("Setting {} (parameter {}) to {})", GetParFancyName(parNo),  parNo, parVal);
  if (pca) PCAObj->TransferToPCA();
}

// ************************************************
// Propose a step for the set of systematics parameters this covariance class holds
void ParameterHandlerBase::ProposeStep() {
// ************************************************
  // Make the random numbers for the step proposal
  Randomize();
  CorrelateSteps();

  // KS: According to Dr Wallace we update using previous not proposed step
  // this way we do special proposal after adaptive after.
  // This way we can shortcut and skip rest of proposal
  if(!doSpecialStepProposal) return;

  SpecialStepProposal();
}

// ************************************************
void ParameterHandlerBase::SpecialStepProposal() {
// ************************************************
  /// @warning KS: Following Asher comment we do "Step->Circular Bounds->Flip"

  // HW It should now automatically set dcp to be with [-pi, pi]
  for (size_t i = 0; i < CircularBoundsIndex.size(); ++i) {
    const int index = CircularBoundsIndex[i];
    if(!IsParameterFixed(index))
      CircularParBounds(index, CircularBoundsValues[i].first, CircularBoundsValues[i].second);
  }

  // Okay now we've done the standard steps, we can add in our nice flips hierarchy flip first
  for (size_t i = 0; i < FlipParameterIndex.size(); ++i) {
    const int index = FlipParameterIndex[i];
    if(!IsParameterFixed(index))
      FlipParameterValue(FlipParameterIndex[i], FlipParameterPoint[i]);
  }
}

// ************************************************
// "Randomize" the parameters in the covariance class for the proposed step
// Used the proposal kernel and the current parameter value to set proposed step
// Also get a new random number for the randParams
void ParameterHandlerBase::Randomize() _noexcept_ {
// ************************************************
  if (!pca) {
    //KS: By multithreading here we gain at least factor 2 with 8 threads with ND only fit
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < _fNumPar; ++i) {
      // If parameter isn't fixed
      if (!IsParameterFixed(i) > 0.0) {
        randParams[i] = random_number[M3::GetThreadIndex()]->Gaus(0, 1);
        // If parameter IS fixed
      } else {
        randParams[i] = 0.0;
      }
    } // end for
  // If we're in the PCA basis we instead throw parameters there (only _fNumParPCA parameter)
  } else {
    // Scale the random parameters by the sqrt of eigen values for the throw
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < PCAObj->GetNumberPCAedParameters(); ++i)
    {
      // If parameter IS fixed or out of bounds
      if (PCAObj->IsParameterFixedPCA(i)) {
        randParams[i] = 0.0;
      } else {
        randParams[i] = random_number[M3::GetThreadIndex()]->Gaus(0,1);
      }
    }
  }
}

// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void ParameterHandlerBase::CorrelateSteps() _noexcept_ {
// ************************************************
  //KS: Using custom function compared to ROOT one with 8 threads we have almost factor 2 performance increase, by replacing TMatrix with just double we increase it even more
  M3::MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

  // If not doing PCA
  if (!pca) {
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < _fNumPar; ++i) {
      if (!IsParameterFixed(i) > 0.) {
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*_fGlobalStepScale*_fIndivStepScale[i];
      }
    }
    // If doing PCA throw uncorrelated in PCA basis (orthogonal basis by definition)
  } else {
    PCAObj->CorrelateSteps(_fIndivStepScale, _fGlobalStepScale, randParams, corr_throw);
  }
}
// ********************************************
// Update so that current step becomes the previously proposed step
void ParameterHandlerBase::AcceptStep() _noexcept_ {
// ********************************************
  if (!pca) {
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < _fNumPar; ++i) {
      // Update state so that current state is proposed state
      _fCurrVal[i] = _fPropVal[i];
    }
  } else {
    PCAObj->AcceptStep();
  }

  if (AdaptiveHandler) {
    AdaptiveHandler->IncrementAcceptedSteps();
  }
}

// *************************************
//HW: This method is a tad hacky but modular arithmetic gives me a headache.
void ParameterHandlerBase::CircularParBounds(const int index, const double LowBound, const double UpBound) {
// *************************************
  if(_fPropVal[index] > UpBound) {
    _fPropVal[index] = LowBound + std::fmod(_fPropVal[index] - UpBound, UpBound - LowBound);
  } else if (_fPropVal[index] < LowBound) {
    _fPropVal[index] = UpBound - std::fmod(LowBound - _fPropVal[index], UpBound - LowBound);
  }
}

// *************************************
void ParameterHandlerBase::FlipParameterValue(const int index, const double FlipPoint) {
// *************************************
  if(random_number[0]->Uniform() < 0.5) {
    _fPropVal[index] = 2 * FlipPoint - _fPropVal[index];
  }
}

// ********************************************
// Function to print the prior values
void ParameterHandlerBase::PrintNominal() const {
// ********************************************
  MACH3LOG_INFO("Prior values for {} ParameterHandler:", GetName());
  for (int i = 0; i < _fNumPar; i++) {
    MACH3LOG_INFO("    {}   {} ", GetParFancyName(i), GetParInit(i));
  }
}

// ********************************************
// Function to print the prior, current and proposed values
void ParameterHandlerBase::PrintNominalCurrProp() const {
// ********************************************
  MACH3LOG_INFO("Printing parameters for {}", GetName());
  // Dump out the PCA parameters too
  if (pca) {
    PCAObj->Print();
  }
  MACH3LOG_INFO("{:<30} {:<10} {:<10} {:<10}", "Name", "Prior", "Current", "Proposed");
  for (int i = 0; i < _fNumPar; ++i) {
    MACH3LOG_INFO("{:<30} {:<10.2f} {:<10.2f} {:<10.2f}", GetParFancyName(i), _fPreFitValue[i], _fCurrVal[i], _fPropVal[i]);
  }
}

// ********************************************
// Get the likelihood in the case where we want to include priors on the parameters
// _fFlatPrior stores if we want to evaluate the likelihood for the given parameter
//                    true = don't evaluate likelihood (so run without a prior)
//                    false = evaluate likelihood (so run with a prior)
double ParameterHandlerBase::CalcLikelihood() const _noexcept_ {
// ********************************************
  double logL = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:logL)
  #endif
  for(int i = 0; i < _fNumPar; ++i) {
    if(_fFlatPrior[i]){
      //HW: Flat prior, no need to calculate anything
      continue;
    }
    // KS: Precalculate Diff once per "i" without doing this for every "j"
    const double Diff = _fPropVal[i] - _fPreFitValue[i];
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (int j = 0; j <= i; ++j) {
      if (!_fFlatPrior[j]) {
        //KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
        double scale = (i != j) ? 1. : 0.5;
        logL += scale * Diff * (_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];
      }
    }
  }
  return logL;
}

// ********************************************
int ParameterHandlerBase::CheckBounds() const _noexcept_ {
// ********************************************
  int NOutside = 0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:NOutside)
  #endif
  for (int i = 0; i < _fNumPar; ++i){
    if(_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]){
      NOutside++;
    }
  }
  return NOutside;
}

// ********************************************
double ParameterHandlerBase::GetLikelihood() {
// ********************************************
  // Default behaviour is to reject negative values + do std llh calculation
  const int NOutside = CheckBounds();

  if(NOutside > 0) return NOutside*M3::_LARGE_LOGL_;

  return CalcLikelihood();
}

// ********************************************
// Sets the proposed parameters to the prior values
void ParameterHandlerBase::SetParameters(const std::vector<double>& pars) {
// ********************************************
  // If empty, set the proposed to prior
  if (pars.empty()) {
    // For xsec this means setting to the prior (because prior is the prior)
    for (int i = 0; i < _fNumPar; i++) {
      _fPropVal[i] = _fPreFitValue[i];
    }
    // If not empty, set the parameters to the specified
  } else {
    if (pars.size() != size_t(_fNumPar)) {
      MACH3LOG_ERROR("Parameter arrays of incompatible size! Not changing parameters! {} has size {} but was expecting {}", matrixName, pars.size(), _fNumPar);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    int parsSize = int(pars.size());
    for (int i = 0; i < parsSize; i++) {
      //Make sure that you are actually passing a number to set the parameter to
      if(std::isnan(pars[i])) {
        MACH3LOG_ERROR("Trying to set parameter value to a nan for parameter {} in matrix {}. This will not go well!", GetParName(i), matrixName);
        throw MaCh3Exception(__FILE__ , __LINE__ );
      } else {
        _fPropVal[i] = pars[i];
      }
    }
  }
  // And if pca make the transfer
  if (pca) {
    PCAObj->TransferToPCA();
    PCAObj->TransferToParam();
  }
}

// ********************************************
void ParameterHandlerBase::SetBranches(TTree &tree, bool SaveProposal) {
// ********************************************
  // loop over parameters and set a branch
  for (int i = 0; i < _fNumPar; ++i) {
    tree.Branch(_fNames[i].c_str(), &_fCurrVal[i], Form("%s/D", _fNames[i].c_str()));
  }
  // When running PCA, also save PCA parameters
  if (pca) {
    PCAObj->SetBranches(tree, SaveProposal, _fNames);
  }
  if(SaveProposal)
  {
    // loop over parameters and set a branch
    for (int i = 0; i < _fNumPar; ++i) {
      tree.Branch(Form("%s_Prop", _fNames[i].c_str()), &_fPropVal[i], Form("%s_Prop/D", _fNames[i].c_str()));
    }
  }
  if(use_adaptive && AdaptiveHandler->GetUseRobbinsMonro()){
    tree.Branch(Form("GlobalStepScale_%s", GetName().c_str()), &_fGlobalStepScale, Form("GlobalStepScale_%s/D", GetName().c_str()));
  }
}

// ********************************************
void ParameterHandlerBase::SetStepScale(const double scale, const bool verbose) {
// ********************************************
  if(scale <= 0) {
    MACH3LOG_ERROR("You are trying so set StepScale to 0 or negative this will not work");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if(verbose){
    MACH3LOG_INFO("{} setStepScale() = {}", GetName(), scale);
    const double SuggestedScale = 2.38/std::sqrt(_fNumPar);
    if(std::fabs(scale - SuggestedScale)/SuggestedScale > 1) {
      MACH3LOG_WARN("Defined Global StepScale is {}, while suggested suggested {}", scale, SuggestedScale);
    }
  }
  _fGlobalStepScale = scale;
}

// ********************************************
int ParameterHandlerBase::GetParIndex(const std::string& name) const {
// ********************************************
  int Index = M3::_BAD_INT_;
  for (int i = 0; i <_fNumPar; ++i) {
    if(name == _fFancyNames[i]) {
      Index = i;
      break;
    }
  }
  return Index;
}

// ********************************************
void ParameterHandlerBase::SetFixAllParameters() {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  for (int i = 0; i < _fNumPar; ++i)
    if(!IsParameterFixed(i)) ToggleFixParameter(i);
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const int i) {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  if(!IsParameterFixed(i)) ToggleFixParameter(i);
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const std::string& name) {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  if(!IsParameterFixed(name)) ToggleFixParameter(name);
}

// ********************************************
void ParameterHandlerBase::SetFreeAllParameters() {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  for (int i = 0; i < _fNumPar; ++i)
    if(IsParameterFixed(i)) ToggleFixParameter(i);
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const int i) {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  if(IsParameterFixed(i)) ToggleFixParameter(i);
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const std::string& name) {
// ********************************************
  // Check if the parameter is fixed and if not, toggle fix it
  if(IsParameterFixed(name)) ToggleFixParameter(name);
}

// ********************************************
void ParameterHandlerBase::ToggleFixAllParameters() {
// ********************************************
  // toggle fix/free all parameters
  if(!pca) for (int i = 0; i < _fNumPar; i++) ToggleFixParameter(i);
  else PCAObj->ToggleFixAllParameters(_fNames);
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const int i) {
// ********************************************
  if(!pca) {
    if (i > _fNumPar) {
      MACH3LOG_ERROR("Can't {} for parameter {} because size of covariance ={}", __func__, i, _fNumPar);
      MACH3LOG_ERROR("Fix this in your config file please!");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    } else {
      _fError[i] *= -1.0;
      if(IsParameterFixed(i)) MACH3LOG_INFO("Setting {}(parameter {}) to fixed at {}", GetParFancyName(i), i, _fCurrVal[i]);
      else MACH3LOG_INFO("Setting {}(parameter {}) free", GetParFancyName(i), i);
    }
    if( (_fCurrVal[i] > _fUpBound[i] || _fCurrVal[i] < _fLowBound[i]) && IsParameterFixed(i) ) {
      MACH3LOG_ERROR("Parameter {} (index {}) is fixed at {}, which is outside of its bounds [{}, {}]", GetParFancyName(i), i, _fCurrVal[i], _fLowBound[i], _fUpBound[i]);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  } else {
    PCAObj->ToggleFixParameter(i, _fNames);
  }
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const std::string& name) {
// ********************************************
  const int Index = GetParIndex(name);
  if(Index != M3::_BAD_INT_) {
    ToggleFixParameter(Index);
    return;
  }

  MACH3LOG_WARN("I couldn't find parameter with name {}, therefore will not fix it", name);
}

// ********************************************
bool ParameterHandlerBase::IsParameterFixed(const std::string& name) const {
// ********************************************
  const int Index = GetParIndex(name);
  if(Index != M3::_BAD_INT_) {
    return IsParameterFixed(Index);
  }

  MACH3LOG_WARN("I couldn't find parameter with name {}, therefore don't know if it fixed", name);
  return false;
}

// ********************************************
void ParameterHandlerBase::SetFlatPrior(const int i, const bool eL) {
// ********************************************
  if (i > _fNumPar) {
    MACH3LOG_INFO("Can't {} for Cov={}/Param={} because size of Covariance = {}", __func__, GetName(), i, _fNumPar);
    MACH3LOG_ERROR("Fix this in your config file please!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  } else {
    if(eL){
      MACH3LOG_INFO("Setting {} (parameter {}) to flat prior", GetParName(i), i);
    }
    else{
      // HW :: This is useful
      MACH3LOG_INFO("Setting {} (parameter {}) to non-flat prior", GetParName(i), i);
    }
    _fFlatPrior[i] = eL;
  }
}

// ********************************************
void ParameterHandlerBase::SetIndivStepScale(const std::vector<double>& stepscale) {
// ********************************************
  if (int(stepscale.size()) != _fNumPar)
  {
    MACH3LOG_WARN("Stepscale vector not equal to number of parameters. Quitting..");
    MACH3LOG_WARN("Size of argument vector: {}", stepscale.size());
    MACH3LOG_WARN("Expected size: {}", _fNumPar);
    return;
  }

  for (int iParam = 0 ; iParam < _fNumPar; iParam++) {
    _fIndivStepScale[iParam] = stepscale[iParam];
  }
  PrintIndivStepScale();
}

// ********************************************
void ParameterHandlerBase::PrintIndivStepScale() const {
// ********************************************
  MACH3LOG_INFO("============================================================");
  MACH3LOG_INFO("{:<{}} | {:<11}", "Parameter:", PrintLength, "Step scale:");
  for (int iParam = 0; iParam < _fNumPar; iParam++) {
    MACH3LOG_INFO("{:<{}} | {:<11}", _fFancyNames[iParam].c_str(), PrintLength, _fIndivStepScale[iParam]);
  }
  MACH3LOG_INFO("============================================================");
}

// ********************************************
//Makes sure that matrix is positive-definite by adding a small number to on-diagonal elements
void ParameterHandlerBase::MakePosDef(TMatrixDSym *cov) {
// ********************************************
  if(cov == nullptr){
    cov = &*covMatrix;
    MACH3LOG_WARN("Passed nullptr to cov matrix in {}", matrixName);
  }

  M3::MakeMatrixPosDef(cov);
}

// ********************************************
void ParameterHandlerBase::ResetIndivStepScale() {
// ********************************************
  std::vector<double> stepScales(_fNumPar);
  for (int i = 0; i <_fNumPar; i++) {
    stepScales[i] = 1.;
  }
  _fGlobalStepScale = 1.0;
  SetIndivStepScale(stepScales);
}

void ParameterHandlerBase::SetIndivStepScaleForSkippedAdaptParams() {
  if (!param_skip_adapt_flags.size()) {
    MACH3LOG_ERROR("Parameter skip adapt flags not set, cannot set individual step scales for skipped parameters.");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  // HH: Cancel the effect of global step scale change for parameters that are not adapting
  for (int i = 0; i <_fNumPar; i++) {
    if (param_skip_adapt_flags[i]) {
      _fIndivStepScale[i] = _fIndivStepScaleInitial[i] * _fGlobalStepScaleInitial / _fGlobalStepScale;
    }
  }
  MACH3LOG_INFO("Updating individual step scales for non-adapting parameters to cancel global step scale change.");
  MACH3LOG_INFO("Global step scale initial: {}, current: {}", _fGlobalStepScaleInitial, _fGlobalStepScale);
  PrintIndivStepScale();
}

// ********************************************
// HW: Code for throwing from separate throw matrix, needs to be set after init to ensure pos-def
void ParameterHandlerBase::SetThrowMatrix(TMatrixDSym *cov){
// ********************************************

   if (cov == nullptr) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to {}", __func__);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if (covMatrix->GetNrows() != cov->GetNrows()) {
    MACH3LOG_ERROR("Matrix given for throw Matrix is not the same size as the covariance matrix stored in object!");
    MACH3LOG_ERROR("Stored covariance matrix size: {}", covMatrix->GetNrows());
    MACH3LOG_ERROR("Given matrix size: {}", cov->GetNrows());
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  throwMatrix = static_cast<TMatrixDSym*>(cov->Clone());
  if(use_adaptive && AdaptiveHandler->AdaptionUpdate()) MakeClosestPosDef(throwMatrix);
  else MakePosDef(throwMatrix);

  auto throwMatrix_CholDecomp = M3::GetCholeskyDecomposedMatrix(*throwMatrix, matrixName);

  //KS: ROOT has bad memory management, using standard double means we can decrease most operation by factor 2 simply due to cache hits
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  for (int i = 0; i < _fNumPar; ++i)
  {
    for (int j = 0; j < _fNumPar; ++j)
    {
      throwMatrixCholDecomp[i][j] = throwMatrix_CholDecomp[i][j];
    }
  }
}

void ParameterHandlerBase::SetSubThrowMatrix(int first_index, int last_index,
                                             TMatrixDSym const &subcov) {

  if ((last_index - first_index) >= subcov.GetNrows()) {
    MACH3LOG_ERROR("Trying to SetSubThrowMatrix into range: ({},{}) with a "
                   "submatrix with only {} rows {}",
                   first_index, last_index, subcov.GetNrows(), __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TMatrixDSym *current_ThrowMatrix =
      static_cast<TMatrixDSym *>(throwMatrix->Clone());
  for (int i = first_index; i <= last_index; ++i) {
    for (int j = first_index; j <= last_index; ++j) {
      current_ThrowMatrix->operator()(i, j) =
          subcov(i - first_index, j - first_index);
    }
  }

  SetThrowMatrix(current_ThrowMatrix);
  delete current_ThrowMatrix;
}

// ********************************************
void ParameterHandlerBase::UpdateThrowMatrix(TMatrixDSym *cov){
// ********************************************
  delete throwMatrix;
  throwMatrix = nullptr;
  SetThrowMatrix(cov);
}

// ********************************************
// HW : Here be adaption
void ParameterHandlerBase::InitialiseAdaption(const YAML::Node& adapt_manager){
// ********************************************
  if(PCAObj){
    MACH3LOG_ERROR("PCA has been enabled and now trying to enable Adaption. Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if(AdaptiveHandler){
    MACH3LOG_ERROR("Adaptive Handler has already been initialise can't do it again so skipping.");
    return;
  }
  AdaptiveHandler = std::make_unique<adaptive_mcmc::AdaptiveMCMCHandler>();

  // HH: Backing up _fIndivStepScale and _fGlobalStepScale before adaption
  _fIndivStepScaleInitial = _fIndivStepScale;
  _fGlobalStepScaleInitial = _fGlobalStepScale;

  // HH: adding these here because they will be used to set the individual step scales for non-adapting parameters
  std::vector<std::string> params_to_skip = GetFromManager<std::vector<std::string>>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ParametersToSkip"], {});
  // Build a list of skip flags
  param_skip_adapt_flags.resize(_fNumPar, false);
  for (int i = 0; i <_fNumPar; ++i) {
    for (const auto& name : params_to_skip) {
      if(name == _fFancyNames[i]) {
        param_skip_adapt_flags[i] = true;
        break;
      }
    }
  }
  // HH: Loop over correlations to check if any skipped parameter is correlated with adapted one
  // We don't want to change one parameter while keeping the other fixed as this would
  // lead to weird penalty terms in the prior after adapting
  double max_correlation = 0.01; // Define a threshold for significant correlation above which we throw an error
  for (int i = 0; i < _fNumPar; ++i) {
    for (int j = 0; j <= i; ++j) {
      // The symmetry should have been checked during the Init phase
      if(param_skip_adapt_flags[i] && !param_skip_adapt_flags[j]) {
        double corr = (*covMatrix)(i,j)/std::sqrt((*covMatrix)(i,i)*(*covMatrix)(j,j));
        if(std::fabs(corr) > max_correlation) {
          MACH3LOG_ERROR("Correlation between skipped parameter {} ({}) and non-skipped parameter {} ({}) is {:.6e}, above the allowed threshold of {:.6e}.",
                         i, _fFancyNames[i], j, _fFancyNames[j], corr, max_correlation);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }
  // Now we read the general settings [these SHOULD be common across all matrices!]
  bool success = AdaptiveHandler->InitFromConfig(adapt_manager, matrixName,
    &_fFancyNames, &_fCurrVal, &_fError,
    &param_skip_adapt_flags, throwMatrix, _fGlobalStepScaleInitial
  );
  if (success) {
    AdaptiveHandler->Print();
  }
  else {
    MACH3LOG_INFO("Not using adaptive MCMC for {}. Checking external matrix options...", matrixName);
  }

  // HH: Adjusting the external matrix reading logic such that you can not do adaptive
  // and still read an external matrix
  // Logic:
  // if read external matrix:
  //    set throw matrix regardless of adaptive or not
  // else:
  //    if adaptive:
  //       create new adaptive matrix from scratch
  //    else:
  //       do nothing
  if(GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["UseExternalMatrix"], false, __FILE__ , __LINE__)) {
    // Finally, we accept that we want to read the matrix from a file!
    auto external_file_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMatrixFileName"], "", __FILE__ , __LINE__);
    auto external_matrix_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMatrixName"], "", __FILE__ , __LINE__);
    auto external_mean_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMeansName"], "", __FILE__ , __LINE__);

    AdaptiveHandler->SetThrowMatrixFromFile(external_file_name, external_matrix_name, external_mean_name, use_adaptive);
    SetThrowMatrix(AdaptiveHandler->GetAdaptiveCovariance());

    ResetIndivStepScale();
    // HH: Set individual step scales for non-adapting parameters to the default individual step scales
    // global step scale should be 1 so no need to adjust for that
    SetIndivStepScaleForSkippedAdaptParams();

    MACH3LOG_INFO("Successfully Set External Throw Matrix Stored in {}", external_file_name);
  } else {
    MACH3LOG_INFO("Not using external matrix for {}", matrixName);
    if (!success) return; // Not adaptive either so nothing to do
    MACH3LOG_INFO("Initialising adaption from scratch");
    // If we don't have a covariance matrix to start from for adaptive tune we need to make one!
    use_adaptive = true;
    AdaptiveHandler->CheckMatrixValidityForAdaption(GetCovMatrix());
    AdaptiveHandler->CreateNewAdaptiveCovariance();
    return;
  }
}

// ********************************************
// Truely adaptive MCMC!
void ParameterHandlerBase::UpdateAdaptiveCovariance(){
// ********************************************
  // Updates adaptive matrix
  // First we update the total means

  // Skip this if we're at a large number of steps
  if(AdaptiveHandler->SkipAdaption()) {
    AdaptiveHandler->IncrementNSteps();
    return;
  }

  /// Need to adjust the scale every step
  if(AdaptiveHandler->GetUseRobbinsMonro()){
    bool verbose=false;
    #ifdef DEBUG
    verbose=true;
    #endif
    AdaptiveHandler->UpdateRobbinsMonroScale();
    SetStepScale(AdaptiveHandler->GetAdaptionScale(), verbose);
    SetIndivStepScaleForSkippedAdaptParams();
  }

  // Call main adaption function
  AdaptiveHandler->UpdateAdaptiveCovariance();

  // Set scales to 1 * optimal scale
  if(AdaptiveHandler->IndivStepScaleAdapt()) {
    ResetIndivStepScale();
    SetStepScale(AdaptiveHandler->GetAdaptionScale());
    SetIndivStepScaleForSkippedAdaptParams();
  }

  if(AdaptiveHandler->UpdateMatrixAdapt()) {
    TMatrixDSym* update_matrix = static_cast<TMatrixDSym*>(AdaptiveHandler->GetAdaptiveCovariance()->Clone());
    UpdateThrowMatrix(update_matrix); //Now we update and continue!
    //Also Save the adaptive to file
    AdaptiveHandler->SaveAdaptiveToFile(AdaptiveHandler->GetOutFileName(), GetName());
  }

  AdaptiveHandler->IncrementNSteps();
}

// ********************************************
//HW: Finds closest possible positive definite matrix in Frobenius Norm ||.||_frob
// Where ||X||_frob=sqrt[sum_ij(x_ij^2)] (basically just turns an n,n matrix into vector in n^2 space
// then does Euclidean norm)
void ParameterHandlerBase::MakeClosestPosDef(TMatrixDSym *cov) {
// ********************************************
  // Want to get cov' = (cov_sym+cov_polar)/2
  // cov_sym=(cov+cov^T)/2
  // cov_polar-> SVD cov to cov=USV^T then cov_polar=VSV^T

  //Get frob norm of cov
  //  Double_t cov_norm=cov->E2Norm();

  TMatrixDSym* cov_trans = cov;
  cov_trans->T();
  TMatrixDSym cov_sym = 0.5*(*cov+*cov_trans); //If cov is symmetric does nothing, otherwise just ensures symmetry

  //Do SVD to get polar form
  TDecompSVD cov_sym_svd=TDecompSVD(cov_sym);
  if(!cov_sym_svd.Decompose()){
    MACH3LOG_WARN("Cannot do SVD on input matrix, trying MakePosDef() first!");
    MakePosDef(&cov_sym);
  }

  TMatrixD cov_sym_v = cov_sym_svd.GetV();
  TMatrixD cov_sym_vt = cov_sym_v;
  cov_sym_vt.T();
  //SVD returns as vector (grrr) so need to get into matrix form for multiplying!
  TVectorD cov_sym_sigvect = cov_sym_svd.GetSig();

  const Int_t nCols = cov_sym_v.GetNcols(); //square so only need rows hence lack of cols
  TMatrixDSym cov_sym_sig(nCols);
  TMatrixDDiag cov_sym_sig_diag(cov_sym_sig);
  cov_sym_sig_diag=cov_sym_sigvect;

  //Can finally get H=VSV
  TMatrixDSym cov_sym_polar = cov_sym_sig.SimilarityT(cov_sym_vt);//V*S*V^T (this took forver to find!)

  //Now we can construct closest approximater Ahat=0.5*(B+H)
  TMatrixDSym cov_closest_approx  = 0.5*(cov_sym+cov_sym_polar);//Not fully sure why this is even needed since symmetric B -> U=V
  //Get norm of transformed
  //  Double_t approx_norm=cov_closest_approx.E2Norm();
  //MACH3LOG_INFO("Initial Norm: {:.6f} | Norm after transformation: {:.6f} | Ratio: {:.6f}", cov_norm, approx_norm, cov_norm / approx_norm);

  *cov = cov_closest_approx;
  //Now can just add a makeposdef!
  MakePosDef(cov);
}

// ********************************************
// KS: Convert covariance matrix to correlation matrix and return TH2D which can be used for fancy plotting
TH2D* ParameterHandlerBase::GetCorrelationMatrix() {
// ********************************************
  TH2D* hMatrix = new TH2D(GetName().c_str(), GetName().c_str(), _fNumPar, 0.0, _fNumPar, _fNumPar, 0.0, _fNumPar);
  hMatrix->SetDirectory(nullptr);
  for(int i = 0; i < _fNumPar; i++)
  {
    hMatrix->SetBinContent(i+1, i+1, 1.);
    hMatrix->GetXaxis()->SetBinLabel(i+1, GetParFancyName(i).c_str());
    hMatrix->GetYaxis()->SetBinLabel(i+1, GetParFancyName(i).c_str());
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < _fNumPar; i++)
  {
    for(int j = 0; j <= i; j++)
    {
      const double Corr = (*covMatrix)(i,j) / ( GetDiagonalError(i) * GetDiagonalError(j));
      hMatrix->SetBinContent(i+1, j+1, Corr);
      hMatrix->SetBinContent(j+1, i+1, Corr);
    }
  }
  return hMatrix;
}

// ********************************************
// KS: After step scale, prefit etc. value were modified save this modified config.
void ParameterHandlerBase::SaveUpdatedMatrixConfig() {
// ********************************************
  if (!_fYAMLDoc)
  {
    MACH3LOG_CRITICAL("Yaml node hasn't been initialised for matrix {}, something is not right", matrixName);
    MACH3LOG_CRITICAL("I am not throwing error but should be investigated");
    return;
  }

  YAML::Node copyNode = _fYAMLDoc;
  int i = 0;

  for (YAML::Node param : copyNode["Systematics"])
  {
    //KS: Feel free to update it, if you need updated prefit value etc
    param["Systematic"]["StepScale"]["MCMC"] = MaCh3Utils::FormatDouble(_fIndivStepScale[i], 4);
    i++;
  }
  // Save the modified node to a file
  std::ofstream fout("Modified_Matrix.yaml");
  fout << copyNode;
  fout.close();
}

// ********************************************
bool ParameterHandlerBase::AppliesToSample(const int SystIndex, const std::string& SampleName) const {
// ********************************************
  // Empty means apply to all
  if (_fSampleNames[SystIndex].size() == 0) return true;

  // Make a copy and to lower case to not be case sensitive
  std::string SampleNameCopy = SampleName;
  std::transform(SampleNameCopy.begin(), SampleNameCopy.end(), SampleNameCopy.begin(), ::tolower);

  // Check for unsupported wildcards in SampleNameCopy
  if (SampleNameCopy.find('*') != std::string::npos) {
    MACH3LOG_ERROR("Wildcards ('*') are not supported in sample name: '{}'", SampleName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  bool Applies = false;

  for (size_t i = 0; i < _fSampleNames[SystIndex].size(); i++) {
    // Convert to low case to not be case sensitive
    std::string pattern = _fSampleNames[SystIndex][i];
    std::transform(pattern.begin(), pattern.end(), pattern.begin(), ::tolower);

    // Replace '*' in the pattern with '.*' for regex matching
    std::string regexPattern = "^" + std::regex_replace(pattern, std::regex("\\*"), ".*") + "$";
    try {
      std::regex regex(regexPattern);
      if (std::regex_match(SampleNameCopy, regex)) {
        Applies = true;
        break;
      }
    } catch (const std::regex_error& e) {
      // Handle regex error (for invalid patterns)
      MACH3LOG_ERROR("Regex error: {}", e.what());
    }
  }
  return Applies;
}

// ********************************************
// Set proposed parameter values vector to be base on tune values
void ParameterHandlerBase::SetTune(const std::string& TuneName) {
// ********************************************
  if(Tunes == nullptr) {
    MACH3LOG_ERROR("Tunes haven't been initialised, which are being loaded from YAML, have you used some deprecated constructor");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  auto Values = Tunes->GetTune(TuneName);

  SetParameters(Values);
}

// *************************************
/// @brief Matches branches in a TTree to parameters in a systematic handler.
///
/// @param PosteriorFile Pointer to the ROOT TTree from MaCh3 fit.
/// @param Systematic Pointer to the systematic parameter handler.
/// @param[out] BranchValues Vector to store the values of the branches (resized inside).
/// @param[out] BranchNames Vector to store the names of the branches (resized inside).
///
/// @throws MaCh3Exception if any parameter branch is uninitialized.
void ParameterHandlerBase::MatchMaCh3OutputBranches(TTree *PosteriorFile,
                              std::vector<double>& BranchValues,
                              std::vector<std::string>& BranchNames) {
// *************************************
  BranchValues.resize(GetNumParams());
  BranchNames.resize(GetNumParams());

  for (int i = 0; i < GetNumParams(); ++i) {
    BranchNames[i] = GetParName(i);
    if (!PosteriorFile->GetBranch(BranchNames[i].c_str())) {
      MACH3LOG_ERROR("Branch '{}' does not exist in the TTree!", BranchNames[i]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    PosteriorFile->SetBranchStatus(BranchNames[i].c_str(), true);
    PosteriorFile->SetBranchAddress(BranchNames[i].c_str(), &BranchValues[i]);
  }
}

