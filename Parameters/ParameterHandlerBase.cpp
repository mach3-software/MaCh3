#include "Parameters/ParameterHandlerBase.h"
#include "regex"

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string name, std::string file, double threshold, int FirstPCA, int LastPCA) : inputFile(file), pca(false),
eigen_threshold(threshold), FirstPCAdpar(FirstPCA), LastPCAdpar(LastPCA) {
// ********************************************
  MACH3LOG_DEBUG("Constructing instance of ParameterHandler");
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
  if (pca) ConstructPCA();
}
// ********************************************
ParameterHandlerBase::ParameterHandlerBase(const std::vector<std::string>& YAMLFile, std::string name, double threshold, int FirstPCA, int LastPCA) : inputFile(YAMLFile[0].c_str()), matrixName(name), pca(true), eigen_threshold(threshold), FirstPCAdpar(FirstPCA), LastPCAdpar(LastPCA) {
// ********************************************
  MACH3LOG_INFO("Constructing instance of ParameterHandler using");
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
  if (pca) ConstructPCA();
}

// ********************************************
//Destructor
ParameterHandlerBase::~ParameterHandlerBase(){
// ********************************************
  delete[] randParams;
  delete[] corr_throw;

  if (covMatrix != nullptr) delete covMatrix;
  if (invCovMatrix != nullptr) delete invCovMatrix;
  if (throwMatrix_CholDecomp != nullptr) delete throwMatrix_CholDecomp;
  if (throwMatrix != nullptr) delete throwMatrix;
  for(int i = 0; i < _fNumPar; i++) {
    delete[] throwMatrixCholDecomp[i];
  }
  delete[] throwMatrixCholDecomp;
}

// ********************************************
void ParameterHandlerBase::ConstructPCA() {
// ********************************************
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
void ParameterHandlerBase::Init(std::string name, std::string file) {
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
  _fYAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    YAML::Node YAMLDocTemp = M3OpenConfig(YAMLFile[i]);
    for (const auto& item : YAMLDocTemp["Systematics"]) {
      _fYAMLDoc["Systematics"].push_back(item);
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
    _fGenerated[i] = Get<double>(param["Systematic"]["ParameterValues"]["Generated"], __FILE__ , __LINE__);
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

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for(const auto &file : YAMLFile){
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total", _fNumPar);
  MACH3LOG_INFO("----------------");
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
  _fGenerated = std::vector<double>(SizeVec);
  _fPreFitValue = std::vector<double>(SizeVec);
  _fError = std::vector<double>(SizeVec);
  _fCurrVal = std::vector<double>(SizeVec);
  _fPropVal = std::vector<double>(SizeVec);
  _fLowBound = std::vector<double>(SizeVec);
  _fUpBound = std::vector<double>(SizeVec);
  _fFlatPrior = std::vector<bool>(SizeVec);
  _fIndivStepScale = std::vector<double>(SizeVec);
  _fSampleNames = std::vector<std::vector<std::string>>(_fNumPar);

  corr_throw = new double[SizeVec]();
  // set random parameter vector (for correlated steps)
  randParams = new double[SizeVec];

  // Set the defaults to true
  for(int i = 0; i < SizeVec; i++) {
    _fGenerated.at(i) = 1.;
    _fPreFitValue.at(i) = 1.;
    _fError.at(i) = 1.;
    _fCurrVal.at(i) = 0.;
    _fPropVal.at(i) = 0.;
    _fLowBound.at(i) = -999.99;
    _fUpBound.at(i) = 999.99;
    _fFlatPrior.at(i) = false;
    _fIndivStepScale.at(i) = 1.;
    corr_throw[i] = 0.0;
  }

  _fGlobalStepScale = 1.0;
}

// ********************************************
// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void ParameterHandlerBase::SetPar(int i , double val) {
// ********************************************
  MACH3LOG_INFO("Over-riding {}: ", GetParName(i));
  MACH3LOG_INFO("_fPropVal ({}), _fCurrVal ({}), _fPreFitValue ({}) to ({})", _fPropVal[i], _fCurrVal[i], _fPreFitValue[i], val);

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

  // KS: We use PCA very rarely on top PCA functionality isn't implemented for this function.
  // Use __builtin_expect to give compiler a hint which option is more likely, which should help
  // with better optimisation. This isn't critical but more to have example
  if (__builtin_expect(!pca, 1)) {
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

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
        const double corr_throw_single = MatrixVectorMultiSingle(throwMatrixCholDecomp, randParams, _fNumPar, i);
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
    MACH3LOG_CRITICAL("Hold on, you are trying to run Prior Predictive Code with PCA, which is wrong");
    MACH3LOG_CRITICAL("Sorry I have to kill you, I mean your job");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
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
}

// *************************************
// Set a single parameter
void ParameterHandlerBase::SetSingleParameter(const int parNo, const double parVal) {
// *************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  MACH3LOG_DEBUG("Setting {} (parameter {}) to {})", GetParName(parNo),  parNo, parVal);
  if (pca) PCAObj->TransferToPCA();
}
// ********************************************
void ParameterHandlerBase::SetParCurrProp(const int parNo, const double parVal) {
// ********************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  MACH3LOG_DEBUG("Setting {} (parameter {}) to {})", GetParName(parNo),  parNo, parVal);
  if (pca) PCAObj->TransferToPCA();
}

// ************************************************
// Propose a step for the set of systematics parameters this covariance class holds
void ParameterHandlerBase::ProposeStep() {
// ************************************************
  // Make the random numbers for the step proposal
  Randomize();
  CorrelateSteps();
  if(use_adaptive) UpdateAdaptiveCovariance();
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
      if (_fError[i] > 0.0) {
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
    for (int i = 0; i < _fNumPar; ++i)
    {
      if (PCAObj->IsParameterFixedPCA(i) > 0. && i < PCAObj->GetNumberPCAedParameters())
      {
        randParams[i] = random_number[M3::GetThreadIndex()]->Gaus(0,1);
      } else { // If parameter IS fixed or out of bounds
        randParams[i] = 0.0;
      }
    }
  }
}

// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void ParameterHandlerBase::CorrelateSteps() _noexcept_ {
// ************************************************
  //KS: Using custom function compared to ROOT one with 8 threads we have almost factor 2 performance increase, by replacing TMatrix with just double we increase it even more
  MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

  // If not doing PCA
  if (!pca) {
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < _fNumPar; ++i) {
      if (_fError[i] > 0.) {
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
}

// ********************************************
// Throw the proposed parameter by mag sigma
// Should really just have the user specify this throw by having argument double
void ParameterHandlerBase::ThrowParProp(const double mag) {
// ********************************************
  Randomize();
  if (!pca) {
    // Make the correlated throw
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);
    // Number of sigmas we throw
    for (int i = 0; i < _fNumPar; i++) {
      if (_fError[i] > 0.)
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*mag;
    }
  } else {
    PCAObj->ThrowParProp(mag, randParams);
  }
}
// ********************************************
// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void ParameterHandlerBase::ThrowParCurr(const double mag) {
// ********************************************
  Randomize();
  if (!pca) {
    // Get the correlated throw vector
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);
    // The number of sigmas to throw
    // Should probably have this as a default parameter input to the function instead
    for (int i = 0; i < _fNumPar; i++) {
      if (_fError[i] > 0.){
        _fCurrVal[i] = corr_throw[i]*mag;
      }
    }
  } else {
    PCAObj->ThrowParCurr(mag, randParams);
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
// fParEvalLikelihood stores if we want to evaluate the likelihood for the given parameter
//                    true = evaluate likelihood (so run with a prior)
//                    false = don't evaluate likelihood (so run without a prior)
double ParameterHandlerBase::CalcLikelihood() const _noexcept_ {
// ********************************************
  double logL = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:logL)
  #endif
  for(int i = 0; i < _fNumPar; ++i){
    if(_fFlatPrior[i]){
      //HW: Flat prior, no need to calculate anything
      continue;
    }

    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (int j = 0; j <= i; ++j) {
      if (!_fFlatPrior[j]) {
        //KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
        double scale = (i != j) ? 1. : 0.5;
        logL += scale * (_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];
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
}

// ********************************************
void ParameterHandlerBase::SetStepScale(const double scale) {
// ********************************************
  if(scale <= 0) {
    MACH3LOG_ERROR("You are trying so set StepScale to 0 or negative this will not work");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("{} setStepScale() = {}", GetName(), scale);
  const double SuggestedScale = 2.38*2.38/_fNumPar;
  if(std::fabs(scale - SuggestedScale)/SuggestedScale > 1) {
    MACH3LOG_WARN("Defined Global StepScale is {}, while suggested suggested {}", scale, SuggestedScale);
  }

  _fGlobalStepScale = scale;
}

// ********************************************
void ParameterHandlerBase::ToggleFixAllParameters() {
// ********************************************
  // fix or unfix all parameters by multiplying by -1
  if(!pca) {
    for (int i = 0; i < _fNumPar; i++) _fError[i] *= -1.0;
  } else{
    PCAObj->ToggleFixAllParameters();
  }
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
      MACH3LOG_INFO("Setting {}(parameter {}) to fixed at {}", GetParFancyName(i), i, _fCurrVal[i]);
    }
  } else {
    PCAObj->ToggleFixParameter(i, _fNames);
  }
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const std::string& name) {
// ********************************************
  for (int i = 0; i <_fNumPar; ++i) {
    if(name == _fFancyNames[i]) {
      ToggleFixParameter(i);
      return;
    }
  }
  MACH3LOG_WARN("I couldn't find parameter with name {}, therefore will not fix it", name);
}

// ********************************************
bool ParameterHandlerBase::IsParameterFixed(const std::string& name) const {
// ********************************************
  for (int i = 0; i <_fNumPar; ++i) {
    if(name == _fFancyNames[i]) {
      return IsParameterFixed(i);
    }
  }
  MACH3LOG_WARN("I couldn't find parameter with name {}, therefore will not fix it", name);
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
//KS: Custom function to perform multiplication of matrix and vector with multithreading
void ParameterHandlerBase::MatrixVectorMulti(double* _restrict_ VecMulti, double** _restrict_ matrix, const double* _restrict_ vector, const int n) const {
// ********************************************
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; ++i)
  {
    double result = 0.0;
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (int j = 0; j < n; ++j)
    {
      result += matrix[i][j]*vector[j];
    }
    VecMulti[i] = result;
  }
}

// ********************************************
double ParameterHandlerBase::MatrixVectorMultiSingle(double** _restrict_ matrix, const double* _restrict_ vector, const int Length, const int i) const {
// ********************************************
  double Element = 0.0;
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int j = 0; j < Length; ++j) {
    Element += matrix[i][j]*vector[j];
  }
  return Element;
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

  //DB Save original warning state and then increase it in this function to suppress 'matrix not positive definite' messages
  //Means we no longer need to overload
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  //DB Loop 1000 times adding 1e-9 which tops out at 1e-6 shift on the diagonal before throwing error
  int MaxAttempts = 1e5;
  int iAttempt = 0;
  bool CanDecomp = false;
  TDecompChol chdcmp;
  
  for (iAttempt = 0; iAttempt < MaxAttempts; iAttempt++) {
    chdcmp = TDecompChol(*cov);
    if (chdcmp.Decompose()) {
      CanDecomp = true;
      break;
    } else {
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif
      for (int iVar = 0 ; iVar < _fNumPar; iVar++) {
        (*cov)(iVar,iVar) += pow(10, -9);
      }
    }
  }

  if (!CanDecomp) {
    MACH3LOG_ERROR("Tried {} times to shift diagonal but still can not decompose the matrix", MaxAttempts);
    MACH3LOG_ERROR("This indicates that something is wrong with the input matrix");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if(!use_adaptive || AdaptiveHandler->total_steps < 2) {
    MACH3LOG_INFO("Had to shift diagonal {} time(s) to allow the covariance matrix to be decomposed", iAttempt);
  }
  //DB Resetting warning level
  gErrorIgnoreLevel = originalErrorWarning;
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
  
  TDecompChol TDecompChol_throwMatrix(*throwMatrix);
  
  if(!TDecompChol_throwMatrix.Decompose()) {
    MACH3LOG_ERROR("Cholesky decomposition failed for {} trying to make positive definite", matrixName);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  throwMatrix_CholDecomp = new TMatrixD(TDecompChol_throwMatrix.GetU());
  throwMatrix_CholDecomp->T();

  //KS: ROOT has bad memory management, using standard double means we can decrease most operation by factor 2 simply due to cache hits
  #ifdef MULTITHREAD
  #pragma omp parallel for collapse(2)
  #endif
  for (int i = 0; i < _fNumPar; ++i)
  {
    for (int j = 0; j < _fNumPar; ++j)
    {
      throwMatrixCholDecomp[i][j] = (*throwMatrix_CholDecomp)(i,j);
    }
  }
}

// ********************************************
void ParameterHandlerBase::UpdateThrowMatrix(TMatrixDSym *cov){
// ********************************************
  delete throwMatrix;
  throwMatrix = nullptr;
  delete throwMatrix_CholDecomp;
  throwMatrix_CholDecomp = nullptr;
  SetThrowMatrix(cov);
}

// ********************************************
// HW : Here be adaption
void ParameterHandlerBase::InitialiseAdaption(const YAML::Node& adapt_manager){
// ********************************************
  AdaptiveHandler = std::make_unique<adaptive_mcmc::AdaptiveMCMCHandler>();
  // Now we read the general settings [these SHOULD be common across all matrices!]
  bool success = AdaptiveHandler->InitFromConfig(adapt_manager, matrixName, &_fCurrVal, &_fError);
  if(!success) return;
  AdaptiveHandler->Print();

  // Next let"s check for external matrices
  // We"re going to grab this info from the YAML manager
  if(!GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["UseExternalMatrix"], false, __FILE__ , __LINE__)) {
    MACH3LOG_WARN("Not using external matrix for {}, initialising adaption from scratch", matrixName);
    // If we don't have a covariance matrix to start from for adaptive tune we need to make one!
    use_adaptive = true;
    AdaptiveHandler->CreateNewAdaptiveCovariance();
    return;
  }

  // Finally, we accept that we want to read the matrix from a file!
  auto external_file_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMatrixFileName"], "", __FILE__ , __LINE__);
  auto external_matrix_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMatrixName"], "", __FILE__ , __LINE__);
  auto external_mean_name = GetFromManager<std::string>(adapt_manager["AdaptionOptions"]["Covariance"][matrixName]["ExternalMeansName"], "", __FILE__ , __LINE__);

  AdaptiveHandler->SetThrowMatrixFromFile(external_file_name, external_matrix_name, external_mean_name, use_adaptive);
  SetThrowMatrix(AdaptiveHandler->adaptive_covariance);

  MACH3LOG_INFO("Successfully Set External Throw Matrix Stored in {}", external_file_name);
}

// ********************************************
// Truely adaptive MCMC!
void ParameterHandlerBase::UpdateAdaptiveCovariance(){
// ********************************************
  // Updates adaptive matrix
  // First we update the total means

  // Skip this if we're at a large number of steps
  if(AdaptiveHandler->SkipAdaption()) {
    AdaptiveHandler->total_steps++;
    return;
  }

  // Call main adaption function
  AdaptiveHandler->UpdateAdaptiveCovariance();

  //This is likely going to be the slow bit!
  if(AdaptiveHandler->IndivStepScaleAdapt()) {
    ResetIndivStepScale();
  }

  if(AdaptiveHandler->UpdateMatrixAdapt()) {
    TMatrixDSym* update_matrix = static_cast<TMatrixDSym*>(AdaptiveHandler->adaptive_covariance->Clone());
    UpdateThrowMatrix(update_matrix); //Now we update and continue!
    //Also Save the adaptive to file
    AdaptiveHandler->SaveAdaptiveToFile(AdaptiveHandler->output_file_name,GetName());
  }

  AdaptiveHandler->total_steps++;
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
    MACH3LOG_ERROR("Cannot do SVD on input matrix!");
    throw MaCh3Exception(__FILE__ , __LINE__ );
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
std::vector<double> ParameterHandlerBase::GetNominalArray() {
// ********************************************
  std::vector<double> prior(_fNumPar);
  for (int i = 0; i < _fNumPar; ++i)
  {
    prior[i] = _fPreFitValue[i];
  }
  return prior;
}

// ********************************************
// KS: Convert covariance matrix to correlation matrix and return TH2D which can be used for fancy plotting
TH2D* ParameterHandlerBase::GetCorrelationMatrix() {
// ********************************************
  TH2D* hMatrix = new TH2D(GetName().c_str(), GetName().c_str(), _fNumPar, 0.0, _fNumPar, _fNumPar, 0.0, _fNumPar);

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
