#include "covarianceBase.h"
// ********************************************
covarianceBase::covarianceBase(const char *name, const char *file) : inputFile(std::string(file)), pca(false) {
// ********************************************
#ifdef MULTITHREAD
  int nThreads = omp_get_max_threads();
#else
  int nThreads = 1;
#endif
  //KS: set Random numbers for each thread so each thread has differnt seed
  //or for one thread if without MULTITHREAD 
  random_number = new TRandom3*[(const int)nThreads]();
  for (int iThread=0;iThread<nThreads;iThread++) {
    random_number[iThread] = new TRandom3(0);
  }

  std::cout << "Constructing instance of covarianceBase" << std::endl;
  init(name, file);
  FirstPCAdpar = -999;
  LastPCAdpar = -999;

  PrintLength = 35;
}

covarianceBase::covarianceBase(const char *YAMLFile) : inputFile(std::string(YAMLFile)), pca(false) {
// ********************************************
#ifdef MULTITHREAD
  int nThreads = omp_get_max_threads();
#else
  int nThreads = 1;
#endif
  //KS: set Random numbers for each thread so each thread has differnt seed
  //or for one thread if without MULTITHREAD 
  random_number = new TRandom3*[(const int)nThreads]();
  for (int iThread=0;iThread<nThreads;iThread++) {
    random_number[iThread] = new TRandom3(0);
  }

  std::cout << "Constructing instance of covarianceBase using " << YAMLFile << " as an input" << std::endl;
  init(YAMLFile);
  FirstPCAdpar = -999;
  LastPCAdpar = -999;

  PrintLength = 35;
}

covarianceBase::covarianceBase(const char *name, const char *file, int seed) : inputFile(std::string(file)), pca(false) {
#ifdef MULTITHREAD
  int nThreads = omp_get_max_threads();
  if(seed != 0)
  {
     std::cerr<<"You have set seed to "<<seed<<std::endl;
     std::cerr<<"And you are running with MULTITHREAD"<<std::endl;
     std::cerr<<"TRandom for each thread will have same seed"<<std::endl;
     std::cerr<<"This is fine if this was your intetnion"<<std::endl;
  }
  #else
  int nThreads = 1;
#endif
   //KS: set Random numbers for each thread so each thread has differnt seed
  //or for one thread if without MULTITHREAD 
  random_number = new TRandom3*[(const int)nThreads]();
  for (int iThread=0;iThread<nThreads;iThread++) {
    random_number[iThread] = new TRandom3(seed);
  }

  init(name, file);
  FirstPCAdpar = -999;
  LastPCAdpar = -999;

  PrintLength = 35;
}

covarianceBase::covarianceBase(const char *name, const char *file, int seed, double threshold, int firstpcapar, int lastpcapar) : inputFile(std::string(file)), pca(true), eigen_threshold(threshold), FirstPCAdpar(firstpcapar), LastPCAdpar(lastpcapar) {

  if (threshold < 0 || threshold >= 1) {
    std::cerr << "*** NOTE: " << name << " " << file << std::endl;
    std::cerr << "    Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work!" << std::endl;
    std::cerr << "    Please specify a number between 0 and 1" << std::endl;
    std::cerr << "    You specified: " << threshold << std::endl;
    std::cerr << "    Am instead calling the usual non-PCA constructor..." << std::endl;
    pca = false;
  }
#ifdef MULTITHREAD
  int nThreads = omp_get_max_threads();
  if(seed != 0)
  {
     std::cerr<<"You have set seed to "<<seed<<std::endl;
     std::cerr<<"And you are running with MULTITHREAD"<<std::endl;
     std::cerr<<"TRandom for each thread will have same seed"<<std::endl;
     std::cerr<<"This is fine if this was your intetnion"<<std::endl;
  }
#else
  int nThreads = 1;
#endif
  //KS: set Random numbers for each thread so each thread has differnt seed
  //or for one thread if without MULTITHREAD 
  random_number = new TRandom3*[(const int)nThreads]();
  for (int iThread=0;iThread<nThreads;iThread++) {
    random_number[iThread] = new TRandom3(seed);
  }
  std::cout << "Constructing instance of covarianceBase" << std::endl;
  init(name, file);
  // Call the innocent helper function
  if (pca) ConstructPCA();

   PrintLength = 35;
}

// ********************************************
//Desctructor
covarianceBase::~covarianceBase(){
// ********************************************
  /*delete[] fParInit;
  delete[] fParSigma;
  delete[] fParCurr;
  delete[] fParProp;
  delete[] fParLoLimit;
  delete[] fParHiLimit;
  delete[] fIndivStepScale;
  delete[] fParEvalLikelihood;
  */

  _fPreFitValue.clear();
  _fError.clear();
  _fCurrVal.clear();
  _fPropVal.clear();
  _fLowBound.clear();
  _fUpBound.clear();
  _fIndivStepScale.clear();
  _fFlatPrior.clear();

  delete[] randParams;

  if (covMatrix != NULL) delete covMatrix;
  if (invCovMatrix != NULL) delete invCovMatrix;
  if (throwMatrix_CholDecomp != NULL) delete throwMatrix_CholDecomp;

  for(int i = 0; i < size; i++) 
  {
    delete[] InvertCovMatrix[i];
    delete[] throwMatrixCholDecomp[i];
  }
  delete[] InvertCovMatrix;
  delete[] throwMatrixCholDecomp;
  
#ifdef MULTITHREAD
  int nThreads = omp_get_max_threads();
#else
  int nThreads = 1;
#endif
  for (int iThread=0;iThread < nThreads; iThread++)  delete random_number[iThread];
  delete[] random_number;
}

void covarianceBase::ConstructPCA() {

  // Check that covariance matrix exists
  if (covMatrix == NULL) {
    std::cerr << "Covariance matrix for " << matrixName << " has not yet been set" << std::endl;
    std::cerr << "Can not construct PCA until it is set" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  //Check whether first and last pcadpar are set and if not just PCA everything
  if(FirstPCAdpar == -999 || LastPCAdpar == -999){
    if(FirstPCAdpar == -999 && LastPCAdpar == -999){
      FirstPCAdpar = 0;
      LastPCAdpar = covMatrix->GetNrows()-1;
    }
    else{
      std::cerr << "You must either leave FirstPCAdpar and LastPCAdpar at -999 or set them both to something"<<std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
  }
  if(FirstPCAdpar>covMatrix->GetNrows()-1||LastPCAdpar>covMatrix->GetNrows()-1){
    std::cerr <<"FirstPCAdpar and LastPCAdpar are higher than the number of parameters"<<std::endl;
    std::cerr <<"first: "<<FirstPCAdpar<<" last: "<<LastPCAdpar<<std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  if(FirstPCAdpar<0||LastPCAdpar<0){
    std::cerr <<"FirstPCAdpar and LastPCAdpar are less than 0 but not default -999"<<std::endl;
    std::cerr <<"first: "<<FirstPCAdpar<<" last: "<<LastPCAdpar<<std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  std::cout<<"PCAing parameters "<<FirstPCAdpar<<" through "<<LastPCAdpar<<" inclusive"<<std::endl;

  int numunpcadpars = covMatrix->GetNrows()-(LastPCAdpar-FirstPCAdpar+1);

  TMatrixDSym submat(covMatrix->GetSub(FirstPCAdpar,LastPCAdpar,FirstPCAdpar,LastPCAdpar));

  // Calculate how many eigen values this threshold corresponds to
  TMatrixDSymEigen eigen(submat);
  eigen_values.ResizeTo(eigen.GetEigenValues());
  eigen_vectors.ResizeTo(eigen.GetEigenVectors());
  eigen_values = eigen.GetEigenValues();
  eigen_vectors = eigen.GetEigenVectors();
  double sum = 0;
  // Loop over eigen values and sum them up
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    sum += eigen_values(i);
  }
  nKeptPCApars=eigen_values.GetNrows();
  // Now go through again and see how many eigen values correspond to threshold
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    // Get the relative size of the eigen value
    double sig = eigen_values(i)/sum;
    // Check against the threshold
    if (sig < eigen_threshold) {
      nKeptPCApars = i;
      break;
    }
  }
  npars = numunpcadpars+nKeptPCApars;
  std::cout << "Threshold of " << eigen_threshold << " on eigen values relative sum of eigen value (" << sum << ") generates " << nKeptPCApars << " eigen vectors, plus we have "<<numunpcadpars<<" unpcad pars, for a total of "<<npars << std::endl;

  //DB Create array of correct size so eigen_values can be used in CorrelateSteps
  eigen_values_master = std::vector<double>(npars,1.0);
  for (int i = FirstPCAdpar; i < FirstPCAdpar+nKeptPCApars; ++i) {eigen_values_master[i] = eigen_values(i-FirstPCAdpar);}

  // Now construct the transfer matrices
  //These matrices will be as big as number of unPCAd pars plus number of eigenvalues kept
  TransferMat.ResizeTo(covMatrix->GetNrows(), npars);
  TransferMatT.ResizeTo(covMatrix->GetNrows(), npars);

  // Get a subset of the eigen vector matrix
  TMatrixD temp(eigen_vectors.GetSub(0, eigen_vectors.GetNrows()-1, 0, nKeptPCApars-1));
  
  //Make transfer matrix which is two blocks of identity with a block of the PCA transfer matrix in between
  TMatrixD temp2;
  temp2.ResizeTo(covMatrix->GetNrows(), npars);

  //First set the whole thing to 0
  for(int iRow=0;iRow<covMatrix->GetNrows();iRow++){
    for(int iCol=0;iCol<npars;iCol++){
      temp2[iRow][iCol]=0;
    }
  }
  //Set the first identity block
  if(FirstPCAdpar!=0){
    for(int iRow=0;iRow<FirstPCAdpar;iRow++){
      temp2[iRow][iRow]=1;
    }
  }

  //Set the transfer matrix block for the PCAd pars
  temp2.SetSub(FirstPCAdpar,FirstPCAdpar,temp);

  //Set the second identity block
  if(LastPCAdpar!=covMatrix->GetNrows()-1){
    for(int iRow=0;iRow<(covMatrix->GetNrows()-1)-LastPCAdpar;iRow++){
      temp2[LastPCAdpar+1+iRow][FirstPCAdpar+nKeptPCApars+iRow]=1;
    }
  }
   
  TransferMat = temp2;
  // Copy the contents
  TransferMatT = TransferMat;
  // And then transpose
  TransferMatT.T();

  // Make a note that we have now done PCA
  pca = true;

  // Resize the size of the random parameters for step proposal
  randParams = new double[size];

  // Make the PCA parameter arrays
  fParCurr_PCA.ResizeTo(npars);
  fParProp_PCA.ResizeTo(npars);
  
  //KS: make easy map so we could easily find un-decomposed parameters
  isDecomposed_PCA.resize(npars);
  fParSigma_PCA.resize(npars);
  for (int i = 0; i < npars; ++i)
  {
      fParSigma_PCA[i] = 1;
      isDecomposed_PCA[i] = -1;
  }
  for (int i = 0; i < FirstPCAdpar; ++i) isDecomposed_PCA[i] = i;
  
  for (int i = FirstPCAdpar+nKeptPCApars+1; i < npars; ++i) isDecomposed_PCA[i] = i+(size-npars);

  //KS: Let's dump all usefull matrices to properly validate PCA
#ifdef DEBUG_PCA
DebugPCA(sum, temp, submat);
#endif
}

void covarianceBase::init(const char *name, const char *file)
{
  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  TFile *infile = new TFile(file, "READ");
  if (infile->IsZombie()) {
    std::cerr << "ERROR: Could not open input covariance ROOT file " << file << " !!!" << std::endl;
    std::cerr << "Was about to retrieve matrix with name " << name << std::endl;
    throw;
  }

  // Should put in a 
  TMatrixDSym *covMatrix = (TMatrixDSym*)(infile->Get(name));
  if (covMatrix == NULL) {
    std::cerr << "Could not find covariance matrix name " << name << " in file " << file << std::endl;
    std::cerr << "Are you really sure " << name << " exists in the file?" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__  << std::endl;
    throw;
  } 

  // Not using adaptive by default
  use_adaptive=false;
  
  // Set the covariance matrix
  size = covMatrix->GetNrows();
  _fNumPar = size;
    
  InvertCovMatrix = new double*[size]();
  throwMatrixCholDecomp = new double*[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++)
  {
    InvertCovMatrix[i] = new double[size]();
    throwMatrixCholDecomp[i] = new double[size]();
    for (int j = 0; j < size; j++)
    {
        InvertCovMatrix[i][j] = 0.;
        throwMatrixCholDecomp[i][j] = 0.;
    }
  }

  setName(name);
  size = covMatrix->GetNrows();
  MakePosDef(covMatrix);
  setCovMatrix(covMatrix);

  if (size <= 0) {
    std::cerr << "Covariance matrix " << getName() << " has " << size << " entries!" << std::endl;
    throw;
  }
  npars = size;

  _fNames = std::vector<std::string>(size);
  _fPreFitValue = std::vector<double>(size);
  _fError = std::vector<double>(size);
  _fCurrVal = std::vector<double>(size);
  _fPropVal = std::vector<double>(size);
  _fLowBound = std::vector<double>(size);
  _fUpBound = std::vector<double>(size);
  _fFlatPrior = std::vector<bool>(size);
  _fIndivStepScale = std::vector<double>(size);

  corr_throw = new Double_t[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++) {
	_fPreFitValue.at(i) = 1.;
	_fError.at(i) = 0.;
	_fCurrVal.at(i) = 0.;
	_fPropVal.at(i) = 0.;
	_fLowBound.at(i) = -999.99;
	_fUpBound.at(i) = 999.99;
	_fFlatPrior.at(i) = true;
	_fIndivStepScale.at(i) = 1.;
	corr_throw[i] = 0.0;
  }

  // Set the logLs to very large so next step is accepted
  currLogL = __LARGE_LOGL__;
  propLogL = __LARGE_LOGL__;

  // Set random parameter vector (for correlated steps)
  randParams = new double[size];

  infile->Close();

  _fGlobalStepScale = 1.0;

  std::cout << "Created covariance matrix named: " << getName() << std::endl;
  std::cout << "from file: " << file << std::endl;

  delete infile;
}

// ETA
// An init function for the YAML constructor
// All you really need from the YAML file is the number of Systematics
// Then get all the info from the YAML file in the covarianceXsec::ParseYAML function
void covarianceBase::init(const char *YAMLFile)
{
  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  /*TFile *infile = new TFile(file, "READ");
  if (infile->IsZombie()) {
    std::cerr << "ERROR: Could not open input covariance ROOT file " << file << " !!!" << std::endl;
    std::cerr << "Was about to retrieve matrix with name " << name << std::endl;
    throw;
  }

  // Should put in a 
  TMatrixDSym *covMatrix = (TMatrixDSym*)(infile->Get(name));
  if (covMatrix == NULL) {
    std::cerr << "Could not find covariance matrix name " << name << " in file " << file << std::endl;
    std::cerr << "Are you really sure " << name << " exists in the file?" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__  << std::endl;
    throw;
  }
  */

  _fYAMLDoc = YAML::LoadFile(YAMLFile);

  // Not using adaptive by default
  use_adaptive=false;
  
  // Set the covariance matrix
  _fNumPar = _fYAMLDoc["Systematics"].size();
  size = _fNumPar;
    
  _fCovMatrix = new TMatrixDSym(size);
  InvertCovMatrix = new double*[size]();
  throwMatrixCholDecomp = new double*[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++)
  {
    InvertCovMatrix[i] = new double[size]();
    throwMatrixCholDecomp[i] = new double[size]();
    for (int j = 0; j < size; j++)
    {
        InvertCovMatrix[i][j] = 0.;
        throwMatrixCholDecomp[i][j] = 0.;
    }
  }

  setName(YAMLFile);
  //size = covMatrix->GetNrows();
  //MakePosDef(covMatrix);
  //setCovMatrix(covMatrix);

  if (size <= 0) {
    std::cerr << "Covariance matrix " << getName() << " has " << size << " entries!" << std::endl;
    throw;
  }
  npars = size;

  _fNames = std::vector<std::string>(size);
  _fPreFitValue = std::vector<double>(size);
  _fError = std::vector<double>(size);
  _fCurrVal = std::vector<double>(size);
  _fPropVal = std::vector<double>(size);
  _fLowBound = std::vector<double>(size);
  _fUpBound = std::vector<double>(size);
  _fFlatPrior = std::vector<bool>(size);
  _fIndivStepScale = std::vector<double>(size);

  corr_throw = new Double_t[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++) {
	_fPreFitValue.at(i) = 1.;
	_fError.at(i) = 0.;
	_fCurrVal.at(i) = 0.;
	_fPropVal.at(i) = 0.;
	_fLowBound.at(i) = -999.99;
	_fUpBound.at(i) = 999.99;
	_fFlatPrior.at(i) = true;
	_fIndivStepScale.at(i) = 1.;
	corr_throw[i] = 0.0;
  }

  // Set the logLs to very large so next step is accepted
  currLogL = __LARGE_LOGL__;
  propLogL = __LARGE_LOGL__;

  // Set random parameter vector (for correlated steps)
  randParams = new double[size];

  _fGlobalStepScale = 1.0;

  std::cout << "Created covariance matrix named: " << getName() << std::endl;
  std::cout << "from file: " << YAMLFile << std::endl;
}

void covarianceBase::init(TMatrixDSym* covMat) {

  size = covMat->GetNrows();
  InvertCovMatrix = new double*[size]();
  throwMatrixCholDecomp = new double*[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++) 
  {
    InvertCovMatrix[i] = new double[size]();
    throwMatrixCholDecomp[i] = new double[size]();
    for (int j = 0; j < size; j++) 
    {
        InvertCovMatrix[i][j] = 0.;  
        throwMatrixCholDecomp[i][j] = 0.;
    }
  }
  
  setCovMatrix(covMat);

  _fNames = std::vector<std::string>(size);
  _fPreFitValue = std::vector<double>(size);
  _fError = std::vector<double>(size);
  _fCurrVal = std::vector<double>(size);
  _fPropVal = std::vector<double>(size);
  _fLowBound = std::vector<double>(size);
  _fUpBound = std::vector<double>(size);
  _fFlatPrior = std::vector<bool>(size);
  _fIndivStepScale = std::vector<double>(size);

  corr_throw = new Double_t[size]();
  // Set the defaults to true
  for(int i = 0; i < size; i++) {
	_fPreFitValue.at(i) = 1.;
	_fError.at(i) = 0.;
	_fCurrVal.at(i) = 0.;
	_fPropVal.at(i) = 0.;
	_fLowBound.at(i) = -999.99;
	_fUpBound.at(i) = 999.99;
	_fFlatPrior.at(i) = true;
	_fIndivStepScale.at(i) = 1.;
	corr_throw[i] = 0.0;
  }

  // dont need these 2 i tihnk
  currLogL = __LARGE_LOGL__;
  propLogL = __LARGE_LOGL__;

  // set random parameter vector (for correlated steps)
  randParams = new double[size];

  _fGlobalStepScale = 1.0;

  std::cout << "Created covariance matrix named: " << getName() << std::endl;
}

// Set the covariance matrix for this class
void covarianceBase::setCovMatrix(TMatrixDSym *cov) {
  if (cov == NULL) {
    std::cerr << "Could not find covariance matrix you provided to setCovMatrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  covMatrix = cov;
  invCovMatrix = (TMatrixDSym*)cov->Clone();
  invCovMatrix->Invert();
  //KS: ROOT has bad memory managment, using standard double means we can decrease most operation by factor 2 simply due to cache hits
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; ++j)
    {
        InvertCovMatrix[i][j] = (*invCovMatrix)(i,j);
    }
  }

  setThrowMatrix(cov);
}

// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void covarianceBase::setPar(int i , double val) {

  std::cout << "Over-riding " << getParName(i) << ": " << std::endl;
  std::cout << "_fPropVal (" << _fPropVal[i];
  std::cout << "), _fCurrVal (" << _fCurrVal[i];
  std::cout << "), _fPreFitValue (" << _fPreFitValue[i]; 
  std::cout << ") to " << val << std::endl;

  _fPropVal[i] = val;
  _fCurrVal[i] = val;
  _fPreFitValue[i] = val;

  // Transfer the parameter values to the PCA basis
  if (pca) TransferToPCA();
}

// Transfer a parameter variation in the parameter basis to the eigen basis
void covarianceBase::TransferToPCA() {
  if (!pca) {
    std::cerr << "Can not transfer to PCA if PCA isn't enabled" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  // Make the temporary vectors
  TVectorD fParCurr_vec(size);
  TVectorD fParProp_vec(size);
  for (int i = 0; i < size; ++i) {
    fParCurr_vec(i) = _fCurrVal[i];
    fParProp_vec(i) = _fPropVal[i];
  }

  fParCurr_PCA = TransferMatT*fParCurr_vec;
  fParProp_PCA = TransferMatT*fParProp_vec;
}

// Transfer a parameter variation in the eigen basis to the parameter basis
void covarianceBase::TransferToParam() {
  if (!pca) {
    std::cerr << "Can not transfer to PCA if PCA isn't enabled" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Make the temporary vectors
  TVectorD fParProp_vec = TransferMat*fParProp_PCA;
  TVectorD fParCurr_vec = TransferMat*fParCurr_PCA;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < size; ++i) {
    _fPropVal[i] = fParProp_vec(i);
    _fCurrVal[i] = fParCurr_vec(i);
  }
}

const std::vector<double> covarianceBase::getProposed() const {
  std::vector<double> props;
  for (int i = 0; i < size; ++i) props.push_back(_fPropVal[i]);
  return props;
}

// Throw nominal values
void covarianceBase::throwNominal(bool nomValues, int seed) {

  TVectorD* vec = new TVectorD(size);
  for (int i = 0; i < size; i++) {
    (*vec)(i) = 1.0;
  }

  ThrowParms* nom_throws = new ThrowParms(*vec, (*covMatrix));
  nom_throws->SetSeed(seed);
  std::vector<double> nominal = getNominalArray();
  nominal.clear();
  nominal.resize(size);

  // If we want to put the nominals somewhere else than user specified
  // Don't fully understand this though: won't we have to reweight the MC somehow?
  // nominal[i] is used in GetLikelihood() as the penalty term, so we're essentially setting a random parameter penalty term?
  if (!nomValues)
  {
    bool throw_again = true;

    while(throw_again == true)
    {
      throw_again = false;
      std::cout << "- setting " << getName() << " nominal values to random throws." << std::endl;
      nom_throws->ThrowSet(nominal);

      for (int i = 0; i < size; i++)
      {
      // if parameter is fixed, dont throw
        if (_fError[i] < 0) {
          nominal[i] = 1.0;
          continue;
        }

        if (nominal[i] < 0) {
          nominal[i] = 0.0;
          throw_again = true;
        }
      }
    }
  } else {

    // If we want nominal values, set all entries to 1 (defined as nominal in MaCh3)
    for (int i = 0; i < int(nominal.size()); i++) {
      nominal[i] = 1.0;
    }
  }

  delete nom_throws;
  delete vec;
}

// *************************************
// Throw the parameters according to the covariance matrix
// This shoulnd't be used in MCMC code ase it can break Detailed Balance;
void covarianceBase::throwParameters() {
// *************************************

  // First draw new randParams
  randomize();

  if (!pca) {
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, size);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i) {
      // Check if parameter is fixed first: if so don't randomly throw
      if (isParameterFixed(i)) continue;

      _fPropVal[i] = _fPreFitValue[i] + corr_throw[i];
      int throws = 0;
      // Try again if we the initial parameter proposal falls outside of the range of the parameter
      while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
#ifdef MULTITHREAD
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0, 1);
#else
        randParams[i] = random_number[0]->Gaus(0,1);
#endif
        const double corr_throw_single = MatrixVectorMultiSingle(throwMatrixCholDecomp, randParams, size, i);
        _fPropVal[i] = _fPreFitValue[i] + corr_throw_single;
        if (throws > 10000) 
        {
          //KS: Since we are mulithreading there is danger that those messages
		  //will be all over the place, small price to pay for faster code
          std::cerr << "Tried " << throws << " times to throw parameter "; 
		  std::cerr << i << " but failed" << std::endl;
          std::cerr << "Matrix: " << matrixName << std::endl;
          std::cerr << "Param:  " << _fNames[i] << std::endl;
          std::cerr << "Setting _fPropVal:  " << _fPropVal[i] <<" to "<< _fPreFitValue[i]<<std::endl;
          std::cerr << "I live at " << __FILE__ << ":" << __LINE__ << std::endl;
          _fPropVal[i] = _fPreFitValue[i];
          //throw;
        }
        throws++;
      }
      _fCurrVal[i] = _fPropVal[i];
    }
  }
  else
  {      
      std::cerr << "Hold on, you are trying to run Prior Predicitve Code with PCA, which is wrong" << std::endl;
      std::cerr << "Sorry I have to kill you, I mean your job" << std::endl;
      std::cerr << "I live at " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
}

// *************************************
// Throw each parameter within their 1 sigma range
// Used to start the chain in different states
void covarianceBase::RandomConfiguration() {
  // *************************************
  // Have the 1 sigma for each parameter in each covariance class, sweet!
  // Don't want to change the nominal array because that's what determines our likelihood
  // Want to change the fParProp, fParCurr, fParInit
  // fParInit and the others will already be set
  for (int i = 0; i < size; ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (isParameterFixed(i)) continue;
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    const double paramrange = _fUpBound[i] - _fLowBound[i];
    const double sigma = sqrt((*covMatrix)(i,i));
    double throwrange = sigma;
    if (paramrange < sigma) throwrange = paramrange;

    _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
    // Try again if we the initial parameter proposal falls outside of the range of the parameter
    // Really only relevant for the xsec parameters; the flux and ND280 have -999 and 999 set to the limits!
    int throws = 0;
    while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
      if (throws > 1000) {
        std::cerr << "Tried " << throws << " times to throw parameter " << i << " but failed" << std::endl;
        std::cerr << "Matrix: " << matrixName << std::endl;
        std::cerr << "Param:  " << _fNames[i].c_str() << std::endl;
        throw;
      }
      _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
      throws++;
    }
    std::cout << "Setting current step in " << matrixName << " param " << i << " = " << _fPropVal[i] << " from " << _fCurrVal[i] << std::endl;
    _fCurrVal[i] = _fPropVal[i];
  }
  if (pca) TransferToPCA();

}

// *************************************
// Set a single parameter
void covarianceBase::setSingleParameter(const int parNo, const double parVal) {
  // *************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  std::cout << "Setting " << getParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;

  if (pca) TransferToPCA();
}

void covarianceBase::setParCurrProp(const int parNo, const double parVal) {
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  std::cout << "Setting " << getParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;
  if (pca) TransferToPCA();
}

// Propose a step for the set of systematics parameters this covariance class holds
void covarianceBase::proposeStep() {
  // Make the random numbers for the step proposal
  randomize();
  CorrelateSteps();
  if(use_adaptive && total_steps<upper_adapt) updateAdaptiveCovariance();
}

// ************************************************
// "Randomize" the parameters in the covariance class for the proposed step
// Used the proposal kernel and the current parameter value to set proposed step
// Also get a new random number for the randParams
void covarianceBase::randomize() {
// ************************************************
  if (!pca) {
//KS: By multithreading here we gain at least factor 2 with 8 threads with ND only fit      
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; i++) {
      // If parameter isn't fixed
      if (_fError[i] > 0.0) {
#ifdef MULTITHREAD
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0, 1);
#else
        randParams[i] = random_number[0]->Gaus(0, 1);
#endif 
        // If parameter IS fixed
      } else {
        randParams[i] = 0.0;
      }
    } // end for
  // If we're in the PCA basis we instead throw parameters there (only npars parameter)
  } else {
    // Scale the random parameters by the sqrt of eigen values for the throw
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; i++)
    {
      if (fParSigma_PCA[i] > 0. && i < npars)
      {
#ifdef MULTITHREAD        
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0,1);
#else        
        randParams[i] = random_number[0]->Gaus(0,1);
#endif
      } else { // If parameter IS fixed or out od bounds
        randParams[i] = 0.0;
      }
    }
  }
}


// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void covarianceBase::CorrelateSteps() {
// ************************************************
  //KS: Using custom fucntion comapred to ROOT one with 8 threads we have almost factor 2 performance increase, by replacing TMatrix with just double we increase it even more
  MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, size);

  // If not doing PCA
  if (!pca) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; i++) {
      if (_fError[i] > 0.) {
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*_fGlobalStepScale*_fIndivStepScale[i];
      }
    }
    // If doing PCA throw uncorrelated in PCA basis (orthogonal basis by definition)
  } else { 
    // Throw around the current step
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < npars; i++) 
    {
      if (fParSigma_PCA[i] > 0.) 
      {
        double IndStepScale = 1;
        //KS: If undecomposed parameter apply individual step scale and cholesky for better acceptance rate
        if(isDecomposed_PCA[i] >= 0)
        {
            IndStepScale *= _fIndivStepScale[isDecomposed_PCA[i]];
            IndStepScale *= corr_throw[isDecomposed_PCA[i]];
        }
        //If decomposed apply only random number
        else
        {
         IndStepScale *= randParams[i];
         //KS: All PCA-ed parameters have the same step scale
         IndStepScale *= _fIndivStepScale[FirstPCAdpar];
        }
        fParProp_PCA(i) = fParCurr_PCA(i)+_fGlobalStepScale*IndStepScale*eigen_values_master[i];
      }
    }
    // Then update the parameter basis
    TransferToParam();
  }
}

// Update so that current step becomes the previously proposed step
void covarianceBase::acceptStep() {
  if (!pca) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; i++) {
      // Update state so that current state is proposed state
      _fCurrVal[i] = _fPropVal[i];
    }
  } else {
  // Update the book-keeping for the output
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < npars; i++) {
      fParCurr_PCA(i) = fParProp_PCA(i);
    }
    // Then update the parameter basis
    TransferToParam();
  }
}

// Throw the proposed parameter by mag sigma
// Should really just have the user specify this throw by having argument double
void covarianceBase::throwParProp(const double mag) {
  randomize();
  if (!pca) {
    // Make the correlated throw
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, size);
    // Number of sigmas we throw
    for (int i = 0; i < size; i++) {
      if (_fError[i] > 0.)
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*mag;
    }
  } else {
    for (int i = 0; i < size; i++) {
      fParProp_PCA(i) = fParCurr_PCA(i)+mag*randParams[i];
    }
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD proposed = TransferMat*fParProp_PCA;
    for (int i = 0; i < size; ++i) {
      if (fParSigma_PCA[i] > 0.) {
        _fPropVal[i] = proposed(i);
      }
    }
  }
}

// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void covarianceBase::throwParCurr(const double mag)
{
  randomize();
  if (!pca) {
    // Get the correlated throw vector
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, size);
    // The number of sigmas to throw
    // Should probably have this as a default parameter input to the function instead
    for (int i = 0; i < size; i++) {
      if (_fError[i] > 0.){
        _fCurrVal[i] = corr_throw[i]*mag;
      }
    }
  } else {
    for (int i = 0; i < size; i++) {
      fParProp_PCA(i) = mag*randParams[i];
    }
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD current = TransferMat*fParCurr_PCA;
    for (int i = 0; i < size; ++i) {
      if (fParSigma_PCA[i] > 0.) {
        _fCurrVal[i] = current(i);
      }
    }
  }
}

// Function to print the nominal values
void covarianceBase::printNominal() {
  std::cout << "Prior values for " << getName() << " covarianceBase: " << std::endl;
  for (int i = 0; i < size; i++) {
    std::cout << "    " << getParName(i) << "   " << getParInit(i) << "\n";
  }
  std::cout << std::endl;
}

// Function to print the nominal, current and proposed values
void covarianceBase::printNominalCurrProp() {

  std::cout << "Printing parameters for " << getName() << std::endl;
  // Dump out the PCA parameters too
  if (pca) {
    std::cout << "PCA:" << "\n";
    for (int i = 0; i < npars; ++i) {
      std::cout << std::setw(PrintLength) << std::left << "PCA " << i << " Current: " << fParCurr_PCA(i) << " Proposed: " << fParProp_PCA(i) << std::endl;
    }
  }
  std::cout << std::setw(PrintLength) << std::left << "Name" << std::setw(PrintLength) << "Prior" << std::setw(PrintLength) << "Current" << std::setw(35) << "Proposed" << "\n";
  for (int i = 0; i < size; ++i) {
    std::cout << std::setw(PrintLength) << std::left << getParName(i) << std::setw(PrintLength) << _fPreFitValue[i] << std::setw(PrintLength) << _fCurrVal[i] << std::setw(PrintLength) << _fPropVal[i] << "\n";
  }
   //KS: "\n" is faster performance wise, keep std::endl at the end to flush just in case, also looks pretty
  std::cout << std::endl;
}

// Get the likelihood in the case where we want to include priors on the parameters
// fParEvalLikelihood stores if we want to evaluate the likelihood for the given parameter
//                    true = evaluate likelihood (so run with a prior)
//                    false = don't evaluate likelihood (so run without a prior)
double covarianceBase::CalcLikelihood() {
  double logL = 0.0;
  //TStopwatch clock;
  ///clock.Start();
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:logL)
#endif
  for(int i = 0; i < size; i++){
    for (int j = 0; j <= i; ++j) {
      if (_fFlatPrior[i] && _fFlatPrior[j]) {
        //KS: Since matrix is symetric we can calcaute non daigonal elements only once and multiply by 2, can bring up to factor speed decrease.   
        int scale = 1;
        if(i != j) scale = 2;
        logL += scale * 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];
      }
    }
  }
  //clock.Stop();
  //std::cout << __FILE__ << "::GetLikelihood took " << clock.RealTime() << "s" << std::endl;

  return logL;
}

int covarianceBase::CheckBounds(){
  int NOutside=0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:NOutside)
  #endif
  for (int i = 0; i < _fNumPar; i++){
      //if(_fPropVal[i] > xsec_param_ub_a[i] || _fPropVal[i] < xsec_param_lb_a[i]){
      if(_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]){
		std::cout << "_fPropVal at param " << i << "( " << getParName(i) << ") is out of bounds, param is at " << _fPropVal[i] << " and UB is " << _fUpBound[i] << " and LB is " << _fLowBound[i] << std::endl;
        NOutside++;
      }
  }
  return NOutside;
}

double covarianceBase::GetLikelihood(){
  // Checkbounds and calclikelihood are virtual
  // Default behaviour is to reject negative values + do std llh calculation
  const int NOutside = CheckBounds();
  
  if(NOutside>0){

	//ETA - don't need to check here. CheckBounds actually checks _fLowBound
	// and _fUpBound rather than 0 all the time
	//for (int i = 0; i < size; i++) {
	  // If the proposed step is negative, set a incredibly unlikely likelihood (but don't return yet for openMP!)
	  //if (_fPropVal[i] < 0) {std::cout << "_fPropVal at " << i << " is " << _fPropVal[i] << std::endl;}
	//}
	std::cout << "Parameters outside of bounds!" << std::endl;
	std::cout << "NOutside is " << NOutside << std::endl;
    return NOutside*__LARGE_LOGL__;
  }

  return CalcLikelihood();
}


void covarianceBase::printPars() {

  std::cout << "Number of pars: " << size << std::endl;
  std::cout << "current " << matrixName << " parameters:" << std::endl;
  for(int i = 0; i < size; i++) {
    std::cout << std::fixed << std::setprecision(5) << _fNames[i].c_str() << " current: \t" << _fCurrVal[i] << "   \tproposed: \t" << _fPropVal[i] << std::endl;
  }

  return;
}

// Sets the proposed parameters to the nominal values
void covarianceBase::setParameters(std::vector<double> pars) {

  // If empty, set the proposed to nominal
  if (pars.empty()) {
    // For xsec this means setting to the prior (because nominal is the prior)
    for (int i = 0; i < size; i++) {
      _fPropVal[i] = _fPreFitValue[i];
    }
    // If not empty, set the parameters to the specified
  } else {

	if (pars.size() != size_t(size)) {
      std::cerr << "Warning: parameter arrays of incompatible size! Not changing parameters! " << matrixName << " has size " << pars.size() << " but was expecting " << size << std::endl;
      throw;
    }

    unsigned int parsSize = pars.size();
    for (unsigned int i = 0; i < parsSize; i++) {
	  //Make sure that you are actually passing a number to set the parameter to
	  if(isnan(pars[i])) {
		std::cerr << "Error: trying to set parameter value to a nan for parameter " << getParName(i) << " in matrix " << matrixName << ". This will not go well!" << std::endl;
		throw;
	  } else {
		_fPropVal[i] = pars[i];
	  }
    }
  }

  // And if pca make the transfer
  if (pca) {
    TransferToPCA();
    TransferToParam();
  }

  return;
}

void covarianceBase::setBranches(TTree &tree) {
  // loop over parameters and set a branch
  for (int i = 0; i < size; ++i) {
    tree.Branch(_fNames[i].c_str(), &_fCurrVal[i], Form("%s/D", _fNames[i].c_str()));
  }
  // When running PCA, also save PCA parameters
  if (pca) {
    for (int i = 0; i < npars; ++i) { 
      tree.Branch(Form("%s_PCA", _fNames[i].c_str()), (double*)&(fParCurr_PCA.GetMatrixArray()[i]), Form("%s_PCA/D", _fNames[i].c_str()));
    }
  }
}

void covarianceBase::setStepScale(double scale) {
    if(scale == 0)
    {
        std::cerr << "You are trying so set StepScale to 0 this will not work"<< std::endl;
        throw;   
    }
  std::cout << getName() << " setStepScale() = " << scale << std::endl;
  _fGlobalStepScale = scale;
}

void covarianceBase::toggleFixAllParameters() {
  // fix or unfix all parameters by multiplying by -1
  if(!pca)
  {
    for (int i = 0; i < size; i++) _fError[i] *= -1.0;
  } else{
     for (int i = 0; i < npars; i++) fParSigma_PCA[i] *= -1.0;
  }
  return;
}

void covarianceBase::toggleFixParameter(const int i) {
  if(!pca) {
	if (i > size) {
	  std::cerr << "Can't toggleFixParameter for parameter " << i << " because size of covariance =" << size << std::endl;
	  std::cerr << "Fix this in your config file please!" << std::endl;
	  throw;
	} else {
	  _fError[i] *= -1.0;
	  std::cout << "Setting " << getParName(i) << "(parameter " << i << ") to fixed at " << _fCurrVal[i] << std::endl;
	} 
  } else {
	int isDecom = -1;
	for (int im = 0; im < npars; ++im) { 
	  if(isDecomposed_PCA[im] == i) {isDecom = im;}
	}
	if(isDecom < 0) {
	  std::cerr << "Parameter " << getParName(i) << " is PCA decomposed can't fix this" << std::endl;
	  //throw; 
	} else {
	  fParSigma_PCA[isDecom] *= -1.0;
	  std::cout << "Setting un-decomposed " << getParName(i) << "(parameter " << i <<"/"<< isDecom<< " in PCA base) to fixed at " << _fCurrVal[i] << std::endl;
	}
  }
  
  return;
}

void covarianceBase::setEvalLikelihood(int i, bool eL) {

  if (i > size) {
    std::cerr << "Can't setEvalLikelihood for " << getName() << "_" << i << " because size of covarianceXsec2015 = " << size << std::endl;
    std::cerr << "Fix this in your config file please!" << std::endl;
    throw;
  } else {
    std::cout << "Setting " << getParName(i) << " (parameter " << i << ") to flat prior? " << eL << std::endl;
    _fFlatPrior[i] = eL;
  }
}

// Multi-threaded matrix multiplication
TMatrixD covarianceBase::MatrixMult(TMatrixD A, TMatrixD B) {
  double *C_mon = MatrixMult(A.GetMatrixArray(), B.GetMatrixArray(), A.GetNcols());
  TMatrixD C;
  C.Use(A.GetNcols(), A.GetNrows(), C_mon);
  return C;
}

double** covarianceBase::MatrixMult(double **A, double **B, int n) {
  // First make into monolithic array
  double *A_mon = new double[n*n];
  double *B_mon = new double[n*n];

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A_mon[i*n+j] = A[i][j];
      B_mon[i*n+j] = B[i][j];
    }
  }
  // Now call the monolithic calculator
  double *C_mon = MatrixMult(A_mon, B_mon, n);
  delete A_mon;
  delete B_mon;

  // Return the double pointer
  double **C = new double*[n];
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    C[i] = new double[n];
    for (int j = 0; j < n; ++j) {
      C[i][j] = C_mon[i*n+j];
    }
  }
  delete C_mon;

  return C;
}

double* covarianceBase::MatrixMult(double *A, double *B, int n) {
  // First transpose to increse cache hits
  double *BT = new double[n*n];
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      BT[j*n+i] = B[i*n+j];
    }
  }

  // Now multiply
  double *C = new double[n*n];
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      for (int k = 0; k < n; k++) {
        sum += A[i*n+k]*BT[j*n+k];
      }
      C[i*n+j] = sum;
    }
  }
  delete BT;

  return C;
}

//KS: Custom function to perform multiplication of matrix and vector with mulithreadeing
void covarianceBase::MatrixVectorMulti(double* VecMulti, double** matrix, const double* vector, const int n) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++) {
	VecMulti[i] = 0.0;
	for (int j = 0; j < n; j++)
	{
	  VecMulti[i] += matrix[i][j]*vector[j];
	}
  }
  return;
}

double covarianceBase::MatrixVectorMultiSingle(double** matrix, const double* vector, const int Length, const int i)
{
  double Element = 0.0;
  for (int j = 0; j < Length; j++) {
	Element += matrix[i][j]*vector[j];
  }
  return Element;
}

void covarianceBase::setIndivStepScale(std::vector<double> stepscale) {
  if ((int)stepscale.size()!=size) {
    std::cout << "Stepscale vector not equal to number of parameters. Qutting.." << std::endl;
    std::cout << "Size of argument vector:" << stepscale.size() << std::endl;
    std::cout << "Expected size:" << size << std::endl;
    return;
  }

  for (int iParam=0;iParam<size;iParam++) {
    _fIndivStepScale[iParam] = stepscale[iParam];
  }

  printIndivStepScale();

  return;
}


void covarianceBase::printIndivStepScale() {
  std::cout << "============================================================" << std::endl;
  std::cout << std::setw(PrintLength) << "Parameter:" << " | " << std::setw(11) << "Step scale:" << std::endl;
  for (int iParam=0;iParam<size;iParam++) {
    std::cout << std::setw(PrintLength) << _fNames[iParam].c_str() << " | " << std::setw(11) << _fIndivStepScale[iParam] << std::endl;
  }
  std::cout << "============================================================" << std::endl;
}

//Makes sure that matrix is positive-definite (so no error is thrown when
//throwNominal() is called) by adding a small number to on-diagonal elements
void covarianceBase::MakePosDef(TMatrixDSym *cov) {
  //DB Save original warning state and then increase it in this function to suppress 'matrix not positive definite' messages
  //Means we no longer need to overload
  if(cov==NULL){
    cov = &*covMatrix;
  }

  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  //DB Loop 1000 times adding 1e-9 which tops out at 1e-6 shift on the diagonal before throwing error
  int MaxAttempts = 1e5;
  int iAttempt = 0;
  bool CanDecomp = false;
  TDecompChol chdcmp;
  
  for (iAttempt=0;iAttempt<MaxAttempts;iAttempt++) {
    chdcmp = TDecompChol(*cov);
    if (chdcmp.Decompose()) {
      CanDecomp = true;
      break;
    } else {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
      for (int iVar=0;iVar<size;iVar++) {
        (*cov)(iVar,iVar) += pow(10,-9);
      }
    }
  }

  if (!CanDecomp) {
    std::cerr << "Tried " << MaxAttempts << " times to shift diagonal but still can not decompose the matrix" << std::endl;
    std::cerr << "This indicates that something is wrong with the input matrix" << std::endl;
    throw;
  }
  if(total_steps<2) {
    std::cout << "Had to shift diagonal " << iAttempt << " time(s) to allow the covariance matrix to be decomposed" << std::endl;
  }
  //DB Reseting warning level
  gErrorIgnoreLevel = originalErrorWarning;

  return;
}

void covarianceBase::resetIndivStepScale() {
  std::vector<double> stepScales(size);
  for (int i=0;i<size;i++) {
    stepScales[i] = 1.;
  }
  _fGlobalStepScale=1.0;
  setIndivStepScale(stepScales);
}

//KS: Convert TMatrix to TH2D, mostly usefull for making fancy plots
TH2D* TMatrixIntoTH2D(const TMatrix &Matrix, std::string title)       
{
  TH2D* hMatrix = new TH2D(title.c_str(), title.c_str(), Matrix.GetNrows(), 0.0, Matrix.GetNrows(), Matrix.GetNcols(), 0.0, Matrix.GetNcols());
  for(int i = 0; i < Matrix.GetNrows(); i++)
  {
    for(int j = 0; j < Matrix.GetNcols(); j++)
    {
      //KS: +1 becasue there is offset in histogram realtive to TMatrix
      hMatrix->SetBinContent(i+1,j+1, (Matrix)(i,j));
    }
  }
  
  return hMatrix;
}

//KS: Let's dump all usefull matrices to properly validate PCA
#ifdef DEBUG_PCA
void covarianceBase::DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat)
{
  TFile *PCA_Debug = new TFile("Debug_PCA.root", "RECREATE");
  PCA_Debug->cd();
  
  bool PlotText = true;
  //KS: If we have more than 200 plot becomes unreadable :(
  if(size > 200) PlotText = false;
  
  TH1D* heigen_values = new TH1D("eigen_values", "Eigen Values", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_cumulative = new TH1D("heigen_cumulative", "heigen_cumulative", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_frac = new TH1D("heigen_fractional", "heigen_fractional", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  heigen_values->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative = 0;
  for(int i = 0; i < eigen_values.GetNrows(); i++) 
  {
      heigen_values->SetBinContent(i+1, (eigen_values)(i));
      heigen_cumulative->SetBinContent(i+1, (eigen_values)(i)/sum + Cumulative);
      heigen_frac->SetBinContent(i+1, (eigen_values)(i)/sum);
      Cumulative += (eigen_values)(i)/sum;
  }
  heigen_values->Write("heigen_values");
  eigen_values.Write("eigen_values");
  heigen_cumulative->Write("heigen_values_cumulative");
  heigen_frac->Write("heigen_values_frac");

  TH2D* heigen_vectors = TMatrixIntoTH2D(eigen_vectors, "eigen_vectors");  
  heigen_vectors->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors->Write("heigen_vectors");
  eigen_vectors.Write("eigen_vectors");
    
  TH2D* SubsetPCA = TMatrixIntoTH2D(temp, "SubsetPCA");  
  SubsetPCA->GetXaxis()->SetTitle("Parameter in Normal Base");
  SubsetPCA->GetYaxis()->SetTitle("Parameter in Decomposed Base");

  SubsetPCA->Write("hSubsetPCA");
  temp.Write("SubsetPCA");
  TH2D* hTransferMat = TMatrixIntoTH2D(TransferMat, "hTransferMat"); 
  hTransferMat->GetXaxis()->SetTitle("Parameter in Normal Base");
  hTransferMat->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  TH2D* hTransferMatT = TMatrixIntoTH2D(TransferMatT, "hTransferMatT"); 

  hTransferMatT->GetXaxis()->SetTitle("Parameter in Decomposed Base");
  hTransferMatT->GetYaxis()->SetTitle("Parameter in Normal Base");

  hTransferMat->Write("hTransferMat");
  TransferMat.Write("TransferMat");
  hTransferMatT->Write("hTransferMatT");
  TransferMatT.Write("TransferMatT");
  
  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 1024, 1024);
  c1->SetBottomMargin(0.1);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetGrid();
  
  gStyle->SetPaintTextFormat("4.1f"); 
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // Make pretty correlation colors (red to blue)
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
  
  double maxz = 0;
  double minz = 0;

  c1->Print("Debug_PCA.pdf[");
  TLine *EigenLine = new TLine(nKeptPCApars, 0, nKeptPCApars, heigen_cumulative->GetMaximum());
  EigenLine->SetLineColor(kPink);
  EigenLine->SetLineWidth(2);
  EigenLine->SetLineStyle(kSolid);
    
  TText* text = new TText(0.5, 0.5, Form("Threshold = %g", eigen_threshold));
  text->SetTextFont (43);
  text->SetTextSize (40);
  
  heigen_values->SetLineColor(kRed);
  heigen_values->SetLineWidth(2);
  heigen_cumulative->SetLineColor(kGreen);
  heigen_cumulative->SetLineWidth(2);
  heigen_frac->SetLineColor(kBlue);
  heigen_frac->SetLineWidth(2);

  c1->SetLogy();
  heigen_values->SetMaximum(heigen_cumulative->GetMaximum()+heigen_cumulative->GetMaximum()*0.4);
  heigen_values->Draw();
  heigen_frac->Draw("SAME");
  heigen_cumulative->Draw("SAME");
  EigenLine->Draw("Same");
  text->DrawTextNDC(0.42, 0.84,Form("Threshold = %g", eigen_threshold));
  
  TLegend *leg = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg->SetTextSize(0.04);
  leg->AddEntry(heigen_values, "Absolute", "l");
  leg->AddEntry(heigen_frac, "Fractional", "l");
  leg->AddEntry(heigen_cumulative, "Cumulative", "l");

  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("Same");
   
  c1->Print("Debug_PCA.pdf");
  c1->SetRightMargin(0.15);
  c1->SetLogy(0);
  delete EigenLine;
  delete leg;
  delete text;
  delete heigen_values;
  delete heigen_frac;
  delete heigen_cumulative;

  heigen_vectors->SetMarkerSize(0.2);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");
  
  TLine *Eigen_Line = new TLine(0, nKeptPCApars, LastPCAdpar-FirstPCAdpar, nKeptPCApars);
  Eigen_Line->SetLineColor(kGreen);
  Eigen_Line->SetLineWidth(2);
  Eigen_Line->SetLineStyle(kDotted);
  Eigen_Line->Draw("SAME");
  c1->Print("Debug_PCA.pdf");
  delete Eigen_Line;
  
  SubsetPCA->SetMarkerSize(0.2);
  minz = SubsetPCA->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) SubsetPCA->Draw("COLZ TEXT");
  else SubsetPCA->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete SubsetPCA;
  
  hTransferMat->SetMarkerSize(0.15);
  minz = hTransferMat->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMat->Draw("COLZ TEXT");
  else hTransferMat->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete hTransferMat;

  hTransferMatT->SetMarkerSize(0.15);
  minz = hTransferMatT->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMatT->Draw("COLZ TEXT");
  else hTransferMatT->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete hTransferMatT;

  
//KS: Crosscheck againts Eigen library
#if DEBUG_PCA == 2
  Eigen::MatrixXd Submat_Eigen(submat.GetNrows(), submat.GetNcols());
  
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for(int i = 0; i < submat.GetNrows(); i++)
  {
    for(int j = 0; j < submat.GetNcols(); j++)
    {
      Submat_Eigen(i,j) = (submat)(i,j);
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> EigenSolver;
  EigenSolver.compute(Submat_Eigen);
  Eigen::VectorXd eigen_val = EigenSolver.eigenvalues().real();
  Eigen::MatrixXd eigen_vect = EigenSolver.eigenvectors().real();
  std::vector<std::tuple<double, Eigen::VectorXd>> eigen_vectors_and_values; 
  double Sum_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++)
  {
    std::tuple<double, Eigen::VectorXd> vec_and_val(eigen_val[i], eigen_vect.row(i));
    eigen_vectors_and_values.push_back(vec_and_val);
    Sum_Eigen += eigen_val[i];
  }
  std::sort(eigen_vectors_and_values.begin(), eigen_vectors_and_values.end(), 
    [&](const std::tuple<double, Eigen::VectorXd>& a, const std::tuple<double, Eigen::VectorXd>& b) -> bool
  { return std::get<0>(a) > std::get<0>(b); } );
  int index = 0;
  for(auto const vect : eigen_vectors_and_values)
  {
        eigen_val(index) = std::get<0>(vect);
        eigen_vect.row(index) = std::get<1>(vect);
        index++;
  }
  TH1D* heigen_values_Eigen = new TH1D("eig_values", "Eigen Values", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_cumulative_Eigen = new TH1D("eig_cumulative", "heigen_cumulative", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_frac_Eigen = new TH1D("eig_fractional", "heigen_fractional", eigen_val.size(), 0.0, eigen_val.size());
  heigen_values_Eigen->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values_Eigen->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++) 
  {
      heigen_values_Eigen->SetBinContent(i+1, eigen_val(i));
      heigen_cumulative_Eigen->SetBinContent(i+1, eigen_val(i)/sum + Cumulative_Eigen);
      heigen_frac_Eigen->SetBinContent(i+1, eigen_val(i)/sum);
      Cumulative_Eigen += eigen_val(i)/sum;
  }
  heigen_values_Eigen->SetLineColor(kRed);
  heigen_values_Eigen->SetLineWidth(2);
  heigen_cumulative_Eigen->SetLineColor(kGreen);
  heigen_cumulative_Eigen->SetLineWidth(2);
  heigen_frac_Eigen->SetLineColor(kBlue);
  heigen_frac_Eigen->SetLineWidth(2);
  
  c1->SetLogy();
  heigen_values_Eigen->SetMaximum(heigen_cumulative_Eigen->GetMaximum()+heigen_cumulative_Eigen->GetMaximum()*0.4);
  heigen_values_Eigen->Draw();
  heigen_cumulative_Eigen->Draw("SAME");
  heigen_frac_Eigen->Draw("SAME");
  
  TLegend *leg_Eigen = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg_Eigen->SetTextSize(0.04);
  leg_Eigen->AddEntry(heigen_values_Eigen, "Absolute", "l");
  leg_Eigen->AddEntry(heigen_frac_Eigen, "Fractional", "l");
  leg_Eigen->AddEntry(heigen_cumulative_Eigen, "Cumulative", "l");

  leg_Eigen->SetLineColor(0);
  leg_Eigen->SetLineStyle(0);
  leg_Eigen->SetFillColor(0);
  leg_Eigen->SetFillStyle(0);
  leg_Eigen->Draw("Same");
  
  c1->Print( "Debug_PCA.pdf");
  c1->SetLogy(0);
  delete heigen_values_Eigen;
  delete heigen_cumulative_Eigen;
  delete heigen_frac_Eigen;
  delete leg_Eigen;
  
  TH2D* heigen_vectors_Eigen = new TH2D("Eigen_Vectors", "Eigen_Vectors", eigen_val.size(), 0.0, eigen_val.size(), eigen_val.size(), 0.0, eigen_val.size());
  
  for(int i = 0; i < eigen_val.size(); i++)
  {
    for(int j = 0; j < eigen_val.size(); j++)
    {
      //KS: +1 becasue there is offset in histogram realtive to TMatrix
      heigen_vectors_Eigen->SetBinContent(i+1,j+1, eigen_vect(i,j));
    }
  }
  heigen_vectors_Eigen->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors_Eigen->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors_Eigen->SetMarkerSize(0.15);
  minz = heigen_vectors_Eigen->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  
  if(PlotText) heigen_vectors_Eigen->Draw("COLZ TEXT");
  else heigen_vectors_Eigen->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");

  heigen_vectors->SetTitle("ROOT minus Eigen");
  heigen_vectors->Add(heigen_vectors_Eigen, -1.);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete heigen_vectors_Eigen;
  
#endif 
  delete heigen_vectors;

  c1->Print( "Debug_PCA.pdf]");
  delete c1;
  PCA_Debug->Close();
  delete PCA_Debug;
}
#endif


// HI : Code for throwing from separate throw matrix, needs to be set after init to ensure pos-def
void covarianceBase::setThrowMatrix(TMatrixDSym *cov){
   if (cov == NULL) {
    std::cerr << "Could not find covariance matrix you provided to setThrowMatrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  if (covMatrix->GetNrows() != cov->GetNrows()) {
    std::cerr << "Matrix given for throw Matrix is not the same size as the covariance matrix stored in object!" << std::endl;
    std::cerr << "Stored covariance matrix size:" << covMatrix->GetNrows() << std::endl;
    std::cerr << "Given matrix size:" << cov->GetNrows() << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  throwMatrix = (TMatrixDSym*)cov->Clone();
  if(total_steps<=lower_adapt) makeClosestPosDef(throwMatrix);
  else MakePosDef(throwMatrix);
  
  TDecompChol TDecompChol_throwMatrix(*throwMatrix);
  
  if(!TDecompChol_throwMatrix.Decompose()) {
    std::cerr << "Cholesky decomposition failed for " << matrixName << " trying to make positive definite" << std::endl;
    
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  throwMatrix_CholDecomp = new TMatrixD(TDecompChol_throwMatrix.GetU());
  throwMatrix_CholDecomp->T();

  //KS: ROOT has bad memory managment, using standard double means we can decrease most operation by factor 2 simply due to cache hits
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; ++j)
    {
        throwMatrixCholDecomp[i][j] = (*throwMatrix_CholDecomp)(i,j);
    }
  }
}

void covarianceBase::updateThrowMatrix(TMatrixDSym *cov){
  delete throwMatrix;
  throwMatrix = NULL;
  delete throwMatrix_CholDecomp;
  throwMatrix_CholDecomp = NULL;
  setThrowMatrix(cov);
}

// The setter
void covarianceBase::useSeparateThrowMatrix(TString throwMatrixFileName, TString throwMatrixName, TString meansVectorName){
  // Firstly let's check if the file exists
  TFile* throwMatrixFile = new TFile(throwMatrixFileName);
  resetIndivStepScale();
  if(throwMatrixFile->IsZombie()) {
    std::cerr<<"ERROR : Couldn't find throw Matrix file : "<<throwMatrixFileName<<std::endl;
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  } //We're done for now

  TMatrixDSym* tmp_throwMatrix = (TMatrixDSym*)throwMatrixFile->Get(throwMatrixName);
  TVectorD* tmp_meansvec = (TVectorD*)throwMatrixFile->Get(meansVectorName);

  if(!tmp_throwMatrix){
    std::cerr<<"ERROR : Couldn't find throw matrix "<<throwMatrixName<<" in "<<throwMatrixFileName<<std::endl;
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  if(!tmp_meansvec){
    std::cerr<<"ERROR : Couldn't find means vector "<<meansVectorName<<" in "<<throwMatrixFileName<<std::endl;
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  adaptiveCovariance = (TMatrixDSym*)tmp_throwMatrix->Clone();
  for(int iMean=0; iMean<tmp_meansvec->GetNrows(); iMean++){
      par_means.push_back((*tmp_meansvec)(iMean));
      updateThrowMatrix(adaptiveCovariance);
  }
  delete tmp_throwMatrix;
  delete tmp_meansvec;
  throwMatrixFile->Close();
  delete throwMatrixFile; // Just in case
}


void covarianceBase::useSeparateThrowMatrix(){
    initialiseNewAdaptiveChain(); 
    updateThrowMatrix(covMatrix);
}

// Truely adaptive MCMC!
void covarianceBase::updateAdaptiveCovariance(){
  // https://projecteuclid.org/journals/bernoulli/volume-7/issue-2/An-adaptive-Metropolis-algorithm/bj/1080222083.full
  // Updates adaptive matrix
  // First we update the total means

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for(int iRow=0; iRow<size; iRow++){
    par_means_prev[iRow]=par_means[iRow];
    par_means[iRow]=(_fCurrVal[iRow]+par_means[iRow]*total_steps)/(total_steps+1);  
  }

  //Now we update the covariances using cov(x,y)=E(xy)-E(x)E(y)
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for(int iRow=0; iRow<size; iRow++){
    for(int iCol=0; iCol<=iRow; iCol++){

      double cov_val = (*adaptiveCovariance)(iRow, iCol)*size/5.6644;
      cov_val += par_means_prev[iRow]*par_means_prev[iCol]; //First we remove the current means
      cov_val = (cov_val*total_steps+_fCurrVal[iRow]*_fCurrVal[iCol])/(total_steps+1); //Now get mean(iRow*iCol)
      cov_val -= par_means[iCol]*par_means[iRow];
      cov_val*=5.6644/size;
      (*adaptiveCovariance)(iRow, iCol) = cov_val;
      (*adaptiveCovariance)(iCol, iRow) = cov_val;
    }
  }
  //This is likely going to be the slow bit!
  total_steps+=1;
  if(total_steps==lower_adapt) 
  {
    resetIndivStepScale();
  }

  if(total_steps>=lower_adapt) {
    updateThrowMatrix(adaptiveCovariance); //Now we update and continue!
  }
  // if(total_steps%1000==0 && total_steps>1000){
  //   suboptimality_vals.push_back(calculatesuboptimality(covSqrt, adaptiveCovariance));
  // }
}

void covarianceBase::initialiseNewAdaptiveChain(){
  // If we don't have a covariance matrix to start from for adaptive tune we need to make one!
  adaptiveCovariance = new TMatrixDSym(size);
  //  par_means.reserve(size);
  for(int i=0; i<size; i++){
    par_means.push_back(0);
    par_means_prev.push_back(0);
    for(int j=0; j<=i; j++){
      (*adaptiveCovariance)(i,j)=0.0; // Just make an empty matrix
      (*adaptiveCovariance)(j,i)=0.0; // Just make an empty matrix
    }
  }
}

void covarianceBase::saveAdaptiveToFile(TString outFileName, TString systematicName){
  TFile* outFile = new TFile(outFileName, "UPDATE");
  if(outFile->IsZombie()){
    std::cerr<<"ERROR : Couldn't find "<<outFileName<<std::endl;
    throw;
  }
  TVectorD* outMeanVec = new TVectorD((int)par_means.size());
  for(int i=0; i<(int)par_means.size(); i++){
    (*outMeanVec)(i)=par_means[i];
  }
  outFile->cd();
  adaptiveCovariance->Write(systematicName+"_postfit_matrix");
  outMeanVec->Write(systematicName+"_mean_vec");
  outFile->Close();
  delete outFile;

}

//HI Finds closest possible positive definite matrix in Frobenius Norm ||.||_frob
// Where ||X||_frob=sqrt[sum_ij(x_ij^2)] (basically just turns an n,n matrix into vector in n^2 space
// then does Euclidean norm)
void covarianceBase::makeClosestPosDef(TMatrixDSym *cov)
{
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
    std::cerr<<"Cannot do SVD on input matrix!"<<std::endl;
    throw;
  }
  
  TMatrixD cov_sym_v = (TMatrixD)cov_sym_svd.GetV();
  TMatrixD cov_sym_vt = cov_sym_v;
  cov_sym_vt.T();
  //SVD returns as vector (grrr) so need to get into matrix form for multiplying!
  TVectorD cov_sym_sigvect = (TVectorD)cov_sym_svd.GetSig();
  const Int_t nCols = cov_sym_v.GetNcols(); //square so only need rows hence lack of cols
  TMatrixDSym cov_sym_sig(nCols);
  TMatrixDDiag cov_sym_sig_diag(cov_sym_sig);
  cov_sym_sig_diag=cov_sym_sigvect;

  
  //Can finally get H=VSV
  TMatrixDSym cov_sym_polar =cov_sym_sig.SimilarityT(cov_sym_vt);//V*S*V^T (this took forver to find!)
  
  //Now we can construct closest approximater Ahat=0.5*(B+H)
  TMatrixDSym cov_closest_approx  = 0.5*(cov_sym+cov_sym_polar);//Not fully sure why this is even needed since symmetric B -> U=V
  //Get norm of transformed
//  Double_t approx_norm=cov_closest_approx.E2Norm();

  //std::cout<<"Initial Norm : "<<cov_norm<<" | Norm after transformation : "<<approx_norm<<" | Ratio : "<<cov_norm/approx_norm<<std::endl;
  
  *cov = cov_closest_approx;
  //Now can just add a makeposdef!
  MakePosDef(cov);
}

std::vector<double> covarianceBase::getNominalArray()
{
 std::vector<double> nominal;
  for (int i = 0; i < size; i++)
  {
    nominal.push_back(_fPreFitValue[i]);
  }

 return nominal;

}

