#include "covarianceBase.h"


covarianceBase::covarianceBase(const char *name, const char *file) : inputFile(std::string(file)), pca(false) {
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
  
  init(name, file);
  firstpcadpar=-999;
  lastpcadpar=-999;
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
  firstpcadpar=-999;
  lastpcadpar=-999;
}


covarianceBase::covarianceBase(const char *name, const char *file, int seed, double threshold, int firstpcapar, int lastpcapar) : inputFile(std::string(file)), pca(true), eigen_threshold(threshold), firstpcadpar(firstpcapar), lastpcadpar(lastpcapar) {

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
  init(name, file);

  // Call the innocent helper function
  if (pca) ConstructPCA();
}


// ********************************************
//Desctructor
covarianceBase::~covarianceBase(){
// ********************************************
  delete[] fParInit;
  delete[] fParSigma;
  delete[] fParCurr;
  delete[] fParProp;
  delete[] fParLoLimit;
  delete[] fParHiLimit;
  delete[] fIndivStepScale;
  delete[] fParEvalLikelihood;
  
  if (covMatrix != NULL) delete covMatrix;
  if (invCovMatrix != NULL) delete invCovMatrix;
  if (chel != NULL) delete chel;
    fPropKernel = new TF1*[size]();
  for(int i = 0; i < size; i++) 
  {
    delete  fPropKernel[i];
  }
  delete[] fPropKernel;
  
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
  if(firstpcadpar==-999||lastpcadpar==-999){
    if(firstpcadpar==-999&&lastpcadpar==-999){
      firstpcadpar=0;
      lastpcadpar=covMatrix->GetNrows()-1;
    }
    else{
      std::cerr << "You must either leave firstpcadpar and lastpcadpar at -999 or set them both to something"<<std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
  }
  if(firstpcadpar>covMatrix->GetNrows()-1||lastpcadpar>covMatrix->GetNrows()-1){
    std::cerr <<"firstpcadpar and lastpcadpar are higher than the number of parameters"<<std::endl;
    std::cerr <<"first: "<<firstpcadpar<<" last: "<<lastpcadpar<<std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  if(firstpcadpar<0||lastpcadpar<0){
    std::cerr <<"firstpcadpar and lastpcadpar are less than 0 but not default -999"<<std::endl;
    std::cerr <<"first: "<<firstpcadpar<<" last: "<<lastpcadpar<<std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  std::cout<<"PCAing parameters "<<firstpcadpar<<" through "<<lastpcadpar<<" inclusive"<<std::endl;

  int numunpcadpars=covMatrix->GetNrows()-(lastpcadpar-firstpcadpar+1);

  TMatrixDSym submat(covMatrix->GetSub(firstpcadpar,lastpcadpar,firstpcadpar,lastpcadpar));

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
  nkeptpcapars=eigen_values.GetNrows();
  // Now go through again and see how many eigen values correspond to threshold
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    // Get the relative size of the eigen value
    double sig = eigen_values(i)/sum;
    // Check against the threshold
    if (sig < eigen_threshold) {
      nkeptpcapars = i;
      break;
    }
  }
  npars=numunpcadpars+nkeptpcapars;
  std::cout << "Threshold of " << eigen_threshold << " on eigen values relative sum of eigen value (" << sum << ") generates " << nkeptpcapars << " eigen vectors, plus we have "<<numunpcadpars<<" unpcad pars, for a total of "<<npars << std::endl;

  //DB Create array of correct size so eigen_values can be used in CorrelateSteps
  eigen_values_master = std::vector<double>(npars,1.0);
  for (int i=firstpcadpar;i<npars;i++) {eigen_values_master[i] = eigen_values(i-firstpcadpar);}

  // Now construct the transfer matrices
  //These matrices will be as big as number of unPCAd pars plus number of eigenvalues kept
  TransferMat.ResizeTo(covMatrix->GetNrows(), npars);
  TransferMatT.ResizeTo(covMatrix->GetNrows(), npars);


  // Get a subset of the eigen vector matrix
  TMatrixD temp(eigen_vectors.GetSub(0, eigen_vectors.GetNrows()-1, 0, nkeptpcapars-1));
  
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
  if(firstpcadpar!=0){
    for(int iRow=0;iRow<firstpcadpar;iRow++){
      temp2[iRow][iRow]=1;
    }
  }


  //Set the transfer matrix block for the PCAd pars
  temp2.SetSub(firstpcadpar,firstpcadpar,temp);


  //Set the second identity block
  if(lastpcadpar!=covMatrix->GetNrows()-1){
    for(int iRow=0;iRow<(covMatrix->GetNrows()-1)-lastpcadpar;iRow++){
      temp2[lastpcadpar+1+iRow][firstpcadpar+nkeptpcapars+iRow]=1;
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
  randParams.ResizeTo(npars);

  // Make the PCA parameter arrays
  fParCurr_PCA.ResizeTo(npars);
  fParProp_PCA.ResizeTo(npars);
  
}

void covarianceBase::init(const char *name, const char *file) {

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

  // Set the covariance matrix
  setCovMatrix(covMatrix);
  setName(name);
  size = covMatrix->GetNrows();
  if (size <= 0) {
    std::cerr << "Covariance matrix " << getName() << " has " << size << " entries!" << std::endl;
    throw;
  }
  npars = size;

  // Set the nominal values to 1
  for (int i = 0; i < size; i++) {
    nominal.push_back(1.0);
  }

  fPropKernel = new TF1*[size]();
  fParNames = new Char_t*[size]();
  fParInit = new Double_t[size]();
  fParSigma = new Double_t[size]();
  fParCurr = new Double_t[size]();
  fParProp = new Double_t[size]();
  fParLoLimit = new Double_t[size]();
  fParHiLimit = new Double_t[size]();
  fIndivStepScale = new Double_t[size];
  fParEvalLikelihood =  new bool[size]();

  // Set the defaults to true
  for(int i = 0; i < size; i++) {
    fParInit[i] = 0;
    fParSigma[i] = 0;
    fParCurr[i] = 0;
    fParProp[i] = 0;
    fParLoLimit[i] = -999.99;
    fParHiLimit[i] = 999.99;
    fParEvalLikelihood[i] = true;
    fIndivStepScale[i] = 1.;
  }

  // Set the logLs to very large so next step is accepted
  currLogL = __LARGE_LOGL__;
  propLogL = __LARGE_LOGL__;

  // Set random parameter vector (for correlated steps)
  randParams.ResizeTo(size);

  infile->Close();

  fStepScale = 1.0;

  std::cout << "Created covariance matrix named: " << getName() << std::endl;
  std::cout << "from file: " << file << std::endl;

  delete infile;
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
}

// Set all the covariance matrix parameters to a user-defined value
void covarianceBase::setPar(int i , double val) {

  std::cout << "Over-riding " << getParName(i) << ": " << std::endl;
  std::cout << "fParProp (" << fParProp[i];
  std::cout << "), fParCurr (" << fParCurr[i];
  std::cout << "), nominal (" << nominal[i];
  std::cout << "), fParInit (" << fParInit[i]; 
  std::cout << ") to " << val << std::endl;

  fParProp[i] = val;
  fParCurr[i] = val;
  fParInit[i] = val;
  nominal[i]  = val;

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
    fParCurr_vec(i) = fParCurr[i];
    fParProp_vec(i) = fParProp[i];
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
    fParProp[i] = fParProp_vec(i);
    fParCurr[i] = fParCurr_vec(i);
  }
}

const std::vector<double> covarianceBase::getProposed() const {
  std::vector<double> props;
  for (int i = 0; i < size; ++i) props.push_back(fParProp[i]);
  return props;
}

// Throw nominal values
void covarianceBase::throwNominal(bool nomValues, int seed) {

  TDecompChol chdcmp(*covMatrix);

  if (!chdcmp.Decompose()) {
    std::cerr << "Cholesky decomposition failed for " << matrixName << std::endl;
    throw;
  }

  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);

  TVectorD* vec = new TVectorD(size);
  for (int i = 0; i < size; i++) {
    (*vec)(i) = 1.0;
  }

  ThrowParms* nom_throws = new ThrowParms(*vec, (*covMatrix));
  nom_throws->SetSeed(seed);
  nominal.clear();
  nominal.resize(size);

  // If we want to put the nominals somewhere else than user specified
  // Don't fully understand this though: won't we have to reweight the MC somehow?
  // nominal[i] is used in GetLikelihood() as the penalty term, so we're essentially setting a random parameter penalty term?
  if (!nomValues) {
    bool throw_again = true;

    while(throw_again == true) {

      throw_again = false;
      std::cout << "- setting " << getName() << " nominal values to random throws." << std::endl;
      nom_throws->ThrowSet(nominal);

      for (int i = 0; i < size; i++) {

      // if parameter is fixed, dont throw
        if (fParSigma[i] < 0) {
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

// This looks like a custom Cholesky decomposition function
// Someone trying to reinvent the ROOT wheel?
void covarianceBase::CholeskyDecomp(int npars, TMatrixD &chel_mat) {

  // Loop over the size of chel_mat
  if (npars != chel_mat.GetNrows()) {
    std::cerr << "Number of parameters passed to CholeskyDecomp != size of matrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Loop over the entries columns
  for (int i = 0; i < npars; i++) {

    // Loop over the entries rows
    for (int j = 0; j < npars; j++) {

      // If diagonal element
      if (i == j) {
        chel_mat(i,i) = (*covMatrix)(i,i);
        // For a given entry, loop over the entries up to i
        for (int k = 0; k <= i-1; k++) {
          chel_mat(i,i) = chel_mat(i,i) - pow(chel_mat(i,k),2);
        }
        chel_mat(i,i) = sqrt(chel_mat(i,i));

      // If lower half
      } else if (j < i) {

        chel_mat(i,j) = (*covMatrix)(i,j);
        for (int k = 0; k <= j - 1; k++) {
          chel_mat(i,j) = chel_mat(i,j)-chel_mat(i,k)*chel_mat(j,k);
        }
        chel_mat(i,j) = chel_mat(i,j)/chel_mat(j,j);
      } else {
        chel_mat(i,j) = 0.;
      }
    }
  }

}

// Update the Cholesky matrix with post-fit
void covarianceBase::UpdateCholeskyMatrix() {

}

// *************************************
// Make the proposal kernel from Gaussians and the errors we have on the parameters
void covarianceBase::genPropKernels() {
// *************************************
  Char_t aname[1024];

  for (int i = 0; i < size; i++) {
    fPropKernel[i] = NULL;
    if (fParSigma[i] > 0.0) {
      // This is the error on the ith parameter
      double e = TMath::Sqrt((*covMatrix)(i,i));
      sprintf(aname, "%s_%03d",getName(), i);
      fPropKernel[i] = new TF1(aname, "gaus", -e*10.0, e*10.0);
      fPropKernel[i]->SetParameters(1. / (e*TMath::Sqrt(2*3.14159265)), 0., e);
    }
  }
}

// *************************************
// Throw the parameters according to the covariance matrix
void covarianceBase::throwParameters() {
// *************************************

  // First draw new randParams
  randomize();

  if (!pca) {
    // Then make the TVectorD
    TVectorD CorrThrows = (*chel)*randParams;

    for (int i = 0; i < size; ++i) {
      // Check if parameter is fixed first: if so don't randomly throw
      if (isParameterFixed(i)) continue;

      fParProp[i] = fParInit[i] + CorrThrows(i);
      int throws = 0;
      // Try again if we the initial parameter proposal falls outside of the range of the parameter
      while (fParProp[i] > fParHiLimit[i] || fParProp[i] < fParLoLimit[i]) {
        randParams(i) = random_number[0]->Gaus(0,1);
   
        fParProp[i] = fParInit[i] + ((*chel)*randParams)(i);
        if (throws > 10000) 
        {
          std::cerr << "Tried " << throws << " times to throw parameter " << i << " but failed" << std::endl;
          std::cerr << "Matrix: " << matrixName << std::endl;
          std::cerr << "Param:  " << fParNames[i] << std::endl;
          std::cerr << "Setting fParProp:  " << fParProp[i] <<" to "<<fParInit[i]<<std::endl;
          std::cerr << "I live at " << __FILE__ << ":" << __LINE__ << std::endl;
          fParProp[i] = fParInit[i];
          //throw;
        }
        throws++;
      }
      fParCurr[i] = fParProp[i];
    }
  }
  else
  {      
      std::cerr << "Hold on, you are trying to run Prior/Posterior Predicitve Code with PCA, which is wrong" << std::endl;
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
    double paramrange = fParHiLimit[i] - fParLoLimit[i];
    double sigma = sqrt((*covMatrix)(i,i));
    double throwrange = sigma;
    if (paramrange < sigma) throwrange = paramrange;

    fParProp[i] = fParInit[i] + random_number[0]->Gaus(0, 1)*throwrange;
    // Try again if we the initial parameter proposal falls outside of the range of the parameter
    // Really only relevant for the xsec parameters; the flux and ND280 have -999 and 999 set to the limits!
    double oldval = fParProp[i];
    int throws = 0;
    while (fParProp[i] > fParHiLimit[i] || fParProp[i] < fParLoLimit[i]) {
      if (throws > 1000) {
        std::cerr << "Tried " << throws << " times to throw parameter " << i << " but failed" << std::endl;
        std::cerr << "Matrix: " << matrixName << std::endl;
        std::cerr << "Param:  " << fParNames[i] << std::endl;
        throw;
      }
      fParProp[i] = oldval;
      fParProp[i] = fParInit[i] + random_number[0]->Gaus(0, 1)*throwrange;
      throws++;
    }
    std::cout << "Setting current step in " << matrixName << " param " << i << " = " << fParProp[i] << " from " << fParCurr[i] << std::endl;
    fParCurr[i] = fParProp[i];
  }
  if (pca) TransferToPCA();

}

// *************************************
// Set a single parameter
void covarianceBase::setSingleParameter(int parNo, double parVal) {
  // *************************************
  fParProp[parNo] = parVal;
  fParCurr[parNo] = parVal;
  std::cout << "Setting " << getParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;

  if (pca) TransferToPCA();
}

void covarianceBase::setParCurrProp(int parNo, double parVal) {
  fParProp[parNo] = parVal;
  fParCurr[parNo] = parVal;
  std::cout << "Setting " << getParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;
  if (pca) TransferToPCA();
}

void covarianceBase::setPropFunct(int i, TF1 *func) {

  if (i >= size) {
    std::cerr << "Parameter index out of range!" << std::endl;
    throw;
  }

  fPropKernel[i] = func;
  fParSigma[i] = func->GetParameter(1); // this maybe bad if the function isnt gaussian...
  std::cout << "Set " << getName() << " parameter " << i << "'s proposal function to " << func->GetName() << std::endl;
}

// Generate the posterior histograms (without filling them!)
// CWRET this is probably deprecated by the MCMC class
// Used to get filled in acceptStep()
void covarianceBase::genPosteriorHists() {
  Char_t name[1024];
  for (int i = 0; i < size; i++) {
    posterior[i] = NULL;
    sprintf(name, "post_%s_%03i",getName(), i);
    posterior[i] = new TH1D(name, name, 100, -5, 5);
  }
}

// Propose a step for the set of systematics parameters this covariance class holds
void covarianceBase::proposeStep() {
  // Make the random numbers for the step proposal
  randomize();
  CorrelateSteps();
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
      if (fParSigma[i] > 0.0) {
#ifdef MULTITHREAD
        //fParProp[i] = fParCurr[i] + fPropKernel[i]->GetRandom(); //this is overwritten in CorrelateSteps()
        randParams(i) = random_number[omp_get_thread_num()]->Gaus(0, 1);
#else
        //fParProp[i] = fParCurr[i] + fPropKernel[i]->GetRandom(); //this is overwritten in CorrelateSteps()
        randParams(i) = random_number[0]->Gaus(0, 1);
#endif 
        // If parameter IS fixed
      } else {
        randParams(i) = 0.0;
      }
    } // end for
  // If we're in the PCA basis we instead throw parameters there (only npars parameter)
  } else {
    // Scale the random parameters by the sqrt of eigen values for the throw
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < npars; i++) {
#ifdef MULTITHREAD        
      randParams(i) = random_number[omp_get_thread_num()]->Gaus(0,1);
#else        
      randParams(i) = random_number[0]->Gaus(0,1);
#endif
    }
  }
}


// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void covarianceBase::CorrelateSteps() {
// ************************************************

  // If not doing PCA
  if (!pca) {
    TVectorD corr_throw = (*chel)*randParams;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < size; i++) {
      if (fParSigma[i] > 0.) {
        fParProp[i] = fParCurr[i] + corr_throw(i)*fStepScale*fIndivStepScale[i];
      }
    }
    // If doing PCA throw uncorrelated in PCA basis (orthogonal basis by definition)
  } else { 

    // Throw around the current step
    //fParProp_PCA = fStepScale*randParams;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < npars; i++) {
      fParProp_PCA(i) = fParCurr_PCA(i)+fStepScale*fIndivStepScale[i]*randParams(i)*eigen_values_master[i];
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
      fParCurr[i] = fParProp[i];
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
void covarianceBase::throwParProp() {
  randomize();
  double mag = 1.;
  if (!pca) {
    // Make the correlated throw
    TVectorD corr_throw = (*chel)*randParams;
    // Number of sigmas we throw
    for (int i = 0; i < size; i++) {
      if (fParSigma[i] > 0.)
        fParProp[i] = fParCurr[i] + corr_throw(i)*mag;
    }
  } else {
    fParProp_PCA = fParCurr_PCA+mag*randParams;
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD proposed = TransferMat*fParProp_PCA;
    for (int i = 0; i < size; ++i) {
      fParProp[i] = proposed(i);
    }
  }
}

// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void covarianceBase::throwParCurr() {
  randomize();
  double mag = 1.;
  if (!pca) {
    // Get the correlated throw vector
    TVectorD corr_throw = (*chel)*randParams;
    // The number of sigmas to throw
    // Should probably have this as a default parameter input to the function instead
    for (int i = 0; i < size; i++) {
      if (fParSigma[i] > 0.) {
        fParCurr[i] = corr_throw(i)*mag;
      }
    }
  } else {
    fParCurr_PCA = mag*randParams;
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD current = TransferMat*fParCurr_PCA;
    for (int i = 0; i < size; ++i) {
      fParCurr[i] = current(i);
    }
  }

}

// Function to print the nominal values
void covarianceBase::printNominal() {
  std::cout << "Nominal values for " << getName() << " covarianceBase: " << std::endl;
  for (int i = 0; i < size; i++) {
    std::cout << "    " << getParName(i) << "   " << getNominal(i) << std::endl;
  }
  std::cout << std::endl;
}

// Function to print the nominal, current and proposed values
void covarianceBase::printNominalCurrProp() {


std::cout << "Printing parameters for " << getName() << std::endl;
  // Dump out the PCA parameters too
  if (pca) {
    std::cout << "PCA:" << std::endl;
    for (int i = 0; i < npars; ++i) {
      std::cout << std::setw(35) << std::left << "PCA " << i << " Current: " << fParCurr_PCA(i) << " Proposed: " << fParProp_PCA(i) << std::endl;
    }
  }
  std::cout << std::setw(35) << std::left << "Name" << std::setw(35) << "Nominal" << std::setw(35) << "Init" << std::setw(35) << "Current" << std::setw(35) << "Proposed" << std::endl;
  for (size_t i = 0; i < nominal.size(); ++i) {
    std::cout << std::setw(35) << std::left << getParName(i) << std::setw(35) << nominal[i] << std::setw(35) << fParInit[i] << std::setw(35) << fParCurr[i] << std::setw(35) << fParProp[i] << std::endl;
  }
}

// Get the likelihood in the case where we want to include priors on the parameters
// fParEvalLikelihood stores if we want to evaluate the likelihood for the given parameter
//                    true = evaluate likelihood (so run with a prior)
//                    false = don't evaluate likelihood (so run without a prior)
double covarianceBase::getLikelihood() {

  double logL = 0.0;
  //TStopwatch clock;
  ///clock.Start();
  // Profiling shows that invCovMatrix almost always cache misses, change this to an array instead
  // Probably some strange structure of the TMatrixTSym

#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:logL)
#endif
  for (int i = 0; i < size; i++) {

    // If the proposed step is negative, set a incredibly unlikely likelihood (but don't return yet for openMP!)
    if (fParProp[i] < 0) {
      logL += __LARGE_LOGL__;
      std::cout << "Param " << i << " in " << matrixName << " is " << fParProp[i] << std::endl;
      std::cout << "Rejecting!" << std::endl;
    }

    for (int j = 0; j <= i; ++j) {
      if (fParEvalLikelihood[i] && fParEvalLikelihood[j]) {
        //KS: Since matrix is symetric we can calcaute non daigonal elements only once and multiply by 2, can bring up to factor speed decrease.   
        int scale = 1;
        if(i != j) scale = 2;
        logL += scale * 0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*(*invCovMatrix)(i,j);
      }
    }
  }

  //clock.Stop();
  //std::cout << __FILE__ << "::getLikelihood took " << clock.RealTime() << "s" << std::endl;

  return logL;
}


void covarianceBase::printPars() {

  std::cout << "Number of pars: " << size << std::endl;
  std::cout << "current " << matrixName << " parameters:" << std::endl;
  for(int i = 0; i < size; i++) {
    std::cout << std::fixed << std::setprecision(5) << fParNames[i] << " current: \t" << fParCurr[i] << "   \tproposed: \t" << fParProp[i] << std::endl;
  }

}

// Sets the proposed parameters to the nominal values
void covarianceBase::setParameters(std::vector<double> pars) {

  // If empty, set the proposed to nominal
  if (pars.empty()) {
    // For xsec this means setting to the prior (because nominal is the prior)
    for (int i = 0; i < size; i++) {
      fParProp[i] = nominal[i];
    }
    // If not empty, set the parameters to the specified
  } else {

    if (pars.size() != size_t(size)) {
      std::cerr << "Warning: parameter arrays of incompatible size! Not changing parameters! " << matrixName << " has size " << pars.size() << " but was expecting " << size << std::endl;
      throw;
    }

    unsigned int parsSize = pars.size();
    for (unsigned int i = 0; i < parsSize; i++) {
      fParProp[i] = pars[i];
    }
  }

  // And if pca make the transfer
  if (pca) {
    TransferToPCA();
    TransferToParam();
  }
}

void covarianceBase::setBranches(TTree &tree) {
  // loop over parameters and set a branch
  for (int i = 0; i < size; ++i) {
    tree.Branch(fParNames[i], (double*)&fParCurr[i], Form("%s/D", fParNames[i]));
  }
  // When running PCA, also save PCA parameters
  if (pca) {
    for (int i = 0; i < npars; ++i) { 
      tree.Branch(Form("%s_PCA", fParNames[i]), (double*)&(fParCurr_PCA.GetMatrixArray()[i]), Form("%s_PCA/D", fParNames[i]));
    }
  }
}

void covarianceBase::setStepScale(double scale) {
  std::cout << getName() << " setStepScale() = " << scale << std::endl;
  fStepScale = scale;
}

void covarianceBase::toggleFixAllParameters() {
  // fix or unfix all parameters by multiplying by -1
  for (int i = 0; i < size; i++) fParSigma[i] *= -1.0;
}

void covarianceBase::toggleFixParameter(int i) {

  if (i > size) {
    std::cerr << "Can't toggleFixParameter for parameter " << i << " because size of covariance =" << size << std::endl;
    std::cerr << "Fix this in your config file please!" << std::endl;
    throw;
  } else {
    fParSigma[i] *= -1.0;
    std::cout << "Setting " << getParName(i) << "(parameter " << i << ") to fixed at " << fParCurr[i] << std::endl;
  }
}

void covarianceBase::setEvalLikelihood(int i, bool eL) {

  if (i > size) {
    std::cerr << "Can't setEvalLikelihood for " << getName() << "_" << i << " because size of covarianceXsec2015 = " << size << std::endl;
    std::cerr << "Fix this in your config file please!" << std::endl;
    throw;
  } else {
    std::cout << "Setting " << getParName(i) << " (parameter " << i << ") to flat prior? " << eL << std::endl;
    fParEvalLikelihood[i] = eL;
  }
}


void covarianceBase::SetBANFFCov() {
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

#if MULTITHREAD
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
#if MULTITHREAD
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
#if MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      BT[j*n+i] = B[i*n+j];
    }
  }

  // Now multiply
  double *C = new double[n*n];
#if MULTITHREAD
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

void covarianceBase::setIndivStepScale(std::vector<double> stepscale) {
  if ((int)stepscale.size()!=size) {
    std::cout << "Stepscale vector not equal to number of parameters. Qutting.." << std::endl;
    std::cout << "Size of argument vector:" << stepscale.size() << std::endl;
    std::cout << "Expected size:" << size << std::endl;
    return;
  }

  for (int iParam=0;iParam<size;iParam++) {
    fIndivStepScale[iParam] = stepscale[iParam];
  }

  printIndivStepScale();

  return;
}


void covarianceBase::printIndivStepScale() {
  std::cout << "============================================================" << std::endl;
  std::cout << std::setw(35) << "Parameter:" << " | " << std::setw(11) << "Step scale:" << std::endl;
  for (int iParam=0;iParam<size;iParam++) {
    std::cout << std::setw(35) << fParNames[iParam] << " | " << std::setw(11) << fIndivStepScale[iParam] << std::endl;
  }
  std::cout << "============================================================" << std::endl;
}

void covarianceBase::MakePosDef() //Makes sure that matrix is positive-definite (so no error is thrown when throwNominal() is called) by adding a small number to all on-diagonal elements           
{
  //DB Save original warning state and then increase it in this function to suppress 'matrix not positive definite' messages
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  //DB Loop 1000 times adding 1e-9 which tops out at 1e-6 shift on the diagonal before throwing error
  int MaxAttempts = 1000;
  int iAttempt = 0;
  bool CanDecomp = false;
  TDecompChol chdcmp;
  
  for (iAttempt=0;iAttempt<MaxAttempts;iAttempt++) {
    chdcmp = TDecompChol(*covMatrix);
    if (chdcmp.Decompose()) {
      CanDecomp = true;
      break;
    } else {
      for (int iVar=0;iVar<size;iVar++) {
	(*covMatrix)(iVar,iVar) += pow(10,-9);
      }
    }
  }

  if (!CanDecomp) {
    std::cerr << "Tried " << MaxAttempts << " times to shift diagonal but still can not decompose the matrix" << std::endl;
    std::cerr << "This indicates that something is wrong with the input matrix" << std::endl;
    throw;
  }

  std::cout << "Had to shift diagonal " << iAttempt << " time(s) to allow the covariance matrix to be decomposed" << std::endl;
  //DB Reseting warning level
  gErrorIgnoreLevel = originalErrorWarning;

  // TH: Allows LLH calculator to pick up correct matrix
  setCovMatrix(covMatrix);
}

