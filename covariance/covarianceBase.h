#ifndef _covarianceBase_h_
#define _covarianceBase_h_

// ROOT includes
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TMath.h"
#include "math.h"
#include "TDecompChol.h"
#include "TStopwatch.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDEigen.h"
#include "TDecompSVD.h"

// MaCh3 includes
#include "samplePDF/Structs.h"
#include "throwParms/ThrowParms.h"

// Don't forget yaml!
#include "yaml-cpp/yaml.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

#ifndef __LARGE_LOGL__
#define __LARGE_LOGL__ 1234567890.0
#endif

//#define DEBUG_PCA
#ifdef DEBUG_PCA
//KS: When debuging we produce some fancy plots, but we don't need it during normal work flow
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"
#include "TText.h"
#include "TLegend.h"

#if DEBUG_PCA == 2
#include "Eigen/Eigenvalues"
#endif

#endif


class covarianceBase {
 public:
  // The constructors
  covarianceBase(){};
  //ETA - construcotr for a YAML file
  covarianceBase(const char *YAMLFile);
  //"Usual" constructors from root file
  covarianceBase(const char *name, const char *file);
  covarianceBase(const char *name, const char *file, int seed);
  // For Eigen Value decomp
  covarianceBase(const char *name, const char *file, int seed, double threshold,int FirstPCAdpar, int LastPCAdpar);
  virtual ~covarianceBase();
  
  // Setters
  // need virtual for eRPA parameter over-ride
  // ETA - maybe need to add checks to index on the setters? i.e. if( i > _fPropVal.size()){throw;}
  void setCovMatrix(TMatrixDSym *cov);
  void setName(const char *name) { matrixName = name; }
  virtual void setParName(int i, char *name) { _fNames.at(i) = std::string(name); }
  void setSingleParameter(const int parNo, const double parVal);
  void setPar(const int i, const double val);
  void setParCurrProp(int i, double val);
  void setParProp(int i, double val) {
    std::cout << "Setting " << getParName(i) << " to " << val << std::endl; 
    _fPropVal[i] = val;
    if (pca) TransferToPCA();
  };
  void setParameters(std::vector<double> pars = std::vector<double>());    
  virtual void setEvalLikelihood(int i, bool eL);
  
  // set branches for output file
  void setBranches(TTree &tree);
  void setStepScale(double scale);

  //DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  void setIndivStepScale(int ParameterIndex, double StepScale){
	_fIndivStepScale.at(ParameterIndex) = StepScale;
  };

  void setIndivStepScale(std::vector<double> stepscale);
  //KS: In case someone really want to change this
  void setPrintLength(unsigned int PriLen) { PrintLength = PriLen; };

  // set a custom proposal function
  //DEPRECATED
  void setPropFunct(int i, TF1 *func) {};

  // Throwers
  void throwParProp(const double mag = 1.);
  void throwParCurr(const double mag = 1.);
  void throwParameters();
  void throwNominal(bool nomValues=false, int seed = 0);
  // Randomly throw the parameters in their 1 sigma range
  void RandomConfiguration();
  
  // Getters
  TMatrixDSym *getCovMatrix() { return covMatrix; };
  TMatrixDSym *getInvCovMatrix() { return invCovMatrix; };
  bool getEvalLikelihood(const int i) { return _fFlatPrior[i]; };

  virtual int CheckBounds();
  double CalcLikelihood();
  virtual double GetLikelihood();

  const char *getName() { return matrixName; };
  // ETA - Why is this virtual?
  virtual const char* getParName(const int i) const {
    return _fNames[i].c_str();
  };
  std::string const getInputFile() const { return inputFile; };

  // Get diagonal error for ith parameter
  const double getDiagonalError(const int i) { 
    return sqrt((*covMatrix)(i,i));
  }


  // Adaptive Step Tuning Stuff
  void resetIndivStepScale();
  void useSeparateThrowMatrix(TString throwMatrixName, TString throwMatrixFile, TString parameterMeansName="");
  void useSeparateThrowMatrix();
  void saveAdaptiveToFile(TString outFileName, TString systematicName);

  void setThrowMatrix(TMatrixDSym *cov);
  void updateThrowMatrix(TMatrixDSym *cov);
  void setNumberOfSteps(const int nsteps){
    total_steps=nsteps;
    if(total_steps>=lower_adapt)resetIndivStepScale();
  }
  // Set thresholds for MCMC steps
  void setAdaptiveThresholds(int low_threshold=10000, int up_threshold=1000000){lower_adapt=low_threshold; upper_adapt=up_threshold;}

  TMatrixDSym *getThrowMatrix(){return throwMatrix;}
  TMatrixD *getThrowMatrix_CholDecomp(){return throwMatrix_CholDecomp;}
  std::vector<double> getParameterMeans(){return par_means;}

  //========
  //DB Pointer return
  //ETA - This might be a bit squiffy? If the vector gots moved from say a
  //push_back then the pointer is no longer valid... maybe need a better 
  //way to deal with this? It was fine before when the return was to an 
  //element of a new array. There must be a clever C++ way to be careful
  //========
  const double* retPointer(int iParam) {return &(_fPropVal.data()[iParam]);}

  virtual std::vector<double> getNominalArray();
  const std::vector<double> getProposed() const;
  const double getParProp(const int i) {
    return _fPropVal[i]; 
  };
  const double getParCurr(const int i) {
    return _fCurrVal[i];
  };
  const double getParInit(const int i) {
    return _fPreFitValue[i];
  };

  virtual const double getNominal(const int i) {
    return getParInit(i);
  };

  double GetUpperBound(const int i){
    return _fUpBound[i];
  }

  double GetLowerBound(const int i){
    return _fLowBound[i];
  }

  const double getParProp_PCA(const int i) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return fParProp_PCA(i);
  };
  const double getParCurr_PCA(const int i) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return fParCurr_PCA(i);
  };

  bool isParameterFixedPCA(const int i) {
    if (fParSigma_PCA[i] < 0) {
      return true;
    } else {
      return false;
    }
  }

  const TMatrixD getTransferMatrix() {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return TransferMat;
  }

  const TMatrixD getEigenVectors() {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return eigen_vectors;
  }

  const TVectorD getEigenValues() {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return eigen_values;
  }

  inline const std::vector<double> getEigenValuesMaster() {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return eigen_values_master;
  }

  void setParProp_PCA(const int i, const double value) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    fParProp_PCA(i) = value;
    // And then transfer back to the parameter basis
    TransferToParam();
  }

  void setParCurr_PCA(const int i, const double value) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    fParCurr_PCA(i) = value;
    // And then transfer back to the parameter basis
    TransferToParam();
  }

  inline void setParameters_PCA(std::vector<double> pars)
  {
    if (!pca)
    {
      std::cerr<<" PCA disabled"<<std::endl;
      throw;
    }
    if (pars.size() != size_t(npars)) {
      std::cerr << "Warning: parameter arrays of incompatible size! Not changing parameters! " << matrixName << " has size " << pars.size() << " but was expecting " << size << std::endl;
      throw;
    }
    unsigned int parsSize = pars.size();
    for (unsigned int i = 0; i < parsSize; i++)
    {
      fParProp_PCA(i) = pars[i];
    }
    //KS: Transfer to normal base
    TransferToParam();
  }

  int getSize() { return size; };
  int getNpars() { 
    if (pca) return npars;
    else return size;
  }

  // Printers
  void printNominal();
  void printNominalCurrProp();
  void printPars();
  void printIndivStepScale();

  // Steppers
  virtual void proposeStep(); // generate a new proposed state
  void acceptStep(); // accepted this step

  // fix parameters at nominal values
  void toggleFixAllParameters();
  virtual void toggleFixParameter(const int i);
  bool isParameterFixed(const int i) {
    if (_fError[i] < 0) {
      return true;
    } else {
      return false;
    }
  }
  void ConstructPCA();
#ifdef DEBUG_PCA
  void DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat);
#endif
  // is PCA, can use to query e.g. LLH scans
  bool IsPCA() { return pca; };

  double* MatrixMult(double*, double*, int);
  double** MatrixMult(double**, double**, int);
  TMatrixD MatrixMult(TMatrixD, TMatrixD);
  inline void MatrixVectorMulti(double* VecMulti, double** matrix, const double* vector, const int n);
  inline double MatrixVectorMultiSingle(double** matrix, const double* vector, const int Length, const int i);

  //Turn on/off true adaptive MCMC
  //Also set thresholds for use (having a lower threshold gives us some data to adapt from!)
  void enableAdaptiveMCMC(bool enable=true){
    use_adaptive=enable;
    total_steps=0; //Set these to default values
    lower_adapt=10000;
    upper_adapt=10000000;
  }
  void updateAdaptiveCovariance();

 protected:
  void init(const char *name, const char *file);
  //YAML init
  void init(const char *YAMLFile);
  void init(TMatrixDSym* covMat);

  void MakePosDef(TMatrixDSym *cov = NULL);
  void makeClosestPosDef(TMatrixDSym *cov);
  void TransferToPCA();
  void TransferToParam();

  //Handy function to return 1 for any systs
  const double* ReturnUnity(){return &Unity;}

  // The input root file we read in
  const std::string inputFile;

  int size;
  // Name of cov matrix
  const char *matrixName;
  // The covariance matrix
  TMatrixDSym *covMatrix;
  // The inverse covariance matrix
  TMatrixDSym *invCovMatrix;
  //KS: Same as above but much faster as TMatrixDSym cache miss
  double **InvertCovMatrix;
    
  //KS: set Random numbers for each thread so each thread has differnt seed
  TRandom3 **random_number;

  // For Cholesky decomposed parameter throw
  double* randParams;
  //  TMatrixD *chel;
  double _fGlobalStepScale;

  //KS: This is used when printing parameters, sometimes we have super long parmaeters name, we want to flexibly adjust couts
  unsigned int PrintLength;

  Double_t currLogL;
  Double_t propLogL;
  Double_t *corr_throw;

  // state info (now mostly vectors)
  //ETA - duplication of some of these
  //ideally these should all be private and we have setters be protected 
  //setters and public getters
  std::vector<std::string> _fNames;
  int _fNumPar;
  YAML::Node _fYAMLDoc;
  std::vector<double> _fPreFitValue;
  std::vector<double> _fCurrVal;
  std::vector<double> _fPropVal;
  std::vector<double> _fGenerated;
  std::vector<double> _fError;
  std::vector<double> _fLowBound;
  std::vector<double> _fUpBound;
  std::vector<int> _fDetID;
  std::vector<std::string> _fDetString;
  std::vector<std::string> _fParamType;
  std::vector<double> _fIndivStepScale;
  std::vector<bool> _fFlatPrior;
  TMatrixDSym *_fCovMatrix;
  //TMatrixT<double> *_fCovMatrix;

  //Some "usual" variables. Don't think we really need the ND/FD split
  std::vector<int> _fNormModes;
  std::vector<std::string> _fNDSplineNames;
  std::vector<std::string> _fFDSplineNames;
  std::vector<std::vector<int>> _fFDSplineModes;

  //Information to be able to apply generic cuts
  std::vector<std::vector<std::string>> _fKinematicPars;
  std::vector<std::vector<std::vector<double>>> _fKinematicBounds;

  //Unity for null systs to point back to
  const double Unity = 1.0;

  // PCA
  bool pca;
  double eigen_threshold;
  int npars;
  int FirstPCAdpar;
  int LastPCAdpar;
  int nKeptPCApars;
  TVectorD eigen_values;
  TMatrixD eigen_vectors;
  std::vector<double> eigen_values_master;
  TMatrixD TransferMat;
  TMatrixD TransferMatT;
  // Also save current and proposed parameter in PCA basis
  TVectorD fParProp_PCA;
  TVectorD fParCurr_PCA;
  std::vector<double> fParSigma_PCA;
  std::vector<int> isDecomposed_PCA;

  TMatrixDSym* throwMatrix;
  TMatrixD* throwMatrix_CholDecomp;
  //Same as above but much faster as TMatrixDSym cache miss
  double **throwMatrixCholDecomp;

  void randomize();
  void CorrelateSteps();

  // Truely Adaptive Stuff
  void initialiseNewAdaptiveChain();
  // TMatrixD* getMatrixSqrt(TMatrixDSym* inputMatrix);
  // double calculateSubmodality(TMatrixD* sqrtVectorCov, TMatrixDSym* throwCov);
  bool use_adaptive;
  int total_steps;
  int lower_adapt, upper_adapt; //Thresholds for when to turn on/off adaptive MCMC
  // TMatrixD* covSqrt;
  std::vector<double> par_means;
  std::vector<double> par_means_prev;
  TMatrixDSym* adaptiveCovariance;
};

TH2D* TMatrixIntoTH2D(const TMatrix &Matrix, std::string title);
#endif
