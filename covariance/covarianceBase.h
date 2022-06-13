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

// MaCh3 includes
#include "samplePDF/Structs.h"
#include "throwParms/ThrowParms.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

#ifndef __LARGE_LOGL__
#define __LARGE_LOGL__ 1234567890
#endif

class covarianceBase {
 public:
  // The constructors
  covarianceBase();
  covarianceBase(const char *name, const char *file);
  covarianceBase(const char *name, const char *file, int seed);
  // For Eigen Value decomp
  covarianceBase(const char *name, const char *file, int seed, double threshold,int firstpcapar, int lastpcapar);
  virtual ~covarianceBase() {};
  
  // Setters
  // need virtual for eRPA parameter over-ride
  void setCovMatrix(TMatrixDSym *cov);
  void setName(const char *name) { matrixName = name; }
  virtual void setParName(int i, char *name) { fParNames[i] = name; }
  void setSingleParameter(int parNo, double parVal);
  void setPar(int i, double val);
  void setParCurrProp(int i, double val);
  void setParProp(int i, double val) { 
    std::cout << "Setting " << getParName(i) << " to " << val << std::endl; 
    fParProp[i] = val;
    if (pca) TransferToPCA();
  };
  void setParameters(std::vector<double> pars = std::vector<double>());    
  virtual void setEvalLikelihood(int i, bool eL);
  void setDoStep(int i, bool doStep) { fParDoStep[i] = doStep; }
  
  // set branches for output file
  void setBranches(TTree &tree);
  void setStepScale(double scale);

  //DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  void setIndivStepScale(std::vector<double> stepscale);

  // set a custom proposal function
  void setPropFunct(int i, TF1 *func);

  virtual void SetBANFFCov();

  // Throwers
  void throwParProp();
  void throwParCurr();
  void throwParameters();
  virtual void throwNominal(bool nomValues=true, int seed = 0);
  // Randomly throw the parameters in their 1 sigma range
  void RandomConfiguration();
  
  // Getters
  TMatrixDSym *getCovMatrix() { return covMatrix; };
  TMatrixDSym *getInvCovMatrix() { return invCovMatrix; };
  bool getEvalLikelihood(int i) { return fParEvalLikelihood[i]; };
  bool getDoStep(int i) { return fParDoStep[i]; };

  virtual double getLikelihood();

  const char *getName() { return matrixName; };
  virtual const char* getParName(int i) const { 
    return fParNames[i];
  };
  std::string const getInputFile() const { return inputFile; };

  // Get diagonal error for ith parameter
  const double getDiagonalError(const int i) { 
    return sqrt((*covMatrix)(i,i));
  }

  //========
  //DB Pointer return
  //========
  const double* retPointer(int iParam) {return &fParProp[iParam];}

  const std::vector<double> & getNominalArray() const { return nominal; };
  const std::vector<double> getProposed() const;
  const double getNominal(int i) { 
    return nominal[i]; };
  const double getParProp(int i) { 
    return fParProp[i]; 
  };
  const double getParCurr(int i) { 
    return fParCurr[i];
  };
  const double getParInit(int i) { 
    return fParInit[i];
  };
  const double getParProp_PCA(int i) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return fParProp_PCA(i);
  };
  const double getParCurr_PCA(int i) {
    if (!pca) {
      std::cerr << "Am not running in PCA mode" << std::endl;
      throw;
    }
    return fParCurr_PCA(i);
  };

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

  const int getSize() { return size; };
  const int getNpars() { 
    if (pca) return npars;
    else return size;
  }
  TF1 *getPropFunct(int i) { return fPropKernel[i]; };

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
  virtual void toggleFixParameter(int i);
  bool isParameterFixed(int i) {
    if (fParSigma[i] < 0) {
      return true;
    } else {
      return false;
    }
  }
  void ConstructPCA();

  // is PCA, can use to query e.g. LLH scans
  bool IsPCA() { return pca; };

  double* MatrixMult(double*, double*, int);
  double** MatrixMult(double**, double**, int);
  TMatrixD MatrixMult(TMatrixD, TMatrixD);

 protected:
  void init(const char *name, const char *file);
  void init(TMatrixDSym* covMat);

  void MakePosDef();

  void TransferToPCA();
  void TransferToParam();

  void UpdateCholeskyMatrix();

  // The input root file we read in
  const std::string inputFile;

  int size;
  // Name of cov matrix
  const char *matrixName;
  // The covariance matrix
  TMatrixDSym *covMatrix;
  // The inverse covariance matrix
  TMatrixDSym *invCovMatrix;
  // The nominal
  std::vector<double> nominal;
    
  // Random numbers
  TRandom3 *random_number;
  // For Cholesky decomposed parameter throw
  TVectorD randParams;
  TMatrixD *chel;
  double fStepScale;

  // proposal kernels
  TF1 **fPropKernel;

  // state info (arrays of each parameter)
  Char_t  **fParNames;
  Double_t *fParInit;
  Double_t *fParCurr;
  Double_t *fParProp;
  Double_t *fParSigma;
  Double_t *fParLoLimit;
  Double_t *fParHiLimit;
  Double_t *fIndivStepScale;
  bool     *fParEvalLikelihood;
  bool     *fParDoStep;    
  Double_t currLogL;
  Double_t propLogL;

  // posteriors (binned for now)
  void genPosteriorHists();
  TH1D **posterior;

  // PCA
  bool pca;
  double eigen_threshold;
  int npars;
  int firstpcadpar;
  int lastpcadpar;
  int nkeptpcapars;
  TVectorD eigen_values;
  TMatrixD eigen_vectors;
  std::vector<double> eigen_values_master;
  TMatrixD TransferMat;
  TMatrixD TransferMatT;
  // Also save current and proposed parameter in PCA basis
  TVectorD fParProp_PCA;
  TVectorD fParCurr_PCA;

  void CorrelateSteps();
  void genPropKernels();
  void CholeskyDecomp(int npars, TMatrixD &chel_mat);
  void randomize(); 
};

#endif
