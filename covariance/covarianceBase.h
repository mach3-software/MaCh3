#pragma once

// MaCh3 includes
#include "samplePDF/Structs.h"
#include "covariance/CovarianceUtils.h"
#include "covariance/ThrowParms.h"

// Don't forget yaml!
#include "yaml-cpp/yaml.h"
#include "manager/MaCh3Logger.h"

#ifndef __LARGE_LOGL__
#define __LARGE_LOGL__ 1234567890.0
#endif

//#define DEBUG_PCA 1
#ifdef DEBUG_PCA
//KS: When debugging we produce some fancy plots, but we don't need it during normal work flow
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
  //ETA - constructor for a YAML file
  covarianceBase(std::vector<std::string> YAMLFile, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
  //"Usual" constructors from root file
  covarianceBase(const char *name, const char *file);
  covarianceBase(const char *name, const char *file, int seed);
  // For Eigen Value decomp
  covarianceBase(const char *name, const char *file, int seed, double threshold,int FirstPCAdpar, int LastPCAdpar);
  virtual ~covarianceBase();
  
  // Setters
  // ETA - maybe need to add checks to index on the setters? i.e. if( i > _fPropVal.size()){throw;}
  void setCovMatrix(TMatrixDSym *cov);
  void setName(const char *name) { matrixName = name; }
  void setParName(int i, char *name) { _fNames.at(i) = std::string(name); }
  void setSingleParameter(const int parNo, const double parVal);
  void setPar(const int i, const double val);
  void setParCurrProp(int i, double val);
  void setParProp(int i, double val) {
    _fPropVal[i] = val;
    if (pca) TransferToPCA();
  };
  void setParameters(std::vector<double> pars = std::vector<double>());    
  void setEvalLikelihood(int i, bool eL);
  
  // set branches for output file
  void setBranches(TTree &tree);
  void setStepScale(double scale);
  //DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  void setIndivStepScale(int ParameterIndex, double StepScale){ _fIndivStepScale.at(ParameterIndex) = StepScale; };
  void setIndivStepScale(std::vector<double> stepscale);
  //KS: In case someone really want to change this
  inline void setPrintLength(const unsigned int PriLen) { PrintLength = PriLen; };

  // Throwers
  void throwParProp(const double mag = 1.);
  void throwParCurr(const double mag = 1.);
  void throwParameters();
  void throwNominal(bool nomValues = false, int seed = 0);
  // Randomly throw the parameters in their 1 sigma range
  void RandomConfiguration();
  
  //LLH Related
  virtual int CheckBounds();
  double CalcLikelihood();
  virtual double GetLikelihood();

  // Getters
  TMatrixDSym *getCovMatrix() { return covMatrix; };
  TMatrixDSym *getInvCovMatrix() { return invCovMatrix; };
  bool getEvalLikelihood(const int i) { return _fFlatPrior[i]; };

  const char *getName() { return matrixName; };
  std::string GetParName(const int i) {return _fNames[i];};
  const char* GetParName(const int i) const { return _fNames[i].c_str(); };
  std::string GetParFancyName(const int i) {return _fFancyNames[i];};
  const char* GetParFancyName(const int i) const { return _fFancyNames[i].c_str(); };
  std::string const getInputFile() const { return inputFile; };

  // Get diagonal error for ith parameter
  inline double getDiagonalError(const int i) { return std::sqrt((*covMatrix)(i,i)); }

  // Adaptive Step Tuning Stuff
  void resetIndivStepScale();
  void useSeparateThrowMatrix(TString throwMatrixName, TString throwMatrixFile, TString parameterMeansName="");
  void useSeparateThrowMatrix();
  void saveAdaptiveToFile(TString outFileName, TString systematicName);

  void setThrowMatrix(TMatrixDSym *cov);
  void updateThrowMatrix(TMatrixDSym *cov);
  void setNumberOfSteps(const int nsteps){
    total_steps = nsteps;
    if(total_steps >= lower_adapt)resetIndivStepScale();
  }
  // Set thresholds for MCMC steps
  inline void setAdaptiveThresholds(const int low_threshold = 10000, const int up_threshold = 1000000){lower_adapt=low_threshold; upper_adapt = up_threshold;}

  inline TMatrixDSym *getThrowMatrix(){return throwMatrix;}
  inline TMatrixD *getThrowMatrix_CholDecomp(){return throwMatrix_CholDecomp;}
  inline std::vector<double> getParameterMeans(){return par_means;}
  TH2D* GetCorrelationMatrix();

  // What parameter Gets reweighted by what amount according to MCMC
  inline double calcReWeight(const int bin) {return _fPropVal[bin];}
  //========
  //DB Pointer return
  //ETA - This might be a bit squiffy? If the vector gots moved from say a
  //push_back then the pointer is no longer valid... maybe need a better 
  //way to deal with this? It was fine before when the return was to an 
  //element of a new array. There must be a clever C++ way to be careful
  //========
  const double* retPointer(int iParam) {return &(_fPropVal.data()[iParam]);}

  //Some Getters
  int    GetNumParams()               {return _fNumPar;}
  virtual std::vector<double> getNominalArray();
  const std::vector<double>& getPreFitValues(){return _fPreFitValue;}
  const std::vector<double>& getGeneratedValues(){return _fGenerated;}
  const std::vector<double> getProposed() const;
  inline double getParProp(const int i) { return _fPropVal[i]; }
  inline double getParCurr(const int i) { return _fCurrVal[i]; }
  inline double getParInit(const int i) { return _fPreFitValue[i]; }

  virtual double getNominal(const int i) { return getParInit(i); }
  inline double GetGenerated(const int i) { return _fGenerated[i];}
  inline double GetUpperBound(const int i){ return _fUpBound[i];}
  inline double GetLowerBound(const int i){ return _fLowBound[i]; }

  double getParProp_PCA(const int i) {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    return fParProp_PCA(i);
  };
  
  double getParCurr_PCA(const int i) {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
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
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    return TransferMat;
  }

  const TMatrixD getEigenVectors() {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    return eigen_vectors;
  }

  const TVectorD getEigenValues() {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    return eigen_values;
  }

  inline const std::vector<double> getEigenValuesMaster() {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    return eigen_values_master;
  }

  void setParProp_PCA(const int i, const double value) {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw;
    }
    fParProp_PCA(i) = value;
    // And then transfer back to the parameter basis
    TransferToParam();
  }

  void setParCurr_PCA(const int i, const double value) {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
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
      MACH3LOG_ERROR("Am not running in PCA mode");
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

  inline int getSize() { return size; };
  inline int getNpars() {
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
  void toggleFixParameter(const int i);
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

  inline void MatrixVectorMulti(double* VecMulti, double** matrix, const double* vector, const int n);
  inline double MatrixVectorMultiSingle(double** matrix, const double* vector, const int Length, const int i);

  //Turn on/off true adaptive MCMC
  //Also set thresholds for use (having a lower threshold gives us some data to adapt from!)
  void enableAdaptiveMCMC(bool enable = true){
    use_adaptive = enable;
    total_steps = 0; //Set these to default values
    lower_adapt = 10000;
    upper_adapt = 10000000;
  }
  void updateAdaptiveCovariance();

 protected:
  void init(const char *name, const char *file);
  //YAML init
  void init(std::vector<std::string> YAMLFile);
  void init(TMatrixDSym* covMat);
  void ReserveMemory(const int size);

  void randomize();
  void CorrelateSteps();

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
    
  //KS: set Random numbers for each thread so each thread has different seed
  TRandom3 **random_number;

  // For Cholesky decomposed parameter throw
  double* randParams;
  double* corr_throw;
  double _fGlobalStepScale;

  //KS: This is used when printing parameters, sometimes we have super long parameters name, we want to flexibly adjust couts
  unsigned int PrintLength;

  // state info (now mostly vectors)
  //ETA - duplication of some of these
  //ideally these should all be private and we have setters be protected 
  //setters and public getters
  //_fNames is set automatically in the covariance class to be something like xsec_i
  //this is currently to make things compatible with the Diagnostic tools
  std::vector<std::string> _fNames;
  std::vector<std::string> _fFancyNames;
  int _fNumPar;
  YAML::Node _fYAMLDoc;
  std::vector<double> _fPreFitValue;
  std::vector<double> _fCurrVal;
  std::vector<double> _fPropVal;
  std::vector<double> _fGenerated;
  std::vector<double> _fError;
  std::vector<double> _fLowBound;
  std::vector<double> _fUpBound;
  std::vector<double> _fIndivStepScale;
  std::vector<bool> _fFlatPrior;

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

  //Adaptive MCMC
  TMatrixDSym* throwMatrix;
  TMatrixD* throwMatrix_CholDecomp;
  //Same as above but much faster as TMatrixDSym cache miss
  double **throwMatrixCholDecomp;

  // Truly Adaptive Stuff
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
