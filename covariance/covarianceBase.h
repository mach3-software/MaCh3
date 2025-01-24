#pragma once

// MaCh3 includes
#include "manager/manager.h"
#include "samplePDF/Structs.h"
#include "covariance/CovarianceUtils.h"
#include "covariance/AdaptiveMCMCHandler.h"
#include "covariance/PCAHandler.h"

/// Large Likelihood is used it parameter go out of physical boundary, this indicates in MCMC that such step should eb removed
constexpr static const double _LARGE_LOGL_ = 1234567890.0;

/// @brief Base class responsible for handling of systematic error parameters. Capable of using PCA or using adaptive throw matrix
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class covarianceBase {
 public:
  /// @brief ETA - constructor for a YAML file
  /// @param YAMLFile A vector of strings representing the YAML files used for initialisation of matrix
  /// @param name Matrix name
  /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  covarianceBase(const std::vector<std::string>& YAMLFile, std::string name, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief "Usual" constructors from root file
  /// @param name Matrix name
  /// @param file Path to matrix root file
  covarianceBase(std::string name, std::string file, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);

  /// @brief Destructor
  virtual ~covarianceBase();
  
  /// \defgroup Setters
  /// Group of functions to set various parameters, names, and values.

  // ETA - maybe need to add checks to index on the setters? i.e. if( i > _fPropVal.size()){throw;}
  /// @brief Set covariance matrix
  /// @param cov Covariance matrix which we set and will be used later for evaluation of penalty term
  /// \ingroup Setters
  void setCovMatrix(TMatrixDSym *cov);
  /// @brief Set matrix name
  void setName(std::string name) { matrixName = name; }
  /// @brief change parameter name
  /// @param i Parameter index
  /// @param name new name which will be set
  /// \ingroup Setters
  void setParName(int i, std::string name) { _fNames.at(i) = name; }
  /// \ingroup Setters
  void setSingleParameter(const int parNo, const double parVal);
  /// @brief Set all the covariance matrix parameters to a user-defined value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// \ingroup Setters
  void setPar(const int i, const double val);
  /// @brief Set current parameter value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// \ingroup Setters
  void setParCurrProp(const int i, const double val);
  /// @brief Set proposed parameter value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// \ingroup Setters
  void setParProp(const int i, const double val) {
    _fPropVal[i] = val;
    if (pca) TransferToPCA();
  }
  /// @brief Set parameter values using vector, it has to have same size as covariance class
  /// @param pars Vector holding new values for every parameter
  /// \ingroup Setters
  void setParameters(const std::vector<double>& pars = {});
  /// @brief Set if parameter should have flat prior or not
  /// @param i Parameter index
  /// @param eL bool telling if it will be flat or not
  /// \ingroup Setters
  void setFlatPrior(const int i, const bool eL);
  
  /// @brief Set random value useful for debugging/CI
  /// @param i Parameter index
  /// @param rand New value for random number
  /// \ingroup Setters
  void SetRandomThrow(const int i, const double rand) { randParams[i] = rand;}
  /// @brief Get random value useful for debugging/CI
  /// @param i Parameter index
  double GetRandomThrow(const int i) { return randParams[i];}

  /// @brief set branches for output file
  /// @param tree Tree to which we will save branches
  /// @param SaveProposal Normally we only save parameter after is accepted, for debugging purpose it is helpful to see also proposed values. That's what this variable controls
  /// \ingroup Setters
  void SetBranches(TTree &tree, bool SaveProposal = false);
  /// @brief Set global step scale for covariance object
  /// @param scale Value of global step scale
  /// @cite luengo2020survey
  /// \ingroup Setters
  void setStepScale(const double scale);
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  /// @param ParameterIndex Parameter Index
  /// @param StepScale Value of individual step scale
  /// \ingroup Setters
  void setIndivStepScale(const int ParameterIndex, const double StepScale){ _fIndivStepScale.at(ParameterIndex) = StepScale; }
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  /// @param stepscale Vector of individual step scale, should have same
  /// \ingroup Setters
  void setIndivStepScale(const std::vector<double>& stepscale);
  /// @brief KS: In case someone really want to change this
  /// \ingroup Setters
  inline void setPrintLength(const unsigned int PriLen) { PrintLength = PriLen; }

  /// @brief KS: After step scale, prefit etc. value were modified save this modified config.
  void SaveUpdatedMatrixConfig();

  /// @brief Throw the proposed parameter by mag sigma. Should really just have the user specify this throw by having argument double
  void throwParProp(const double mag = 1.);

  /// @brief Helper function to throw the current parameter by mag sigma. Can study bias in MCMC with this; put different starting parameters
  void throwParCurr(const double mag = 1.);
  /// @brief Throw the parameters according to the covariance matrix. This shouldn't be used in MCMC code ase it can break Detailed Balance;
  void throwParameters();
  /// @brief Randomly throw the parameters in their 1 sigma range
  void RandomConfiguration();
  
  /// @brief Check if parameters were proposed outside physical boundary
  int CheckBounds() const;
  /// @brief Calc penalty term based on inverted covariance matrix
  double CalcLikelihood() _noexcept_;
  /// @brief Return CalcLikelihood if some params were thrown out of boundary return _LARGE_LOGL_
  virtual double GetLikelihood();

  /// @brief Return covariance matrix
  TMatrixDSym *getCovMatrix() { return covMatrix; }
  /// @brief Return inverted covariance matrix
  TMatrixDSym *getInvCovMatrix() { return invCovMatrix; }
  /// @brief Return inverted covariance matrix
  double GetInvCovMatrix(const int i, const int j) { return InvertCovMatrix[i][j]; }

  /// @brief Return correlated throws
  /// @param i Parameter index
  double GetCorrThrows(const int i) { return corr_throw[i]; }

  /// @brief Get if param has flat prior or not
  /// @param i Parameter index
  inline bool getFlatPrior(const int i) { return _fFlatPrior[i]; }

  /// @brief Get name of covariance
  /// \ingroup Setters
  std::string getName() const { return matrixName; }
  /// @brief Get name of covariance
  /// @param i Parameter index
  /// \ingroup Setters
  std::string GetParName(const int i) const {return _fNames[i];}
  /// @brief Get fancy name of the Parameter
  /// @param i Parameter index
  /// \ingroup Setters
  std::string GetParFancyName(const int i) const {return _fFancyNames[i];}
  /// @brief Get name of input file
  /// \ingroup Setters
  std::string getInputFile() const { return inputFile; }

  /// @brief Get diagonal error for ith parameter
  /// @param i Parameter index
  /// \ingroup Setters
  inline double getDiagonalError(const int i) const { return std::sqrt((*covMatrix)(i,i)); }
  /// @brief Get the error for the ith parameter
  /// @param i Parameter index
  /// \ingroup Setters
  inline double GetError(const int i) const {return _fError[i];}

  /// @brief Adaptive Step Tuning Stuff
  void resetIndivStepScale();

  /// @brief Initialise adaptive MCMC
  /// @param adapt_manager Node having from which we load all adaptation options
  void initialiseAdaption(const YAML::Node& adapt_manager);
  /// @brief Save adaptive throw matrix to file
  void saveAdaptiveToFile(const TString& outFileName, const TString& systematicName) {
    AdaptiveHandler.SaveAdaptiveToFile(outFileName, systematicName); }

  /// @brief Do we adapt or not
  bool getDoAdaption(){return use_adaptive;}
  /// @brief Use new throw matrix, used in adaptive MCMC
  void setThrowMatrix(TMatrixDSym *cov);
  void updateThrowMatrix(TMatrixDSym *cov);
  /// @brief Set number of MCMC step, when running adaptive MCMC it is updated with given frequency. We need number of steps to determine frequency.
  inline void setNumberOfSteps(const int nsteps) {
    AdaptiveHandler.total_steps = nsteps;
    if(AdaptiveHandler.AdaptionUpdate()) resetIndivStepScale();
  }

  /// @brief Get matrix used for step proposal
  inline TMatrixDSym *getThrowMatrix(){return throwMatrix;}
  /// @brief Get matrix used for step proposal
  double GetThrowMatrix(int i, int j) { return throwMatrixCholDecomp[i][j];}
  /// @brief Get the Cholesky decomposition of the throw matrix
  inline TMatrixD *getThrowMatrix_CholDecomp(){return throwMatrix_CholDecomp;}
  /// @brief Get the parameter means used in the adaptive handler
  inline std::vector<double> getParameterMeans(){return AdaptiveHandler.par_means;}
  /// @brief KS: Convert covariance matrix to correlation matrix and return TH2D which can be used for fancy plotting
  /// @details This function converts the covariance matrix to a correlation matrix and
  ///          returns a TH2D object, which can be used for advanced plotting purposes.
  /// @return A pointer to a TH2D object representing the correlation matrix
  TH2D* GetCorrelationMatrix();

  /// @brief DB Pointer return to param position
  ///
  /// @param iParam The index of the parameter in the vector.
  /// @return A pointer to the parameter value at the specified index.
  ///
  /// @warning ETA - This might be a bit squiffy? If the vector gots moved from say a
  /// push_back then the pointer is no longer valid... maybe need a better
  /// way to deal with this? It was fine before when the return was to an
  /// element of a new array. There must be a clever C++ way to be careful
  inline const double* retPointer(const int iParam) {return &(_fPropVal.data()[iParam]);}

  /// @brief Get a reference to the proposed parameter values
  /// Can be useful if you want to track these without having to copy values using getProposed()
  inline const std::vector<double> &getParPropVec() {return _fPropVal;}

  //Some Getters
  /// @brief Get total number of parameters
  inline int  GetNumParams() const {return _fNumPar;}
  /// @brief Get the prior array for parameters.
  virtual std::vector<double> getNominalArray();
  /// @brief Get the pre-fit values of the parameters.
  std::vector<double> getPreFitValues() const {return _fPreFitValue;}
  /// @brief Get the generated values of the parameters.
  std::vector<double> getGeneratedValues() const {return _fGenerated;}
  /// @brief Get vector of all proposed parameter values
  std::vector<double> getProposed() const;
  /// @brief Get proposed parameter value
  /// @param i Parameter index
  inline double getParProp(const int i) const { return _fPropVal[i]; }
  /// @brief Get current parameter value
  /// @param i Parameter index
  inline double getParCurr(const int i) const { return _fCurrVal[i]; }
  /// @brief Get prior parameter value
  /// @param i Parameter index
  inline double getParInit(const int i) const { return _fPreFitValue[i]; }
  /// @brief Return generated value, although is virtual so class inheriting might actual get prior not generated.
  /// @param i Parameter index
  virtual double getNominal(const int i) { return getParInit(i); }
  /// @brief Return generated value for a given parameter
  /// @param i Parameter index
  inline double GetGenerated(const int i) const { return _fGenerated[i];}
  /// @brief Get upper parameter bound in which it is physically valid
  /// @param i Parameter index
  inline double GetUpperBound(const int i) const { return _fUpBound[i];}
  /// @brief Get lower parameter bound in which it is physically valid
  /// @param i Parameter index
  inline double GetLowerBound(const int i) const { return _fLowBound[i]; }
  /// @brief Get individual step scale for selected parameter
  /// @param ParameterIndex Parameter index
  inline double GetIndivStepScale(const int ParameterIndex) const {return _fIndivStepScale.at(ParameterIndex); }
  /// @brief Get global step scale for covariance object
  inline double GetGlobalStepScale() const {return _fGlobalStepScale; }
  /// @brief Get current parameter value using PCA
  /// @param i Parameter index
  inline double getParProp_PCA(const int i)  {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return fParProp_PCA(i);
  }
  
  /// @brief Get current parameter value using PCA
  /// @param i Parameter index
  inline double getParCurr_PCA(const int i) {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return fParCurr_PCA(i);
  }

  /// @brief Is parameter fixed in PCA base or not
  /// @param i Parameter index
  inline bool isParameterFixedPCA(const int i) {
    if (fParSigma_PCA[i] < 0) { return true;  }
    else                      { return false; }
  }
  /// @brief Get transfer matrix allowing to go from PCA base to normal base
  inline const TMatrixD getTransferMatrix() {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return PCAObj.TransferMat;
  }
  /// @brief Get eigen vectors of covariance matrix, only works with PCA
  inline const TMatrixD getEigenVectors() {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return PCAObj.eigen_vectors;
  }
  /// @brief Get eigen values for all parameters, if you want for decomposed only parameters use getEigenValuesMaster
  inline const TVectorD getEigenValues() {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return PCAObj.eigen_values;
  }
  /// @brief Get eigen value of only decomposed parameters, if you want for all parameters use getEigenValues
  inline const std::vector<double> getEigenValuesMaster() {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    return PCAObj.eigen_values_master;
  }
  /// @brief Set proposed value for parameter in PCA base
  /// @param i Parameter index
  /// @param value new value
  inline void setParProp_PCA(const int i, const double value) {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    fParProp_PCA(i) = value;
    // And then transfer back to the parameter basis
    TransferToParam();
  }
  /// @brief Set current value for parameter in PCA base
  /// @param i Parameter index
  /// @param value new value
  inline void setParCurr_PCA(const int i, const double value) {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    fParCurr_PCA(i) = value;
    // And then transfer back to the parameter basis
    TransferToParam();
  }

  /// @brief Set values for PCA parameters in PCA base
  /// @param pars vector with new values of PCA params
  inline void setParameters_PCA(const std::vector<double> &pars) {
    if (!pca) { MACH3LOG_ERROR("Am not running in PCA mode"); throw MaCh3Exception(__FILE__ , __LINE__ ); }
    if (int(pars.size()) != _fNumParPCA) {
      MACH3LOG_ERROR("Warning: parameter arrays of incompatible size! Not changing parameters! {} has size {} but was expecting {}", matrixName, pars.size(), _fNumPar);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    int parsSize = int(pars.size());
    for (int i = 0; i < parsSize; i++) {
      fParProp_PCA(i) = pars[i];
    }
    //KS: Transfer to normal base
    TransferToParam();
  }

  /// @brief Get number of params which will be different depending if using Eigen decomposition or not
  inline int getNpars() {
    if (pca) return _fNumParPCA;
    else return _fNumPar;
  }

  /// @brief Print prior value for every parameter
  void printNominal() const;
  /// @brief Print prior, current and proposed value for each parameter
  void printNominalCurrProp() const;
  /// @warning only for backward compatibility
  /// @todo remove it
  void printPars() const {printNominalCurrProp();};
  /// @brief Print step scale for each parameter
  void printIndivStepScale() const;

  /// @brief Generate a new proposed state
  virtual void proposeStep();
  /// @brief "Randomize" the parameters in the covariance class for the proposed step. Used the proposal kernel and the current parameter value to set proposed step
  void randomize() _noexcept_;
  /// @brief Use Cholesky throw matrix for better step proposal
  void CorrelateSteps() _noexcept_;

  /// @brief Accepted this step
  void acceptStep() _noexcept_;

  /// @brief fix parameters at prior values
  void toggleFixAllParameters();
  /// @brief fix parameter at prior values
  /// @param i Parameter index
  void toggleFixParameter(const int i);
  /// @brief Fix parameter at prior values
  /// @param name Name of parameter you want to fix
  void toggleFixParameter(const std::string& name);

  /// @brief Is parameter fixed or not
  /// @param i Parameter index
  bool isParameterFixed(const int i) const {
    if (_fError[i] < 0) { return true; }
    else                { return false; }
  }
  /// @brief Is parameter fixed or not
  /// @param name Name of parameter you want to check if is fixed
  bool isParameterFixed(const std::string& name) const;

  /// @brief CW: Calculate eigen values, prepare transition matrices and remove param based on defined threshold
  /// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/03.-Eigen-Decomposition-%E2%80%90-PCA).
  void ConstructPCA();

  /// @brief is PCA, can use to query e.g. LLH scans
  inline bool IsPCA() { return pca; }

  /// @brief KS: Custom function to perform multiplication of matrix and vector with multithreading
  /// @param VecMulti Output Vector, VecMulti = matrix x vector
  /// @param matrix This matrix is used for multiplication VecMulti = matrix x vector
  /// @param VecMulti This vector is used for multiplication VecMulti = matrix x vector
  /// @param n this is size of matrix and vector, we assume matrix is symmetric
  inline void MatrixVectorMulti(double* _restrict_ VecMulti, double** _restrict_ matrix, const double* _restrict_ vector, const int n) const;
  /// @brief KS: Custom function to perform multiplication of matrix and single element which is thread safe
  inline double MatrixVectorMultiSingle(double** _restrict_ matrix, const double* _restrict_ vector, const int Length, const int i) const;

  /// @brief Getter to return a copy of the YAML node
  YAML::Node GetConfig() const { return _fYAMLDoc; }

protected:
  /// @brief Initialisation of the class using matrix from root file
  void init(std::string name, std::string file);
  /// @brief Initialisation of the class using config
  /// @param YAMLFile A vector of strings representing the YAML files used for initialisation of matrix
  void init(const std::vector<std::string>& YAMLFile);
  /// @brief Initialise vectors with parameters information
  /// @param size integer telling size to which we will resize all vectors/allocate memory
  void ReserveMemory(const int size);

  /// @brief Make matrix positive definite by adding small values to diagonal, necessary for inverting matrix
  /// @param cov Matrix which we evaluate Positive Definitiveness
  void MakePosDef(TMatrixDSym *cov = nullptr);

  /// @brief HW: Finds closest possible positive definite matrix in Frobenius Norm ||.||_frob Where ||X||_frob=sqrt[sum_ij(x_ij^2)] (basically just turns an n,n matrix into vector in n^2 space then does Euclidean norm)
  void makeClosestPosDef(TMatrixDSym *cov);
  /// @brief Transfer param values from normal base to PCA base
  void TransferToPCA();
  /// @brief Transfer param values from PCA base to normal base
  void TransferToParam();

  /// @brief sets throw matrix from a file
  /// @param matrix_file_name name of file matrix lives in
  /// @param matrix_name name of matrix in file
  /// @param means_name name of means vec in file
  void setThrowMatrixFromFile(const std::string& matrix_file_name, const std::string& matrix_name, const std::string& means_name);

  /// @brief Method to update adaptive MCMC
  /// @cite haario2001adaptive
  void updateAdaptiveCovariance();

  /// @brief Check if parameter is affecting given det ID
  /// @param SystIndex number of parameter
  /// @param DetID The Detector ID used to filter parameters.
  bool AppliesToDetID(const int SystIndex, const std::string& DetID) const;

  /// The input root file we read in
  const std::string inputFile;

  /// Name of cov matrix
  std::string matrixName;
  /// The covariance matrix
  TMatrixDSym *covMatrix;
  /// The inverse covariance matrix
  TMatrixDSym *invCovMatrix;
  /// KS: Same as above but much faster as TMatrixDSym cache miss
  std::vector<std::vector<double>> InvertCovMatrix;
    
  /// KS: Set Random numbers for each thread so each thread has different seed
  std::vector<std::unique_ptr<TRandom3>> random_number;

  /// Random number taken from gaussian around prior error used for corr_throw
  double* randParams;
  /// Result of multiplication of Cholesky matrix and randParams
  double* corr_throw;
  /// Global step scale applied to all params in this class
  double _fGlobalStepScale;

  /// KS: This is used when printing parameters, sometimes we have super long parameters name, we want to flexibly adjust couts
  int PrintLength;

  /// ETA _fNames is set automatically in the covariance class to be something like xsec_i, this is currently to make things compatible with the Diagnostic tools
  std::vector<std::string> _fNames;
  /// Fancy name for example rather than xsec_0 it is MAQE, useful for human reading
  std::vector<std::string> _fFancyNames;
  /// Stores config describing systematics
  YAML::Node _fYAMLDoc;
  /// Number of systematic parameters
  int _fNumPar;
  /// Parameter value dictated by the prior model. Based on it penalty term is calculated
  std::vector<double> _fPreFitValue;
  /// Current value of the parameter
  std::vector<double> _fCurrVal;
  /// Proposed value of the parameter
  std::vector<double> _fPropVal;
  /// Generated value of the parameter
  std::vector<double> _fGenerated;
  /// Prior error on the parameter
  std::vector<double> _fError;
  /// Lowest physical bound, parameter will not be able to go beyond it
  std::vector<double> _fLowBound;
  /// Upper physical bound, parameter will not be able to go beyond it
  std::vector<double> _fUpBound;
  /// Individual step scale used by MCMC algorithm
  std::vector<double> _fIndivStepScale;
  /// Whether to apply flat prior or not
  std::vector<bool> _fFlatPrior;
  /// Tells to which samples object param should be applied
  std::vector<std::vector<std::string>> _fDetID;

  /// perform PCA or not
  bool pca;
  /// CW: Threshold based on which we remove parameters in eigen base
  double eigen_threshold;
  /// Number of parameters in PCA base
  int _fNumParPCA;
  /// Index of the first param that is being decomposed
  int FirstPCAdpar;
  /// Index of the last param that is being decomposed
  int LastPCAdpar;

  /// Prefit value for PCA params
  std::vector<double> _fPreFitValue_PCA;
  /// CW: Current parameter value in PCA base
  TVectorD fParProp_PCA;
  /// CW: Proposed parameter value in PCA base
  TVectorD fParCurr_PCA;
  /// Tells if parameter is fixed in PCA base or not
  std::vector<double> fParSigma_PCA;
  /// If param is decomposed this will return -1, if not this will return enumerator to param in normal base. This way we can map stuff like step scale etc between normal base and undecomposed param in eigen base.
  std::vector<int> isDecomposed_PCA;

  /// Matrix which we use for step proposal before Cholesky decomposition (not actually used for step proposal)
  TMatrixDSym* throwMatrix;
  /// Matrix which we use for step proposal after Cholesky decomposition
  TMatrixD* throwMatrix_CholDecomp;
  /// Throw matrix that is being used in the fit, much faster as TMatrixDSym cache miss
  double** throwMatrixCholDecomp;

  /// Are we using AMCMC?
  bool use_adaptive;

  /// Struct containing information about PCA
  PCAHandler PCAObj;
  /// Struct containing information about adaption
  adaptive_mcmc::AdaptiveMCMCHandler AdaptiveHandler;
};
