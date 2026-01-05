#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"
#include "Parameters/AdaptiveMCMCHandler.h"
#include "Parameters/PCAHandler.h"
#include "Parameters/ParameterTunes.h"

/// @brief Base class responsible for handling of systematic error parameters. Capable of using PCA or using adaptive throw matrix
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class ParameterHandlerBase {
 public:
  /// @brief ETA - constructor for a YAML file
  /// @param YAMLFile A vector of strings representing the YAML files used for initialisation of matrix
  /// @param name Matrix name
  /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  ParameterHandlerBase(const std::vector<std::string>& YAMLFile, std::string name, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief "Usual" constructors from root file
  /// @param name Matrix name
  /// @param file Path to matrix root file
  ParameterHandlerBase(std::string name, std::string file, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);

  /// @brief Destructor
  virtual ~ParameterHandlerBase();
  
  /// @defgroup ParameterHandlerSetters Parameter Handler Setters
  /// Group of functions to set various parameters, names, and values.

  /// @defgroup ParameterHandlerGetters Parameter Handler Getters
  /// Group of functions to get various parameters, names, and values.

  // ETA - maybe need to add checks to index on the setters? i.e. if( i > _fPropVal.size()){throw;}
  /// @brief Set covariance matrix
  /// @param cov Covariance matrix which we set and will be used later for evaluation of penalty term
  /// @ingroup ParameterHandlerSetters
  void SetCovMatrix(TMatrixDSym *cov);
  /// @brief Set matrix name
  /// @ingroup ParameterHandlerSetters
  void SetName(const std::string& name) { matrixName = name; }
  /// @brief change parameter name
  /// @param i Parameter index
  /// @param name new name which will be set
  /// @ingroup ParameterHandlerSetters
  void SetParName(const int i, const std::string& name) { _fNames.at(i) = name; }
  /// @brief Set value of single param to a given value
  /// @ingroup ParameterHandlerSetters
  void SetSingleParameter(const int parNo, const double parVal);
  /// @brief Set all the covariance matrix parameters to a user-defined value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// @ingroup ParameterHandlerSetters
  void SetPar(const int i, const double val);
  /// @brief Set current parameter value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// @ingroup ParameterHandlerSetters
  void SetParCurrProp(const int i, const double val);
  /// @brief Set proposed parameter value
  /// @param i Parameter index
  /// @param val new value which will be set
  /// @ingroup ParameterHandlerSetters
  void SetParProp(const int i, const double val) {
    _fPropVal[i] = val;
    if (pca) PCAObj->TransferToPCA();
  }
  /// @brief Set parameter values using vector, it has to have same size as covariance class
  /// @param pars Vector holding new values for every parameter
  /// @ingroup ParameterHandlerSetters
  void SetParameters(const std::vector<double>& pars = {});
  /// @brief Set if parameter should have flat prior or not
  /// @param i Parameter index
  /// @param eL bool telling if it will be flat or not
  /// @ingroup ParameterHandlerSetters
  void SetFlatPrior(const int i, const bool eL);
  
  /// @brief Set random value useful for debugging/CI
  /// @param i Parameter index
  /// @param rand New value for random number
  /// @ingroup ParameterHandlerSetters
  void SetRandomThrow(const int i, const double rand) { randParams[i] = rand;}
  /// @brief Get random value useful for debugging/CI
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetRandomThrow(const int i) const { return randParams[i];}

  /// @brief set branches for output file
  /// @param tree Tree to which we will save branches
  /// @param SaveProposal Normally we only save parameter after is accepted, for debugging purpose it is helpful to see also proposed values. That's what this variable controls
  /// @ingroup ParameterHandlerSetters
  void SetBranches(TTree &tree, const bool SaveProposal = false);
  /// @brief Set global step scale for covariance object
  /// @param scale Value of global step scale
  /// @param verbose Print that we've changed scale + use warnings [default: true]
  /// @cite luengo2020survey
  /// @ingroup ParameterHandlerSetters
  void SetStepScale(const double scale, const bool verbose=true);
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  /// @param ParameterIndex Parameter Index
  /// @param StepScale Value of individual step scale
  /// @ingroup ParameterHandlerSetters
  void SetIndivStepScale(const int ParameterIndex, const double StepScale){ _fIndivStepScale.at(ParameterIndex) = StepScale; }
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from execs and inside covariance constructors)
  /// @param stepscale Vector of individual step scale, should have same
  /// @ingroup ParameterHandlerSetters
  void SetIndivStepScale(const std::vector<double>& stepscale);
  /// @brief KS: In case someone really want to change this
  /// @ingroup ParameterHandlerSetters
  inline void SetPrintLength(const unsigned int PriLen) { PrintLength = PriLen; }

  /// @brief KS: After step scale, prefit etc. value were modified save this modified config.
  void SaveUpdatedMatrixConfig();

  /// @brief Throw the parameters according to the covariance matrix. This shouldn't be used in MCMC code ase it can break Detailed Balance;
  void ThrowParameters();
  /// @brief Randomly throw the parameters in their 1 sigma range
  void RandomConfiguration();
  
  /// @brief Check if parameters were proposed outside physical boundary
  int CheckBounds() const _noexcept_;
  /// @brief Calc penalty term based on inverted covariance matrix
  ///
  /// @details
  /// The log-likelihood is computed as:
  /// \f[
  ///   \log \mathcal{L} = \frac{1}{2} \sum_{i}^{\textrm{pars}} \sum_{j}^{\textrm{pars}} \Delta \vec{p}_i \left( V^{-1} \right)_{i,j} \Delta \vec{p}_j
  /// \f]
  /// where:
  /// - \f$\Delta \vec{p}_i = \theta_i - \theta_{i,0}\f$ is the difference between the current and pre-fit parameter values,
  /// - \f$V^{-1}\f$ is the inverted covariance matrix.
  ///
  /// @note
  /// - If `_fFlatPrior[i]` is `true`, the parameter is excluded from the calculation.
  double CalcLikelihood() const _noexcept_;
  /// @brief Return CalcLikelihood if some params were thrown out of boundary return _LARGE_LOGL_
  /// @ingroup ParameterHandlerGetters
  virtual double GetLikelihood();

  /// @brief Return covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *GetCovMatrix() const { return covMatrix; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *GetInvCovMatrix() const { return invCovMatrix; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  double GetInvCovMatrix(const int i, const int j) const { return InvertCovMatrix[i][j]; }

  /// @brief Return correlated throws
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetCorrThrows(const int i) const { return corr_throw[i]; }

  /// @brief Get if param has flat prior or not
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline bool GetFlatPrior(const int i) const { return _fFlatPrior[i]; }

  /// @brief Get name of covariance
  /// @ingroup ParameterHandlerGetters
  std::string GetName() const { return matrixName; }
  /// @brief Get name of parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParName(const int i) const {return _fNames[i];}

  /// @brief Get index based on name
  /// @ingroup ParameterHandlerGetters
  int GetParIndex(const std::string& name) const;

  /// @brief Get fancy name of the Parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParFancyName(const int i) const {return _fFancyNames[i];}
  /// @brief Get name of input file
  /// @ingroup ParameterHandlerGetters
  std::string GetInputFile() const { return inputFile; }

  /// @brief Get diagonal error for ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetDiagonalError(const int i) const { return std::sqrt((*covMatrix)(i,i)); }
  /// @brief Get the error for the ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetError(const int i) const {return _fError[i];}

  /// @brief Adaptive Step Tuning Stuff
  void ResetIndivStepScale();

  /// @brief Set individual step scale for parameters which are skipped during adaption to initial values
  void SetIndivStepScaleForSkippedAdaptParams();

  /// @brief Initialise adaptive MCMC
  /// @param adapt_manager Node having from which we load all adaptation options
  void InitialiseAdaption(const YAML::Node& adapt_manager);
  /// @brief Save adaptive throw matrix to file
  void SaveAdaptiveToFile(const std::string& outFileName, const std::string& systematicName) {
    AdaptiveHandler->SaveAdaptiveToFile(outFileName, systematicName); }

  /// @brief Do we adapt or not
  /// @ingroup ParameterHandlerGetters
  bool GetDoAdaption() const {return use_adaptive;}
  /// @brief Use new throw matrix, used in adaptive MCMC
  /// @ingroup ParameterHandlerSetters
  void SetThrowMatrix(TMatrixDSym *cov);
  /// @brief Replaces old throw matrix with new one
  void UpdateThrowMatrix(TMatrixDSym *cov);
  /// @brief Set number of MCMC step, when running adaptive MCMC it is updated with given frequency. We need number of steps to determine frequency.
   /// @ingroup ParameterHandlerSetters
  inline void SetNumberOfSteps(const int nsteps) {
    AdaptiveHandler->SetTotalSteps(nsteps);
    if(AdaptiveHandler->AdaptionUpdate()) ResetIndivStepScale();
  }

  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  inline TMatrixDSym *GetThrowMatrix() const {return throwMatrix;}
  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  double GetThrowMatrix(const int i, const int j) const { return throwMatrixCholDecomp[i][j];}

  /// @brief KS: Convert covariance matrix to correlation matrix and return TH2D which can be used for fancy plotting
  /// @details This function converts the covariance matrix to a correlation matrix and
  ///          returns a TH2D object, which can be used for advanced plotting purposes.
  /// @return A pointer to a TH2D object representing the correlation matrix
  /// @ingroup ParameterHandlerGetters
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
  inline const double* RetPointer(const int iParam) {return &(_fPropVal.data()[iParam]);}

  /// @brief Get a reference to the proposed parameter values
  /// Can be useful if you want to track these without having to copy values using getProposed()
  inline const std::vector<double> &GetParPropVec() {return _fPropVal;}

  /// @brief Get total number of parameters
  /// @ingroup ParameterHandlerGetters
  inline int  GetNumParams() const {return _fNumPar;}
  /// @brief Get the pre-fit values of the parameters.
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetPreFitValues() const {return _fPreFitValue;}
  /// @brief Get vector of all proposed parameter values
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetProposed() const;
  /// @brief Get proposed parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParProp(const int i) const { return _fPropVal[i]; }
  /// @brief Get current parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParCurr(const int i) const { return _fCurrVal[i]; }
  /// @brief Get vector of current parameter values
  /// @ingroup ParameterHandlerGetters
  inline const std::vector<double> &GetParCurrVec() const { return _fCurrVal; }

  /// @brief Get prior parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParInit(const int i) const { return _fPreFitValue[i]; }
  /// @brief Get upper parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetUpperBound(const int i) const { return _fUpBound[i];}
  /// @brief Get lower parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetLowerBound(const int i) const { return _fLowBound[i]; }
  /// @brief Get individual step scale for selected parameter
  /// @param ParameterIndex Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetIndivStepScale(const int ParameterIndex) const {return _fIndivStepScale.at(ParameterIndex); }
  /// @brief Get global step scale for covariance object
  /// @ingroup ParameterHandlerGetters
  inline double GetGlobalStepScale() const {return _fGlobalStepScale; }

  /// @brief Get number of params which will be different depending if using Eigen decomposition or not
  /// @ingroup ParameterHandlerGetters
  inline int GetNParameters() const {
    if (pca) return PCAObj->GetNumberPCAedParameters();
    else return _fNumPar;
  }

  /// @brief Print prior value for every parameter
  void PrintNominal() const;
  /// @brief Print prior, current and proposed value for each parameter
  void PrintNominalCurrProp() const;
  /// @warning only for backward compatibility
  /// @todo remove it
  void PrintParameters() const {PrintNominalCurrProp();};
  /// @brief Print step scale for each parameter
  void PrintIndivStepScale() const;

  /// @brief Generate a new proposed state
  virtual void ProposeStep();
  /// @brief "Randomize" the parameters in the covariance class for the proposed step. Used the proposal kernel and the current parameter value to set proposed step
  void Randomize() _noexcept_;
  /// @brief Use Cholesky throw matrix for better step proposal
  void CorrelateSteps() _noexcept_;
  /// @brief Method to update adaptive MCMC
  /// @cite haario2001adaptive
  void UpdateAdaptiveCovariance();

  /// @brief Accepted this step
  void AcceptStep() _noexcept_;

  /// @brief Set all parameters to be fixed at prior values
  void SetFixAllParameters();
  /// @brief Set parameter to be fixed at prior value
  /// @param i Parameter index
  void SetFixParameter(const int i);
  /// @brief Set parameter to be fixed at prior value
  /// @param name Name of the parameter to be fixed
  void SetFixParameter(const std::string& name);

  /// @brief Set all parameters to be treated as free
  void SetFreeAllParameters();
  /// @brief Set parameter to be treated as free
  /// @param i Parameter index
  void SetFreeParameter(const int i);
  /// @brief Set parameter to be treated as free
  /// @param name Name of the parameter to be treated as free
  void SetFreeParameter(const std::string& name);

  /// @brief Is parameter fixed or not
  /// @param i Parameter index
  bool IsParameterFixed(const int i) const {
    if (_fError[i] < 0) { return true; }
    else                { return false; }
  }
  /// @brief Is parameter fixed or not
  /// @param name Name of parameter you want to check if is fixed
  bool IsParameterFixed(const std::string& name) const;

  /// @brief CW: Calculate eigen values, prepare transition matrices and remove param based on defined threshold
  /// @param eigen_threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  /// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/03.-Eigen-Decomposition-%E2%80%90-PCA).
  void ConstructPCA(const double eigen_threshold, int FirstPCAdpar, int LastPCAdpar);

  /// @brief is PCA, can use to query e.g. LLH scans
  inline bool IsPCA() const { return pca; }

  /// @brief Getter to return a copy of the YAML node
  /// @ingroup ParameterHandlerGetters
  YAML::Node GetConfig() const { return _fYAMLDoc; }

  /// @brief Get pointer for AdaptiveHandler
  /// @ingroup ParameterHandlerGetters
  inline adaptive_mcmc::AdaptiveMCMCHandler* GetAdaptiveHandler() const  {
    if (!use_adaptive) {
      MACH3LOG_ERROR("Am not running in Adaptive mode");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    return AdaptiveHandler.get();
  }

  /// @brief KS: Set proposed parameter values vector to be base on tune values, for example set proposed values to be of generated or maybe PostND
  /// @ingroup ParameterHandlerSetters
  void SetTune(const std::string& TuneName);
  
  /// @brief Get pointer for PCAHandler
  inline PCAHandler* GetPCAHandler() const {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    return PCAObj.get();
  }

  /// @brief Matches branches in a TTree to parameters in a systematic handler.
  ///
  /// @param PosteriorFile Pointer to the ROOT TTree from MaCh3 fit.
  /// @param[out] BranchValues Vector to store the values of the branches (resized inside).
  /// @param[out] BranchNames Vector to store the names of the branches (resized inside).
  ///
  /// @throws MaCh3Exception if any parameter branch is uninitialized.
  void MatchMaCh3OutputBranches(TTree *PosteriorFile,
                                std::vector<double>& BranchValues,
                                std::vector<std::string>& BranchNames);
 protected:
  /// @brief Toggle fixing parameters at prior values
  void ToggleFixAllParameters();
  /// @brief Toggle fixing parameter at prior values
  /// @param i Parameter index
  void ToggleFixParameter(const int i);
  /// @brief Toggle fixing parameter at prior values
  /// @param name Name of parameter you want to fix
  void ToggleFixParameter(const std::string& name);

  /// @brief Initialisation of the class using matrix from root file
  void Init(const std::string& name, const std::string& file);
  /// @brief Initialisation of the class using config
  /// @param YAMLFile A vector of strings representing the YAML files used for initialisation of matrix
  void Init(const std::vector<std::string>& YAMLFile);
  /// @brief Initialise vectors with parameters information
  /// @param size integer telling size to which we will resize all vectors/allocate memory
  void ReserveMemory(const int size);

  /// @brief Make matrix positive definite by adding small values to diagonal, necessary for inverting matrix
  /// @param cov Matrix which we evaluate Positive Definitiveness
  void MakePosDef(TMatrixDSym *cov = nullptr);

  /// @brief HW: Finds closest possible positive definite matrix in Frobenius Norm ||.||_frob Where ||X||_frob=sqrt[sum_ij(x_ij^2)] (basically just turns an n,n matrix into vector in n^2 space then does Euclidean norm)
  void MakeClosestPosDef(TMatrixDSym *cov);

  /// @brief sets throw matrix from a file
  /// @param matrix_file_name name of file matrix lives in
  /// @param matrix_name name of matrix in file
  /// @param means_name name of means vec in file
  void SetThrowMatrixFromFile(const std::string& matrix_file_name, const std::string& matrix_name, const std::string& means_name);

  /// @brief Check if parameter is affecting given sample name
  /// @param SystIndex number of parameter
  /// @param SampleName The Sample name used to filter parameters.
  bool AppliesToSample(const int SystIndex, const std::string& SampleName) const;

  /// @brief KS: Flip parameter around given value, for example mass ordering around 0
  /// @param index parameter index you want to flip
  /// @param FlipPoint Value around which flipping is done
  void FlipParameterValue(const int index, const double FlipPoint);

  /// @brief HW :: This method is a tad hacky but modular arithmetic gives me a headache.
  /// @author Henry Wallace
  void CircularParBounds(const int i, const double LowBound, const double UpBound);

  /// @brief Enable special proposal
  void EnableSpecialProposal(const YAML::Node& param, const int Index);

  /// @brief Perform Special Step Proposal
  /// @warning KS: Following Asher comment we do "Step->Circular Bounds->Flip"
  void SpecialStepProposal();

  /// Check if any of special step proposal were enabled
  bool doSpecialStepProposal;

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
  std::vector<std::vector<std::string>> _fSampleNames;
  
  /// Backup of _fIndivStepScale for parameters which are skipped during adaption
  std::vector<double> _fIndivStepScaleInitial;

  /// Backup of _fGlobalStepScale for parameters which are skipped during adaption
  double _fGlobalStepScaleInitial;

  /// Flags telling if parameter should be skipped during adaption
  std::vector<bool> param_skip_adapt_flags;

  /// Matrix which we use for step proposal before Cholesky decomposition (not actually used for step proposal)
  TMatrixDSym* throwMatrix;
  /// Throw matrix that is being used in the fit, much faster as TMatrixDSym cache miss
  double** throwMatrixCholDecomp;

  /// perform PCA or not
  bool pca;
  /// Are we using AMCMC?
  bool use_adaptive;

  /// Struct containing information about PCA
  std::unique_ptr<PCAHandler> PCAObj;
  /// Struct containing information about adaption
  std::unique_ptr<adaptive_mcmc::AdaptiveMCMCHandler> AdaptiveHandler;
  /// Struct containing information about adaption
  std::unique_ptr<ParameterTunes> Tunes;

  /// Indices of parameters with flip symmetry
  std::vector<int>    FlipParameterIndex;
  /// Central points around which parameters are flipped
  std::vector<double> FlipParameterPoint;
  /// Indices of parameters with circular bounds
  std::vector<int>    CircularBoundsIndex;
  /// Circular bounds for each parameter (lower, upper)
  std::vector<std::pair<double,double>> CircularBoundsValues;
};
