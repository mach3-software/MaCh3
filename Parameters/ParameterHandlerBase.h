#pragma once

// MaCh3 includes
#include "Manager/Manager.h"

#include "Parameters/ParameterHandlerUtils.h"
#include "Parameters/ParameterList.h"
#include "Parameters/StepProposer.h"

/// @brief Base class responsible for handling of systematic error parameters.
/// Capable of using PCA or using adaptive throw matrix
/// @see For more details, visit the
/// [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class ParameterHandlerBase {
public:
  /// @brief ETA - constructor for a YAML file
  /// @param YAMLFile A vector of strings representing the YAML files used for
  /// initialisation of matrix
  /// @param name Matrix name
  /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  ParameterHandlerBase(const std::vector<std::string> &YAMLFile,
                       std::string name, double threshold = -1,
                       int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief "Usual" constructors from root file
  /// @param name Matrix name
  /// @param file Path to matrix root file
  ParameterHandlerBase(std::string name, std::string file,
                       double threshold = -1, int FirstPCAdpar = -999,
                       int LastPCAdpar = -999);

  /// @brief Destructor
  virtual ~ParameterHandlerBase() {}

  /// @defgroup ParameterHandlerSetters Parameter Handler Setters
  /// Group of functions to set various parameters, names, and values.

  /// @defgroup ParameterHandlerGetters Parameter Handler Getters
  /// Group of functions to get various parameters, names, and values.

  // ETA - maybe need to add checks to index on the setters? i.e. if( i >
  // _fPropVal.size()){throw;}
  /// @brief Set covariance matrix
  /// @param cov Covariance matrix which we set and will be used later for
  /// evaluation of penalty term
  /// @ingroup ParameterHandlerSetters
  void SetCovMatrix(TMatrixDSym *cov);
  /// @brief Set matrix name
  /// @ingroup ParameterHandlerSetters
  std::string matrixName;
  void SetName(const std::string &name) { matrixName = name; }
  /// @brief change parameter name
  /// @param i Parameter index
  /// @param name new name which will be set
  /// @ingroup ParameterHandlerSetters
  void SetParName(const int i, const std::string &name) {
    parlist.params.name.at(i) = name;
  }
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
  void SetParProp(const int i, const double val);
  /// @brief Set parameter values using vector, it has to have same size as
  /// covariance class
  /// @param pars Vector holding new values for every parameter
  /// @ingroup ParameterHandlerSetters
  void SetParameters(const std::vector<double> &pars = {});
  /// @brief Set if parameter should have flat prior or not
  /// @param i Parameter index
  /// @param eL bool telling if it will be flat or not
  /// @ingroup ParameterHandlerSetters
  void SetFlatPrior(const int i, const bool eL);

  /// @brief Set random value useful for debugging/CI
  /// @param i Parameter index
  /// @param rand New value for random number
  /// @ingroup ParameterHandlerSetters
  void SetRandomThrow(const int i, const double rand);
  /// @brief Get random value useful for debugging/CI
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetRandomThrow(const int i) const;

  /// @brief set branches for output file
  /// @param tree Tree to which we will save branches
  /// @param SaveProposal Normally we only save parameter after is accepted, for
  /// debugging purpose it is helpful to see also proposed values. That's what
  /// this variable controls
  /// @ingroup ParameterHandlerSetters
  void SetBranches(TTree &tree, const bool SaveProposal = false);
  /// @brief Set global step scale for covariance object
  /// @param scale Value of global step scale
  /// @param verbose Print that we've changed scale + use warnings [default:
  /// true]
  /// @cite luengo2020survey
  /// @ingroup ParameterHandlerSetters
  void SetStepScale(const double scale, const bool verbose = true);
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from
  /// execs and inside covariance constructors)
  /// @param ParameterIndex Parameter Index
  /// @param StepScale Value of individual step scale
  /// @ingroup ParameterHandlerSetters
  void SetIndivStepScale(const int ParameterIndex, const double StepScale);
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from
  /// execs and inside covariance constructors)
  /// @param stepscale Vector of individual step scale, should have same
  /// @ingroup ParameterHandlerSetters
  void SetIndivStepScale(const std::vector<double> &stepscale);
  /// @brief KS: In case someone really want to change this
  /// @ingroup ParameterHandlerSetters
  int PrintLength = 35;
  void SetPrintLength(const unsigned int PriLen) { PrintLength = PriLen; }

  /// @brief KS: After step scale, prefit etc. value were modified save this
  /// modified config.
  void SaveUpdatedMatrixConfig();

  /// @brief Throw the proposed parameter by mag sigma. Should really just have
  /// the user specify this throw by having argument double
  void ThrowParProp(const double mag = 1.);

  /// @brief Helper function to throw the current parameter by mag sigma. Can
  /// study bias in MCMC with this; put different starting parameters
  void ThrowParCurr(const double mag = 1.);
  /// @brief Throw the parameters according to the covariance matrix. This
  /// shouldn't be used in MCMC code ase it can break Detailed Balance;
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
  ///   \log \mathcal{L} = \frac{1}{2} \sum_{i}^{\textrm{pars}}
  ///   \sum_{j}^{\textrm{pars}} \Delta \vec{p}_i \left( V^{-1} \right)_{i,j}
  ///   \Delta \vec{p}_j
  /// \f]
  /// where:
  /// - \f$\Delta \vec{p}_i = \theta_i - \theta_{i,0}\f$ is the difference
  /// between the current and pre-fit parameter values,
  /// - \f$V^{-1}\f$ is the inverted covariance matrix.
  ///
  /// @note
  /// - If `_fFlatPrior[i]` is `true`, the parameter is excluded from the
  /// calculation.
  double CalcLikelihood() _noexcept_;
  /// @brief Return CalcLikelihood if some params were thrown out of boundary
  /// return _LARGE_LOGL_
  /// @ingroup ParameterHandlerGetters
  virtual double GetLikelihood();

  /// @brief Return covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *covMatrix;
  TMatrixDSym *GetCovMatrix() const { return covMatrix; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *invCovMatrix;
  TMatrixDSym *GetInvCovMatrix() const { return invCovMatrix; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  double GetInvCovMatrix(const int i, const int j) const;

  /// @brief Return correlated throws
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetCorrThrows(const int i) const;

  /// @brief Get if param has flat prior or not
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  bool GetFlatPrior(const int i) const;

  /// @brief Get name of covariance
  /// @ingroup ParameterHandlerGetters
  std::string GetName() const { return matrixName; }
  /// @brief Get name of parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParName(const int i) const {
    return parlist.params.name.at(i);
  }

  /// @brief Get index based on name
  /// @ingroup ParameterHandlerGetters
  int GetParIndex(const std::string &name) const;

  /// @brief Get fancy name of the Parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParFancyName(const int i) const {
    return parlist.params.fancy_name.at(i);
    ;
  }
  /// @brief Get name of input file
  /// @ingroup ParameterHandlerGetters
  std::string inputFile;
  std::string GetInputFile() const { return inputFile; }

  /// @brief Get diagonal error for ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetDiagonalError(const int i) const {
    return std::sqrt(parlist.params.covariance(i, i));
  }
  /// @brief Get the error for the ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetError(const int i) const { return parlist.params.error[i]; }

  /// @brief Adaptive Step Tuning Stuff
  void ResetIndivStepScale();

  /// @brief Initialise adaptive MCMC
  /// @param adapt_manager Node having from which we load all adaptation options
  void InitialiseAdaption(const YAML::Node &adapt_manager);
  /// @brief Save adaptive throw matrix to file
  void SaveAdaptiveToFile(const std::string &outFileName,
                          const std::string &systematicName);

  /// @brief Do we adapt or not
  /// @ingroup ParameterHandlerGetters
  bool GetDoAdaption() const { return false; }
  /// @brief Use new throw matrix, used in adaptive MCMC
  /// @ingroup ParameterHandlerSetters
  void SetThrowMatrix(TMatrixDSym *cov);
  /// @brief Replaces old throw matrix with new one
  void UpdateThrowMatrix(TMatrixDSym *cov);
  /// @brief Set number of MCMC step, when running adaptive MCMC it is updated
  /// with given frequency. We need number of steps to determine frequency.
  /// @ingroup ParameterHandlerSetters
  void SetNumberOfSteps(const int nsteps);

  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *throwMatrix;
  TMatrixDSym *GetThrowMatrix() const { return throwMatrix; }
  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  double GetThrowMatrix(const int i, const int j) const;

  /// @brief KS: Convert covariance matrix to correlation matrix and return TH2D
  /// which can be used for fancy plotting
  /// @details This function converts the covariance matrix to a correlation
  /// matrix and
  ///          returns a TH2D object, which can be used for advanced plotting
  ///          purposes.
  /// @return A pointer to a TH2D object representing the correlation matrix
  /// @ingroup ParameterHandlerGetters
  TH2D *GetCorrelationMatrix();

  /// @brief DB Pointer return to param position
  ///
  /// @param iParam The index of the parameter in the vector.
  /// @return A pointer to the parameter value at the specified index.
  ///
  /// @warning ETA - This might be a bit squiffy? If the vector gots moved from
  /// say a push_back then the pointer is no longer valid... maybe need a better
  /// way to deal with this? It was fine before when the return was to an
  /// element of a new array. There must be a clever C++ way to be careful
  const double *RetPointer(const int iParam);

  /// @brief Get a reference to the proposed parameter values
  /// Can be useful if you want to track these without having to copy values
  /// using getProposed()
  mutable std::vector<double> _fPropVal;
  const std::vector<double> &GetParPropVec() const {
    _fPropVal = M3::EigenToStdVector(proposed);
    return _fPropVal;
  }

  /// @brief Get total number of parameters
  /// @ingroup ParameterHandlerGetters
  int GetNumParams() const { return parlist.NumSystematicBasisParameters(); }
  /// @brief Get the pre-fit values of the parameters.
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetPreFitValues() const {
    return M3::EigenToStdVector(parlist.params.prefit);
  }
  /// @brief Get vector of all proposed parameter values
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetProposed() const;
  /// @brief Get proposed parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetParProp(const int i) const { return proposed[i]; }
  /// @brief Get current parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetParCurr(const int i) const { return current[i]; }
  /// @brief Get vector of current parameter values
  /// @ingroup ParameterHandlerGetters
  mutable std::vector<double> _fCurrVal;
  const std::vector<double> &GetParCurrVec() const {
    _fCurrVal = M3::EigenToStdVector(current);
    return _fCurrVal;
  }

  /// @brief Get prior parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetParInit(const int i) const { return parlist.params.prefit[i]; }
  /// @brief Get upper parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetUpperBound(const int i) const { return parlist.params.upbound[i]; }
  /// @brief Get lower parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetLowerBound(const int i) const { return parlist.params.lowbound[i]; }
  /// @brief Get individual step scale for selected parameter
  /// @param ParameterIndex Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetIndivStepScale(const int ParameterIndex) const;
  /// @brief Get global step scale for covariance object
  /// @ingroup ParameterHandlerGetters
  double GetGlobalStepScale() const { return proposer.params.global_scale; }

  /// @brief Get number of params which will be different depending if using
  /// Eigen decomposition or not
  /// @ingroup ParameterHandlerGetters
  int GetNParameters() const { return parlist.NumSystematicBasisParameters(); }

  /// @brief Print prior value for every parameter
  void PrintNominal() const;
  /// @brief Print prior, current and proposed value for each parameter
  void PrintNominalCurrProp() const;
  /// @warning only for backward compatibility
  /// @todo remove it
  void PrintParameters() const { PrintNominalCurrProp(); };
  /// @brief Print step scale for each parameter
  void PrintIndivStepScale() const;

  /// @brief Generate a new proposed state
  virtual void ProposeStep();
  /// @brief "Randomize" the parameters in the covariance class for the proposed
  /// step. Used the proposal kernel and the current parameter value to set
  /// proposed step
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
  void SetFixParameter(const std::string &name);

  /// @brief Set all parameters to be treated as free
  void SetFreeAllParameters();
  /// @brief Set parameter to be treated as free
  /// @param i Parameter index
  void SetFreeParameter(const int i);
  /// @brief Set parameter to be treated as free
  /// @param name Name of the parameter to be treated as free
  void SetFreeParameter(const std::string &name);

  /// @brief Toggle fixing parameters at prior values
  void ToggleFixAllParameters();
  /// @brief Toggle fixing parameter at prior values
  /// @param i Parameter index
  void ToggleFixParameter(const int i);
  /// @brief Toggle fixing parameter at prior values
  /// @param name Name of parameter you want to fix
  void ToggleFixParameter(const std::string &name);
  /// @brief Is parameter fixed or not
  /// @param i Parameter index
  bool IsParameterFixed(const int i) const;
  /// @brief Is parameter fixed or not
  /// @param name Name of parameter you want to check if is fixed
  bool IsParameterFixed(const std::string &name) const;

  /// @brief CW: Calculate eigen values, prepare transition matrices and remove
  /// param based on defined threshold
  /// @param eigen_threshold PCA threshold from 0 to 1. Default is -1 and means
  /// no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  /// @see For more details, visit the
  /// [Wiki](https://github.com/mach3-software/MaCh3/wiki/03.-Eigen-Decomposition-%E2%80%90-PCA).
  void ConstructPCA(const double eigen_threshold, int FirstPCAdpar,
                    int LastPCAdpar);

  /// @brief is PCA, can use to query e.g. LLH scans
  bool IsPCA() const { return parlist.pca.enabled; }

  /// @brief Getter to return a copy of the YAML node
  /// @ingroup ParameterHandlerGetters
  YAML::Node _fYAMLDoc;
  YAML::Node GetConfig() { return _fYAMLDoc; }

  /// @brief KS: Set proposed parameter values vector to be base on tune values,
  /// for example set proposed values to be of generated or maybe PostND
  /// @ingroup ParameterHandlerSetters
  void SetTune(const std::string &TuneName);

  /// @brief Matches branches in a TTree to parameters in a systematic handler.
  ///
  /// @param PosteriorFile Pointer to the ROOT TTree from MaCh3 fit.
  /// @param[out] BranchValues Vector to store the values of the branches
  /// (resized inside).
  /// @param[out] BranchNames Vector to store the names of the branches (resized
  /// inside).
  ///
  /// @throws MaCh3Exception if any parameter branch is uninitialized.
  void MatchMaCh3OutputBranches(TTree *PosteriorFile,
                                std::vector<double> &BranchValues,
                                std::vector<std::string> &BranchNames);

protected:
  /// @brief Check if parameter is affecting given sample name
  /// @param SystIndex number of parameter
  /// @param SampleName The Sample name used to filter parameters.
  bool AppliesToSample(const int SystIndex,
                       const std::string &SampleName) const;

  ParameterList parlist;
  StepProposer proposer;

  mutable Eigen::ArrayXd current;
  mutable Eigen::ArrayXd proposed;
};
