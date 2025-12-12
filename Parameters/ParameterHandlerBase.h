#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Parameters/AdaptiveMCMCHandler.h"
#include "Parameters/PCAHandler.h"
#include "Parameters/ParameterHandlerUtils.h"
#include "Parameters/ParameterTunes.h"

#include "Eigen/Dense"

/// @brief Base class responsible for handling of systematic error parameters.
/// Capable of using PCA or using adaptive throw matrix
/// @see For more details, visit the
/// [Wiki](https://github.com/mach3-software/MaCh3/wiki/02.-Implementation-of-Systematic).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
class ParameterHandlerBase {
public:
  struct ParamInfo {
    std::string name, fancy_name;
    double prefit, error, stepscale;
    std::array<double, 2> bounds;
    bool flatprior, fixed;
    std::vector<std::string> affected_samples;
  };

  void AddParameters(std::vector<ParamInfo> const &params);
  void AddParameter(ParamInfo const &param) {
    AddParameters({
        params,
    });
  }

  /// @brief ETA - constructor for a YAML file
  /// @param YAMLFile A vector of strings representing the YAML files used for
  /// initialisation of matrix
  /// @param name Matrix name
  /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  static ParameterHandlerBase
  MakeFromYAML(const std::vector<std::string> &YAMLFiles, std::string name,
               double threshold = -1, int FirstPCAdpar = -999,
               int LastPCAdpar = -999) {
    return ParameterHandlerBase(YAMLFiles, name, threshold, FirstPCAdpar,
                                LastPCAdpar);
  }
  /// @brief "Usual" constructors from root file
  /// @param name Matrix name
  /// @param file Path to matrix root file
  static ParameterHandlerBase
  MakeFromTMatrix(std::string name, std::string file, double threshold = -1,
                  int FirstPCAdpar = -999, int LastPCAdpar = -999) {
    return ParameterHandlerBase(name, file, threshold, FirstPCAdpar,
                                LastPCAdpar);
  }

  /// @brief Destructor
  virtual ~ParameterHandlerBase() {};

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
  void SetName(const std::string &name) { settings.name = name; }
  /// @brief change parameter name
  /// @param i Parameter index
  /// @param name new name which will be set
  /// @ingroup ParameterHandlerSetters
  void SetParName(const int i, const std::string &name) {
    params.name.at(i) = name;
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
  void SetParProp(const int i, const double val) { params.proposed[i] = val; }
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
  void SetRandomThrow(const int i, const double rand) {
    throws.random_vector[i] = rand;
  }
  /// @brief Get random value useful for debugging/CI
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetRandomThrow(const int i) const { return throws.random_vector[i]; }

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
  void SetIndivStepScale(const int ParameterIndex, const double StepScale) {
    steps.scale[ParameterIndex] = StepScale;
  }
  /// @brief DB Function to set fIndivStepScale from a vector (Can be used from
  /// execs and inside covariance constructors)
  /// @param stepscale Vector of individual step scale, should have same
  /// @ingroup ParameterHandlerSetters
  void SetIndivStepScale(const std::vector<double> &stepscale);
  /// @brief KS: In case someone really want to change this
  /// @ingroup ParameterHandlerSetters
  inline void SetPrintLength(const unsigned int PriLen) {
    settings.PrintLength = PriLen;
  }

  /// @brief KS: After step scale, prefit etc. value were modified save this
  /// modified config.
  void
  SaveUpdatedMatrixConfig(std::string const &filename = "Modified_Matrix.yaml");

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
  double CalcLikelihood() const _noexcept_;
  /// @brief Return CalcLikelihood if some params were thrown out of boundary
  /// return _LARGE_LOGL_
  /// @ingroup ParameterHandlerGetters
  virtual double GetLikelihood();

  /// @brief Return covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *GetCovMatrix() const { return &root_copies.covariance; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  TMatrixDSym *GetInvCovMatrix() const { return &root_copies.inv_covariance; }
  /// @brief Return inverted covariance matrix
  /// @ingroup ParameterHandlerGetters
  double GetInvCovMatrix(const int i, const int j) const {
    return params.inv_covariance(i, j);
  }

  /// @brief Return correlated throws
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  double GetCorrThrows(const int i) const { return corr_throw[i]; }

  /// @brief Get if param has flat prior or not
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline bool GetFlatPrior(const int i) const { return params.flatprior[i]; }

  /// @brief Get name of covariance
  /// @ingroup ParameterHandlerGetters
  std::string GetName() const { return params.matrix_name; }
  /// @brief Get name of parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParName(const int i) const { return params.name[i]; }

  /// @brief Get index based on name
  /// @ingroup ParameterHandlerGetters
  int GetParIndex(const std::string &name) const;

  /// @brief Get fancy name of the Parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  std::string GetParFancyName(const int i) const { return params.fancyname[i]; }
  /// @brief Get name of input file
  /// @ingroup ParameterHandlerGetters
  std::string GetInputFile() const { return config.inputFile; }

  /// @brief Get diagonal error for ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetDiagonalError(const int i) const {
    return std::sqrt(params.covariance(i, i));
  }
  /// @brief Get the error for the ith parameter
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetError(const int i) const { return params.error[i]; }

  /// @brief Adaptive Step Tuning Stuff
  void ResetIndivStepScale();

  /// @brief Initialise adaptive MCMC
  /// @param adapt_manager Node having from which we load all adaptation options
  void InitialiseAdaption(const YAML::Node &adapt_manager);
  /// @brief Save adaptive throw matrix to file
  void SaveAdaptiveToFile(const std::string &outFileName,
                          const std::string &systematicName) {
    AdaptiveHandler.SaveAdaptiveToFile(outFileName, systematicName);
  }

  /// @brief Do we adapt or not
  /// @ingroup ParameterHandlerGetters
  bool GetDoAdaption() const { return settings.use_adaptive; }
  /// @brief Use new throw matrix, used in adaptive MCMC
  /// @ingroup ParameterHandlerSetters
  void SetThrowMatrix(TMatrixDSym *cov);
  /// @brief Replaces old throw matrix with new one
  void UpdateThrowMatrix(TMatrixDSym *cov);
  /// @brief Set number of MCMC step, when running adaptive MCMC it is updated
  /// with given frequency. We need number of steps to determine frequency.
  /// @ingroup ParameterHandlerSetters
  inline void SetNumberOfSteps(const int nsteps) {
    AdaptiveHandler.SetTotalSteps(nsteps);
    if (AdaptiveHandler.AdaptionUpdate()) {
      ResetIndivStepScale();
    }
  }

  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  inline TMatrixDSym *GetThrowMatrix() const { return &root_copies.proposal; }
  /// @brief Get matrix used for step proposal
  /// @ingroup ParameterHandlerGetters
  double GetThrowMatrix(const int i, const int j) const {
    return params.l_proposal(i, j);
  }

  /// @brief KS: Convert covariance matrix to correlation matrix and return TH2D
  /// which can be used for fancy plotting
  /// @details This function converts the covariance matrix to a correlation
  /// matrix and
  ///          returns a TH2D object, which can be used for advanced plotting
  ///          purposes.
  /// @return A pointer to a TH2D object representing the correlation matrix
  /// @ingroup ParameterHandlerGetters
  TH2D *GetCorrelationMatrix();

  /// @brief Get total number of parameters
  /// @ingroup ParameterHandlerGetters
  inline int GetNumParams() const { return params.prefit.size(); }
  /// @brief Get the pre-fit values of the parameters.
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetPreFitValues() const {
    return M3::EigenToStdVector(params.prefit);
  }
  /// @brief Get vector of all proposed parameter values
  /// @ingroup ParameterHandlerGetters
  std::vector<double> GetProposed() const {
    return M3::EigenToStdVector(params.proposed);
  }
  /// @brief Get proposed parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParProp(const int i) const { return params.proposed[i]; }
  /// @brief Get current parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParCurr(const int i) const { return params.current[i]; }

  /// @brief Get prior parameter value
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetParInit(const int i) const { return params.prefit[i]; }
  /// @brief Get upper parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetUpperBound(const int i) const { return params.upbound[i]; }
  /// @brief Get lower parameter bound in which it is physically valid
  /// @param i Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetLowerBound(const int i) const { return params.lowbound[i]; }
  /// @brief Get individual step scale for selected parameter
  /// @param ParameterIndex Parameter index
  /// @ingroup ParameterHandlerGetters
  inline double GetIndivStepScale(const int i) const {
    return params.stepscale[i];
  }
  /// @brief Get global step scale for covariance object
  /// @ingroup ParameterHandlerGetters
  inline double GetGlobalStepScale() const { return settings.GlobalStepScale; }

  /// @brief Get number of params which will be different depending if using
  /// Eigen decomposition or not
  /// @ingroup ParameterHandlerGetters
  inline int GetNParameters() const {
    if (pca)
      return PCAObj->GetNumberPCAedParameters();
    else
      return GetNumParams();
  }

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
  bool IsParameterFixed(const int i) const { return params.fixed[i]; }

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
  inline bool IsPCA() const { return pca; }

  /// @brief Getter to return a copy of the YAML node
  /// @ingroup ParameterHandlerGetters
  YAML::Node GetConfig() const { return config.YAMLDoc; }

  /// @brief Get pointer for AdaptiveHandler
  /// @ingroup ParameterHandlerGetters
  inline adaptive_mcmc::AdaptiveMCMCHandler &GetAdaptiveHandler() const {
    if (!use_adaptive) {
      MACH3LOG_ERROR("Am not running in Adaptive mode");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    return AdaptiveHandler;
  }

  /// @brief KS: Set proposed parameter values vector to be base on tune values,
  /// for example set proposed values to be of generated or maybe PostND
  /// @ingroup ParameterHandlerSetters
  void SetTune(const std::string &TuneName);

  /// @brief Get pointer for PCAHandler
  inline PCAHandler &GetPCAHandler() const {
    if (!pca) {
      MACH3LOG_ERROR("Am not running in PCA mode");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    return PCAObj;
  }

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
  ParameterHandlerBase();

  /// @brief ETA - constructor for a YAML file
  /// @param YAMLFile A vector of strings representing the YAML files used for
  /// initialisation of matrix
  /// @param name Matrix name
  /// @param threshold PCA threshold from 0 to 1. Default is -1 and means no PCA
  /// @param FirstPCAdpar First PCA parameter that will be decomposed.
  /// @param LastPCAdpar First PCA parameter that will be decomposed.
  ParameterHandlerBase(const std::vector<std::string> &YAMLFiles,
                       std::string name, double threshold = -1,
                       int FirstPCAdpar = -999, int LastPCAdpar = -999);
  /// @brief "Usual" constructors from root file
  /// @param name Matrix name
  /// @param file Path to matrix root file
  ParameterHandlerBase(std::string name, std::string file,
                       double threshold = -1, int FirstPCAdpar = -999,
                       int LastPCAdpar = -999);

  /// @brief sets throw matrix from a file
  /// @param matrix_file_name name of file matrix lives in
  /// @param matrix_name name of matrix in file
  /// @param means_name name of means vec in file
  void SetThrowMatrixFromFile(const std::string &matrix_file_name,
                              const std::string &matrix_name,
                              const std::string &means_name);

  /// @brief Check if parameter is affecting given sample name
  /// @param SystIndex number of parameter
  /// @param SampleName The Sample name used to filter parameters.
  bool AppliesToSample(const int SystIndex,
                       const std::string &SampleName) const;

  /// @brief KS: Flip parameter around given value, for example mass ordering
  /// around 0
  /// @param index parameter index you want to flip
  /// @param FlipPoint Value around which flipping is done
  void FlipParameterValue(const int index, const double FlipPoint);

  /// @brief HW :: This method is a tad hacky but modular arithmetic gives me a
  /// headache.
  /// @author Henry Wallace
  void CircularParBounds(const int i, const double LowBound,
                         const double UpBound);

  /// @brief Enable special proposal
  void EnableSpecialProposal(const YAML::Node &param, const int Index);

  /// @brief Perform Special Step Proposal
  /// @warning KS: Following Asher comment we do "Step->Circular Bounds->Flip"
  void SpecialStepProposal();

  struct {
    bool enabled;
    std::array<int, 2> block_indices;
    Eigen::MatrixXd rotation;
  } pca;

  struct {
    Eigen::VectorXd random_vector, values;
  } throws;

  struct {
    Eigen::VectorXd current, proposed, scale;
    double global_scale;
    Eigen::MatrixXd l_proposal;
  } steps;

  struct {
    bool enabled;
    /// Indices of parameters with flip symmetry
    std::vector<int> FlipParameterIndex;
    /// Central points around which parameters are flipped
    std::vector<double> FlipParameterPoint;
    /// Indices of parameters with circular bounds
    std::vector<int> CircularBoundsIndex;
    /// Circular bounds for each parameter (lower, upper)
    std::vector<std::pair<double, double>> CircularBoundsValues;
  } special_proposal;

  struct {
    std::vector<std::string> name, fancy_name;
    Eigen::VectorXd prefit, error, lowbound, upbound;
    Eigen::VectorXi flatprior, fixed;
    std::vector<std::vector<std::string>> samples;

    Eigen::MatrixXd covariance;
    Eigen::MatrixXd inv_covariance;
  } params;

  struct {
    std::ranlux48 e1;
    std::normal_distribution<double> gaus;
  } rng;

  // to retain external interface where possible
  struct {
    TMatrixDSym covariance;
    TMatrixDSym inv_covariance;
    TMatrixDSym proposal;
  } root_copies;

  struct {
    std::string name;

    /// Are we using AMCMC?
    bool use_adaptive;

    /// KS: This is used when printing parameters, sometimes we have super long
    /// parameters name, we want to flexibly adjust couts
    int PrintLength;
  } settings;

  struct {
    /// The input root file we read in
    std::string inputFile;
    /// Stores config describing systematics
    YAML::Node YAMLDoc;
  } config;

  /// Struct containing information about adaption
  adaptive_mcmc::AdaptiveMCMCHandler AdaptiveHandler;
  /// Struct containing information about adaption
  ParameterTunes Tunes;
};
