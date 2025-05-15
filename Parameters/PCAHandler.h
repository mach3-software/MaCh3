#pragma once

// MaCh3 Includes
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerUtils.h"

#ifdef DEBUG
  #define DEBUG_PCA 1
#endif

#ifdef DEBUG_PCA
//KS: When debugging we produce some fancy plots, but we don't need it during normal work flow
#include "TCanvas.h"
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

/// @brief Class responsible for handling Principal Component Analysis (PCA) of covariance matrix
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/03.-Eigen-Decomposition-%E2%80%90-PCA).
/// @author Clarence Wret
class PCAHandler{
 public:
  /// @brief Constructor
  PCAHandler();

  /// @brief Destructor
  virtual ~PCAHandler();

  /// @brief KS: Setup pointers to current and proposed parameter value which we need to convert them to PCA base each step
  /// @param fCurr_Val pointer to current position of parameter
  /// @param fProp_Val pointer to proposed position of parameter
  void SetupPointers(std::vector<double>* fCurr_Val,
                     std::vector<double>* fProp_Val);

  /// @brief CW: Calculate eigen values, prepare transition matrices and remove param based on defined threshold
  /// @param CovMatrix       Symmetric covariance matrix used for eigen decomposition.
  /// @param firstPCAd       Index of the first PCA component to include.
  /// @param lastPCAd        Index of the last PCA component to include.
  /// @param eigen_thresh    Threshold for eigenvalues below which parameters are discarded.
  /// @param _fNumParPCA     Output: number of parameters retained in PCA basis.
  /// @param _fNumPar        Total number of parameters in the original (non-PCA) basis.
  void ConstructPCA(TMatrixDSym * CovMatrix, const int firstPCAd, const int lastPCAd,
                    const double eigen_thresh, int& _fNumParPCA, const int _fNumPar);

  /// @brief Transfer param values from normal base to PCA base
  void TransferToPCA();
  /// @brief Transfer param values from PCA base to normal base
  void TransferToParam();

  /// @brief Accepted this step
  void AcceptStep() _noexcept_;
  /// @brief Use Cholesky throw matrix for better step proposal
  void CorrelateSteps(const std::vector<double>& IndivStepScale,
                      const double GlobalStepScale,
                      const double* _restrict_ randParams,
                      const double* _restrict_ corr_throw) _noexcept_;

  /// @brief KS: Transfer the starting parameters to the PCA basis, you don't want to start with zero..
  void SetInitialParameters();

  #ifdef DEBUG_PCA
  /// @brief KS: Let's dump all useful matrices to properly validate PCA
  void DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat, int NumPar);
  #endif

  /// Prefit value for PCA params
  std::vector<double> _fPreFitValuePCA;
  /// CW: Current parameter value in PCA base
  TVectorD _fParPropPCA;
  /// CW: Proposed parameter value in PCA base
  TVectorD _fParCurrPCA;
  /// Tells if parameter is fixed in PCA base or not
  std::vector<double> _fParSigmaPCA;
  /// If param is decomposed this will return -1, if not this will return enumerator to param in normal base. This way we can map stuff like step scale etc between normal base and undecomposed param in eigen base.
  std::vector<int> isDecomposedPCA;
  /// Number of parameters in PCA base
  int NumParPCA;

  /// Matrix used to converting from PCA base to normal base
  TMatrixD TransferMat;
  /// Matrix used to converting from normal base to PCA base
  TMatrixD TransferMatT;
  /// Eigen value only of particles which are being decomposed
  TVectorD eigen_values;
  /// Eigen vectors only of params which are being decomposed
  TMatrixD eigen_vectors;
  /// Eigen values which have dimension equal to _fNumParPCA, and can be used in CorrelateSteps
  std::vector<double> eigen_values_master;

  /// Total number that remained after applying PCA Threshold
  int nKeptPCApars;
  /// Index of the first param that is being decomposed
  int FirstPCAdpar;
  /// Index of the last param that is being decomposed
  int LastPCAdpar;
  /// CW: Threshold based on which we remove parameters in eigen base
  double eigen_threshold;

  /// Pointer to current value of the parameter
  std::vector<double>* _pCurrVal;
  /// Pointer to proposed value of the parameter
  std::vector<double>* _pPropVal;
};

