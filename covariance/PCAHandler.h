#pragma once

// MaCh3 Includes
#include "manager/manager.h"
#include "covariance/CovarianceUtils.h"

//#define DEBUG_PCA 1
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

class PCAHandler{
 public:

  /// @brief Constructor
  PCAHandler();

  /// @brief Destructor
  ~PCAHandler();

  /// @brief CW: Calculate eigen values, prepare transition matrices and remove param based on defined threshold
  /// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/03.-Eigen-Decomposition-%E2%80%90-PCA).
  void ConstructPCA(TMatrixDSym * covMatrix, const int firstPCAd, const int lastPCAd, const double eigen_thresh, int& _fNumParPCA);

  #ifdef DEBUG_PCA
  /// @brief KS: Let's dump all useful matrices to properly validate PCA
  void DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat, int NumPar);
  #endif

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
};

