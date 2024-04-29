#pragma once

// ROOT includes
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TMath.h"
#include "math.h"
#include "TDecompChol.h"
#include "TStopwatch.h"
#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDEigen.h"
#include "TDecompSVD.h"


#ifdef MULTITHREAD
#include "omp.h"
#endif

namespace MaCh3Utils
{
  /// @brief Handy conversion for matrix to TH2D for plotting
  TH2D* TMatrixIntoTH2D(const TMatrix &Matrix, std::string title);

  /// @brief number of threads which we need for example for TRandom3
  int GetNThreads();
  /// @brief CW: Multi-threaded matrix multiplication
  double* MatrixMult(double*, double*, int);
  /// @brief CW: Multi-threaded matrix multiplication
  double** MatrixMult(double**, double**, int);
  /// @brief CW: Multi-threaded matrix multiplication
  TMatrixD MatrixMult(TMatrixD, TMatrixD);

}
