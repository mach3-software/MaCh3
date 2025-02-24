#pragma once

// MaCh3 includes
#include "SampleHandler/Structs.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TDecompChol.h"
#include "TStopwatch.h"
#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDEigen.h"
#include "TDecompSVD.h"
_MaCh3_Safe_Include_End_ //}

namespace M3
{
  /// @brief CW: Multi-threaded matrix multiplication
  inline double* MatrixMult(double *A, double *B, int n) {
    //CW: First transpose to increse cache hits
    double *BT = new double[n*n];
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        BT[j*n+i] = B[i*n+j];
      }
    }

    // Now multiply
    double *C = new double[n*n];
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double sum = 0;
        for (int k = 0; k < n; k++) {
          sum += A[i*n+k]*BT[j*n+k];
        }
        C[i*n+j] = sum;
      }
    }
    delete BT;

    return C;
  }

  /// @brief CW: Multi-threaded matrix multiplication
  inline double** MatrixMult(double **A, double **B, int n) {
    // First make into monolithic array
    double *A_mon = new double[n*n];
    double *B_mon = new double[n*n];

    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        A_mon[i*n+j] = A[i][j];
        B_mon[i*n+j] = B[i][j];
      }
    }
    //CW: Now call the monolithic calculator
    double *C_mon = MatrixMult(A_mon, B_mon, n);
    delete A_mon;
    delete B_mon;

    // Return the double pointer
    double **C = new double*[n];
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < n; ++i) {
      C[i] = new double[n];
      for (int j = 0; j < n; ++j) {
        C[i][j] = C_mon[i*n+j];
      }
    }
    delete C_mon;

    return C;
  }

  /// @brief CW: Multi-threaded matrix multiplication
  inline TMatrixD MatrixMult(TMatrixD A, TMatrixD B)
  {
    double *C_mon = MatrixMult(A.GetMatrixArray(), B.GetMatrixArray(), A.GetNcols());
    TMatrixD C;
    C.Use(A.GetNcols(), A.GetNrows(), C_mon);
    return C;
  }

}
