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

namespace covutils
{
  TH2D* TMatrixIntoTH2D(const TMatrix &Matrix, std::string title);

  int GetNThreads();
  double* MatrixMult(double*, double*, int);
  double** MatrixMult(double**, double**, int);
  TMatrixD MatrixMult(TMatrixD, TMatrixD);

}
