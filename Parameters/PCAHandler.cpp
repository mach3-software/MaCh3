#include "Parameters/PCAHandler.h"

// ********************************************
PCAHandler::PCAHandler() {
// ********************************************
  _pCurrVal = nullptr;
  _pPropVal = nullptr;
}

// ********************************************
PCAHandler::~PCAHandler() {
// ********************************************
}

// ********************************************
void PCAHandler::SetupPointers(std::vector<double>* fCurr_Val,
                               std::vector<M3::float_t>* fProp_Val) {
// ********************************************
  _pCurrVal = fCurr_Val;
  _pPropVal = fProp_Val;
}

// ********************************************
void PCAHandler::ConstructPCA(TMatrixDSym* CovMatrix, const int firstPCAd, const int lastPCAd,
                              const double eigen_thresh, const int _fNumPar) {
// ********************************************
  FirstPCAdpar = firstPCAd;
  LastPCAdpar = lastPCAd;
  eigen_threshold = eigen_thresh;

  // Check that covariance matrix exists
  if (CovMatrix == nullptr) {
    MACH3LOG_ERROR("Covariance matrix for has not yet been set");
    MACH3LOG_ERROR("Can not construct PCA until it is set");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if(FirstPCAdpar > CovMatrix->GetNrows()-1 || LastPCAdpar>CovMatrix->GetNrows()-1) {
    MACH3LOG_ERROR("FirstPCAdpar and LastPCAdpar are higher than the number of parameters");
    MACH3LOG_ERROR("first: {} last: {}, params: {}", FirstPCAdpar, LastPCAdpar, CovMatrix->GetNrows()-1);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  if(FirstPCAdpar < 0 || LastPCAdpar < 0){
    MACH3LOG_ERROR("FirstPCAdpar and LastPCAdpar are less than 0 but not default -999");
    MACH3LOG_ERROR("first: {} last: {}", FirstPCAdpar, LastPCAdpar);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("PCAing parameters {} through {} inclusive", FirstPCAdpar, LastPCAdpar);
  int NumUnPCAdPars = CovMatrix->GetNrows()-(LastPCAdpar-FirstPCAdpar+1);

  // KS: Make sure we are not doing anything silly with PCA
  SanitisePCA(CovMatrix);

  TMatrixDSym submat(CovMatrix->GetSub(FirstPCAdpar,LastPCAdpar,FirstPCAdpar,LastPCAdpar));

  //CW: Calculate how many eigen values this threshold corresponds to
  TMatrixDSymEigen eigen(submat);
  eigen_values.ResizeTo(eigen.GetEigenValues());
  eigen_vectors.ResizeTo(eigen.GetEigenVectors());
  eigen_values = eigen.GetEigenValues();
  eigen_vectors = eigen.GetEigenVectors();
  double sum = 0;
  // Loop over eigen values and sum them up
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    sum += eigen_values(i);
  }
  nKeptPCApars = eigen_values.GetNrows();
  //CW: Now go through again and see how many eigen values correspond to threshold
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    // Get the relative size of the eigen value
    double sig = eigen_values(i)/sum;
    // Check against the threshold
    if (sig < eigen_threshold) {
      nKeptPCApars = i;
      break;
    }
  }
  NumParPCA = NumUnPCAdPars+nKeptPCApars;
  MACH3LOG_INFO("Threshold of {} on eigen values relative sum of eigen value ({}) generates {} eigen vectors, plus we have {} unpcad pars, for a total of {}", eigen_threshold, sum, nKeptPCApars, NumUnPCAdPars, NumParPCA);
  //DB Create array of correct size so eigen_values can be used in CorrelateSteps
  eigen_values_master = std::vector<double>(NumParPCA, 1.0);
  for (int i = FirstPCAdpar; i < FirstPCAdpar+nKeptPCApars; ++i) {eigen_values_master[i] = eigen_values(i-FirstPCAdpar);}

  // Now construct the transfer matrices
  //These matrices will be as big as number of unPCAd pars plus number of eigenvalues kept
  TransferMat.ResizeTo(CovMatrix->GetNrows(), NumParPCA);
  TransferMatT.ResizeTo(CovMatrix->GetNrows(), NumParPCA);

  // Get a subset of the eigen vector matrix
  TMatrixD temp(eigen_vectors.GetSub(0, eigen_vectors.GetNrows()-1, 0, nKeptPCApars-1));

  //Make transfer matrix which is two blocks of identity with a block of the PCA transfer matrix in between
  TMatrixD temp2;
  temp2.ResizeTo(CovMatrix->GetNrows(), NumParPCA);

  //First set the whole thing to 0
  temp2.Zero();

  //Set the first identity block for non-PCAed params before PCA block, before PCA XRows == YRows
  for(int iRow = 0; iRow < FirstPCAdpar; iRow++) {
    temp2(iRow, iRow) = 1.0;
  }

  //Set the transfer matrix block for the PCAd pars
  temp2.SetSub(FirstPCAdpar,FirstPCAdpar,temp);

  //Set the second identity block starting after PCA block, remember XRows != YRows.
  // XRows -> normal base, YRows, PCA base
  for(int iRow = 0;iRow < (CovMatrix->GetNrows()-1)-LastPCAdpar; iRow++) {
    const int OrigRow = LastPCAdpar + 1 + iRow;
    const int PCARow = FirstPCAdpar + nKeptPCApars + iRow;
    temp2(OrigRow, PCARow) = 1.;
  }

  TransferMat = temp2;
  // Copy the contents
  TransferMatT = TransferMat;
  // And then transpose
  TransferMatT.T();

  // Make the PCA parameter arrays
  _fParCurrPCA.ResizeTo(NumParPCA);
  _fParPropPCA.ResizeTo(NumParPCA);
  _fPreFitValuePCA.resize(NumParPCA);

  //KS: make easy map so we could easily find un-decomposed parameters
  _fErrorPCA.assign(NumParPCA, 1);
  isDecomposedPCA.assign(NumParPCA, -1);
  // First non PCA-ed block, since this is before PCA-ed block we don't need any mapping
  for (int i = 0; i < FirstPCAdpar; ++i) {
    isDecomposedPCA[i] = i;
  }

  // Second non-PCA-ed block, keep in mind this is in PCA base so we cannot use LastPCAdpar
  // we must shift them back into the original parameter index range.
  const int shift = _fNumPar - NumParPCA;
  for (int i = FirstPCAdpar + nKeptPCApars; i < NumParPCA; ++i) {
    isDecomposedPCA[i] = i + shift;
  }

  #ifdef MACH3_DEBUG
  //KS: Let's dump all useful matrices to properly validate PCA
  M3::DebugPCA(sum, temp, eigen_vectors, eigen_values, TransferMat, TransferMatT,
               CovMatrix->GetNrows(), FirstPCAdpar, LastPCAdpar,
               nKeptPCApars, eigen_threshold);
  #endif
}

// ********************************************
// Make sure decomposed matrix isn't correlated with undecomposed
void PCAHandler::SanitisePCA(TMatrixDSym* CovMatrix) {
// ********************************************
  constexpr double correlation_threshold = 1e-6;

  bool found_significant_correlation = false;

  int N = CovMatrix->GetNrows();
  for (int i = FirstPCAdpar; i <= LastPCAdpar; ++i) {
    for (int j = 0; j < N; ++j) {
      // Skip if j is inside the decomposed range (we only want cross-correlations)
      if (j >= FirstPCAdpar && j <= LastPCAdpar) continue;

      double corr_val = (*CovMatrix)(i, j);
      if (std::fabs(corr_val) > correlation_threshold) {
        found_significant_correlation = true;
        MACH3LOG_ERROR("Significant correlation detected between decomposed parameter '{}' "
        "and undecomposed parameter '{}': {:.6e}", i, j, corr_val);
      }
    }
  }

  if (found_significant_correlation) {
    MACH3LOG_ERROR("There are correlations between undecomposed and decomposed part of matrices, this will not work");
    throw MaCh3Exception(__FILE__ , __LINE__);
  }
}

// ********************************************
// Update so that current step becomes the previously proposed step
void PCAHandler::AcceptStep() _noexcept_ {
// ********************************************
  // Update the book-keeping for the output
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < NumParPCA; ++i) {
    _fParCurrPCA(i) = _fParPropPCA(i);
  }
  // Then update the parameter basis
  TransferToParam();
}

// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void PCAHandler::CorrelateSteps(const std::vector<double>& IndivStepScale,
                                const double GlobalStepScale,
                                const double* _restrict_ randParams,
                                const double* _restrict_ corr_throw) _noexcept_ {
// ************************************************
  // Throw around the current step
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < NumParPCA; ++i)
  {
    if (_fErrorPCA[i] > 0.)
    {
      double IndStepScale = 1.;
      //KS: If undecomposed parameter apply individual step scale and Cholesky for better acceptance rate
      if(isDecomposedPCA[i] >= 0)
      {
        IndStepScale *= IndivStepScale[isDecomposedPCA[i]];
        IndStepScale *= corr_throw[isDecomposedPCA[i]];
      }
      //If decomposed apply only random number
      else
      {
        IndStepScale *= randParams[i];
        //KS: All PCA-ed parameters have the same step scale
        IndStepScale *= IndivStepScale[FirstPCAdpar];
      }
      _fParPropPCA(i) = _fParCurrPCA(i)+GlobalStepScale*IndStepScale*eigen_values_master[i];
    }
  }
  // Then update the parameter basis
  TransferToParam();
}

// ********************************************
// Transfer a parameter variation in the parameter basis to the eigen basis
void PCAHandler::TransferToPCA() {
// ********************************************
  // Make the temporary vectors
  TVectorD fParCurr_vec(static_cast<Int_t>(_pCurrVal->size()));
  TVectorD fParProp_vec(static_cast<Int_t>(_pCurrVal->size()));
  for(int i = 0; i < static_cast<int>(_pCurrVal->size()); ++i) {
    fParCurr_vec(i) = (*_pCurrVal)[i];
    fParProp_vec(i) = (*_pPropVal)[i];
  }

  _fParCurrPCA = TransferMatT*fParCurr_vec;
  _fParPropPCA = TransferMatT*fParProp_vec;
}

// ********************************************
void PCAHandler::SetInitialParameters() {
// ********************************************
  TransferToPCA();
  for (int i = 0; i < NumParPCA; ++i) {
    _fPreFitValuePCA[i] = _fParCurrPCA(i);
  }
}

// ********************************************
// Transfer a parameter variation in the eigen basis to the parameter basis
void PCAHandler::TransferToParam() {
// ********************************************
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  // Make the temporary vectors
  TVectorD fParProp_vec = TransferMat*_fParPropPCA;
  TVectorD fParCurr_vec = TransferMat*_fParCurrPCA;
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < static_cast<int>(_pCurrVal->size()); ++i) {
    (*_pPropVal)[i] = static_cast<M3::float_t>(fParProp_vec(i));
    (*_pCurrVal)[i] = fParCurr_vec(i);
  }
  #pragma GCC diagnostic pop
}

// ********************************************
void PCAHandler::Print() const {
// ********************************************
  MACH3LOG_INFO("PCA:");
  for (int i = 0; i < NumParPCA; ++i) {
    MACH3LOG_INFO("PCA {:<2} Current: {:<10.2f} Proposed: {:<10.2f}", i, _fParCurrPCA(i), _fParPropPCA(i));
  }
}

// ********************************************
void PCAHandler::SetBranches(TTree &tree, const bool SaveProposal, const std::vector<std::string>& Names) {
// ********************************************
  for (int i = 0; i < NumParPCA; ++i) {
    tree.Branch(Form("%s_PCA", Names[i].c_str()), &_fParCurrPCA.GetMatrixArray()[i], Form("%s_PCA/D", Names[i].c_str()));
  }

  if(SaveProposal)
  {
    for (int i = 0; i < NumParPCA; ++i) {
      tree.Branch(Form("%s_PCA_Prop", Names[i].c_str()), &_fParPropPCA.GetMatrixArray()[i], Form("%s_PCA_Prop/D", Names[i].c_str()));
    }
  }
}

// ********************************************
void PCAHandler::ToggleFixAllParameters(const std::vector<std::string>& Names) {
// ********************************************
  for (int i = 0; i < NumParPCA; i++) ToggleFixParameter(i, Names);
}

// ********************************************
void PCAHandler::ToggleFixParameter(const int i, const std::vector<std::string>& Names) {
// ********************************************
  int isDecom = -1;
  for (int im = 0; im < NumParPCA; ++im) {
    if(isDecomposedPCA[im] == i) {isDecom = im;}
  }
  if(isDecom < 0) {
    MACH3LOG_ERROR("Parameter {} is PCA decomposed can't fix this", Names[i]);
    //throw MaCh3Exception(__FILE__ , __LINE__ );
  } else {
    _fErrorPCA[isDecom] *= -1.0;
    if(IsParameterFixedPCA(i)) MACH3LOG_INFO("Setting un-decomposed {}(parameter {}/{} in PCA base) to fixed at {}", Names[i], i, isDecom, (*_pCurrVal)[i]);
    else MACH3LOG_INFO("Setting un-decomposed {}(parameter {}/{} in PCA base) free", Names[i], i, isDecom);
  }
}

// ********************************************
void PCAHandler::ThrowParameters(const std::vector<std::unique_ptr<TRandom3>>& random_number,
                                 double** throwMatrixCholDecomp,
                                 double* randParams,
                                 double* corr_throw,
                                 const std::vector<double>& fPreFitValue,
                                 const std::vector<double>& fLowBound,
                                 const std::vector<double>& fUpBound,
                                 const int _fNumPar) {
// ********************************************
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuseless-cast"
  //KS: Do not multithread!
  for (int i = 0; i < NumParPCA; ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (IsParameterFixedPCA(i)) continue;

    if(!IsParameterDecomposed(i))
    {
      (*_pPropVal)[i] = static_cast<M3::float_t>(fPreFitValue[i] + corr_throw[i]);
      int throws = 0;
      // Try again if we the initial parameter proposal falls outside of the range of the parameter
      while ((*_pPropVal)[i] > fUpBound[i] || (*_pPropVal)[i] < fLowBound[i]) {
        randParams[i] = random_number[M3::GetThreadIndex()]->Gaus(0, 1);
        const double corr_throw_single = M3::MatrixVectorMultiSingle(throwMatrixCholDecomp, randParams, _fNumPar, i);
        (*_pPropVal)[i] = static_cast<M3::float_t>(fPreFitValue[i] + corr_throw_single);
        if (throws > 10000)
        {
          //KS: Since we are multithreading there is danger that those messages
          //will be all over the place, small price to pay for faster code
          MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws, i);
          MACH3LOG_WARN("Setting _fPropVal:  {} to {}", (*_pPropVal)[i], fPreFitValue[i]);
          MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
         (*_pPropVal)[i] = static_cast<M3::float_t>(fPreFitValue[i]);
        }
        throws++;
      }
      (*_pCurrVal)[i] = (*_pPropVal)[i];

    } else {
      // KS: We have to multiply by number of parameters in PCA base
      SetParPropPCA(i, GetPreFitValuePCA(i) + randParams[i] * eigen_values_master[i] * (LastPCAdpar - FirstPCAdpar));
      SetParCurrPCA(i, GetParPropPCA(i));
    }
  } // end of parameter loop

  /// @todo KS: We don't check if param is out of bounds. This is more problematic for PCA params.
  for (int i = 0; i < _fNumPar; ++i) {
    auto low = static_cast<M3::float_t>(fLowBound[i]);
    auto up  = static_cast<M3::float_t>(fUpBound[i]);
    (*_pPropVal)[i] = std::max(up, std::min((*_pPropVal)[i], low));
    (*_pCurrVal)[i] = (*_pPropVal)[i];
  }
  #pragma GCC diagnostic pop
}
