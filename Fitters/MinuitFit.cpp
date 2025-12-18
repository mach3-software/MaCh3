#include "MinuitFit.h"

// *******************
// Run the Minuit Fit with all the systematic objects added
MinuitFit::MinuitFit(manager *man) : LikelihoodFit(man) {
  // *******************
  AlgorithmName = "MinuitFit";
  /// @todo KS: Make this in future configurable, for more see:
  /// https://root.cern.ch/doc/master/classROOT_1_1Math_1_1Minimizer.html
  // Minimizer type: determines the underlying implementation.
  // Available types include:
  //   - "Minuit2" (recommended modern option)
  //   - "Minuit"  (legacy)
  //   - "Fumili"
  //   - "GSLMultiMin" (for gradient-free minimization)
  //   - "GSLMultiFit"
  //   - "GSLSimAn" (Simulated Annealing)
  const std::string MinimizerType = "Minuit2";
  // Minimizer algorithm (specific to the selected type).
  // For Minuit2, the following algorithms are available:
  //   - "Migrad"   : gradient-based minimization (default)
  //   - "Simplex"  : Nelder-Mead simplex method (derivative-free)
  //   - "Combined" : combination of Simplex and Migrad
  //   - "Scan"     : parameter grid scan
  const std::string MinimizerAlgo = "Migrad";

  MACH3LOG_INFO("Creating instance of Minimizer with {} and {}", MinimizerType,
                MinimizerAlgo);

  minuit = std::unique_ptr<ROOT::Math::Minimizer>(
      ROOT::Math::Factory::CreateMinimizer(MinimizerType.c_str(),
                                           MinimizerAlgo.c_str()));
}

// *************************
// Destructor: close the logger and output file
MinuitFit::~MinuitFit() {
  // *************************
}

// *******************
// Run the Minuit with all the systematic objects added
void MinuitFit::RunMCMC() {
  // *******************
  PrepareFit();

  // Remove obsolete memory and make other checks before fit starts
  SanitiseInputs();

  // KS: Set SetFunction we will Minimize
  ROOT::Math::Functor fChi2(this, &MinuitFit::CalcChi2PC, NParsPCA);
  minuit->SetFunction(fChi2);

  // KS: add config or something
  minuit->SetPrintLevel(2);
  minuit->SetTolerance(0.01);
  minuit->SetMaxFunctionCalls(
      fitMan->raw()["General"]["Minuit2"]["NSteps"].as<unsigned>());
  minuit->SetMaxIterations(10000);

  MACH3LOG_INFO("Preparing Minuit");
  int ParCounter = 0;

  for (auto &parhandlr : systematics) {

    for (int i = 0; i < parhandlr->GetNumProposalParams(); ++i, ++ParCounter) {
      // KS: Index, name, prior, step scale [different to MCMC],
      minuit->SetVariable(ParCounter, (parhandlr->GetPCParName(i)),
                          parhandlr->GetPCParInit(i),
                          parhandlr->GetPCDiagonalError(i) / 10.0);
      minuit->SetVariableValue(ParCounter, parhandlr->GetPCParInit(i));
      // KS: lower bound, upper bound, if Mirroring enabled then ignore
      if (!fMirroring)
        minuit->SetVariableLimits(ParCounter, parhandlr->GetPCLowerBound(i),
                                  parhandlr->GetPCUpperBound(i));
      if (parhandlr->IsPCParameterFixed(i)) {
        minuit->FixVariable(ParCounter);
      }
    }
  }

  minuit->SetPrintLevel(2);

  MACH3LOG_INFO("Starting MIGRAD");
  minuit->Minimize();

  MACH3LOG_INFO("Starting HESSE");
  minuit->Hesse();
  outputFile->cd();

  TVectorD *MinuitParValue = new TVectorD(NParsPCA);
  TVectorD *MinuitParError = new TVectorD(NParsPCA);
  TMatrixDSym *Postmatrix = new TMatrixDSym(NParsPCA);

  for (int i = 0; i < NParsPCA; ++i) {
    (*MinuitParValue)(i) = 0;
    (*MinuitParError)(i) = 0;
    for (int j = 0; j < NParsPCA; ++j) {
      (*Postmatrix)(i, j) = 0;
      (*Postmatrix)(i, j) = minuit->CovMatrix(i, j);
    }
  }

  ParCounter = 0;
  const double *X = minuit->X();
  const double *err = minuit->Errors();
  for (auto &parhandlr : systematics) {

    // set the parameters in the sampler basis
    std::copy_n(X + ParCounter, parhandlr->GetNumProposalParams(),
                parhandlr->proposer.params.proposed.data());

    // this accepts the step for the proposer and rotates the parameters back to
    // the systematic basis
    parhandlr->AcceptStep();

    for (int i = 0; i < parhandlr->GetNumSystematicParams(); ++i) {
      double ParVal = parhandlr->GetParCurr(i);
      // KS: Basically apply mirroring for parameters out of bounds
      if (fMirroring) { // mirror in the systematic basis as it is where the
                        // bounds are defined
        if (ParVal < parhandlr->GetLowerBound(i)) {
          ParVal = parhandlr->GetLowerBound(i) +
                   (parhandlr->GetLowerBound(i) - ParVal);
        } else if (ParVal > parhandlr->GetUpperBound(i)) {
          ParVal = parhandlr->GetUpperBound(i) -
                   (ParVal - parhandlr->GetUpperBound(i));
        }
        parhandlr->SetParCurrProp(i, ParVal);
      }
    }

    for (int i = 0; i < parhandlr->GetNumProposalParams(); ++i) {

      (*MinuitParValue)(ParCounter + i) =
          parhandlr->proposer.params.proposed[i];
      (*MinuitParError)(ParCounter + i) = err[ParCounter + i];

      // KS: For fixed params HESS will not calculate error so we need to pass
      // prior error
      if (parhandlr->IsPCParameterFixed(i)) {
        (*MinuitParError)(ParCounter) = parhandlr->GetPCDiagonalError(i);
        (*Postmatrix)(ParCounter, ParCounter) =
            (*MinuitParError)(ParCounter) * (*MinuitParError)(ParCounter);
      }
    }

    ParCounter += parhandlr->GetNumProposalParams();
  }

  MinuitParValue->Write("MinuitParValue");
  MinuitParError->Write("MinuitParError");
  Postmatrix->Write("Postmatrix");
  delete MinuitParValue;
  delete MinuitParError;
  delete Postmatrix;
  // Save all the output
  SaveOutput();
}
