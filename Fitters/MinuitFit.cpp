#include "MinuitFit.h"

// *******************
// Run the Minuit Fit with all the systematic objects added
MinuitFit::MinuitFit(manager *man) : LikelihoodFit(man) {
// *******************
  AlgorithmName = "MinuitFit";
  /// @todo KS: Make this in future configurable, for more see: https://root.cern.ch/doc/master/classROOT_1_1Math_1_1Minimizer.html
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

  MACH3LOG_INFO("Creating instance of Minimizer with {} and {}", MinimizerType, MinimizerAlgo);

  minuit = std::unique_ptr<ROOT::Math::Minimizer>(
    ROOT::Math::Factory::CreateMinimizer(MinimizerType.c_str(), MinimizerAlgo.c_str()));
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

  //KS: For none PCA this will be equal to normal parameters
  const int NparsMinuitFull = NPars;
  const int NparsMinuit = NParsPCA;

  //KS: Set SetFunction we will Minimize
  ROOT::Math::Functor fChi2(this, &MinuitFit::CalcChi2, NparsMinuit);
  minuit->SetFunction(fChi2);

  //KS: add config or something
  minuit->SetPrintLevel(2);
  minuit->SetTolerance(0.01);
  minuit->SetMaxFunctionCalls(fitMan->raw()["General"]["Minuit2"]["NSteps"].as<unsigned>());
  minuit->SetMaxIterations(10000);

  MACH3LOG_INFO("Preparing Minuit");
  int ParCounter = 0;

  for (std::vector<ParameterHandlerBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if(!(*it)->IsPCA())
    {
      for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
      {
        //KS: Index, name, prior, step scale [different to MCMC],
        minuit->SetVariable(ParCounter, ((*it)->GetParName(i)), (*it)->GetParInit(i), (*it)->GetDiagonalError(i)/10);
        minuit->SetVariableValue(ParCounter, (*it)->GetParInit(i));
        //KS: lower bound, upper bound, if Mirroring enabled then ignore
        if(!fMirroring) minuit->SetVariableLimits(ParCounter, (*it)->GetLowerBound(i), (*it)->GetUpperBound(i));
        if((*it)->IsParameterFixed(i))
        {
          minuit->FixVariable(ParCounter);
        }
      }
    }
    else
    {
      for(int i = 0; i < (*it)->GetNParameters(); ++i, ++ParCounter)
      {
        minuit->SetVariable(ParCounter, Form("%i_PCA", i), (*it)->GetPCAHandler()->GetParPropPCA(i), (*it)->GetPCAHandler()->GetEigenValuesMaster()[i]/10);
        if((*it)->GetPCAHandler()->IsParameterFixedPCA(i))
        {
          minuit->FixVariable(ParCounter);
        }
      }
    }
  }

  minuit->SetPrintLevel(2);

  MACH3LOG_INFO("Starting MIGRAD");
  minuit->Minimize();

  MACH3LOG_INFO("Starting HESSE");
  minuit->Hesse();
  outputFile->cd();
  
  TVectorD* MinuitParValue = new TVectorD(NparsMinuitFull);
  TVectorD* MinuitParError = new TVectorD(NparsMinuitFull);
  TMatrixDSym* Postmatrix = new TMatrixDSym(NparsMinuitFull);

  for(int i = 0; i < NparsMinuitFull; ++i)
  {
    (*MinuitParValue)(i) = 0;
    (*MinuitParError)(i) = 0;
    for(int j = 0; j < NparsMinuitFull; ++j)
    {
      (*Postmatrix)(i,j) = 0;
      (*Postmatrix)(i,j) = minuit->CovMatrix(i,j);
    }
  }

  ParCounter = 0;
  const double *X = minuit->X();
  const double *err = minuit->Errors();
  for (std::vector<ParameterHandlerBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if(!(*it)->IsPCA())
    {
      for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
      {
        double ParVal = X[ParCounter];
        //KS: Basically apply mirroring for parameters out of bounds
        if(fMirroring)
        {
          if(ParVal < (*it)->GetLowerBound(i))
          {
            ParVal = (*it)->GetLowerBound(i) + ((*it)->GetLowerBound(i) - ParVal);
          }
          else if (ParVal > (*it)->GetUpperBound(i))
          {
            ParVal = (*it)->GetUpperBound(i) - ( ParVal - (*it)->GetUpperBound(i));
          }
        }
        (*MinuitParValue)(ParCounter) = ParVal;
        (*MinuitParError)(ParCounter) = err[ParCounter];
        //KS: For fixed params HESS will not calculate error so we need to pass prior error
        if((*it)->IsParameterFixed(i))
        {
          (*MinuitParError)(ParCounter) = (*it)->GetDiagonalError(i);
          (*Postmatrix)(ParCounter,ParCounter) = (*MinuitParError)(ParCounter) * (*MinuitParError)(ParCounter);
        }
      }
    }
    else
    {
      //KS: We need to convert parameters from PCA to normal base
      TVectorD ParVals((*it)->GetNumParams());
      TVectorD ParVals_PCA((*it)->GetNParameters());

      TVectorD ErrorVals((*it)->GetNumParams());
      TVectorD ErrorVals_PCA((*it)->GetNParameters());

      TMatrixD MatrixVals((*it)->GetNumParams(), (*it)->GetNumParams());
      TMatrixD MatrixVals_PCA((*it)->GetNParameters(), (*it)->GetNParameters());

      //First save them
      //KS: This code is super convoluted as MaCh3 can store separate matrices while Minuit has one matrix. In future this will be simplified, keep it like this for now.
      const int StartVal = ParCounter;
      for(int i = 0; i < (*it)->GetNParameters(); ++i, ++ParCounter)
      {
        ParVals_PCA(i) = X[ParCounter];
        ErrorVals_PCA(i) = err[ParCounter];
        int ParCounterMatrix = StartVal;
        for(int j = 0; j < (*it)->GetNParameters(); ++j, ++ParCounterMatrix)
        {
          MatrixVals_PCA(i,j) = minuit->CovMatrix(ParCounter,ParCounterMatrix);
        }
      }
      ParVals = ((*it)->GetPCAHandler()->GetTransferMatrix())*ParVals_PCA;
      ErrorVals = ((*it)->GetPCAHandler()->GetTransferMatrix())*ErrorVals_PCA;
      MatrixVals.Mult(((*it)->GetPCAHandler()->GetTransferMatrix()),MatrixVals_PCA);

      ParCounter = StartVal;
      //KS: Now after going from PCA to normal let';s save it
      for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
      {
        (*MinuitParValue)(ParCounter) = ParVals(i);
        (*MinuitParError)(ParCounter) = std::fabs(ErrorVals(i));
        int ParCounterMatrix = StartVal;
        for(int j = 0; j < (*it)->GetNumParams(); ++j, ++ParCounterMatrix)
        {
          (*Postmatrix)(ParCounter,ParCounterMatrix)  = MatrixVals(i,j);
        }
        //If fixed take prior
        if((*it)->GetPCAHandler()->IsParameterFixedPCA(i))
        {
          (*MinuitParError)(ParCounter) = (*it)->GetDiagonalError(i);
          (*Postmatrix)(ParCounter,ParCounter) = (*MinuitParError)(ParCounter) * (*MinuitParError)(ParCounter);
        }
      }
    }
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

