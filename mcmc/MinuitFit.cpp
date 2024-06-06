#include "MinuitFit.h"

// *******************
// Run the Markov chain with all the systematic objects added
MinuitFit::MinuitFit(manager *man) : LikelihoodFit(man) {
// *******************

  MACH3LOG_INFO("Creating instance of Minimizer with Minuit2 and Migrad");

  //Other then Migrad available are Simplex,Combined,Scan
  minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
}

// *************************
// Destructor: close the logger and output file
MinuitFit::~MinuitFit() {
// *************************
  if(minuit != NULL) delete minuit;
}


// *******************
// Run the Minuit with all the systematic objects added
void MinuitFit::runMCMC() {
// *******************

  PrepareFit();

  //KS: For none PCA this will be equal to normal parameters
  const int NparsMinuitFull = NPars;
  const int NparsMinuit = NParsPCA;

  //KS: Set SetFunction we will Minimize
  ROOT::Math::Functor fChi2(this, &MinuitFit::CalcChi2, NparsMinuit);
  minuit->SetFunction(fChi2);

  //KS: add config or something
  minuit->SetPrintLevel(2);
  minuit->SetTolerance(0.01);
  minuit->SetMaxFunctionCalls(fitMan->raw()["General.NSteps"].as<double>());
  minuit->SetMaxIterations(10000);

  MACH3LOG_INFO("Preparing Minuit");
  int ParCounter = 0;

  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if(!(*it)->IsPCA())
    {
      for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
      {
        //KS: Index, name, prior, step scale [differrent to MCMC],
        minuit->SetVariable(ParCounter, ((*it)->GetParName(i)), (*it)->getParInit(i), (*it)->getDiagonalError(i)/10);
        minuit->SetVariableValue(ParCounter, (*it)->getParInit(i));
        //KS: lower bound, upper bound, if Mirroring eneabled then ignore
        if(!fMirroring) minuit->SetVariableLimits(ParCounter, (*it)->GetLowerBound(i), (*it)->GetUpperBound(i));
        if((*it)->isParameterFixed(i))
        {
          minuit->FixVariable(ParCounter);
        }
      }
    }
    else
    {
      for(int i = 0; i < (*it)->getNpars(); ++i, ++ParCounter)
      {
        minuit->SetVariable(ParCounter, Form("%i_PCA", i), (*it)->getParProp_PCA(i), (*it)->getEigenValuesMaster()[i]/10);
        if((*it)->isParameterFixedPCA(i))
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
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
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
        //KS: For fixed params HESS will not calcuate error so we need to pass prior error
        if((*it)->isParameterFixed(i))
        {
          (*MinuitParError)(ParCounter) = (*it)->getDiagonalError(i);
          (*Postmatrix)(ParCounter,ParCounter) = (*MinuitParError)(ParCounter) * (*MinuitParError)(ParCounter);
        }
      }
    }
    else
    {
      //KS: We need to convert parameters from PCA to normal base
      TVectorD ParVals((*it)->GetNumParams());
      TVectorD ParVals_PCA((*it)->getNpars());

      TVectorD ErrorVals((*it)->GetNumParams());
      TVectorD ErrorVals_PCA((*it)->getNpars());

      TMatrixD MatrixVals((*it)->GetNumParams(), (*it)->GetNumParams());
      TMatrixD MatrixVals_PCA((*it)->getNpars(), (*it)->getNpars());

      //First save them
      //KS: This code is super convoluted as MaCh3 can store separate matrices while Minuit has one matrix. In future this will be simplified, keep it like this for now.
      const int StartVal = ParCounter;
      for(int i = 0; i < (*it)->getNpars(); ++i, ++ParCounter)
      {
        ParVals_PCA(i) = X[ParCounter];
        ErrorVals_PCA(i) = err[ParCounter];
        int ParCounterMatrix = StartVal;
        for(int j = 0; j < (*it)->getNpars(); ++j, ++ParCounterMatrix)
        {
          MatrixVals_PCA(i,j) = minuit->CovMatrix(ParCounter,ParCounterMatrix);
        }
      }
      ParVals = ((*it)->getTransferMatrix())*ParVals_PCA;
      ErrorVals = ((*it)->getTransferMatrix())*ErrorVals_PCA;
      MatrixVals.Mult(((*it)->getTransferMatrix()),MatrixVals_PCA);

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
        if((*it)->isParameterFixedPCA(i))
        {
          (*MinuitParError)(ParCounter) = (*it)->getDiagonalError(i);
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

