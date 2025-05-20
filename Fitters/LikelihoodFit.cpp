#include "LikelihoodFit.h"

// *******************
// Run the Markov chain with all the systematic objects added
LikelihoodFit::LikelihoodFit(manager *man) : FitterBase(man) {
// *******************
  NPars = 0;
  NParsPCA = 0;
  fMirroring = GetFromManager<bool>(fitMan->raw()["General"]["Fitter"]["Mirroring"], false);
  if(fMirroring) MACH3LOG_INFO("Mirroring enabled");
}


// *************************
// Destructor: close the logger and output file
LikelihoodFit::~LikelihoodFit() {
// *************************

}

// *******************
void LikelihoodFit::PrepareFit() {
// *******************

  // Save the settings into the output file
  SaveSettings();

  // Prepare the output branches
  PrepareOutput();
  
  for (size_t s = 0; s < systematics.size(); ++s) {
    NPars += systematics[s]->GetNumParams();
    NParsPCA += systematics[s]->GetNParameters();
  }

  //KS: If PCA is note enabled NParsPCA == NPars
  MACH3LOG_INFO("Total number of parameters {}", NParsPCA);
}

// *******************
double LikelihoodFit::CalcChi2(const double* x) {
// *******************

  if (step % 10000 == 0)
  {
    MACH3LOG_INFO("Iteration {}", step);
  }

  stepClock->Start();

  int ParCounter = 0;
  double llh = 0;
  for (std::vector<ParameterHandlerBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if(!(*it)->IsPCA())
    {
      std::vector<double> pars;
      const int NumPar = (*it)->GetNumParams();
      //KS: Avoid push back as they are slow
      pars.resize(NumPar);
      for(int i = 0; i < NumPar; ++i, ++ParCounter)
      {
        double ParVal = x[ParCounter];
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
        pars[i] = ParVal;
      }
      (*it)->SetParameters(pars);
    }
    else
    {
      std::vector<double> pars;
      const int NumPar = (*it)->GetNParameters();
      //KS: Avoid push back as they are slow
      pars.resize(NumPar);
      for(int i = 0; i < NumPar; ++i, ++ParCounter)
      {
        double ParVal = x[ParCounter];
        //KS: Basically apply mirroring for parameters out of bounds
        pars[i] = ParVal;
      }
      (*it)->GetPCAHandler()->SetParametersPCA(pars);
    }
    (*it)->AcceptStep();
  }

  // Loop over the systematics and propose the initial step
  int stdIt = 0;
  for (std::vector<ParameterHandlerBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it, ++stdIt)
  {
    //GetLikelihood will return LargeL if out of bounds, for minimizers this is not the problem, while calcLikelihood will return actual likelihood
    syst_llh[stdIt] = (*it)->CalcLikelihood();
    llh += syst_llh[stdIt];
    #ifdef DEBUG
    if (debug) debugFile << "LLH after " << systematics[stdIt]->GetName() << " " << llh << std::endl;
    #endif
  }
  // Could multi-thread this
  // But since sample reweight is multi-threaded it's probably better to do that
  for (size_t i = 0; i < samples.size(); i++)
  {
    samples[i]->Reweight();
  }

  //DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
  for (size_t i = 0; i < samples.size(); i++) {
    // Get the sample likelihoods and add them
    sample_llh[i] = samples[i]->GetLikelihood();
    llh += sample_llh[i];
    #ifdef DEBUG
    if (debug) debugFile << "LLH after sample " << i << " " << llh << std::endl;
    #endif
  }

  // Save the proposed likelihood (class member)
  logLProp = llh;
  logLCurr = llh;
  accProb = 1;

  stepClock->Stop();
  stepTime = stepClock->RealTime();

  // Write step to output tree
  outTree->Fill();

  // Auto save the output
  if (step % auto_save == 0) outTree->AutoSave();
  step++;
  accCount++;

  llh = 2.0*llh;
  return llh;
}

