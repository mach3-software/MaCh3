#include "LikelihoodFit.h"

// *******************
// Run the Markov chain with all the systematic objects added
LikelihoodFit::LikelihoodFit(manager *man) : FitterBase(man) {
// *******************
    NPars = 0;
    NParsPCA = 0;
    fMirroring = true;
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
    NParsPCA += systematics[s]->getNpars();
  }

  if (osc) {
    std::cerr<<" Osc not supported "<<std::endl;
    throw;
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
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
  {
    if(!(*it)->IsPCA())
    {
      std::vector<double> pars;
      const int Size = (*it)->GetNumParams();
      //KS: Avoid push back as they are slow
      pars.resize(Size);
      for(int i = 0; i < Size; ++i, ++ParCounter)
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
      (*it)->setParameters(pars);
    }
    else
    {
      std::vector<double> pars;
      const int Size = (*it)->getNpars();
      //KS: Avoid push back as they are slow
      pars.resize(Size);
      for(int i = 0; i < Size; ++i, ++ParCounter)
      {
        double ParVal = x[ParCounter];
        //KS: Basically apply mirroring for parameters out of bounds
        pars[i] = ParVal;
      }
      (*it)->setParameters_PCA(pars);
    }
    (*it)->acceptStep();
  }

  // Loop over the systematics and propose the initial step
  int stdIt = 0;
  for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it, ++stdIt)
  {
    //GetLikelihood will return LargeL if out of bounds, for minimizers this is not the problem, while calcLikelihood will return actual likelihood
    syst_llh[stdIt] = (*it)->CalcLikelihood();
    llh += syst_llh[stdIt];
    #ifdef DEBUG
    if (debug) debugFile << "LLH after " << systematics[stdIt]->getName() << " " << llh << std::endl;
    #endif
  }
  // Could multi-thread this
  // But since sample reweight is multi-threaded it's probably better to do that
  for (size_t i = 0; i < samples.size(); i++)
  {
    // If we're running with different oscillation parameters for neutrino and anti-neutrino
    if (osc) {
      samples[i]->reweight(osc->getPropPars());
      // If we aren't using any oscillation
      } else {
        double* fake = NULL;
        samples[i]->reweight(fake);
    }
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

