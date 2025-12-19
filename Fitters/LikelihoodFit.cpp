#include "LikelihoodFit.h"

// *******************
// Run the Markov chain with all the systematic objects added
LikelihoodFit::LikelihoodFit(manager *man) : FitterBase(man) {
  // *******************
  fMirroring = GetFromManager<bool>(
      fitMan->raw()["General"]["Fitter"]["Mirroring"], false);
  if (fMirroring)
    MACH3LOG_INFO("Mirroring enabled");
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

  for (auto phandlr : systematics) {
    NPars += phandlr->GetNumSystematicParams();
    NParsPCA += phandlr->GetNumProposalParams();
  }

  MACH3LOG_INFO("Total number of parameters {}", NPars);
}

// parameters come in in the sampler basis
//  *******************

double LikelihoodFit::CalcChi2PC(const double *x) {

  int ParCounter = 0;
  static std::vector<double> syst_par_values(NPars, 0);

  auto par_it = syst_par_values.begin();
  for (auto &parhandlr : systematics) {

    // set the parameters in the sampler basis
    std::copy_n(x + ParCounter, parhandlr->GetNumProposalParams(),
                parhandlr->proposer.proposal_basis.proposed.data());
    ParCounter += parhandlr->GetNumProposalParams();

    // this accepts the step for the proposer and rotates the parameters back to
    // the systematic basis
    parhandlr->AcceptStep();

    std::copy_n(parhandlr->GetParCurrVec().begin(),
                parhandlr->GetNumSystematicParams(), par_it);
    std::advance(par_it, parhandlr->GetNumSystematicParams());
  }

  return CalcChi2(syst_par_values.data());
}

double LikelihoodFit::CalcChi2(const double *x) {
  // *******************
  if (step % 10000 == 0) {
    MACH3LOG_INFO("Iteration {}", step);
  }

  stepClock->Start();

  int ParCounter = 0;
  for (auto &parhandlr : systematics) {
    for (int i = 0; i < parhandlr->GetNumSystematicParams(); ++i) {
      double ParVal = x[ParCounter++];
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
      }
      parhandlr->SetParCurrProp(i, ParVal);
    }
  }

  double llh = 0;

  // Loop over the systematics and propose the initial step
  int stdIt = 0;
  for (auto &parhandlr : systematics) {

    // GetLikelihood will return LargeL if out of bounds, for minimizers this is
    // not the problem, while calcLikelihood will return actual likelihood
    syst_llh[stdIt] = parhandlr->CalcLikelihood();
    llh += syst_llh[stdIt];
#ifdef DEBUG
    if (debug)
      debugFile << "LLH after " << systematics[stdIt]->GetName() << " " << llh
                << std::endl;
#endif
    stdIt++;
  }
  // Could multi-thread this
  // But since sample reweight is multi-threaded it's probably better to do that
  for (size_t i = 0; i < samples.size(); i++) {
    samples[i]->Reweight();
  }

  // DB for atmospheric event by event sample migration, need to fully reweight
  // all samples to allow event passing prior to likelihood evaluation
  for (size_t i = 0; i < samples.size(); i++) {
    // Get the sample likelihoods and add them
    sample_llh[i] = samples[i]->GetLikelihood();
    llh += sample_llh[i];
#ifdef DEBUG
    if (debug)
      debugFile << "LLH after sample " << i << " " << llh << std::endl;
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
  if (step % auto_save == 0)
    outTree->AutoSave();
  step++;
  accCount++;

  llh = 2.0 * llh;
  return llh;
}
