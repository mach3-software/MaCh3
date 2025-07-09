#include "Fitters/DelayedMR2T2.h"

DelayedMR2T2::DelayedMR2T2(manager * const fitMan) : MR2T2(fitMan) {
    rejected_init_step = false;

    initial_scale = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["InitialStepScale"], 1.0);
    decay_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["DecayRate"], 0.1);

    accProbInit = 0;
    rejected_init_step = false;
}

void DelayedMR2T2::DoMCMCStep() {
    // Reset the rejected flag
    rejected_init_step = false;
    // Call the base class method to do the MCMC step
    ScaleSystematics(initial_scale);
    MCMCBase::DoMCMCStep();
}

void ScaleSystematics(std::vector<ParameterHandlerBase*> &systematics, double scale) {
    for (auto &syst : systematics) {
        syst->SetStepScale(scale*syst->GetGlobalStepScale(), false);
    }
}


void DelayedMR2T2::ProposeStep() {
    // Firstly we do initial MR2T2 step
    MR2T2::ProposeStep();

    if (reject)
    {
        accProbInit = 0;
        rejected_init_step = true
    }
    else{
        accProbInit = MR2T2::CheckStep();
        rejected_init_step = AcceptStep(accProbInit);
    }

    if (not rejected_init_step){
        return;
    }


    // Now we do delayed rejcection!
    // First we scale our steps to init*decay*global
    ScaleSystematics(systematics, initial_scale*decay_rate);

    // We now do another MR2T2 proposal
    MR2T2::ProposeStep();
}

double DelayedMR2T::CheckStep()
{
    // We've calculated this already!
    if(not rejected_init_step){
        return accProbInit;
    }

    /// Dealyed MCMC has a more complicated accept probability
    // This uses the TWO step delayed rejection algorithm
    // This is : min(1, p(y_2)q(y_2, y_1)(1-a(y_2, y_1)/p(x)q(x, y_1)(1-a(x,y_1)))

    double numerator = logLProb;
    double denominator = logLCurr;

    for(int i=0; i<static_cast<int>(systematics.size()); ++i){
        auto prop_vals = systematics[i]->GetParPropVec();
        auto curr_vals = systematics[i]->GetParCurrVec();

        // Get the inverse covariance matrix
        auto inv_cov_matrix = systematics[i]->GetThrowMatrixInv();

        // Calculate the log Gaussian for the proposal
        numerator += CalculateLogGaussian(prop_vals, rejected_step_vals[i], inv_cov_matrix);
        // Calculate the log Gaussian for the current values
        denominator += CalculateLogGaussian(curr_vals, prop_vals, inv_cov_matrix);
    }
    denominator += 1 - accProbInit; // Add the acceptance probability of the initial step
    if (anneal) {
        return TMath::Min(1., TMath::Exp((numerator - denominator) / (TMath::Exp(-step / AnnealTemp))));
    }
    else {
        return TMath::Min(1., TMath::Exp(numerator - denominator));
    }
}

double DelayedMR2T2::CalculateLogGaussian(std::vector<double> &vec_1, std::vector<double> &vec_2, std::vector<std::vector<double>> &inv_cov_matrix) const _noexcept_
{
    double logL = 0.0;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+ : logL)
#endif
    for (int i = 0; i < static_cast<int>(vec_1.size()); ++i)
    {
#ifdef MULTITHREAD
#pragma omp simd
#endif
        for (int j = 0; j <= i; ++j)
        {
            // KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
            double scale = (i != j) ? 1. : 0.5;
            logL += scale * (vec_1[i] - vec_2[i]) * (vec_1[j] - vec_2[j]) * inv_cov_matrix[i][j];
        }
    }
    return logL;
}
