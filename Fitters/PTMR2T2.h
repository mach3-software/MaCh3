#pragma once

#include "Fitters/MR2T2.h"
#include "numeric"

/*

https://academic.oup.com/jrsssb/article/84/2/321/7056147#396806364
https://arxiv.org/pdf/1205.1076
*/

/// @brief Parallel Tempered (Adaptive) MR2T2 algorith
/// @details
/// Run several chains in parallel at different temperatures
/// Temperature is a scales the posterior to
/// \f[
/// \pi(x)_{T} = \pi(x)^{1/T}
/// \f]
/// Chains of different temperatures can swap values with probability
/// \f[
/// \alpha = \min\left(1, \frac{\pi(x_i)^{1/T_j} \pi(x_j)^{1/T_i}}{\pi(x_i)^{1/T_i} \pi(x_j)^{1/T_j}}\right)
/// \f]
/// The algorithm aims to get a swap rates of ~25% which is roughly optimal
/// The only "true" chain is the \f$T=1\f$ chain which samples from the true posterior
///
/// @note This implementation is based on MR2T2 but could easily be adapted to other MCMC algorithms

/// @author Henry Wallace

#ifdef MPIENABLED
class PTMR2T2 : public MR2T2
{
  public:
    PTMR2T2(manager *const manager);
    ~PTMR2T2() = default;

    void RunMCMC() override;

  protected:
    void DoStep() override;
    void ProposeStep() override;

    void Swap();
    void SwapStepInformation(int swap_rank);

    void PrepareOutput() override;
    void AcceptStep() override;

    void PrintProgress() override;

    void GlobalAdaptTemperatures();

    // MPI stuff    

    // PTMCMC parameters
    double TempScale;
    double inv_TempScale; // Cache 1/TempScale for faster division
    double round_trip_rate;
    int temp_adapt_step;
    int n_temp_adapts;

    double LogLCurr_replica;
    double sample_LogLCurr_replica;
    double temp_replica;

    bool do_swap;
    int n_up_swap;
    int n_down_swap;

    std::vector<double> transfer_vector_send;
    std::vector<double> transfer_vector_recv;
    std::vector<double> temp_vec;

    // Timing benchmarks
    double time_base_dostep;
    double time_swap;
    double time_global_adapt;

    // need for swap
    double sample_logLProp;
    double sample_logLCurr;

};



#endif