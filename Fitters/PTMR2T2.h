#pragma once

#include "Fitters/MR2T2.h"

/*


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
    void SynchronizeTemperatures();
    void AdaptTemperature();

    // MPI stuff    

    // PTMCMC parameters
    double temperature;
    double LogL_replica;
    double temp_replica;
    double adaption_rate;


    bool do_swap;
    int n_up_swap;
    int n_down_swap;
    double max_temp;
    double min_temp;
    double target_rate;

    std::vector<double> transfer_vector_send;
    std::vector<double> transfer_vector_recv;

    std::vector<double> temp_vec;
};



#endif