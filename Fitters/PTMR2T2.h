#pragma once

#include "Fitters/MR2T2.h"

/*
https://arxiv.org/pdf/1205.1076
*/

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