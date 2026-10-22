#pragma once
/*
HW: Time to try out HMC
*/

#include "Algorithms/MCMCBase.h"

class HamiltonianMC : public MCMCBase {
public:
    HamiltonianMC(Manager *const FitManager);

    virtual ~HamiltonianMC() = default;

protected:
    void DoStep() override;
    void ProposeStep() override;
    
    double AcceptanceProbability() override;
    
    void RunMCMC();

    void UpdateMomentum();
    void UpdateGradU();
    
    int n_leapfrog;
    float delta_t;

    double ke_init;

    std::vector<double> momentum;
    std::vector<double> gradient;
};