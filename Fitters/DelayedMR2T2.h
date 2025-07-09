#pragma once

#include "Fitters/MR2T2.h"

class DelayedMR2T2 : public MR2T2 {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    DelayedMR2T2(manager * const fitMan);

    /// @brief Destructor
    virtual ~DelayedMR2T2() = default;

 protected:
    // Need to reset stuff
    void DoMCMCStep();

    void ProposeStep() override;
    double CheckStep() override;

    double CalculateLogGaussian(std::vector<double> &vec_1, std::vector<double> &vec_2, std::vector<std::vector<double>> &inv_cov_matrix) const _noexcept_;

    bool rejected_init_step;

    double initial_scale;
    double decay_rate;

    double accProbInit; // Rejected Acceptance probability

};
