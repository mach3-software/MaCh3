#pragma once

#include "Fitters/MR2T2.h"

class DelayedMR2T2 : public MR2T2 {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    DelayedMR2T2(manager * const fitMan);

    /// @brief Destructor
    virtual ~DelayedMR2T2() = default;
    inline std::string GetName() const { return "DelayedMR2T2"; };

 protected:

   double AcceptanceProbability() override;

   void DoStep();

   void StoreCurrentStep();

   void ResetSystScale();
   void ScaleSystematics(double scale);

   double CalculateLogGaussian(std::vector<double> &vec_1, std::vector<double> &vec_2, std::vector<std::vector<double>> &inv_cov_matrix) const _noexcept_;

   std::vector<std::vector<double>> current_step_vals;

   double initial_scale;
   double decay_rate;
   int number_of_iterations;

   double MinLogLikelihood; // Max Likelihood


   std::vector<double> step_syst_scale;

};
