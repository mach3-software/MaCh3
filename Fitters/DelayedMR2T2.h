#pragma once

#include "Fitters/MR2T2.h"

/// @brief Implementation of delayed rejection for MR2T2
/// @author Henry Wallace
/// @cite dram2006
///
/// @section intro Quick Introduction
///
/// Delayed rejection is a method to increase the efficiency of MCMC.
/// It is useful if your chain is not mixing well or if adaptive MCMC is struggling.
/// The idea is to propose increasingly small steps centered on the last rejected step.
///
/// @section config Configuration
///
/// To enable delayed rejection, add the following to your config.yaml file:
/// - Set the algorithm as: General::FittingAlgorithm::DelayedMR2T2
///
/// Delayed MCMC will work "as is" but there are additional options
/// which are set in General::MCMC
///
/// @subsection scaling Scaling Options
///
/// @param DecayRate This should be less than 1 and tells your algorithm how quickly
///                  to decrease the step size. By default this is 0.1
/// @param MaxRejections How many steps to delayed-reject before moving on. By default
///                      this is 1 and increasing is not recommended
/// @param InitialScale The initial scale of your proposed steps. I.e. propose step with
///                     InitialScale*scale in parameter handler. This allows for a really
///                     large/small first step to be proposed. Testing finds the default
///                     value of 1 to be most effective with adaptive MCMC
///
/// @subsection triggers Delay Triggers
///
/// @param DelayOnlyOutBounds Here delaying only happens if a step is thrown out of bounds
/// @param DelayProbability In order to mitigate the slowness of delayed MCMC, there is
///                         a probability of proposing a delayed step. By default this is 1
///
/// @note Delayed rejection MCMC can significantly improve chain mixing but may increase
///       computational overhead per step.
class DelayedMR2T2 : public MR2T2 {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    DelayedMR2T2(manager * const fitMan);

    /// @brief Destructor
    virtual ~DelayedMR2T2() = default;

    /// @brief Get name
    inline std::string GetName() const override { return "DelayedMR2T2"; };

 protected:
     /// @brief Step acceptance probability
    double AcceptanceProbability() override;

    /// @brief The MCMC step proposal
    void DoStep() override;

    /// @brief Store information about the current proposed step
    void StoreCurrentStep();

    /// @brief override to add in rejection information to output
    void PrepareOutput();

    /// @brief Reset scale after delay process
    void ResetSystScale();

    /// @brief Set parameter scale to be original_scale* scale
    /// @param scale Multiplicative scale factor
    void ScaleSystematics(const double scale);

    /// @brief if delay has a chance of occurring
    bool ProbabilisticDelay() const;

    /// Store information about the current proposed step
    std::vector<std::vector<double>> current_step_vals;

    /// Stores the initial scale for all parameters
    std::vector<double> start_step_scale;

    /// Scale for (non-delayed) step is #initial_scale * #start_step_scale
    double initial_scale;
    /// How much to decrease the step scale each step
    double decay_rate;
    /// How many rejections allowed?
    int max_rejections;

    /// Minimum (negative) log-likelihood of proposed steps 
    double MinLogLikelihood; // Max Likelihood

    /// Was the step we just accepted delayed?
    bool accepted_delayed;

    /// Delay only if we go out of bounds
    bool delay_on_oob_only;
    /// Can delay with probability instead
    double delay_probability;

    /// Continue flipping parameters on reject
    bool keep_flipping_on_reject;
    /// bool to store the INITIAL flip setting
    std::vector<bool> initial_flip_setting;

};
