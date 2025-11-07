#pragma once

#include "Fitters/MCMCBase.h"

/// @brief MCMC algorithm implementing the Metropolis–Rosenbluth–Rosenbluth–Teller–Teller (MR\f$^2\f$T\f$^2\f$) method.
/// @note MR\f$^2\f$T\f$^2\f$ is also known as Metropolis-Hastings
/// @cite Hastings:1970aa
///
/// @details
/// The acceptance probability for a proposed step is given by
/// \f[
/// \alpha = \min \Biggl( 1, \exp\bigl( \log \mathcal{L}_\text{prop} - \log \mathcal{L}_\text{curr} \bigr) \Biggr),
/// \f]
/// where \f$\log \mathcal{L}_\text{prop}\f$ is the log-likelihood of the proposed state and
/// \f$\log \mathcal{L}_\text{curr}\f$ is the log-likelihood of the current state.
///
/// If annealing is enabled, the acceptance probability is modified as
/// \f[
/// \alpha_\text{anneal} = \min \Biggl( 1, \exp\Bigl( -\frac{\log \mathcal{L}_\text{prop} - \log \mathcal{L}_\text{curr}}
/// {\exp(-\text{step}/T_\text{anneal})} \Bigr) \Biggr),
/// \f]
/// where \f$T_\text{anneal}\f$ is the current annealing temperature and \f$\text{step}\f$ is the iteration index.
///
/// @author Asher Kaboth
class MR2T2 : public MCMCBase {
 public:
    /// @brief Constructor
    /// @param manager A pointer to a manager object, which will handle all settings.
     MR2T2(manager *const manager);

     /// @brief Destructor
     virtual ~MR2T2() = default;
 protected:
    /// @brief The MCMC step proposal and acceptance
    void DoStep() override;
    /// @brief Propose a step
    void ProposeStep() override;
    /// @brief Step acceptance probability
    double AcceptanceProbability() override;
};

/// For backwards compatibility
using mcmc  = MR2T2;
