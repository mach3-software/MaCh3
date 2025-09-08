#pragma once

#include "Fitters/MCMCBase.h"

/// @brief MCMC algorithm implementing the **Metropolis–Rosenbluth–Rosenbluth–Teller–Teller (MR\f$^2\f$T2)** method.
/// @note MR\f$^2\f$T2 is also known as Metropolis-Hastings
/// @author Asher Kaboth
class MR2T2 : public MCMCBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
     MR2T2(manager *const manager) : MCMCBase(manager) {}

     /// @brief Destructor
     virtual ~MR2T2() = default;
     /// @brief Get name of class
     inline std::string GetName() const override { return "MR2T2"; };
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
