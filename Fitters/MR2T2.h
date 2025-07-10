#pragma once

#include "Fitters/MCMCBase.h"

class MR2T2 : public MCMCBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
     MR2T2(manager *const manager) : MCMCBase(manager) {}

     /// @brief Destructor
     virtual ~MR2T2() = default;

     inline std::string GetName() const { return "MR2T2"; };

 protected:
    void DoStep() override;
    void ProposeStep() override;
    double AcceptanceProbability() override;
};