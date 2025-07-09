#pragma once

#include "Fitters/MCMCBase.h"

class MR2T2 : public MCMCBase {
 public:
    /// @brief Constructor
    /// @param fitMan A pointer to a manager object, which will handle all settings.
    MR2T2(manager * const fitMan) : MCMCBase(fitMan) {}

    /// @brief Destructor
    virtual ~MR2T2() = default;

 protected:
    void ProposeStep() override;
    double CheckStep() override;
};