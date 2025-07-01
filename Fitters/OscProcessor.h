#pragma once

// MaCh3 includes
#include "Fitters/MCMCProcessor.h"
#include "Samples/HistogramUtils.h"

/// @author Clarence Wret
/// @author Kamil Skwarczynski
/// @brief This class extends MCMC and allow specialised for Oscillation parameters analysis which require specialised hardcoding
class OscProcessor : public MCMCProcessor {
  public:
    /// @brief Constructs an OscProcessor object with the specified input file and options.
    /// @param InputFile The path to the input file containing MCMC data.
    OscProcessor(const std::string &InputFile);
    /// @brief Destroys the OscProcessor object.
    virtual ~OscProcessor();

    /// @brief Perform Several Jarlskog Plotting
    /// @author Kevin Wood
    /// @note based on makeJarlskog.C
    void PerformJarlskogAnalysis();

  protected:
    /// @brief Read the Osc cov file and get the input central values and errors
    /// Here we allow Jarlskog Shenanigans
    void LoadAdditionalInfo() override;

    /// @brief Perform Jarlskog Plotting
    /// @author Kevin Wood
    /// @note based on drawJarl_dcpPriorComparison.C
    void MakeJarlskogPlot(const std::unique_ptr<TH1D>& jarl,
                          const std::unique_ptr<TH1D>& jarl_flatsindcp,
                          const std::unique_ptr<TH1D>& jarl_NH,
                          const std::unique_ptr<TH1D>& jarl_NH_flatsindcp,
                          const std::unique_ptr<TH1D>& jarl_IH,
                          const std::unique_ptr<TH1D>& jarl_IH_flatsindcp);

    /// @brief Calculate Jarlskog Invariant using oscillation parameters
    /// @param s2th13  Value of \f$ \sin^2\theta_{13} \f$
    /// @param s2th23  Value of \f$ \sin^2\theta_{23} \f$
    /// @param s2th12  Value of \f$ \sin^2\theta_{12} \f$
    /// @param dcp     CP-violating phase \f$ \delta_{\text{CP}} \f$ (in radians)
    /// @return The value of the Jarlskog invariant \f$ J_{\text{CP}} \f$
    /// @cite Jarlskog:1985ht
    double CalcJarlskog(const double s2th13, const double s2th23, const double s2th12, const double dcp) const;
    /// @brief Draw Prior value
    double SamplePriorForParam(const int paramIndex,
                               const std::unique_ptr<TRandom3>& randGen,
                               const std::vector<double>& FlatBounds) const;

    /// Will plot Jarlskog Invariant using information in the chain
    bool PlotJarlskog;

    /// Will plot Jarlskog Invariant using information in the chain
    bool OscEnabled;

    /// Name of the parameter representing \f$\sin^2\theta_{13}\f$.
    std::string Sin2Theta13Name;
    /// Name of the parameter representing \f$\sin^2\theta_{12}\f$.
    std::string Sin2Theta12Name;
    /// Name of the parameter representing \f$\sin^2\theta_{23}\f$.
    std::string Sin2Theta23Name;
    /// Name of the parameter representing \f$\delta_{\mathrm{CP}}\f$ (the CP-violating phase).
    std::string DeltaCPName;
    /// Name of the parameter representing \f$\Delta m^2_{32}\f$ (mass-squared difference).
    std::string DeltaM2_23Name;

    /// Index of \f$\sin^2\theta_{13}\f$ in the parameter list.
    int Sin2Theta13Index;
    /// Index of \f$\sin^2\theta_{12}\f$ in the parameter list.
    int Sin2Theta12Index;
    /// Index of \f$\sin^2\theta_{23}\f$ in the parameter list.
    int Sin2Theta23Index;
    /// Index of \f$\delta_{\mathrm{CP}}\f$ in the parameter list.
    int DeltaCPIndex;
    /// Index of \f$\Delta m^2_{32}\f$ in the parameter list.
    int DeltaM2_23Index;
};
