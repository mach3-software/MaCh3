// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// MaCh3 includes
#include "Fitters/FitterBase.h"
#include "Fitters/MR2T2.h"
#include "Fitters/DelayedMR2T2.h"
#include "Fitters/MinuitFit.h"
#include "Fitters/PSO.h"

namespace py = pybind11;

/// @brief EW: As FitterBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PyFitterBase : public FitterBase {
public:
    /* Inherit the constructors */
    using FitterBase::FitterBase;

    /* Trampoline (need one for each virtual function) */
    void RunMCMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,        /* Return type */
            FitterBase,  /* Parent class */
            "run",       /* Python name*/
            RunMCMC      /* Name of function in C++ (must match Python name) */
        );
    }
};

/// @brief EW: As LikelihoodFit is an abstract base class we have to do some gymnastics to get it to get it into python
class PyLikelihoodFit : public LikelihoodFit {
public:
    /* Inherit the constructors */
    using LikelihoodFit::LikelihoodFit;

    /* Trampoline (need one for each virtual function) */
    void RunMCMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,        /* Return type */
            LikelihoodFit,/* Parent class */
            "run",       /* Python name*/
            RunMCMC      /* Name of function in C++ (must match Python name) */
        );
    }
};

void initFitters(py::module &m) {
    auto m_fitters = m.def_submodule("fitters");
    m_fitters.doc() =
        "This is a Python binding of MaCh3s C++ fitters library.";
    
    
    py::class_<FitterBase, PyFitterBase /* <--- trampoline*/>(m_fitters, "FitterBase")
        .def(py::init<Manager* const>())
        
        .def(
            "run", 
            &FitterBase::RunMCMC,
            "The implementation of the fitter, you should override this with your own desired fitting algorithm"
        )

        .def(
            "get_name", 
            &FitterBase::GetName, 
            " The name of the algorithm, you should override this with something like:: \n"
            "\n"
            "    return 'mySuperCoolAlgoName' \n"
        )

        .def(
            "run_LLH_scan",
            &FitterBase::RunLLHScan,
            "Perform a 1D likelihood scan"
        )

        .def(
            "get_step_scale_from_LLH_scan",
            &FitterBase::GetStepScaleBasedOnLLHScan,
            "LLH scan is good first estimate of step scale, this will get the rough estimates for the step scales based on running an LLH scan"
        )

        .def(
            "run_2d_LLH_scan",
            &FitterBase::Run2DLLHScan,
            " Perform a 2D likelihood scan. \n"
            " :param warning: This operation may take a significant amount of time, especially for complex models."
        )
 
        .def(
            "run_sigma_var",
            &FitterBase::RunSigmaVar,
            " Perform a 2D and 1D sigma var for all samples. \n"
            " :param warning: Code uses TH2Poly"
        )

        .def(
            "drag_race",
            &FitterBase::DragRace,
            " Calculates the required time for each sample or covariance object in a drag race simulation. Inspired by Dan's feature \n"
            " :param NLaps: number of laps, every part of Fitter will be tested with given number of laps and you will get total and average time",
            py::arg("NLaps") = 100
        )

        // stuff for registering other objects with the fitter
        
        .def(
            "add_sample_handler",
            &FitterBase::AddSampleHandler,
            " This function adds a sample handler object to the analysis framework. The sample handler object will be utilized in fitting procedures or likelihood scans. \n"
            " :param sample: A sample handler object derived from SampleHandlerBase. ",
            py::arg("sample")
        )
        
        .def(
            "add_syst_object",
            &FitterBase::AddSystObj,
            " This function adds a Covariance object to the analysis framework. The Covariance object will be utilized in fitting procedures or likelihood scans. \n"
            " :param cov: A Covariance object derived from ParameterHandlerBase. ",
            py::arg("cov")
        )

    ; // End of FitterBase class binding

    py::class_<MR2T2, FitterBase>(m_fitters, "mcmc")
        .def(py::init<Manager *const>())

        .def(
            "set_chain_length",
            &MR2T2::setChainLength,
            "Set how long chain should be.",
            py::arg("length")); // end of MCMC class binding

    py::class_<DelayedMR2T2, FitterBase>(m_fitters, "DelayedMCMC")
        .def(py::init<Manager *const>())

        .def(
            "set_chain_length",
            &MR2T2::setChainLength,
            "Set how long chain should be.",
            py::arg("length")); // end of MCMC class binding

    py::class_<LikelihoodFit, PyLikelihoodFit /* <--- trampoline*/, FitterBase>(m_fitters, "LikelihoodFit")
        .def(py::init<Manager* const>())
        
        .def(
            "caluclate_chi2",
            [](LikelihoodFit &self, const std::vector<double> &parameterVals)
            {
                return self.CalcChi2(parameterVals.data());
            },
            "Get the Chi2 calculation over all included samples and syst objects for the specified parameter_values \n\
            :param parameter_valuse: The location to evaluate the chi2 at.",
            py::arg("parameter_values")
        )

        .def(
            "get_n_params",
            &LikelihoodFit::GetNPars,
            "Get The total number of parameters across all known covariance objects associated with this LikelihoodFit object."
        )
    ; // end of LikelihoodFit class binding

    py::class_<MinuitFit, LikelihoodFit>(m_fitters, "MinuitFit")
        .def(py::init<Manager* const>())
        
    ; // end of MinuitFit class binding

    py::class_<PSO, LikelihoodFit>(m_fitters, "PSO")
        .def(py::init<Manager* const>())

        .def(
            "init",
            &PSO::init,
            "Initialise the fitter"
        )

    ; // end of PSO class binding
}
