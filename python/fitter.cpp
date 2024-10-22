// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// MaCh3 includes
#include "mcmc/FitterBase.h"
#include "mcmc/mcmc.h"
#include "mcmc/MinuitFit.h"
#include "mcmc/PSO.h"

namespace py = pybind11;

// As FitterBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PyFitterBase : public FitterBase {
public:
    /* Inherit the constructors */
    using FitterBase::FitterBase;

    /* Trampoline (need one for each virtual function) */
    void runMCMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,        /* Return type */
            FitterBase,  /* Parent class */
            "run",       /* Python name*/
            runMCMC      /* Name of function in C++ (must match Python name) */
        );
    }

    std::string GetName() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::string, /* Return type */
            FitterBase,  /* Parent class */
            "get_name",  /* Python name*/
            GetName      /* Name of function in C++ (must match Python name) */
        );
    }
};

// As LikelihoodFit is an abstract base class we have to do some gymnastics to get it to get it into python
class PyLikelihoodFit : public LikelihoodFit {
public:
    /* Inherit the constructors */
    using LikelihoodFit::LikelihoodFit;

    /* Trampoline (need one for each virtual function) */
    void runMCMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,        /* Return type */
            LikelihoodFit,/* Parent class */
            "run",       /* Python name*/
            runMCMC      /* Name of function in C++ (must match Python name) */
        );
    }
};

void initFitter(py::module &m){

    auto m_fitter = m.def_submodule("fitter");
    m_fitter.doc() = 
        "This is a Python binding of MaCh3s C++ mcmc library.";
    
    
    py::class_<FitterBase, PyFitterBase /* <--- trampoline*/>(m_fitter, "FitterBase")
        .def(py::init<manager* const>())
        
        .def(
            "run", 
            &FitterBase::runMCMC, 
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
            "add_sample_PDF",
            &FitterBase::addSamplePDF,
            " This function adds a sample PDF object to the analysis framework. The sample PDF object will be utilized in fitting procedures or likelihood scans. \n"
            " :param sample: A sample PDF object derived from samplePDFBase. ",
            py::arg("sample")
        )
        
        .def(
            "add_syst_object",
            &FitterBase::addSystObj,
            " This function adds a Covariance object to the analysis framework. The Covariance object will be utilized in fitting procedures or likelihood scans. \n"
            " :param cov: A Covariance object derived from covarianceBase. ",
            py::arg("cov")
        )

        .def(
            "add_osc_handler",
            py::overload_cast<covarianceOsc *>(&FitterBase::addOscHandler),
            " Adds an oscillation handler for covariance objects. \n"
            " :param oscf: An oscillation handler object for dealing with neutrino oscillation calculations. ",
            py::arg("osc")
        )
        
    ; // End of FitterBase class binding

    py::class_<mcmc, FitterBase>(m_fitter, "MCMC")
        .def(py::init<manager* const>())
        
        .def(
            "set_chain_length", 
            &mcmc::setChainLength, 
            "Set how long chain should be.",
            py::arg("length")
        )

        .def(
            "set_init_step_num", 
            &mcmc::setInitialStepNumber, 
            "Set initial step number, used when starting from another chain.",
            py::arg("step_num")
        )
    ; // end of MCMC class binding

    py::class_<LikelihoodFit, PyLikelihoodFit /* <--- trampoline*/, FitterBase>(m_fitter, "LikelihoodFit")
        .def(py::init<manager* const>())
        
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

    py::class_<MinuitFit, LikelihoodFit>(m_fitter, "MinuitFit")
        .def(py::init<manager* const>())
        
    ; // end of MinuitFit class binding

    py::class_<PSO, LikelihoodFit>(m_fitter, "PSO")
        .def(py::init<manager* const>())

        .def(
            "init",
            &PSO::init,
            "Initialise the fitter"
        )
        
    ; // end of PSO class binding

}