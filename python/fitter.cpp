// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// MaCh3 includes
#include "mcmc/FitterBase.h"
#include "mcmc/mcmc.h"
#include "mcmc/MinuitFit.h"

namespace py = pybind11;

// As FitterBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PyFitterBase : public FitterBase {
public:
    /* Inherit the constructors */
    using FitterBase::FitterBase;

    /* Trampoline (need one for each virtual function) */
    void runMCMC() override {
        PYBIND11_OVERRIDE_PURE(
            void,        /* Return type */
            FitterBase,  /* Parent class */
            runMCMC,     /* Name of function in C++ (must match Python name) */
                         /* Argument(s) */
        );
    }

    std::string GetName() const override {
        PYBIND11_OVERRIDE_PURE(
            std::string, /* Return type */
            FitterBase,  /* Parent class */
            GetName,     /* Name of function in C++ (must match Python name) */
                         /* Argument(s) */
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
            "runMCMC", 
            &FitterBase::runMCMC, 
            "The implementation of the fitter, you should override this with your own desired fitting algorithm"
        )

        .def(
            "GetName", 
            &FitterBase::GetName, 
            " The name of the algorithm, you should override this with something like \n"
            "``` \n"
            "return 'mySuperCoolAlgoName' \n"
            "``` \n"
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
            " :param cov: A pointer to a Covariance object derived from covarianceBase. \n",
            py::arg("cov")
        )

        .def(
            "add_osc_handler",
            py::overload_cast<covarianceOsc *>(&FitterBase::addOscHandler),
            "  Adds an oscillation handler for covariance objects. \n"
            " :param oscf: A pointer to a covarianceOsc object for forward oscillations. \n",
            py::arg("oscf")
        )
        
        .def(
            "add_osc_handler",
            py::overload_cast<covarianceOsc *, covarianceOsc *>(&FitterBase::addOscHandler),
            "  Adds an oscillation handler for covariance objects. \n"
            " :param osca: A pointer to a covarianceOsc object for the first oscillation. \n"
            " :param oscb: A pointer to a covarianceOsc object for the second oscillation. \n",
            py::arg("osca"),
            py::arg("oscb")
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

    
    py::class_<MinuitFit, FitterBase>(m_fitter, "MinuitFit")
        .def(py::init<manager* const>())
        
    ; // end of MinuitFit class binding

}