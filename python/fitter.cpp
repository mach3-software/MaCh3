#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mcmc/FitterBase.h"

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
        "This module contains the various MaCh3 fitter algorithms which are available, \
        as well as the :py:class:`pyMaCh3.fitter.FitterBase` class which you can use to \
        implement your own!";
    
    
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
            "''' \n"
            "return 'mySuperCoolAlgoName' \n"
            "''' \n"
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
            " *warning* This operation may take a significant amount of time, especially for complex models."
        )
 
        .def(
            "run_sigma_var",
            &FitterBase::RunSigmaVar,
            " Perform a 2D and 1D sigma var for all samples. \n"
            " *warning* Code uses TH2Poly"
        )

        .def(
            "drag_race",
            &FitterBase::DragRace,
            " Calculates the required time for each sample or covariance object in a drag race simulation. Inspired by Dan's feature \n"
            " *NLaps* number of laps, every part of Fitter will be tested with given number of laps and you will get total and average time",
            py::arg("NLaps") = 100
        )

        // stuff for registering other objects with the fitter
        
        .def(
            "add_sample_PDF",
            &FitterBase::addSamplePDF,
            " This function adds a sample PDF object to the analysis framework. The sample PDF object will be utilized in fitting procedures or likelihood scans. \n"
            " *sample* A sample PDF object derived from samplePDFBase. "
        )

        ;

    

}