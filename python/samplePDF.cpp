#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "samplePDF/samplePDFBase.h"
#include "samplePDF/samplePDFFDBase.h"

namespace py = pybind11;

// As SamplePDFBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySamplePDFBase : public samplePDFBase {
public:
    /* Inherit the constructors */
    using samplePDFBase::samplePDFBase;

    /* Trampoline (need one for each virtual function) */
    void reweight() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            reweight,      /* Name of function in C++ (must match Python name) */
                           /* Argument(s) */
        );
    }

    double GetLikelihood() override {
        PYBIND11_OVERRIDE_PURE(
            double,        /* Return type */
            samplePDFBase, /* Parent class */
            GetLikelihood, /* Name of function in C++ (must match Python name) */
                           /* Argument(s) */
        );
    }
    
    void fill1DHist() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            fill1DHist,    /* Name of function in C++ (must match Python name) */
                           /* Argument(s) */
        );
    }
    
    void fill2DHist() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            fill2DHist     /* Name of function in C++ (must match Python name) */
                           /* Argument(s) */
        );
    }
};

void initSamplePDF(py::module &m){

    auto m_sample_pdf = m.def_submodule("sample_pdf");
    m_sample_pdf.doc() = 
        "This module deals with sampling from the posterior density function of your particular experimental model at different points, given your data. \
        In order to do this, you will generally need to create a SamplePDF object derived from :py:class:`pyMaCh3.fitter.SamplePDFBase` for each sample of events for your experiment. \
        For some more details on this you can see [the wiki page](https://github.com/mach3-software/MaCh3/wiki/04.-Making-a-samplePDF-experiment-class) on this. The code examples there are written using c++ however the general ideas are the same. \
        Happy sampling!";

    
    py::class_<samplePDFBase, PySamplePDFBase /* <--- trampoline*/>(m_sample_pdf, "SamplePDFBase")
        .def(py::init())
        
        .def(
            "reweight", 
            &samplePDFBase::reweight, 
            "reweight the MC events in this sample. You will need to override this."
        )
        
        .def(
            "get_likelihood", 
            &samplePDFBase::GetLikelihood, 
            "Get the sample likelihood at the current point in your model space. You will need to override this."
        )
        
        .def(
            "fill_1d_hist", 
            &samplePDFBase::fill1DHist, 
            "Do the initial filling of the sample for 1d histogram. You will need to override this."
        )
        
        .def(
            "fill_2d_hist", 
            &samplePDFBase::fill2DHist, 
            "Do the initial filling of the sample for 2d histogram. You will need to override this."
        )
        
        ;
    
}