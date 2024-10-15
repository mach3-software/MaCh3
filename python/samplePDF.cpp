#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "samplePDF/samplePDFBase.h"
#include "samplePDF/samplePDFFDBase.h"
#include "samplePDF/FDMCStruct.h"

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


// As SamplePDFFDBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySamplePDFFDBase : public samplePDFFDBase {
public:
    /* Inherit the constructors */
    using samplePDFFDBase::samplePDFFDBase;

    /* Trampoline (need one for each virtual function) */
    void SetupWeightPointers() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            samplePDFFDBase, /* Parent class */
            "setup_weight_pointers", /*python name*/
            SetupWeightPointers,     /* Name of function in C++ */
                             /* Argument(s) */
        );
    }

    double ReturnKinematicParameter(std::string, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("sample"),
            py::arg("event")            /* Argument(s) */
        );
    }
    double ReturnKinematicParameter(double, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("sample"),
            py::arg("event")            /* Argument(s) */
        );
    }

    const double *ReturnKinematicParameterByReference(std::string, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            ReturnKinematicParameterByReference, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("sample"),
            py::arg("event")            /* Argument(s) */
        );
    }
    const double *ReturnKinematicParameterByReference(double, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            ReturnKinematicParameterByReference, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("sample"),
            py::arg("event")            /* Argument(s) */
        );
    }
    
    std::vector<double> ReturnKinematicParameterBinning(std::string) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::vector<double>,          /* Return type */
            samplePDFFDBase,              /* Parent class */
            "get_event_kinematic_binning",/* python name*/
            ReturnKinematicParameterBinning, /* Name of function in C++ (must match Python name) */
            py::arg("variable")           /* Argument(s) */
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
        In order to do this, you will generally need to create a SamplePDF object derived from :py:class:`pyMaCh3.fitter.SamplePDFFDBase` for each sample of events for your experiment. \
        For some more details on this you can see [the wiki page](https://github.com/mach3-software/MaCh3/wiki/04.-Making-a-samplePDF-experiment-class) on this. The code examples there are written using c++ however the general ideas are the same. \
        Happy sampling!";

    py::class_<samplePDFBase, PySamplePDFBase /* <--- trampoline*/>(m_sample_pdf, "_SamplePDFBase")
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
    ; // End of samplePDFBase binding

    py::class_<samplePDFFDBase, PySamplePDFFDBase /* <--- trampoline*/, samplePDFBase>(m_sample_pdf, "SamplePDFFDBase")
        .def(py::init())
        
        .def(
            "set_xsec_cov", 
            &samplePDFFDBase::SetXsecCov, 
            "Set the cross section covariance matrix object."
        )
        
        .def(
            "set_osc_cov", 
            &samplePDFFDBase::SetOscCov, 
            "Set the oscillation parameter covariance matrix object."
        )
        
    ;

    /* Not sure if this will be needed in future versions of MaCh3 so leaving commented for now
    py::class_<fdmc_base>(m_sample_pdf, "MCstruct")
        .def(py::init())
        
        // Because a lot of the variables in fdmc_base use c style arrays,
        // we need to provide some setter functions to be able to set them using more
        // "pythony" objects, e.g. lists and numpy arrays
        .def(
            "set_event_variable_values", 
            [](fdmc_base &self, int dim, py::array_t<double, py::array::c_style> &array)
            {
                py::buffer_info bufInfo = array.request();

                if ( dim > 2 )
                    throw MaCh3Exception(__FILE__, __LINE__, "Currently only dimensions of 1 or 2 are supported sorry :(");

                if ( bufInfo.ndim != 1 )
                    throw MaCh3Exception(__FILE__, __LINE__, "Number of dimensions in parameter array must be one if setting only one of the event variable arrays!");
                
                if( dim ==1 )
                    self.x_var = array.data();
                    
                else if ( dim == 2)
                    self.y_var = array.data();
            }
        )
    ;
    */
}