#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "samplePDF/samplePDFFDBase.h"

namespace py = pybind11;

// As SamplePDFBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySamplePDFBase : public samplePDFBase {
public:
    /* Inherit the constructors */
    using samplePDFBase::samplePDFBase;

    /* Trampoline (need one for each virtual function) */
    std::string GetSampleName(int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            samplePDFBase, /* Parent class */
            GetSampleName, /* Name of function in C++ (must match Python name) */
            iSample         /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void reweight() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            reweight       /* Name of function in C++ (must match Python name) */
        );
    }

    double GetLikelihood() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,        /* Return type */
            samplePDFBase, /* Parent class */
            "get_likelihood", /* Python name*/
            GetLikelihood  /* Name of function in C++ (must match Python name) */
                           /* Argument(s) */
        );
    }
    
    void fill1DHist() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            "fill_1d_hist", /* Python Name */
            fill1DHist     /* Name of function in C++ (must match Python name) */
        );
    }
    
    void fill2DHist() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,          /* Return type */
            samplePDFBase, /* Parent class */
            "fill_2d_hist",/* Python name*/
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

    /* Trampoline (need one for each virtual function) */
    void SetupSplines() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            samplePDFFDBase, /* Parent class */
            "setup_splines", /*python name*/
            SetupSplines,     /* Name of function in C++ */
                             /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Init() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            samplePDFFDBase, /* Parent class */
            "init",          /*python name*/
            Init,            /* Name of function in C++ */
                             /* Argument(s) */
        );
    }
    
    /* Trampoline (need one for each virtual function) */
    int setupExperimentMC(int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,             /* Return type */
            samplePDFFDBase, /* Parent class */
            "setup_experiment_MC", /*python name*/
            setupExperimentMC,     /* Name of function in C++ */
            py::arg("sample_id")   /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void setupFDMC(int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            samplePDFFDBase, /* Parent class */
            "setup_FD_MC",   /*python name*/
            setupFDMC,       /* Name of function in C++ */
            py::arg("sample_id") /* Argument(s) */
        );
    }

    int ReturnKinematicParameterFromString(std::string) {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,                     /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_kinematic_by_name",    /* python name*/
            ReturnKinematicParameterFromString, /* Name of function in C++ (must match Python name) */
            py::arg("variable_name")
        );
    }
    
    std::string ReturnStringFromKinematicParameter(int) {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::string,                /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_kinematic_name",       /* python name*/
            ReturnStringFromKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable_id")
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

    const double *GetPointerToKinematicParameter(std::string, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("sample"),
            py::arg("event")            /* Argument(s) */
        );
    }
    const double *GetPointerToKinematicParameter(double, int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            samplePDFFDBase,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
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
        "This is a Python binding of MaCh3s C++ based samplePDF library.";

    // Bind the systematic type enum that lets us set different types of systematics
    py::enum_<TestStatistic>(m_sample_pdf, "TestStatistic")
        .value("Poisson", TestStatistic::kPoisson)
        .value("Barlow_Beeston", TestStatistic::kBarlowBeeston)
        .value("Ice_Cube", TestStatistic::kIceCube)
        .value("Pearson", TestStatistic::kPearson)
        .value("Dembinski_Abdelmottele", TestStatistic::kDembinskiAbdelmotteleb)
        .value("N_Test_Statistics", TestStatistic::kNTestStatistics);

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

        .def(
            "set_test_stat",
            &samplePDFBase::SetTestStatistic,
            "Set the test statistic that should be used when calculating likelihoods. \n\
            :param test_stat: The new test statistic to use",
            py::arg("test_stat")
        )

        .def(
            "get_bin_LLH",
            py::overload_cast<double, double, double>(&samplePDFBase::getTestStatLLH, py::const_),
            "Get the LLH for a bin by comparing the data and MC. The result depends on having previously set the test statistic using :py:meth:`pyMaCh3.sample_pdf.SamplePDFFDBase.set_test_stat` \n\
            :param data: The data content of the bin. \n\
            :param mc: The mc content of the bin \n\
            :param w2: The Sum(w_{i}^2) (sum of weights squared) in the bin, which is sigma^2_{MC stats}",
            py::arg("data"), 
            py::arg("mc"), 
            py::arg("w2")
        )
    ; // End of samplePDFBase binding

    py::class_<samplePDFFDBase, PySamplePDFFDBase /* <--- trampoline*/, samplePDFBase>(m_sample_pdf, "SamplePDFFDBase")
        .def(
            py::init<std::string, covarianceXsec*>(),
            "This should never be called directly as samplePDFFDBase is an abstract base class. \n\
            However when creating a derived class, in the __init__() method, you should call the parent constructor i.e. this one by doing:: \n\
            \n\
            \tsuper(<your derived samplePDF class>, self).__init__(*args) \n\
            \n ",
            py::arg("mc_version"),
            py::arg("xsec_cov")
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
