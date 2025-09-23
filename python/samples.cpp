#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "Samples/SampleHandlerFD.h"

namespace py = pybind11;

/// @brief EW: As SampleHandlerBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerBase : public SampleHandlerBase {
public:
    /* Inherit the constructors */
    using SampleHandlerBase::SampleHandlerBase;

    /* Trampoline (need one for each virtual function) */
    std::string GetSampleHandlerName() const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetSampleHandlerName, /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::string GetSampleTitle(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetSampleTitle,       /* Name of function in C++ (must match Python name) */
            iSample               /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    int GetNOscChannels(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            int,                  /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetNOscChannels,      /* Name of function in C++ (must match Python name) */
            iSample               /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Reweight() override {
        PYBIND11_OVERRIDE_PURE(
            void,              /* Return type */
            SampleHandlerBase, /* Parent class */
            Reweight           /* Name of function in C++ (must match Python name) */
        );
    }


    /* Trampoline (need one for each virtual function) */
    double GetSampleLikelihood(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            double,                /* Return type */
            SampleHandlerBase,     /* Parent class */
            GetSampleLikelihood,   /* Name of function in C++ (must match Python name) */
            iSample                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE(
            void,                /* Return type */
            SampleHandlerBase,   /* Parent class */
            CleanMemoryBeforeFit /* Name of function in C++ (must match Python name) */
        );
    }

    double GetLikelihood() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,            /* Return type */
            SampleHandlerBase, /* Parent class */
            "get_likelihood",  /* Python name*/
            GetLikelihood      /* Name of function in C++ (must match Python name) */
                               /* Argument(s) */
        );
    }
};


/// @brief As SampleHandlerFD is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerFD : public SampleHandlerFD {
public:
    /* Inherit the constructors */
    using SampleHandlerFD::SampleHandlerFD;

    /* Trampoline (need one for each virtual function) */
    void SetupWeightPointers() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_weight_pointers", /*python name*/
            SetupWeightPointers,     /* Name of function in C++ */
                             /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            SampleHandlerFD, /* Parent class */
            CleanMemoryBeforeFit       /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupSplines() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_splines", /*python name*/
            SetupSplines,     /* Name of function in C++ */
                             /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Init() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "init",          /*python name*/
            Init,            /* Name of function in C++ */
                             /* Argument(s) */
        );
    }
    
    /* Trampoline (need one for each virtual function) */
    int SetupExperimentMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,             /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_experiment_MC", /*python name*/
            SetupExperimentMC,     /* Name of function in C++ */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupFDMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_FD_MC",   /*python name*/
            SetupFDMC,       /* Name of function in C++ */
        );
    }

    int ReturnKinematicParameterFromString(std::string) {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_kinematic_by_name",    /* python name*/
            ReturnKinematicParameterFromString, /* Name of function in C++ (must match Python name) */
            py::arg("variable_name")
        );
    }
    
    std::string ReturnStringFromKinematicParameter(int) {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::string,                /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_kinematic_name",       /* python name*/
            ReturnStringFromKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable_id")
        );
    }

    double ReturnKinematicParameter(std::string, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }

    double ReturnKinematicParameter(int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }

    const double *GetPointerToKinematicParameter(std::string, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }
    const double *GetPointerToKinematicParameter(double, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }
    
    void RegisterFunctionalParameters() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,
            SampleHandlerFD,
            "register_functional_parameters",
            RegisterFunctionalParameters
        );
    }
};

void initSamples(py::module &m){

    auto m_samples = m.def_submodule("samples");
    m_samples.doc() =
        "This is a Python binding of MaCh3s C++ based samples library.";

    // Bind the systematic type enum that lets us set different types of systematics
    py::enum_<TestStatistic>(m_samples, "TestStatistic")
        .value("Poisson", TestStatistic::kPoisson)
        .value("Barlow_Beeston", TestStatistic::kBarlowBeeston)
        .value("Ice_Cube", TestStatistic::kIceCube)
        .value("Pearson", TestStatistic::kPearson)
        .value("Dembinski_Abdelmottele", TestStatistic::kDembinskiAbdelmotteleb)
        .value("N_Test_Statistics", TestStatistic::kNTestStatistics);

    py::class_<SampleHandlerBase, PySampleHandlerBase /* <--- trampoline*/>(m_samples, "SampleHandlerBase")
        .def(py::init())
        
        .def(
            "reweight", 
            &SampleHandlerBase::Reweight,
            "reweight the MC events in this sample. You will need to override this."
        )
        
        .def(
            "get_likelihood", 
            &SampleHandlerBase::GetLikelihood,
            "Get the sample likelihood at the current point in your model space. You will need to override this."
        )
        
        .def(
            "set_test_stat",
            &SampleHandlerBase::SetTestStatistic,
            "Set the test statistic that should be used when calculating likelihoods. \n\
            :param test_stat: The new test statistic to use",
            py::arg("test_stat")
        )

        .def(
            "get_bin_LLH",
            py::overload_cast<double, double, double>(&SampleHandlerBase::GetTestStatLLH, py::const_),
            "Get the LLH for a bin by comparing the data and MC. The result depends on having previously set the test statistic using :py:meth:`pyMaCh3.samples.SampleHandlerBase.set_test_stat` \n\
            :param data: The data content of the bin. \n\
            :param mc: The mc content of the bin \n\
            :param w2: The Sum(w_{i}^2) (sum of weights squared) in the bin, which is sigma^2_{MC stats}",
            py::arg("data"), 
            py::arg("mc"), 
            py::arg("w2")
        )
    ; // End of SampleHandlerBase binding

    py::class_<SampleHandlerFD, PySampleHandlerFD /* <--- trampoline*/, SampleHandlerBase>(m_samples, "SampleHandlerFD")
        .def(
            py::init<std::string, ParameterHandlerGeneric*>(),
            "This should never be called directly as SampleHandlerFD is an abstract base class. \n\
            However when creating a derived class, in the __init__() method, you should call the parent constructor i.e. this one by doing:: \n\
            \n\
            \tsuper(<your derived SampleHandler class>, self).__init__(*args) \n\
            \n ",
            py::arg("mc_version"),
            py::arg("xsec_cov")
        )
    ;

    /* Not sure if this will be needed in future versions of MaCh3 so leaving commented for now
    py::class_<fdmc_base>(m_samples, "MCstruct")
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
