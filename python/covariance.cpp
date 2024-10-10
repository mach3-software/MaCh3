// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// MaCh3 includes
#include "covariance/covarianceBase.h"


namespace py = pybind11;

// As FitterBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PyCovarianceBase : public covarianceBase {
public:
    /* Inherit the constructors */
    using covarianceBase::covarianceBase;

    /* Trampoline (need one for each virtual function) */
    int CheckBounds() override {
        PYBIND11_OVERRIDE_NAME(
            int,            /* Return type */
            covarianceBase, /* Parent class */
            "check_bounds", /* Name in python*/
            CheckBounds     /* Name of function in C++ (must match Python name) */
        );
    }

    double GetLikelihood() override {
        PYBIND11_OVERRIDE_NAME(
            double,           /* Return type */
            covarianceBase,   /* Parent class */
            "get_likelihood", /* Name in python*/
            GetLikelihood     /* Name of function in C++ (must match Python name) */
        );
    }
    
    double getNominal(int i) override {
        PYBIND11_OVERRIDE_NAME(
            double,          /* Return type */
            covarianceBase,  /* Parent class */
            "get_nominal",   /* Name in python*/
            getNominal,      /* Name of function in C++ (must match Python name) */
            i                /* Arguments*/
        );
    }
    
    void proposeStep() override {
        PYBIND11_OVERRIDE_NAME(
            void,            /* Return type */
            covarianceBase,  /* Parent class */
            "propose_step",  /* Name in python*/
            proposeStep      /* Name of function in C++ (must match Python name) */
        );
    }

    std::vector<double> getNominalArray() override {
        PYBIND11_OVERRIDE_NAME(
            std::vector<double>, /* Return type */
            covarianceBase,      /* Parent class */
            "get_nominal_array", /* Name in python*/
            getNominalArray      /* Name of function in C++ (must match Python name) */
        );
    }
};


void initCovariance(py::module &m){

    auto m_covariance = m.def_submodule("covariance");
    m_covariance.doc() = 
        "This is a Python binding of MaCh3s C++ covariance library.";

        
    py::class_<covarianceBase, PyCovarianceBase /* <--- trampoline*/>(m_covariance, "CovarianceBase")
        .def(
            py::init<const std::vector<std::string>&, const char *, double, int, int>(),
            "Construct a covariance object from a yaml file \n\
            :param yaml_file: The name of the yaml file to initialise from. \n\
            :param name: the name of this covariance object. \n\
            :param threshold: threshold PCA threshold from 0 to 1. Default is -1 and means no PCA. \n\
            :param first_PCA_par: FirstPCAdpar First PCA parameter that will be decomposed. \n\
            :param last_PCA_par: LastPCAdpar First PCA parameter that will be decomposed.",
            py::arg("yaml_file"),
            py::arg("name"),
            py::arg("threshold") = -1.0,
            py::arg("firs_PCA_par") = -999,
            py::arg("last_PCA_par") = -999
        )
        
        .def(
            "calculate_likelihood", 
            &covarianceBase::CalcLikelihood, 
            "Calculate penalty term based on inverted covariance matrix."
        )

        .def(
            "throw_par_prop",
            &covarianceBase::throwParProp,
            "Throw the proposed parameter by magnitude *mag* X sigma. \n\
            :param mag: This value multiplied by the prior value of each parameter will be the width of the distribution that the parameter values are drawn from. ",
            py::arg("mag") = 1.0
        )

    ;
    

}
