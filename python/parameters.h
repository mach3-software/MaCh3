#pragma once

/// @file parameters.h
/// @author Ewan Miller

// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
// MaCh3 includes
#include "Parameters/ParameterHandlerBase.h"
#include "Parameters/ParameterHandlerGeneric.h"

namespace py = pybind11;

/// @brief EW: As ParameterHandlerBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PyParameterHandlerBase : public ParameterHandlerBase {
public:
    /* Inherit the constructors */
    using ParameterHandlerBase::ParameterHandlerBase;

    /* Trampoline (need one for each virtual function) */
    double GetLikelihood() override {
        PYBIND11_OVERRIDE_NAME(
            double,           /* Return type */
            ParameterHandlerBase,   /* Parent class */
            "get_likelihood", /* Name in python*/
            GetLikelihood     /* Name of function in C++ (must match Python name) */
        );
    }

    void ProposeStep() override {
        PYBIND11_OVERRIDE_NAME(
            void,            /* Return type */
            ParameterHandlerBase,  /* Parent class */
            "propose_step",  /* Name in python*/
            ProposeStep      /* Name of function in C++ (must match Python name) */
        );
    }
};


void initParametersModule(py::module &m){

    auto m_parameters = m.def_submodule("parameters");
    m_parameters.doc() =
        "This is a Python binding of MaCh3s C++ parameters library.";

    
    // Bind the systematic type enum that lets us set different types of systematics
    py::enum_<SystType>(m_parameters, "SystematicType")
            .value("Normalisation", SystType::kNorm)
            .value("Spline", SystType::kSpline)
            .value("Functional", SystType::kFunc)
            .value("N_Systematic_Types", SystType::kSystTypes);

        
    py::class_<ParameterHandlerBase, PyParameterHandlerBase /* <--- trampoline*/>(m_parameters, "ParameterHandlerBase")
        .def(
            py::init<const std::vector<std::string>&, const char *, double, int, int>(),
            "Construct a parameters object from a set of yaml files that define the systematic parameters \n\
            :param yaml_files: The name of the yaml file to initialise from. \n\
            :param name: the name of this ParameterHandler object. \n\
            :param threshold: threshold PCA threshold from 0 to 1. Default is -1 and means no PCA. \n\
            :param first_PCA_par: FirstPCAdpar First PCA parameter that will be decomposed. \n\
            :param last_PCA_par: LastPCAdpar First PCA parameter that will be decomposed.",
            py::arg("yaml_files"),
            py::arg("name"),
            py::arg("threshold") = -1.0,
            py::arg("firs_PCA_par") = -999,
            py::arg("last_PCA_par") = -999
        )
        
        .def(
            "calculate_likelihood", 
            &ParameterHandlerBase::CalcLikelihood,
            "Calculate penalty term based on inverted covariance matrix."
        )

        .def(
            "get_internal_par_name",
            [](ParameterHandlerBase &self, int index)
            {
                // do this to disambiguate between the std::string and const char* version of this fn
                std::string ret;
                ret = self.GetParName(index);
                return ret;
            },
            "Get the internally used name of this parameter. \n\
            :param index: The global index of the parameter",
            py::arg("index")
        )
        
        .def(
            "get_fancy_par_name",
            [](ParameterHandlerBase &self, int index)
            {
                // do this to disambiguate between the std::string and const char* version of this fn
                std::string ret;
                ret = self.GetParFancyName(index);
                return ret;
            },
            "Get the name of this parameter. \n\
            :param index: The global index of the parameter",
            py::arg("index")
        )

        .def(
            "get_n_pars",
            &ParameterHandlerBase::GetNParameters,
            "Get the number of parameters that this ParameterHandler object knows about."
        )
        
        .def(
            "propose_step",
            &ParameterHandlerBase::ProposeStep,
            "Propose a step based on the covariances. Also feel free to overwrite if you want something more funky."
        )

        .def(
            "get_proposal_array",
            [](ParameterHandlerBase &self)
            {
                // Get the number of parameters
                size_t n_pars = self.GetNParameters();
                
                // Get pointer to the data
                const double* data_ptr = self.GetParPropVec().data();
                
                // Create a numpy array that copies the data
                // This ensures the numpy array owns its data and won't have lifetime issues
                py::array_t<double> result(n_pars);
                auto buf = result.request();
                double* result_ptr = static_cast<double*>(buf.ptr);
                
                // Copy the data
                std::memcpy(result_ptr, data_ptr, n_pars * sizeof(double));
                
                return result;
            },
            "Get the parameter proposal values as a numpy array. \n\
            This returns a copy of the current proposal values. \n\
            :return: A numpy array containing the proposal values for all parameters."
        )

        .def("set_parameters", 
             [](ParameterHandlerBase& self, py::object pars_obj = py::none()) {
                 if (pars_obj.is_none()) {
                     self.SetParameters();
                 } else {
                     // This handles both numpy arrays and Python lists
                     std::vector<double> pars_vec = pars_obj.cast<std::vector<double>>();
                     self.SetParameters(pars_vec);
                 }
             },
             py::arg("pars") = py::none(),
             R"pbdoc(
                 Set parameter values using array.
                 
                 Parameters
                 ----------
                 pars : numpy.ndarray or list of float, optional
                     Array holding new values for every parameter.
                     Must have same size as the number of parameters in the covariance class.
                     If not provided, parameters are set to their pre-fit values.
                     
                 Examples
                 --------
                 >>> import numpy as np
                 >>> handler.set_parameters(np.array([1.0, 2.0, 3.0]))
                 >>> handler.set_parameters([1.0, 2.0, 3.0])
                 >>> handler.set_parameters()
             )pbdoc")

    ; // End of ParameterHandlerBase binding

    
    py::class_<ParameterHandlerGeneric, ParameterHandlerBase /* <--- trampoline*/>(m_parameters, "ParameterHandlerGeneric")
        .def(
            py::init<const std::vector<std::string>&, const char *, double, int, int>(),
            "Construct a systematic ParameterHandler object from a set of yaml files that define the systematic parameters \n\
            :param yaml_files: The name of the yaml file to initialise from. \n\
            :param name: the name of this ParameterHandler object. \n\
            :param threshold: threshold PCA threshold from 0 to 1. Default is -1 and means no PCA. \n\
            :param first_PCA_par: FirstPCAdpar First PCA parameter that will be decomposed. \n\
            :param last_PCA_par: LastPCAdpar First PCA parameter that will be decomposed.",
            py::arg("yaml_files"),
            py::arg("name") = "xsec_cov",
            py::arg("threshold") = -1.0,
            py::arg("firs_PCA_par") = -999,
            py::arg("last_PCA_par") = -999
        )
        
        .def(
            "get_par_type",
            &ParameterHandlerGeneric::GetParamType,
            "Get what type of systematic this parameters is (see :py:class:`pyMaCh3.m_parameters.SystematicType` for possible types). \n\
            :param index: The global index of the parameter",
            py::arg("index")
        )
        
        .def(
            "get_par_spline_type",
            &ParameterHandlerGeneric::GetParSplineInterpolation,
            "Get what type of spline this parameter is set to use (assuming that it is a spline type parameter). \n\
            :param index: The index of the spline parameter",
            py::arg("index")
        )
        
        .def(
            "get_par_spline_name",
            &ParameterHandlerGeneric::GetParSplineName,
            "Get the name of the spline associated with a spline parameter. This is generally what it is called in input spline files and can in principle be different to the parameters name. \n\
            :param index: The index of the spline parameter",
            py::arg("index")
        )

    ; // End of ParameterHandlerGeneric binding
}
