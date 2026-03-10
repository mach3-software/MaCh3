#pragma once

/// @file pyMaCh3.h
/// @author Ewan Miller

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "python/plotting.h"
#include "python/fitters.h"
#include "python/samples.h"
#include "python/manager.h"
#include "python/parameters.h"
#include "python/splines.h"

class MaCh3PyBinder {

  public:

    virtual ~MaCh3PyBinder() {};

    void initialise(py::module &m) {
        
        initPlotting(m);
        initFitters(m);
        initSamples(m);
        initManager(m);
        initParameters(m);
        initSplines(m);

    }

  private:
    // these can be overwritten by experiments to set up their own pyMaCh3
    virtual void initParametersExperiment(py::module &m_parameters) { (void)m_parameters; }
    virtual void initPlottingExperiment(py::module &m_plotting)     { (void)m_plotting;   }
    virtual void initFittersExperiment(py::module &m_fitters)       { (void)m_fitters;    }
    virtual void initSamplesExperiment(py::module &m_samples)       { (void)m_samples;    }
    virtual void initManagerExperiment(py::module &m_manager)       { (void)m_manager;    }
    virtual void initSplinesExperiment(py::module &m_splines)       { (void)m_splines;    }
    virtual void initModulesExperiment(py::module &m)               { (void)m;            }

    // These functions will be called to set up each of the python modules
    // they call the core init<module name> functions and then the experiment
    // specific ones
    void initPlotting(py::module &m) {

        auto m_plotting = m.def_submodule("plotting");
        m_plotting.doc() = "This is a Python binding of MaCh3s C++ based plotting library.";

        initPlottingModule(m_plotting); // <- defined in python/plotting.cpp
        initPlottingExperiment(m_plotting);
    }
    void initFitters(py::module &m) {

        auto m_fitters = m.def_submodule("fitters");
        m_fitters.doc() =
            "This is a Python binding of MaCh3s C++ fitters library.";

        initFittersModule(m_fitters); // <- defined in python/fitters.cpp
        initFittersExperiment(m_fitters);
    }
    void initSamples(py::module &m) {

        auto m_samples = m.def_submodule("samples");

        m_samples.doc() =
            "This is a Python binding of MaCh3s C++ based samples library.";
            
        initSamplesModule(m_samples); // <- defined in python/samples.cpp
        initSamplesExperiment(m_samples);
    }
    void initManager(py::module &m) {

        auto m_manager = m.def_submodule("manager");
        m_manager.doc() = 
            "This is a Python binding of MaCh3s C++ based manager library.";
            
        initManagerModule(m_manager); // <- defined in python/manager.cpp
        initManagerExperiment(m_manager);
    }
    void initParameters(py::module &m) {

        auto m_parameters = m.def_submodule("parameters");
        m_parameters.doc() =
            "This is a Python binding of MaCh3s C++ parameters library.";

        initParametersModule(m_parameters); // <- defined in python/parameters.cpp
        initParametersExperiment(m_parameters);
    }
    void initSplines(py::module &m) {

        auto m_splines = m.def_submodule("splines");
        m_splines.doc() = 
            "This is a Python binding of MaCh3s C++ based spline library.";

        initSplinesModule(m_splines); // <- defined in python/splines.cpp
        initSplinesExperiment(m_splines);
    }

};

#define MAKE_PYMACH3_MDULE( PYBINDER_CLASS ) \
PYBIND11_MODULE( _pyMaCh3, m ) {             \
    PYBINDER_CLASS().initialise(m);          \
} 
