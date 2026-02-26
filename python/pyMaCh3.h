#pragma once

/// @file pyMaCh3.h
/// @author Ewan Miller

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "plotting.h"
#include "fitters.h"
#include "samples.h"
#include "manager.h"
#include "parameters.h"
#include "splines.h"

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
    virtual void initParametersExperiment(py::module &m) { (void)m; }
    virtual void initPlottingExperiment(py::module &m) { (void)m; }
    virtual void initFittersExperiment(py::module &m) { (void)m; }
    virtual void initSamplesExperiment(py::module &m) { (void)m; }
    virtual void initManagerExperiment(py::module &m) { (void)m; }
    virtual void initSplinesExperiment(py::module &m) { (void)m; }


    // These functions will be called to set up each of the python modules
    // they call the core init<module name> functions and then the experiment
    // specific ones
    void initPlotting(py::module &m) {
        initPlottingModule(m); // <- defined in python/plotting.cpp
        initPlottingExperiment(m);
    }
    void initFitters(py::module &m) {
        initFittersModule(m); // <- defined in python/fitters.cpp
        initFittersExperiment(m);
    }
    void initSamples(py::module &m) {
        initSamplesModule(m); // <- defined in python/samples.cpp
        initSamplesExperiment(m);
    }
    void initManager(py::module &m) {
        initManagerModule(m); // <- defined in python/manager.cpp
        initManagerExperiment(m);
    }
    void initParameters(py::module &m) {
        initParametersModule(m); // <- defined in python/parameters.cpp
        initParametersExperiment(m);
    }
    void initSplines(py::module &m) {
        initSplinesModule(m); // <- defined in python/splines.cpp
        initSplinesExperiment(m);
    }

};

#define MAKE_PYMACH3_MDULE( PYBINDER_CLASS ) \
PYBIND11_MODULE( _pyMaCh3, m ) {             \
    PYBINDER_CLASS().initialise(m);          \
} 
