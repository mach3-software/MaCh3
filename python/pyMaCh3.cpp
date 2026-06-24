#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

/// @file pyMaCh3.cpp
/// @author Ewan Miller

namespace py = pybind11;

void initPlotting(py::module &); // <- defined in python/plotting.cpp
void initFitters(py::module &); // <- defined in python/fitters.cpp
void initSamples(py::module &); // <- defined in python/samples.cpp
void initManager(py::module &); // <- defined in python/manager.cpp
void initParameters(py::module &); // <- defined in python/parameters.cpp
void initSplines(py::module &);  // <- defined in python/splines.cpp

PYBIND11_MODULE( _pyMaCh3, m ) {
    initPlotting(m);
    initFitters(m);
    initSamples(m);
    initManager(m);
    initParameters(m);
    initSplines(m);
}
