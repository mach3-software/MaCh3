#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void initPlotting(py::module &); // <- defined in python/plotting.cpp
void initFitter(py::module &); // <- defined in python/fitter.cpp
void initSamplePDF(py::module &); // <- defined in python/samplePDF.cpp
void initManager(py::module &); // <- defined in python/manager.cpp
void initCovariance(py::module &); // <- defined in python/covariance.cpp
void initSplines(py::module &); // <- defined in python/splines.cpp

PYBIND11_MODULE( _pyMaCh3, m ) {
    initPlotting(m);
    initFitter(m);
    initSamplePDF(m);
    initManager(m);
    initCovariance(m);
    initSplines(m);
}