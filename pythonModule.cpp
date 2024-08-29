#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void initPlotting(py::module &); // <- defined in plotting/plottingUtils/pythonPlottingModule.cpp

PYBIND11_MODULE( pyMach3, m ) {
    initPlotting(m);
}