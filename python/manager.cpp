#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "manager/manager.h"

namespace py = pybind11;

void initManager(py::module &m){

    auto m_manager = m.def_submodule("manager");
    m_manager.doc() = 
        "This module really just exists to hold the :py:class:`pyMaCh3.manager.Manager` class. The reason for this slightly weird structure is to mimic the structure of the c++ version of Mach3. \
        You can read more about the manager and config files on [the wiki page](https://github.com/mach3-software/MaCh3/wiki/01.-Manager-and-config-handling). \
        Happy managing!";

    
    py::class_<manager>(m_manager, "Manager")
        .def(
            py::init<std::string const &>(),
            "create a Manager object with a specified *config_file*",
            py::arg("config_file")
        )
        
        .def(
            "print", 
            &manager::Print, 
            "Print currently used config."
        )
    ;

};
