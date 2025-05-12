// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// MaCh3 includes
#include "Manager/Manager.h"
// yaml includes
#include "yaml-cpp/yaml.h"

namespace py = pybind11;

void initManager(py::module &m){

    auto m_manager = m.def_submodule("manager");
    m_manager.doc() = 
        "This is a Python binding of MaCh3s C++ based manager library.";

    // Bind some of the cpp-yaml library 
    // shamelessly stolen from stackoverflow: https://stackoverflow.com/questions/62347521/using-pybind11-to-wrap-yaml-cpp-iterator
    py::enum_<YAML::NodeType::value>(m_manager, "NodeType")
        .value("Undefined", YAML::NodeType::Undefined)
        .value("Null", YAML::NodeType::Null)
        .value("Scalar", YAML::NodeType::Scalar)
        .value("Sequence", YAML::NodeType::Sequence)
        .value("Map", YAML::NodeType::Map);

    py::class_<YAML::Node>(m_manager, "YamlNode")
        .def(py::init<const std::string &>())

        .def("data",
            [](const YAML::Node node){
                if ( node.Type() != YAML::NodeType::Scalar )
                {
                    throw MaCh3Exception(__FILE__, __LINE__, "Attempting to access the data of non-scalar yaml node. This is undefined.");
                }
                return node.Scalar();
            }, 
            "Access the data stored in the node. This is only valid if the node is a 'scalar' type, i.e. it is a leaf of the yaml tree structure.")

        .def("__getitem__",
            [](const YAML::Node node, const std::string& key){
              return node[key];
            })

        .def("__getitem__",
            [](const YAML::Node node, const int& key){
              if ( node.Type() != YAML::NodeType::Sequence)
              {
                throw MaCh3Exception(__FILE__, __LINE__, "Trying to access a non sequence yaml node with integer index");
              }
              return node[key];
            })

        .def("__iter__",
            [](const YAML::Node &node) {
              return py::make_iterator(node.begin(), node.end());},
             py::keep_alive<0, 1>())

        .def("__str__",
             [](const YAML::Node& node) {
               YAML::Emitter out;
               out << node;
               return std::string(out.c_str());
             })

        .def("type", &YAML::Node::Type)
        
        .def("__len__", &YAML::Node::size)
        ;

    py::class_<YAML::detail::iterator_value, YAML::Node>(m_manager, "_YamlDetailIteratorValue")
        .def(py::init<>())
        .def("first", [](YAML::detail::iterator_value& val) { return val.first;})
        .def("second", [](YAML::detail::iterator_value& val) { return val.second;})
        ;

    m.def("load_file", &YAML::LoadFile, "");

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
        
        .def(
            "raw", 
            &manager::raw,
            "Get the raw yaml config."
        )

        .def(
            "get_test_stat",
            &manager::GetMCStatLLH,
            "Get the test statistic that was specified in the config file."
        )
    ;

}
