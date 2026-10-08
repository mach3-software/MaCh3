#pragma once

/// @file plotting.h
/// @author Ewan Miller

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <string>

#include "Utils/Plotting/PlottingUtils.h"
#include "Utils/Plotting/PlottingManager.h"
#include "Utils/Plotting/StyleManager.h"
#include "Utils/Plotting/InputManager.h"

namespace py = pybind11;


void initPlottingModule(py::module &m_plotting){
    py::class_<M3::Plotting::PlottingManager>(m_plotting, "PlottingManager")
        .def(
            py::init<const std::string &>(),
            "construct a PlottingManager using the specified config file",
            py::arg("config_file_name")
        ) 
        
        .def(
            py::init(),
            "default constructor, will initialise the PlottingManager with the default plotting config"
        ) 
        
        .def(
            "initialise", 
            &M3::Plotting::PlottingManager::initialise,
            "initalise this PlottingManager"
        )
        
        .def(
            "usage", 
            &M3::Plotting::PlottingManager::usage,
            "Print a usage message for the current executable"
        )
        
        .def(
            "parse_inputs", 
            &M3::Plotting::PlottingManager::parseInputsVec,
            "Parse command line variables",
            py::arg("arguments")
        )
        
        .def(
            "set_exec", 
            &M3::Plotting::PlottingManager::setExec,
            "Set the name of the current executable, which will be used when getting executable specific options from the plotting config file",
            py::arg("exec_name")
        )
        
        .def(
            "get_file_name", 
            &M3::Plotting::PlottingManager::getFileName,
            "Get the path to a particular file",
            py::arg("input_file_id")
        )
        
        .def(
            "get_file_label", 
            &M3::Plotting::PlottingManager::getFileLabel,
            "Get the specified label of a particular input file",
            py::arg("input_file_id")
        )
        
        .def(
            "get_draw_options", 
            &M3::Plotting::PlottingManager::getDrawOptions,
            "Get any additional root drawing options specified by the user"
        )
        
        .def(
            "get_output_name", 
            py::overload_cast<const std::string &>(&M3::Plotting::PlottingManager::getOutputName),
            "Get the output name specified by the user, can specify an additional *suffix* to append to the file name but before the file extension",
            py::arg("suffix") = ""
        )
        
        .def(
            "get_file_names", 
            &M3::Plotting::PlottingManager::getFileNames,
            "Get the list of all file names"
        )
        
        .def(
            "get_file_labels", 
            &M3::Plotting::PlottingManager::getFileLabels,
            "Get the list of all file labels"
        )
        
        .def(
            "get_n_files", 
            &M3::Plotting::PlottingManager::getNFiles,
            "Get the number of specified files"
        )
        
        .def(
            "get_split_by_sample", 
            &M3::Plotting::PlottingManager::getSplitBySample,
            "Get whether or not the user has set the 'split by sample' (-s) option"
        )
        
        .def(
            "get_plot_ratios", 
            &M3::Plotting::PlottingManager::getPlotRatios,
            "Get whether or not the user specified the 'plot ratios' (-r) option"
        )
        
        .def(
            "get_draw_grid", 
            &M3::Plotting::PlottingManager::getDrawGrid,
            "Get wheter or not the user has specified the 'draw grid' (-g) option"
        )
        
        .def(
            "style", 
            &M3::Plotting::PlottingManager::style, py::return_value_policy::reference,
            "Get the StyleManager associated with this PlottingManager"
        )
        
        .def(
            "input", 
            &M3::Plotting::PlottingManager::input, py::return_value_policy::reference,
            "Get the InputManager associated with this PlottingManager"
        )

        // EM: I can't figure out how to add the getOption methods as theres not really a way to deduce the return type from yaml so leaving them out for now :/
        //     I think one solution would be to extend the PlottingManager class inside of python (add a pyPlottingManager (or something like that) that derives 
        //     from this one and add the functions to that) and then fold that python code into the module somehow. But that is currently beyond my pybinding abilities
        ;

    py::class_<M3::Plotting::InputManager>(m_plotting, "InputManager")
        .def(
            "print", 
            &M3::Plotting::InputManager::print,
            "Print a summary of everything this manager knows"
            )
        
        .def(
            "get_llh_scan", 
            &M3::Plotting::InputManager::getLLHScan,
            "Get the LLH scan for a particular parameter from a particular file",
            py::arg("input_file_id"), 
            py::arg("param_name"), 
            py::arg("LLH_type") = "total"
            )
        
        .def(
            "get_llh_scan_by_sample", 
            &M3::Plotting::InputManager::getSampleSpecificLLHScan,
            "Get the LLH scan for a particular parameter from a particular file for a particular sample",
            py::arg("input_file_id"), 
            py::arg("param"), 
            py::arg("sample")
            )
        
        .def(
            "get_enabled_llh", 
            &M3::Plotting::InputManager::getEnabledLLH,
            "Get whether a particular file has LLH scans for a particular parameter", 
            py::arg("input_file_id"), 
            py::arg("param"),
            py::arg("LLH_type") = "total"
            )
        
        .def(
            "get_enabled_llh_by_sample", 
            &M3::Plotting::InputManager::getEnabledLLHBySample,
            "Get whether a particular file has LLH scans for a particular parameter for a particular sample",
            py::arg("input_file_id"), 
            py::arg("param"), 
            py::arg("sample")
            )
        
        .def(
            "get_post_fit_error", 
            &M3::Plotting::InputManager::getPostFitError,
            "Get the post fit error for a parameter from a particular file",
            py::arg("input_file_id"), 
            py::arg("param"),
            py::arg("error_type") = ""
            )
        
        .def(
            "get_post_fit_value", 
            &M3::Plotting::InputManager::getPostFitValue,
            "Get the post fit value for a parameter from a particular file",
            py::arg("input_file_id"), 
            py::arg("param"),
            py::arg("error_type") = ""
            )
        
        .def(
            "get_known_parameters", 
            &M3::Plotting::InputManager::getKnownParameters,
            "Get all the parameters that this manager knows about. Useful for iterating over"
            )
        
        .def(
            "get_known_samples", 
            &M3::Plotting::InputManager::getKnownSamples,
            "Get all the samples that this manager knows about. Useful for iterating over"
            )
        
        .def(
            "get_tagged_parameters", 
            &M3::Plotting::InputManager::getTaggedParameters,
            "Get all the parameters whose tags match some specified list",
            py::arg("tags"),
            py::arg("check_type") = "all"
            )
        
        .def(
            "get_tagged_samples", 
            &M3::Plotting::InputManager::getTaggedSamples,
            "Get all the samples whose tags match some specified list",
            py::arg("tags"),
            py::arg("check_type") = "all"
            )
        
        .def(
            "get_n_input_files", 
            &M3::Plotting::InputManager::getNInputFiles,
            "Get the number of input files registered with this manager"
            )
        
        .def(
            "get_known_llh_parameters", 
            &M3::Plotting::InputManager::getKnownLLHParameters,
            "Get all the parameters that a file has LLH scans for",
            py::arg("file_id")
            )
        
        .def(
            "get_known_llh_samples", 
            &M3::Plotting::InputManager::getKnownLLHSamples,
            "Get all the samples that a file has individual LLH scans for",
            py::arg("file_id")
            )
        
        .def(
            "get_known_post_fit_parameters", 
            &M3::Plotting::InputManager::getKnownPostFitParameters,
            "Get all the parameters that a file has post fit values and errors for",
            py::arg("file_id")
            )
        
        .def(
            "get_known_MCMC_parameters", 
            &M3::Plotting::InputManager::getKnownMCMCParameters,
            "Get all the parameters that a file has MCMC chain entries for",
            py::arg("file_id")
            )
        
        .def(
            "get_known_1d_posterior_parameters", 
            &M3::Plotting::InputManager::getKnown1dPosteriorParameters,
            "Get all the parameters that a file has processed 1d posteriors for",
            py::arg("file_id")
            )
        
        .def(
            "get_MCMC_entry", 
            &M3::Plotting::InputManager::getMCMCentry,
            "Load up a particular step in the MCMC chain for a particular input file",
            py::arg("file_id"),
            py::arg("step")
            )
        
        .def(
            "get_MCMC_value", 
            &M3::Plotting::InputManager::getMCMCvalue,
            "Get the value of a particular parameter for the current entry (set by set_MCMC_entry) in the chain for a particular file",
            py::arg("file_id"),
            py::arg("param")
            )
        
        .def(
            "get_n_MCMC_entries", 
            &M3::Plotting::InputManager::getnMCMCentries,
            "Get the number of entries in the MCMC chain in a particular file"
            )
        
        .def(
            "get_1d_posterior",
            &M3::Plotting::InputManager::get1dPosterior,
            "Get the 1d posterior for a particular parameter from a particular file",
            py::arg("file_id"),
            py::arg("param")
            )
        ;

    py::class_<M3::Plotting::StyleManager>(m_plotting, "StyleManager")
        .def(
            "prettify_parameter_name", 
            &M3::Plotting::StyleManager::prettifyParamName,
            "Convert internally used parameter name to a nice pretty name that can be used in plots",
            py::arg("param")
        )
        
        .def(
            "prettify_sample_name", 
            &M3::Plotting::StyleManager::prettifyParamName,
            "Convert internally used sample name to a nice pretty name that can be used in plots",
            py::arg("sample")
        )
        ;
}
