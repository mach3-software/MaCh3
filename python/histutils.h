#pragma once

#include "TH1.h"
#include "TH2.h"
#include "TH2Poly.h"

#include <pybind11/numpy.h>

/// @file histutils.h
/// @brief HW: Histogram utilities for converting ROOT histograms to numpy arrays for use in Python.
/// @author Henry Wallace
/// @author Ewan Miller

namespace py = pybind11;

using HistTuple = std::tuple<py::array_t<M3::float_t>, py::array_t<M3::float_t>, py::array_t<M3::float_t>>;

void FillEdgesPointer(py::array_t<M3::float_t> edges_buf, TAxis* axis, int nbins) {
    auto edges_ptr = static_cast<M3::float_t*>(edges_buf.request().ptr);
    for (int i = 0; i <= nbins; ++i) {
        edges_ptr[i] = axis->GetBinLowEdge(i + 1);
    }
    edges_ptr[nbins] = axis->GetBinLowEdge(nbins + 1) + axis->GetBinWidth(nbins + 1);
}

HistTuple THNToNumpy(std::unique_ptr<TH1>& hist){
    int nbinsX = hist->GetNbinsX();
    
    
    py::array_t<M3::float_t> edgesX(nbinsX + 1);
    FillEdgesPointer(edgesX, hist->GetXaxis(), nbinsX);
    
    py::array_t<M3::float_t> edgesY;
    
    // If we need to get the contents as a 2D array, we need to check if the histogram is actually 2D
    int nbinsY = hist->GetNbinsY();
    if(nbinsY > 1){
        edgesY = py::array_t<M3::float_t>(nbinsY + 1);
        FillEdgesPointer(edgesY, hist->GetYaxis(), nbinsY);
    }

    // Now we fill the contents
    py::array_t<M3::float_t> contents({nbinsY, nbinsX});
    auto contents_buf = contents.request();
    M3::float_t* contents_ptr = static_cast<M3::float_t*>(contents_buf.ptr);

    for(int iy=0; iy<nbinsY; ++iy){
        for(int ix=0; ix<nbinsX; ++ix){
            contents_ptr[iy * nbinsX + ix] = hist->GetBinContent(ix + 1, iy + 1);
        }
    }

    return std::make_tuple(contents, edgesX, edgesY);
}

inline py::tuple HistToNumpy(std::unique_ptr<TH1>& hist)
{
    if(dynamic_cast<TH2Poly*>(hist.get())){
        throw std::runtime_error("TH2Poly is not supported for conversion to numpy arrays");
    }

    return py::cast(THNToNumpy(hist));
}
