#pragma once

#include "samplePDFFDBase.h"

#include "TH1.h"

#include <memory>
#include <string>
#include <vector>

// LP - These methods take quite a few arguments so that they can exist as util
// methods and not clutter samplePDFFDBase

// LP - if generic binning is 1D, return a nice histogram, if its ND, return
// the global bin histogram.
// arguments:
//            htitle - the title of the constructed TH1. N.B. this can be used
//            to set axis labels as per usual ROOT title rules:
//                       "<plot title>[;<xlabel>[;<ylabel>[;<zlabel>]]]"
//            divide_by_ND_hypervolume - if true the area under the bin
//            corresponds to the count
std::unique_ptr<TH1> GetGenericBinningTH1(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle = "",
                                          bool divide_by_ND_hypervolume = true);

// LP - if generic binning is 2D, return a nice histogram, otherwise throw as
// there is nothing useful we can do here
std::unique_ptr<TH2> GetGenericBinningTH2(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle = "",
                                          bool divide_by_ND_hypervolume = true);

// LP - if generic binning is 3D, return a nice histogram, otherwise throw as
// there is nothing useful we can do here
std::unique_ptr<TH3> GetGenericBinningTH3(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle = "",
                                          bool divide_by_ND_hypervolume = true);

// LP - if generic binning is >1D, return a 1D slice from the slice_definition
// slice definition is a vector specifying a value along each axis to slice
// along, one entry in the vector should be set to kSliceAx to signify the axis
// to show
std::unique_ptr<TH1> GetGenericBinningTH1Slice(
    samplePDFFDBase &sample, std::vector<double> slice_definition,
    std::string const &hname, std::string const &htitle = "",
    bool divide_by_ND_hypervolume = true);

// LP - if generic binning is >2D, return a 2D slice from the slice_definition
// slice definition is a vector specifying a value along each axis to slice
// along, two entries in the vector should be set to kSliceAx to signify the
// axis to show
std::unique_ptr<TH2> GetGenericBinningTH2Slice(
    samplePDFFDBase &sample, std::vector<double> slice_definition,
    std::string const &hname, std::string const &htitle = "",
    bool divide_by_ND_hypervolume = true);

std::vector<std::unique_ptr<TH1>> GetGenericBinningTH1Slices(
    samplePDFFDBase &sample, size_t slice_ax, std::string const &hname_pat,
    std::string const &htitle = "", bool divide_by_ND_hypervolume = true);

std::vector<std::unique_ptr<TH2>> GetGenericBinningTH2Slices(
    samplePDFFDBase &sample, std::array<size_t, 2> slice_axes,
    std::string const &hname_pat, std::string const &htitle = "",
    bool divide_by_ND_hypervolume = true);

constexpr double kSliceAx = 0xdeadbeef;