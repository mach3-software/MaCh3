#pragma once

// MaCh3 includes
#include "manager/MaCh3Logger.h"
#include "manager/MaCh3Exception.h"

// C++
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>

// ROOT
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "THStack.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

namespace MaCh3Plotting {
/// @defgroup Utils Plotting Utility Functions
/// @{

/// @brief This handy little function lets you interpret a TGraph as a TH1D.
/// @param graph The graph you want to convert.
/// @param newName The new name you want to give to the histogram. If not specified, will just use
/// the name of the graph.
/// @param newTitle The new title you want to give to the histogram. If not specified, will just use
/// the title of the graph.
TH1D TGraphToTH1D(TGraph graph, std::string newName = "", std::string newTitle = "");


/// @brief This handy little function lets you interpret a TGraph as a vector containing the same data.
/// @param graph The graph you want to convert.
/// @return A vector of vectors containing the data from the initial graph. The first vector is the x axis, the 2nd the y axis
std::vector<std::vector<float>> TGraphToVector(TGraph graph);

/// @}
} // namespace MaCh3Plotting
