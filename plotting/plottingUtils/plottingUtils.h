#ifndef PLOTTING_UTILS_H
#define PLOTTING_UTILS_H 1

// C++
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <map>

// ROOT
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TKey.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TKey.h"
#include "TH2.h"
#include "TH2Poly.h"

namespace MaCh3Plotting{
    /// \defgroup Utils Plotting Utility Functions
    /// @{

    /// @brief This handy little function lets you interpret a TGraph as a TH1D.
    /// @param graph The graph you want to convert.
    /// @param newName The new name you want to give to the histogram. If not specified, will just use the name of the graph.
    /// @param newTitle The new title you want to give to the histogram. If not specified, will just use the title of the graph.
    TH1D TGraphToTH1D(TGraph graph, std::string newName = "", std::string newTitle = "");


    /// @}
}

#endif
