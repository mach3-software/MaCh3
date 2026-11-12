#pragma once

// MaCh3 includes
#include "Manager/Manager.h"
#include "Manager/Random.h"
//KS: Joy of forward declaration https://gieseanw.wordpress.com/2018/02/25/the-joys-of-forward-declarations-results-from-the-real-world/
class TObject;

/// @file HistogramUtils.h
/// @author Will Parker
/// @author Kamil Skwarczynski

/// @brief WP: Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly* poly);

/// @brief WP: Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly);

/// @brief WP: Poly Projectors
TH1D* PolyProjectionX(TObject* poly, const std::string& TempName, const std::vector<double>& xbins, const bool computeErrors = false);
/// @brief WP: Poly Projectors
TH1D* PolyProjectionY(TObject* poly, const std::string& TempName, const std::vector<double>& ybins, const bool computeErrors = false);

/// @brief KS: Convert TH2D to TH2Poly
TH2D* ConvertTH2PolyToTH2D(TH2Poly *poly, TH2D *TH2Dhist);
/// @brief KS: Convert TH2Poly to TH2D
TH2Poly* ConvertTH2DtoTH2Poly(TH2D *TH2Dhist);

/// @brief WP: Helper to Normalise histograms
TH2Poly* NormalisePoly(TH2Poly* Histogram);

/// @brief Helper to Normalise histograms
/// @param Histogram hist which we normalise
void NormaliseTH2Poly(TH2Poly* Histogram);

/// @brief WP: Helper to scale th2poly analogous to th2d scale with option "width"
TH2Poly* PolyScaleWidth(TH2Poly *Histogram, double scale);

/// @brief WP: Helper to calc integral of th2poly analogous to th2d integra; with option "width"
double PolyIntegralWidth(TH2Poly *Histogram);

/// @brief Helper to make ratio histograms
template<class HistType> HistType* RatioHists(HistType* NumHist, HistType* DenomHist);

/// @brief Helper to make ratio of TH2Polys
TH2Poly* RatioPolys(TH2Poly* NumPoly, TH2Poly* DenomPoly);

/// @brief WP: Helper function to create TH2Poly histogram with uniform binning
/// @param name This will be title of output histogram
/// @param BinArray_x Bin edges for X axis
/// @param BinArray_y Bin edges for Y axis
TH2Poly* MakePolyHist(const std::string& name, const std::vector<double>& BinArray_x, const std::vector<double>& BinArray_y);

/// @brief KS: ROOT changes something with binning when moving from ROOT 5 to ROOT 6. If you open ROOT5 produced file with ROOT6 you will be missing 9 last bins
/// @brief However if you use ROOT6 and have ROOT6 file exactly the same code will work. Something have changed with how TH2Poly bins are stored in TFile
/// @param file ROOT file that we will make version checks
void CheckTH2PolyFileVersion(TFile *file);

/// @brief KS: Remove fitted TF1 from hist to make comparison easier
void RemoveFitter(TH1D* hist, const std::string& name);

/// @brief Make Poisson fluctuation of TH1D hist using default fast method
void MakeFluctuatedHistogramStandard(TH1D *FluctHist, TH1D* PolyHist);
/// @brief Make Poisson fluctuation of TH1D hist using slow method which is only for cross-check
void MakeFluctuatedHistogramAlternative(TH1D *FluctHist, TH1D* PolyHist);

/// @brief Make Poisson fluctuation of TH2D hist using default fast method
void MakeFluctuatedHistogramStandard(TH2D *FluctHist, TH2D* PolyHist);
/// @brief Make Poisson fluctuation of TH2D hist
void MakeFluctuatedHistogramAlternative(TH2D *FluctHist, TH2D* PolyHist);

/// @brief Make Poisson fluctuation of TH2Poly hist using default fast method
void MakeFluctuatedHistogramStandard(TH2Poly *FluctHist, TH2Poly* PolyHist);
/// @brief Make Poisson fluctuation of TH2Poly hist using slow method which is only for cross-check
void MakeFluctuatedHistogramAlternative(TH2Poly *FluctHist, TH2Poly* PolyHist);

/// @brief KS: ROOT developers were too lazy do develop getRanom2 for TH2Poly, this implementation is based on [link](https://root.cern.ch/doc/master/classTH2.html#a883f419e1f6899f9c4255b458d2afe2e)
int GetRandomPoly2(const TH2Poly* PolyHist);

/// @brief KS: Fill Violin histogram with entry from a toy
/// @param violin hist that will be filled
/// @param hist_1d refence hist from which we take entries to be filled
void FastViolinFill(TH2D* violin, TH1D* hist_1d);

/// @brief Converts a vector of pointers from a derived type to a base type.
/// @tparam Derived The derived class type.
/// @tparam Base The base class type.
/// @param inputVec A `std::vector` of pointers to `Derived` objects.
/// @return A `std::vector` of pointers to `Base` objects.
template <typename Derived, typename Base>
std::vector<std::unique_ptr<Base>>
CastVector(std::vector<std::unique_ptr<Derived>>&& inputVec) {
  std::vector<std::unique_ptr<Base>> outputVec;
  outputVec.reserve(inputVec.size());
  for (auto& ptr : inputVec) {
    outputVec.push_back(std::unique_ptr<Base>(std::move(ptr)));
  }
  return outputVec;
}

/// @brief @brief Build a 1D posterior-predictive summary from a violin spectrum.
/// @param Spectra Input TH2D violin/density histogram.
/// @param name    Name of the output TH1D.
std::unique_ptr<TH1D> MakeSummaryFromSpectra(const TH2D* Spectra, const std::string& name);

namespace M3 {
/// @brief KS: Creates a copy of a ROOT-like object and wraps it in a smart pointer.
///
/// @tparam ObjectType The type of the object to clone for example TH1D or TH2Poly.
/// @param obj Pointer to the object to clone.
/// @param name Optional argument allowing to set new name of cloned object
/// @return std::unique_ptr<ObjectType> Owning pointer to the cloned object.
template <typename ObjectType>
std::unique_ptr<ObjectType> Clone(const ObjectType* obj, const std::string& name = "") {
  std::string cloneName = name.empty() ? obj->GetName() : name;

  std::unique_ptr<ObjectType> Hist(static_cast<ObjectType*>(obj->Clone(cloneName.c_str())));
  // Disable ROOT memory management because it causes lot of headache especially as smart pointers are much smarter
  Hist->SetDirectory(nullptr);

  return Hist;
}

/// @brief Scale histogram to get divided by bin width
void ScaleHistogram(TH1* Sample_Hist, const double scale);

/// @brief KS: Helper function check if data and MC binning matches
void CheckBinningMatch(const TH1D* Hist1, const TH1D* Hist2, const std::string& File, const int Line);
/// @brief KS: Helper function check if data and MC binning matches
void CheckBinningMatch(const TH2D* Hist1, const TH2D* Hist2, const std::string& File, const int Line);
/// @brief KS: Helper function check if data and MC binning matches
void CheckBinningMatch(TH2Poly* Hist1, TH2Poly* Hist2, const std::string& File, const int Line);
} //end M3
