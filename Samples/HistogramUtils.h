#pragma once

// MaCh3 includes
#include "Samples/SampleStructs.h"
#include "Parameters/ParameterStructs.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TObjString.h"
#include "TRandom3.h"
_MaCh3_Safe_Include_End_ //}

/// @file HistogramUtils.h
/// @author Will Parker
/// @author Kamil Skwarczynski

/// @brief WP: Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly* poly);

/// @brief WP: Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly);

/// @brief WP: Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, const std::vector<double>& xbins, const bool computeErrors = false);
/// @brief WP: Poly Projectors
TH1D* PolyProjectionY(TObject* poly, std::string TempName, const std::vector<double>& ybins, const bool computeErrors = false);

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
void MakeFluctuatedHistogramStandard(TH1D *FluctHist, TH1D* PolyHist, TRandom3* rand);
/// @brief Make Poisson fluctuation of TH1D hist using slow method which is only for cross-check
void MakeFluctuatedHistogramAlternative(TH1D *FluctHist, TH1D* PolyHist, TRandom3* rand);

/// @brief Make Poisson fluctuation of TH2D hist using default fast method
void MakeFluctuatedHistogramStandard(TH2D *FluctHist, TH2D* PolyHist, TRandom3* rand);
/// @brief Make Poisson fluctuation of TH2D hist
void MakeFluctuatedHistogramAlternative(TH2D *FluctHist, TH2D* PolyHist, TRandom3* rand);

/// @brief Make Poisson fluctuation of TH2Poly hist using default fast method
void MakeFluctuatedHistogramStandard(TH2Poly *FluctHist, TH2Poly* PolyHist, TRandom3* rand);
/// @brief Make Poisson fluctuation of TH2Poly hist using slow method which is only for cross-check
void MakeFluctuatedHistogramAlternative(TH2Poly *FluctHist, TH2Poly* PolyHist, TRandom3* rand);

/// @brief KS: ROOT developers were too lazy do develop getRanom2 for TH2Poly, this implementation is based on [link](https://root.cern.ch/doc/master/classTH2.html#a883f419e1f6899f9c4255b458d2afe2e)
int GetRandomPoly2(const TH2Poly* PolyHist, TRandom3* rand);

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
std::vector<Base*> CastVector(const std::vector<Derived*>& inputVec) {
  std::vector<Base*> outputVec;
  // Reserve space for efficiency
  outputVec.reserve(inputVec.size());
  for (auto* ptr : inputVec) {
    outputVec.push_back(static_cast<Base*>(ptr));
  }
  return outputVec;
}

/// @brief Helper to check if files exist or not
inline std::string file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    MACH3LOG_ERROR("File {} does not exist", filename);
    MACH3LOG_ERROR("Please try again");;
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return filename;
}

/// @brief DB Get the Cherenkov momentum threshold in MeV
/// @param PDG PDG code of the particle for which the Cherenkov threshold is requested.
double returnCherenkovThresholdMomentum(const int PDG);

/// @brief Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
/// @param PLep Shifted outgoing lepton momentum (MeV/c).
/// @param PUpd Updated lepton angle or related kinematic quantity used in the recalculation.
/// @param EnuTrue True neutrino energy (MeV).
/// @param InitialQ2 Optional initial Q^2 value (used as seed or fallback, default = 0).
double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2 = 0.0);

/// @brief Recalculate Enu after Eb shift. Takes in shifted lepton momentum, lepton angle, and binding energy change, and if nu/anu
/// @param PLep Shifted outgoing lepton momentum (MeV/c).
/// @param cosTheta Cosine of the lepton scattering angle.
/// @param EB Binding energy shift applied to the interaction (MeV).
/// @param neutrino True if the interaction is neutrino, false if antineutrino.
double CalculateEnu(double PLep, double cosTheta, double EB, bool neutrino);

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

/// @brief Opens a ROOT file with the given name and mode.
///
/// This function wraps ROOT’s `TFile` constructor and checks whether the file was opened
/// successfully and give some useful debugging information
///
/// @param Name The name or path of the file to open.
/// @param Type The file open mode (e.g., "READ", "RECREATE", "UPDATE").
/// @param File The name of the file where the exception occurred.
/// @param Line The line number where the exception occurred.
/// @return Pointer to the opened ROOT file.
TFile* Open(const std::string& Name, const std::string& Type, const std::string& File, const int Line);

/// @brief Scale histogram to get divided by bin width
void ScaleHistogram(TH1* Sample_Hist, const double scale);

} //end M3
