#pragma once
#include <vector>
#include <string>

#include "Samples/SampleStructs.h"

/// @brief Stores info about each MC event used during reweighting routine
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Kamil Skwarczynski
/// @warning Try to no add more variables here as it will impact RAM usage
/// @ingroup SamplesAndParameters
class EventInfo {
public:
  /// @brief Default constructor.
  EventInfo(){}

  /// @brief Copy constructor (deleted to prevent copying).
  EventInfo(EventInfo const &other) = delete;

  /// @brief Move constructor (defaulted to allow moving).
  EventInfo(EventInfo &&other) = default;

  /// @brief Copy assignment operator (deleted).
  EventInfo& operator=(EventInfo const &other) = delete;

  /// @brief Move assignment operator (deleted).
  EventInfo& operator=(EventInfo &&other) = delete;

  /// @brief default destructor
  ~EventInfo(){}

  ///
  void AddWeightPointer(const M3::float_t* Ptr, const std::string& WeightName="") {
    total_weight_pointers.push_back(Ptr);

#ifdef MACH3_DEBUG
    if (WeightName == "") {
      throw;
    }
    
    weight_names.push_back(WeightName);
#endif
    (void)WeightName;
  }

  ///
  M3::float_t CalcWeightTotal() const _noexcept_ {
    M3::float_t TotalWeight = 1.0;
    const int TotalWeights = static_cast<int>(total_weight_pointers.size());

    //DB Loop over stored pointers
#ifdef MULTITHREAD
#pragma omp simd reduction(*:TotalWeight)
#endif
    for (int iWeight = 0; iWeight < TotalWeights; ++iWeight) {
      TotalWeight *= *(total_weight_pointers[iWeight]);
    }
    
    return TotalWeight;    
  }

  ///
  std::string ReturnDumpedEventWeightString() {
    std::string ReturnString;
    
#ifdef MACH3_DEBUG
    for (int iWeight = 0; iWeight < TotalWeights; ++iWeight) {
      ReturnString += weight_names[iWeight] + ":" + Form("%3.4f, ",*(total_weight_pointers[iWeight]));
    }
    ReturnString += "total_weight:" + Form("%3.4f",CalcWeightTotal());
#else
    throw;
#endif

    return ReturnString;
  }

  /// The x_var and y_vars and beyond that you're binning in
  std::vector<const double*> KinVar;

  /// starting bins for each dimensions allowing to perform quick lookup
  std::vector<int> NomBin;

  /// Nominal sample to which event is associated
  int NominalSample = M3::_BAD_INT_;

  /// PDG of neutrino after oscillation
  int nupdg  = M3::_BAD_INT_;

  /// PDG of neutrino before oscillation
  int nupdgUnosc = M3::_BAD_INT_;

  /// Is event NC or not
  bool isNC = false;

  /// Pointer to true Neutrino Energy
  double enu_true = M3::_BAD_DOUBLE_;

  /// Pointer to true cosine zenith
  double coszenith_true = M3::_BAD_DOUBLE_;

private:
  /// Pointers to weights like oscillation spline, normalisation etc
  std::vector<const M3::float_t*> total_weight_pointers;

#ifdef MACH3_DEBUG
  /// Names of weights used to calculate total weight
  std::vector<const std:string> weight_names;
#endif
};
