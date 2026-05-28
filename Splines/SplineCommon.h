#pragma once

/// @file SplineCommon.h
/// @brief Contains definitions for spline coefficients and structure used in both CPU and GPU code.
///
/// @details This file includes macros and enums for defining spline coefficients.
/// It is designed to be compatible with older CUDA versions, so be cautious
/// when adding new features or including other headers.
///
/// @warning KS: Please add stuff here with super caution. This header is being added to gpuSplineUtils.cu. Right now we support most of CUDA even super old.
/// If you add some header with fancy templates it will not compile for older CUDA.
/// This header is a way to use common macros or Enum in CPU and GPU code. For more sophisticated structs please use SplineStructs.h
///
/// @author Clarence Wret
/// @author Kamil Skwarczynski

/// KS: We store coefficients {y,b,c,d} in one array one by one, this is only to define it once rather then insert "4" all over the code
constexpr int _nCoeff_ = 4;
/// KS: For TF1 we store at most 5 coefficients, we could make it more flexible but for now define it here to make future changes easier to track
constexpr int _nTF1Coeff_ = 2;

/// HW: Coefficients for grabbing items from manycoeff_arr (rather than having y=manycoeffarray[index+0])
enum SplineSegmentCoeffs
{
  kCoeffY = 0, ///< Coefficient Y
  kCoeffB = 1, ///< Coefficient B
  kCoeffC = 2, ///< Coefficient C
  kCoeffD = 3  ///< Coefficient D
};

// *******************
/// @brief Flat representation of a spline index entry
struct SplineIndex {
// *******************
  /// @brief destructor
  virtual ~SplineIndex() = default;

  /// @brief Compare spline lookup indices (excluding value)
  /// @todo Add this to BinnedSplineHandler::GetSplineIndex
  bool operator==(const SplineIndex& idx) const
  {
    return iSample  == idx.iSample  &&
           iOscChan == idx.iOscChan &&
           iSyst    == idx.iSyst    &&
           iMode    == idx.iMode    &&
           iVar1    == idx.iVar1    &&
           iVar2    == idx.iVar2    &&
           iVar3    == idx.iVar3;
  }

  /// Index into the flattened spline weight vector
  int value    = 0;
  /// Sample index
  int iSample  = 0;
  /// Oscillation channel index
  int iOscChan = 0;
  /// Systematic parameter index
  int iSyst    = 0;
  /// Mode index within a systematic
  int iMode    = 0;
  /// First kinematic bin index
  int iVar1    = 0;
  /// Second kinematic bin index
  int iVar2    = 0;
  /// Third kinematic bin index
  int iVar3    = 0;

  #ifndef __CUDACC__
  // Include ClassDef macro for ROOT dictionary generation, but only in C++ code
  ClassDef(SplineIndex, 1);
  #endif
};

// *******************
/// @brief KS: Struct storing information for spline monolith
/// @details This structure holds the X coefficients, other spline coefficients,
/// the number of knots per spline, and the number of points per spline on the CPU.
struct SplineMonoStruct {
// *******************
  /// @brief destructor
  virtual ~SplineMonoStruct() = default;

  /// KS: CPU arrays to hold X coefficient
  std::vector<float> coeff_x;

  /// CPU arrays to hold other coefficients
  std::vector<float> coeff_many;

  /// KS: CPU Number of knots per spline
  std::vector<unsigned int> nKnots_arr;

  /// CW: CPU array with the number of points per spline (not per spline point!)
  std::vector<short int> paramNo_arr;

  #ifndef __CUDACC__
  // Include ClassDef macro for ROOT dictionary generation, but only in C++ code
  ClassDef(SplineMonoStruct, 1);
  #endif
};

