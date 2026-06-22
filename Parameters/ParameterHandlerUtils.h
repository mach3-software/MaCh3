#pragma once

// MaCh3 includes
#include "Parameters/ParameterStructs.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TDecompChol.h"
#include "TStopwatch.h"
#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDEigen.h"
#include "TDecompSVD.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"
#include "TText.h"
#include "TLegend.h"
_MaCh3_Safe_Include_End_ //}


/// @file ParameterHandlerUtils.h
/// @author Clarence Wret
/// @author Dan Barrow
/// @author Kamil Skwarczynski

namespace M3
{
/// @brief Matches a string against a simple wildcard Pattern using regex. Is not case sensitive
/// @param Text    Input string to test.
/// @param Pattern Wildcard pattern to match against.
inline bool RegexMatch(std::string Text, std::string Pattern) {
  // Make a copy and to lower case to not be case sensitive
  std::transform(Text.begin(), Text.end(), Text.begin(), ::tolower);

  // Convert to low case to not be case sensitive
  std::transform(Pattern.begin(), Pattern.end(), Pattern.begin(), ::tolower);
  try {
    // Replace '*' in the Pattern with '.*' for regex matching
    std::string RegexPattern = "^" + std::regex_replace(Pattern, std::regex("\\*"), ".*") + "$";
    std::regex Regex(RegexPattern);
    return std::regex_match(Text, Regex);
  }
  catch (const std::regex_error& e) {
    MACH3LOG_ERROR("Regex error: {}", e.what());
    return false;
  }
}

/// @brief Matches a string against a simple wildcard Pattern using regex. Is not case sensitive
/// @param Text    Input string to test.
/// @param Patterns Collection wildcard patterns to match against.
inline bool RegexMatch(std::string Text, const std::vector<std::string>& Patterns) {
  for (size_t i = 0; i < Patterns.size(); i++) {
    if (M3::RegexMatch(Text, Patterns[i])) {
      return true;
    }
  }
  return false;
}

/// @brief CW: Multi-threaded matrix multiplication
inline double* MatrixMult(double *A, double *B, int n) {
  //CW: First transpose to increse cache hits
  double *BT = new double[n*n];
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      BT[j*n+i] = B[i*n+j];
    }
  }

  // Now multiply
  double *C = new double[n*n];
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      for (int k = 0; k < n; k++) {
        sum += A[i*n+k]*BT[j*n+k];
      }
      C[i*n+j] = sum;
    }
  }
  delete BT;

  return C;
}

/// @brief CW: Multi-threaded matrix multiplication
inline double** MatrixMult(double **A, double **B, int n) {
  // First make into monolithic array
  double *A_mon = new double[n*n];
  double *B_mon = new double[n*n];

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A_mon[i*n+j] = A[i][j];
      B_mon[i*n+j] = B[i][j];
    }
  }
  //CW: Now call the monolithic calculator
  double *C_mon = MatrixMult(A_mon, B_mon, n);
  delete A_mon;
  delete B_mon;

  // Return the double pointer
  double **C = new double*[n];
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; ++i) {
    C[i] = new double[n];
    for (int j = 0; j < n; ++j) {
      C[i][j] = C_mon[i*n+j];
    }
  }
  delete C_mon;

  return C;
}

/// @brief CW: Multi-threaded matrix multiplication
inline TMatrixD MatrixMult(TMatrixD A, TMatrixD B)
{
  double *C_mon = MatrixMult(A.GetMatrixArray(), B.GetMatrixArray(), A.GetNcols());
  TMatrixD C;
  C.Use(A.GetNcols(), A.GetNrows(), C_mon);
  return C;
}

// ********************************************
/// @brief KS: Custom function to perform multiplication of matrix and vector with multithreading
/// @param VecMulti Output Vector, VecMulti = matrix x vector
/// @param matrix This matrix is used for multiplication VecMulti = matrix x vector
/// @param vector This vector is used for multiplication VecMulti = matrix x vector
/// @param n this is size of matrix and vector, we assume matrix is symmetric
inline void MatrixVectorMulti(double* _restrict_ VecMulti, double** _restrict_ matrix, const double* _restrict_ vector, const int n) {
// ********************************************
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; ++i)
  {
    double result = 0.0;
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (int j = 0; j < n; ++j)
    {
      result += matrix[i][j]*vector[j];
    }
    VecMulti[i] = result;
  }
}

// ********************************************
/// @brief KS: Custom function to perform multiplication of matrix and single element which is thread safe
/// @param matrix This matrix is used for multiplication VecMulti = matrix x vector
/// @param vector This vector is used for multiplication VecMulti = matrix x vector
/// @param Length this is size of matrix and vector, we assume matrix is symmetric
/// @param i Element of matrix that we want to multiply
inline double MatrixVectorMultiSingle(double** _restrict_ matrix, const double* _restrict_ vector, const int Length, const int i) {
// ********************************************
  double Element = 0.0;
  #ifdef MULTITHREAD
  #pragma omp simd
  #endif
  for (int j = 0; j < Length; ++j) {
    Element += matrix[i][j]*vector[j];
  }
  return Element;
}

// *************************************
/// @brief KS: Yaml emitter has problem and drops "", if you have special signs in you like * then there is problem. This bit hacky code adds these ""
/// @param yamlStr The YAML string to be processed (modified in-place).
inline void FixSampleNamesQuotes(std::string& yamlStr) {
// *************************************
  std::stringstream input(yamlStr);
  std::string line;
  std::string fixedYaml;
  std::regex sampleNamesRegex(R"(SampleNames:\s*\[([^\]]+)\])");

  while (std::getline(input, line)) {
    std::smatch match;
    if (std::regex_search(line, match, sampleNamesRegex)) {
      std::string contents = match[1];  // inside the brackets
      std::stringstream ss(contents);
      std::string item;
      std::vector<std::string> quotedItems;

      while (std::getline(ss, item, ',')) {
        item = std::regex_replace(item, std::regex(R"(^\s+|\s+$)"), ""); // trim
        quotedItems.push_back("\"" + item + "\"");
      }

      std::string replacement = "SampleNames: [" + fmt::format("{}", fmt::join(quotedItems, ", ")) + "]";
      line = std::regex_replace(line, sampleNamesRegex, replacement);
    }
    fixedYaml += line + "\n";
  }

  yamlStr = fixedYaml;
}

// *************************************
/// @brief KS: Add Tune values to YAML covariance matrix
/// @param root The root YAML node to be updated.
/// @param Values The values to add for the specified tune.
/// @param Tune The name of the tune (e.g., "PostFit").
/// @param FancyNames Optional list of fancy names to match systematics (must match Values size if provided).
inline void AddTuneValues(YAML::Node& root,
                          const std::vector<double>& Values,
                          const std::string& Tune,
                          const std::vector<std::string>& FancyNames = {}) {
// *************************************
  YAML::Node NodeCopy = YAML::Clone(root);
  YAML::Node systematics = NodeCopy["Systematics"];

  if (!systematics || !systematics.IsSequence()) {
    MACH3LOG_ERROR("'Systematics' node is missing or not a sequence in the YAML copy");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (!FancyNames.empty() && FancyNames.size() != Values.size()) {
    MACH3LOG_ERROR("Mismatch in sizes: FancyNames has {}, but Values has {}", FancyNames.size(), Values.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (FancyNames.empty() && systematics.size() != Values.size()) {
    MACH3LOG_ERROR("Mismatch in sizes: Values has {}, but YAML 'Systematics' has {} entries",
                   Values.size(), systematics.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (!FancyNames.empty()) {
    for (std::size_t i = 0; i < FancyNames.size(); ++i) {
      bool matched = false;
      for (std::size_t j = 0; j < systematics.size(); ++j) {
        YAML::Node systematicNode = systematics[j]["Systematic"];
        if (!systematicNode) continue;
        auto nameNode = systematicNode["Names"];
        if (!nameNode || !nameNode["FancyName"]) continue;
        if (nameNode["FancyName"].as<std::string>() == FancyNames[i]) {
          if (!systematicNode["ParameterValues"]) {
            MACH3LOG_ERROR("Missing 'ParameterValues' for matched FancyName '{}'", FancyNames[i]);
            throw MaCh3Exception(__FILE__, __LINE__);
          }
          systematicNode["ParameterValues"][Tune] = M3::Utils::FormatDouble(Values[i], 4);
          matched = true;
          break;
        }
      }
      if (!matched) {
        MACH3LOG_ERROR("Could not find a matching FancyName '{}' in the systematics", FancyNames[i]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  } else {
    for (std::size_t i = 0; i < systematics.size(); ++i) {
      YAML::Node systematicNode = systematics[i]["Systematic"];
      if (!systematicNode || !systematicNode["ParameterValues"]) {
        MACH3LOG_ERROR("Missing 'Systematic' or 'ParameterValues' entry at index {}", i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      systematicNode["ParameterValues"][Tune] = M3::Utils::FormatDouble(Values[i], 4);
    }
  }

  // Convert updated copy to string
  std::string YAMLString = YAMLtoSTRING(NodeCopy);
  FixSampleNamesQuotes(YAMLString);
  // Write to output file
  std::string OutName = "UpdatedMatrixWithTune" + Tune + ".yaml";
  std::ofstream outFile(OutName);
  if (!outFile) {
    MACH3LOG_ERROR("Failed to open file for writing: {}", OutName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  outFile << YAMLString;
  outFile.close();
}

// *************************************
/// @brief KS: Replace correlation matrix and tune values in YAML covariance matrix
/// @param root The root YAML node to be updated.
/// @param Values The new values for each systematic.
/// @param Errors The new errors for each systematic.
/// @param Correlation The new correlation matrix (must be square and match Values size).
/// @param OutYAMLName The output filename for the updated YAML.
/// @param FancyNames Optional list of fancy names to match systematics (must match Values size if provided).
inline void MakeCorrelationMatrix(YAML::Node& root,
                                  const std::vector<double>& Values,
                                  const std::vector<double>& Errors,
                                  const std::vector<std::vector<double>>& Correlation,
                                  const std::string& OutYAMLName,
                                  const std::vector<std::string>& FancyNames = {}) {
// *************************************
  if (Values.size() != Errors.size() || Values.size() != Correlation.size()) {
    MACH3LOG_ERROR("Size mismatch between Values, Errors, and Correlation matrix");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (const auto& row : Correlation) {
    if (row.size() != Correlation.size()) {
      MACH3LOG_ERROR("Correlation matrix is not square");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  YAML::Node NodeCopy = YAML::Clone(root);
  YAML::Node systematics = NodeCopy["Systematics"];

  if (!systematics || !systematics.IsSequence()) {
    MACH3LOG_ERROR("'Systematics' node is missing or not a sequence");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (!FancyNames.empty() && FancyNames.size() != Values.size()) {
    MACH3LOG_ERROR("FancyNames size ({}) does not match Values size ({})", FancyNames.size(), Values.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Map from FancyName to Systematic node
  std::unordered_map<std::string, YAML::Node> nameToNode;
  for (std::size_t i = 0; i < systematics.size(); ++i) {
    YAML::Node syst = systematics[i]["Systematic"];
    if (!syst || !syst["Names"] || !syst["Names"]["FancyName"]) continue;
    std::string name = syst["Names"]["FancyName"].as<std::string>();
    nameToNode[name] = syst;
  }

  if (!FancyNames.empty()) {
    for (std::size_t i = 0; i < FancyNames.size(); ++i) {
      const std::string& name_i = FancyNames[i];
      auto it_i = nameToNode.find(name_i);
      if (it_i == nameToNode.end()) {
        MACH3LOG_ERROR("Could not find FancyName '{}' in YAML", name_i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      YAML::Node& syst_i = it_i->second;

      syst_i["ParameterValues"]["PreFitValue"] = M3::Utils::FormatDouble(Values[i], 4);
      syst_i["Error"] = M3::Utils::FormatDouble(Errors[i], 4);

      YAML::Node correlationsNode = YAML::Node(YAML::NodeType::Sequence);
      for (std::size_t j = 0; j < FancyNames.size(); ++j) {
        if (i == j) continue;
        // KS: Skip if value close to 0
        if (std::abs(Correlation[i][j]) < 1e-8) continue;
        YAML::Node singleEntry;
        singleEntry[FancyNames[j]] = M3::Utils::FormatDouble(Correlation[i][j], 4);
        correlationsNode.push_back(singleEntry);
      }
      syst_i["Correlations"] = correlationsNode;
    }
  } else {
    if (systematics.size() != Values.size()) {
      MACH3LOG_ERROR("Mismatch in sizes: Values has {}, but YAML 'Systematics' has {} entries",
                     Values.size(), systematics.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (std::size_t i = 0; i < systematics.size(); ++i) {
      YAML::Node syst = systematics[i]["Systematic"];
      if (!syst) {
        MACH3LOG_ERROR("Missing 'Systematic' node at index {}", i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      syst["ParameterValues"]["PreFitValue"] = M3::Utils::FormatDouble(Values[i], 4);
      syst["Error"] = M3::Utils::FormatDouble(Errors[i], 4);

      YAML::Node correlationsNode = YAML::Node(YAML::NodeType::Sequence);
      for (std::size_t j = 0; j < Correlation[i].size(); ++j) {
        if (i == j) continue;
        // KS: Skip if value close to 0
        if (std::abs(Correlation[i][j]) < 1e-8) continue;
        YAML::Node singleEntry;
        const std::string& otherName = systematics[j]["Systematic"]["Names"]["FancyName"].as<std::string>();
        singleEntry[otherName] = M3::Utils::FormatDouble(Correlation[i][j], 4);
        correlationsNode.push_back(singleEntry);
      }
      syst["Correlations"] = correlationsNode;
    }
  }

  // Convert and write
  std::string YAMLString = YAMLtoSTRING(NodeCopy);
  FixSampleNamesQuotes(YAMLString);
  std::ofstream outFile(OutYAMLName);
  if (!outFile) {
    MACH3LOG_ERROR("Failed to open file for writing: {}", OutYAMLName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  outFile << YAMLString;
  outFile.close();
}

// *************************************
/// @brief KS: We store configuration macros inside the chain.
/// In the past, multiple configs were stored, which required error-prone hardcoding like "Config_xsec_cov".
/// Therefore, this code maintains backward compatibility by checking the number of macros present and
/// using a hardcoded name only if necessary.
inline TMacro* GetConfigMacroFromChain(TDirectory* CovarianceFolder) {
// *************************************
  if (!CovarianceFolder) {
    MACH3LOG_ERROR("Null TDirectory passed to {}", __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TMacro* foundMacro = nullptr;
  int macroCount = 0;

  TIter next(CovarianceFolder->GetListOfKeys());
  TKey* key;
  while ((key = dynamic_cast<TKey*>(next()))) {
    if (std::string(key->GetClassName()) == "TMacro") {
      ++macroCount;
      if (macroCount == 1) {
        foundMacro = dynamic_cast<TMacro*>(key->ReadObj());
      }
    }
  }

  if (macroCount == 1 && foundMacro) {
    MACH3LOG_INFO("Found single TMacro in directory: using it.");
    return foundMacro;
  } else {
    MACH3LOG_WARN("Found {} TMacro objects. Using hardcoded macro name: Config_xsec_cov.", macroCount);
    TMacro* fallback = CovarianceFolder->Get<TMacro>("Config_xsec_cov");
    if (!fallback) {
      MACH3LOG_WARN("Fallback macro 'Config_xsec_cov' not found in directory.");
    }
    return fallback;
  }
}

// *************************************
/// @brief KS: Retrieve the cross-section covariance matrix from the given TDirectory.
/// Historically, multiple covariance matrices could be stored, requiring fragile hardcoded paths like "CovarianceFolder/xsec_cov".
/// This function maintains backward compatibility by:
/// - Using the single matrix present if only one exists,
/// - Otherwise falling back to the hardcoded path.
/// This avoids error-prone assumptions while supporting both old and new formats.
inline TMatrixDSym* GetCovMatrixFromChain(TDirectory* TempFile) {
// *************************************
  if (!TempFile) {
    MACH3LOG_ERROR("Null TDirectory passed to {}.", __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TMatrixDSym* foundMatrix = nullptr;
  int matrixCount = 0;

  TIter next(TempFile->GetListOfKeys());
  TKey* key;
  while ((key = dynamic_cast<TKey*>(next()))) {
    std::string className = key->GetClassName();
    if (className.find("TMatrix") != std::string::npos) {
      ++matrixCount;
      if (matrixCount == 1) {
        foundMatrix = dynamic_cast<TMatrixDSym*>(key->ReadObj());
      }
    }
  }

  if (matrixCount == 1 && foundMatrix) {
    MACH3LOG_INFO("Found single TMatrixDSym in directory: using it.");
    return foundMatrix;
  } else {
    MACH3LOG_WARN("Found {} TMatrixDSym objects. Using hardcoded path: xsec_cov.", matrixCount);
    TMatrixDSym* fallback = TempFile->Get<TMatrixDSym>("xsec_cov");
    if (!fallback) {
      MACH3LOG_WARN("Fallback matrix 'xsec_cov' not found.");
    }
    return fallback;
  }
}

// *************************************
/// @brief Computes Cholesky decomposition of a symmetric positive definite matrix using custom function which can be even 20 times faster
/// @param matrix Input symmetric positive definite matrix
/// @param matrixName Identifier for error reporting
inline std::vector<std::vector<double>> GetCholeskyDecomposedMatrix(const TMatrixDSym& matrix, const std::string& matrixName) {
// *************************************
  const Int_t n = matrix.GetNrows();
  std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

  for (Int_t j = 0; j < n; ++j) {
    // Compute diagonal element (must be serial)
    double sum_diag = matrix(j, j);
    for (Int_t k = 0; k < j; ++k) {
      sum_diag -= L[j][k] * L[j][k];
    }
    const double tol = 1e-15;
    if (sum_diag <= tol) {
      MACH3LOG_ERROR("Cholesky decomposition failed for {} (non-positive diagonal)", matrixName);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    L[j][j] = std::sqrt(sum_diag);

    // Compute the rest of the column in parallel
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (Int_t i = j + 1; i < n; ++i) {
      double sum = matrix(i, j);
      for (Int_t k = 0; k < j; ++k) {
        sum -= L[i][k] * L[j][k];
      }
      L[i][j] = sum / L[j][j];
    }
  }
  return L;
}

// *************************************
/// @brief Checks if a matrix can be Cholesky decomposed
/// @param matrix Input symmetric matrix to test
inline bool CanDecomposeMatrix(const TMatrixDSym& matrix) {
// *************************************
  TDecompChol chdcmp(matrix);
  return chdcmp.Decompose();
}

// *************************************
/// @brief Makes sure that matrix is positive-definite by adding a small number to on-diagonal elements
inline void MakeMatrixPosDef(TMatrixDSym *cov) {
// *************************************
  //DB Save original warning state and then increase it in this function to suppress 'matrix not positive definite' messages
  //Means we no longer need to overload
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  //DB Loop 1000 times adding 1e-9 which tops out at 1e-6 shift on the diagonal before throwing error
  constexpr int MaxAttempts = 1e5;
  const int matrixSize = cov->GetNrows();
  int iAttempt = 0;
  bool CanDecomp = false;

  for (iAttempt = 0; iAttempt < MaxAttempts; iAttempt++) {
    if (CanDecomposeMatrix(*cov)) {
      CanDecomp = true;
      break;
    } else {
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif
      for (int iVar = 0 ; iVar < matrixSize; iVar++) {
        (*cov)(iVar, iVar) += 1e-9;
      }
    }
  }

  if (!CanDecomp) {
    MACH3LOG_ERROR("Tried {} times to shift diagonal but still can not decompose the matrix", MaxAttempts);
    MACH3LOG_ERROR("This indicates that something is wrong with the input matrix");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  //DB Resetting warning level
  gErrorIgnoreLevel = originalErrorWarning;
}



// ********************************************
/// @brief Dump Matrix to ROOT file, useful when we need to pass matrix info to another fitting group
/// @warning This is mostly used for backward compatibility
inline void DumpParamHandlerToFile(const int _fNumPar,
                                   const std::vector<double>& _fPreFitValue,
                                   const std::vector<double>& _fError,
                                   const std::vector<double>& _fLowBound,
                                   const std::vector<double>& _fUpBound,
                                   const std::vector<double>& _fIndivStepScale,
                                   const std::vector<std::string>& _fFancyNames,
                                   const std::vector<bool>& _fFlatPrior,
                                   const std::vector<SplineParameter>& SplineParams,
                                   TMatrixDSym* covMatrix,
                                   TH2D* CorrMatrix,
                                   const std::string& Name) {
// ********************************************
  TFile* outputFile = new TFile(Name.c_str(), "RECREATE");

  TObjArray* param_names = new TObjArray();
  TObjArray* spline_interpolation = new TObjArray();
  TObjArray* spline_names = new TObjArray();

  TVectorD* param_prior = new TVectorD(_fNumPar);
  TVectorD* flat_prior = new TVectorD(_fNumPar);
  TVectorD* stepscale = new TVectorD(_fNumPar);
  TVectorD* param_lb = new TVectorD(_fNumPar);
  TVectorD* param_ub = new TVectorD(_fNumPar);

  TVectorD* param_knot_weight_lb = new TVectorD(_fNumPar);
  TVectorD* param_knot_weight_ub = new TVectorD(_fNumPar);
  TVectorD* error = new TVectorD(_fNumPar);

  for(int i = 0; i < _fNumPar; ++i)
  {
    TObjString* nameObj = new TObjString(_fFancyNames[i].c_str());
    param_names->AddLast(nameObj);

    TObjString* splineType = new TObjString("TSpline3");
    spline_interpolation->AddLast(splineType);

    TObjString* splineName = new TObjString("");
    spline_names->AddLast(splineName);

    (*param_prior)[i] = _fPreFitValue[i];
    (*flat_prior)[i] = _fFlatPrior[i];
    (*stepscale)[i] = _fIndivStepScale[i];
    (*error)[i] = _fError[i];

    (*param_lb)[i] = _fLowBound[i];
    (*param_ub)[i] = _fUpBound[i];

    //Default values
    (*param_knot_weight_lb)[i] = -9999;
    (*param_knot_weight_ub)[i] = +9999;
  }

  for (size_t SplineIndex = 0; SplineIndex < SplineParams.size(); SplineIndex++ ) {
    auto SystIndex = SplineParams[SplineIndex].index;

    (*param_knot_weight_lb)[SystIndex] = SplineParams.at(SplineIndex)._SplineKnotLowBound;
    (*param_knot_weight_ub)[SystIndex] = SplineParams.at(SplineIndex)._SplineKnotUpBound;

    TObjString* splineType = new TObjString(SplineInterpolation_ToString(SplineParams.at(SplineIndex)._SplineInterpolationType).c_str());
    spline_interpolation->AddAt(splineType, SystIndex);

    TObjString* splineName = new TObjString(SplineParams[SplineIndex]._fSplineNames.c_str());
    spline_names->AddAt(splineName, SystIndex);
  }
  param_names->Write("xsec_param_names", TObject::kSingleKey);
  delete param_names;
  spline_interpolation->Write("xsec_spline_interpolation", TObject::kSingleKey);
  delete spline_interpolation;
  spline_names->Write("xsec_spline_names", TObject::kSingleKey);
  delete spline_names;

  param_prior->Write("xsec_param_prior");
  delete param_prior;
  flat_prior->Write("xsec_flat_prior");
  delete flat_prior;
  stepscale->Write("xsec_stepscale");
  delete stepscale;
  param_lb->Write("xsec_param_lb");
  delete param_lb;
  param_ub->Write("xsec_param_ub");
  delete param_ub;

  param_knot_weight_lb->Write("xsec_param_knot_weight_lb");
  delete param_knot_weight_lb;
  param_knot_weight_ub->Write("xsec_param_knot_weight_ub");
  delete param_knot_weight_ub;
  error->Write("xsec_error");
  delete error;

  covMatrix->Write("xsec_cov");
  CorrMatrix->Write("hcov");

  outputFile->Close();
  delete outputFile;
}

// ********************************************
/// @brief KS: Let's dump all useful matrices to properly validate PCA
//KS: Let's dump all useful matrices to properly validate PCA
inline void DebugPCA(const double sum,
                     const TMatrixD& temp,
                     const TMatrixD& eigen_vectors,
                     const TVectorD& eigen_values,
                     const TMatrixD& TransferMat,
                     const TMatrixD& TransferMatT,
                     const int NumPar,
                     const int FirstPCAdpar,
                     const int LastPCAdpar,
                     const int nKeptPCApars,
                     const double eigen_threshold) {
// ********************************************
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wfloat-conversion"
  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  TDirectory *ogdir = gDirectory;

  TFile *PCA_Debug = new TFile("Debug_PCA.root", "RECREATE");
  PCA_Debug->cd();

  bool PlotText = true;
  //KS: If we have more than 200 plot becomes unreadable :(
  if(NumPar > 200) PlotText = false;

  auto heigen_values     = std::make_unique<TH1D>("eigen_values", "Eigen Values", eigen_values.GetNrows(), 0.0, eigen_values.GetNrows());
  heigen_values->SetDirectory(nullptr);
  auto heigen_cumulative = std::make_unique<TH1D>("heigen_cumulative", "heigen_cumulative", eigen_values.GetNrows(), 0.0, eigen_values.GetNrows());
  heigen_cumulative->SetDirectory(nullptr);
  auto heigen_frac       = std::make_unique<TH1D>("heigen_fractional", "heigen_fractional", eigen_values.GetNrows(), 0.0, eigen_values.GetNrows());
  heigen_frac->SetDirectory(nullptr);
  heigen_values->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative = 0;
  for(int i = 0; i < eigen_values.GetNrows(); i++)
  {
    heigen_values->SetBinContent(i+1, (eigen_values)(i));
    heigen_cumulative->SetBinContent(i+1, (eigen_values)(i)/sum + Cumulative);
    heigen_frac->SetBinContent(i+1, (eigen_values)(i)/sum);
    Cumulative += (eigen_values)(i)/sum;
  }
  heigen_values->Write("heigen_values");
  eigen_values.Write("eigen_values");
  heigen_cumulative->Write("heigen_values_cumulative");
  heigen_frac->Write("heigen_values_frac");

  TH2D* heigen_vectors = new TH2D(eigen_vectors);
  heigen_vectors->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors->Write("heigen_vectors");
  eigen_vectors.Write("eigen_vectors");

  TH2D* SubsetPCA = new TH2D(temp);
  SubsetPCA->GetXaxis()->SetTitle("Parameter in Normal Base");
  SubsetPCA->GetYaxis()->SetTitle("Parameter in Decomposed Base");

  SubsetPCA->Write("hSubsetPCA");
  temp.Write("SubsetPCA");
  TH2D* hTransferMat = new TH2D(TransferMat);
  hTransferMat->GetXaxis()->SetTitle("Parameter in Normal Base");
  hTransferMat->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  TH2D* hTransferMatT = new TH2D(TransferMatT);

  hTransferMatT->GetXaxis()->SetTitle("Parameter in Decomposed Base");
  hTransferMatT->GetYaxis()->SetTitle("Parameter in Normal Base");

  hTransferMat->Write("hTransferMat");
  TransferMat.Write("TransferMat");
  hTransferMatT->Write("hTransferMatT");
  TransferMatT.Write("TransferMatT");

  auto c1 = std::make_unique<TCanvas>("c1", " ", 0, 0, 1024, 1024);
  c1->SetBottomMargin(0.1);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetGrid();

  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // Make pretty correlation colors (red to blue)
  constexpr int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  double maxz = 0;
  double minz = 0;

  c1->Print("Debug_PCA.pdf[");
  auto EigenLine = std::make_unique<TLine>(nKeptPCApars, 0, nKeptPCApars, heigen_cumulative->GetMaximum());
  EigenLine->SetLineColor(kPink);
  EigenLine->SetLineWidth(2);
  EigenLine->SetLineStyle(kSolid);

  auto text = std::make_unique<TText>(0.5, 0.5, Form("Threshold = %g", eigen_threshold));
  text->SetTextFont (43);
  text->SetTextSize (40);

  heigen_values->SetLineColor(kRed);
  heigen_values->SetLineWidth(2);
  heigen_cumulative->SetLineColor(kGreen);
  heigen_cumulative->SetLineWidth(2);
  heigen_frac->SetLineColor(kBlue);
  heigen_frac->SetLineWidth(2);

  c1->SetLogy();
  heigen_values->SetMaximum(heigen_cumulative->GetMaximum()+heigen_cumulative->GetMaximum()*0.4);
  heigen_values->Draw();
  heigen_frac->Draw("SAME");
  heigen_cumulative->Draw("SAME");
  EigenLine->Draw("Same");
  text->DrawTextNDC(0.42, 0.84,Form("Threshold = %g", eigen_threshold));

  auto leg = std::make_unique<TLegend>(0.2, 0.2, 0.6, 0.5);
  leg->SetTextSize(0.04);
  leg->AddEntry(heigen_values.get(), "Absolute", "l");
  leg->AddEntry(heigen_frac.get(), "Fractional", "l");
  leg->AddEntry(heigen_cumulative.get(), "Cumulative", "l");

  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("Same");

  c1->Print("Debug_PCA.pdf");
  c1->SetRightMargin(0.15);
  c1->SetLogy(0);

  heigen_vectors->SetMarkerSize(0.2);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");

  auto Eigen_Line = std::make_unique<TLine>(0, nKeptPCApars, LastPCAdpar - FirstPCAdpar, nKeptPCApars);
  Eigen_Line->SetLineColor(kGreen);
  Eigen_Line->SetLineWidth(2);
  Eigen_Line->SetLineStyle(kDotted);
  Eigen_Line->Draw("SAME");
  c1->Print("Debug_PCA.pdf");

  SubsetPCA->SetMarkerSize(0.2);
  minz = SubsetPCA->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) SubsetPCA->Draw("COLZ TEXT");
  else SubsetPCA->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete SubsetPCA;

  hTransferMat->SetMarkerSize(0.15);
  minz = hTransferMat->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMat->Draw("COLZ TEXT");
  else hTransferMat->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete hTransferMat;

  hTransferMatT->SetMarkerSize(0.15);
  minz = hTransferMatT->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMatT->Draw("COLZ TEXT");
  else hTransferMatT->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete hTransferMatT;

  delete heigen_vectors;

  c1->Print("Debug_PCA.pdf]");
  PCA_Debug->Close();
  delete PCA_Debug;
  gErrorIgnoreLevel = originalErrorWarning;
  ogdir->cd();  // go back to original directory
  #pragma GCC diagnostic pop
}

} // end M3 namespace
