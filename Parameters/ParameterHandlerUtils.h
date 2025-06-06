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
_MaCh3_Safe_Include_End_ //}

namespace M3
{
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
/// @param VecMulti This vector is used for multiplication VecMulti = matrix x vector
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
          systematicNode["ParameterValues"][Tune] = MaCh3Utils::FormatDouble(Values[i], 4);
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
      systematicNode["ParameterValues"][Tune] = MaCh3Utils::FormatDouble(Values[i], 4);
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
inline void MakeCorrelationMatrix(YAML::Node& root,
                                  const std::vector<double>& Values,
                                  const std::vector<double>& Errors,
                                  const std::vector<std::vector<double>>& Correlation,
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

      syst_i["ParameterValues"]["PreFitValue"] = MaCh3Utils::FormatDouble(Values[i], 4);
      syst_i["Error"] = std::round(Errors[i] * 100.0) / 100.0;

      YAML::Node correlationsNode = YAML::Node(YAML::NodeType::Sequence);
      for (std::size_t j = 0; j < FancyNames.size(); ++j) {
        if (i == j) continue;
        YAML::Node singleEntry;
        singleEntry[FancyNames[j]] = MaCh3Utils::FormatDouble(Correlation[i][j], 4);
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

      syst["ParameterValues"]["PreFitValue"] = MaCh3Utils::FormatDouble(Values[i], 4);
      syst["Error"] = std::round(Errors[i] * 100.0) / 100.0;

      YAML::Node correlationsNode = YAML::Node(YAML::NodeType::Sequence);
      for (std::size_t j = 0; j < Correlation[i].size(); ++j) {
        if (i == j) continue;
        YAML::Node singleEntry;
        const std::string& otherName = systematics[j]["Systematic"]["Names"]["FancyName"].as<std::string>();
        singleEntry[otherName] = MaCh3Utils::FormatDouble(Correlation[i][j], 4);
        correlationsNode.push_back(singleEntry);
      }
      syst["Correlations"] = correlationsNode;
    }
  }

  // Convert and write
  std::string YAMLString = YAMLtoSTRING(NodeCopy);
  FixSampleNamesQuotes(YAMLString);
  std::string OutName = "UpdatedCorrelationMatrix.yaml";
  std::ofstream outFile(OutName);
  if (!outFile) {
    MACH3LOG_ERROR("Failed to open file for writing: {}", OutName);
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

}
