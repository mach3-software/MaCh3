/// @file UmbrellaSolver.cpp
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author David Riley

#include <filesystem>
#include "Fitters/Processing/MCMCProcessor.h"
#include "Manager/Manager.h"
#include "Fitters/Algorithms/MulticanonicalMCMCHandler.h"

_MaCh3_Safe_Include_Start_ //{
#include "TSystem.h"
#include "TChain.h"
#include "TSystemDirectory.h"
_MaCh3_Safe_Include_End_ //}

bool debug_mode = false;
/// Structure to hold each window configuration
/// These are individual window settings, function type & params plus filename
struct WindowConfig {
  std::string name;
  M3::BiasFunction umbrellaBiasFunction;
  double center;
  double width;
  std::string input_file;
  double vonMises_kappa;
};

/// Structure to hold umbrella configuration
/// This is the global settings of the umbrella fit eg number of windows and solve settings
struct UmbrellaConfig {
  std::vector<WindowConfig> windows;
  std::string output_file;
  std::string variable_of_interest;
  bool dynamic_files;
  std::string dynamic_pattern;
  int dynamic_n_windows;
  int max_iterations;
  double tolerance;
  int print_frequency;
  bool use_openmp;
};

/// YAML-based config parser using yaml-cpp library
UmbrellaConfig parseYAMLConfig(const std::string &filename) {
  UmbrellaConfig config;

  try {
    YAML::Node yaml_diag_config = YAML::LoadFile(filename);
    YAML::Node yaml_config = yaml_diag_config["UmbrellaSolver"];

    // Parse other configuration
    config.output_file = Get<std::string>(yaml_config["output_file"], __FILE__ , __LINE__);
    config.variable_of_interest = Get<std::string>(yaml_config["variable_of_interest"], __FILE__ , __LINE__);
    config.max_iterations = Get<int>(yaml_config["max_iterations"], __FILE__ , __LINE__);
    config.tolerance = Get<double>(yaml_config["tolerance"], __FILE__ , __LINE__);
    config.print_frequency = GetFromManager<int>(yaml_config["print_frequency"], 0, __FILE__ , __LINE__);
    config.dynamic_files = GetFromManager<bool>(yaml_config["dynamic_files"], false, __FILE__ , __LINE__);
    config.dynamic_pattern = Get<std::string>(yaml_config["dynamic_pattern"], __FILE__ , __LINE__);
    config.dynamic_n_windows = Get<int>(yaml_config["dynamic_n_windows"], __FILE__ , __LINE__);
    config.use_openmp = GetFromManager<bool>(yaml_config["use_openmp"], true, __FILE__ , __LINE__);

    if (!config.dynamic_files) {
      // TODO: these are mostly redundant now that the MaCh3_Config macro is
      // being parsed, work out if its safe to remove them at somepoint!!! Parse
      // windows
      if (yaml_config["windows"]) {
        const YAML::Node &windows = yaml_config["windows"];
        for (size_t i = 0; i < windows.size(); i++) {
          WindowConfig window;
          window.name = Get<std::string>(windows[i]["name"], __FILE__, __LINE__);
          window.center = Get<double>(windows[i]["center"], __FILE__, __LINE__);
          window.width = Get<double>(windows[i]["width"], __FILE__, __LINE__);
          config.windows.push_back(window);
        }
      }

      // Parse input files
      if (yaml_config["input_files"]) {
        const YAML::Node &input_files = yaml_config["input_files"];
        for (size_t i = 0; i < input_files.size() && i < config.windows.size(); i++) {
          config.windows[i].input_file = Get<std::string>(input_files[i], __FILE__, __LINE__);
        }
      }
    } else {
      // set up placeholders to be filled later
      for (int i = 0; i < config.dynamic_n_windows; i++) {
        WindowConfig window;
        window.name = "Window_" + std::to_string(i);
        window.center =
            0.0; // Placeholder, will be updated from MaCh3_Config macro
        window.width =
            1.0; // Placeholder, will be updated from MaCh3_Config macro
        config.windows.push_back(window);
      }
    }
  } catch (const YAML::Exception &e) {
    MACH3LOG_ERROR("Error parsing YAML file {}: {}", filename, e.what());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return config;
}

// Gaussian window function
double gaussianWindow(double x, double center, double width) {
  return exp(-0.5 * pow((x - center) / width, 2)) / (width * sqrt(2 * TMath::Pi()));
}

// von Mises window function (circular analog of Gaussian)
double vonMisesWindow(double x, double center, double kappa) {
  // von Mises PDF: exp(kappa * cos(x - center)) / (2*pi*I0(kappa))
  // For numerical stability, compute in log space when possible
  double I0_kappa;
  if (kappa > 700) {
    // Use asymptotic approximation for large kappa to avoid overflow
    // log(I0(kappa)) ≈ kappa - 0.5*log(2*pi*kappa)
    double log_I0 = kappa - 0.5 * log(2 * TMath::Pi() * kappa);
    return exp(kappa * cos(x - center) - log_I0 - log(2 * TMath::Pi()));
  } else {
    I0_kappa = TMath::BesselI0(kappa);
    return exp(kappa * cos(x - center)) / (2 * TMath::Pi() * I0_kappa);
  }
}

// Generalised Gaussian window allows harsher constraints on the tails of the distribution, better for controlling windows in highly disfavoured regions
double generalisedGaussian2(double x, double mean, double width) {
  constexpr int n = 2; // this controls the tightness of the gaussian fixed at 2
                       // for now due to normalisation
  const double normFactor = 1 / ((M3::UmbrellaGaussianNormFactor) * 2 * std::sqrt(2) * width); // the normalisation is a little ugly (uses gamma functions),
                   // im just going to hardcode them for now
  double likelihood = normFactor * std::exp(-std::pow((std::pow(x - mean, 2) / (2 * std::pow(width, 2))), n));
  return likelihood;
}

/// @todo modify this to instead use atan2 implementation for wrapping
double GetMulticanonicalWeightGenGaussian(double deltacp, double mean, double width) {
  // implementation of the generalised gaussian as a bias function
  // for now with a fixed n = 2 for simplicity

  double g0 = generalisedGaussian2(deltacp, mean, width);
  double g1 = generalisedGaussian2(deltacp, mean - 2 * TMath::Pi(),width); // these two repeats are required for wrapping the gaussian around -+pi
  double g2 = generalisedGaussian2(deltacp, mean + 2 * TMath::Pi(), width);
  double multicanonicalBeta = 1.0;
  return (g0 + g1 + g2) * (multicanonicalBeta);
}

// A sub calculation for the overlap matrix
// Sum of all windows weighted by z values
double summedWindowsWeighted(double x, const std::vector<WindowConfig> &windows, const std::vector<double> &z_values) {
  double sum = 0.0;
  for (size_t k = 0; k < windows.size(); k++) {
    double window_val;
    if (windows[k].umbrellaBiasFunction == M3::BiasFunction::kVonMises) {
      window_val = vonMisesWindow(x, windows[k].center, windows[k].vonMises_kappa);
    } else if (windows[k].umbrellaBiasFunction == M3::BiasFunction::kGaussian) {
      window_val = gaussianWindow(x, windows[k].center, windows[k].width);
    } else if (windows[k].umbrellaBiasFunction == M3::BiasFunction::kGeneralisedGaussian) {
      window_val = GetMulticanonicalWeightGenGaussian(x, windows[k].center, windows[k].width);
    } else {
      MACH3LOG_ERROR("Unrecognised BiasFunction!!");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    sum += window_val / z_values[k];
  }
  return sum;
}

// Precompute all window evaluations once: cache[i][j][s] = window_j evaluated
// at samples[i][s]
// Memory heavy depending on number of steps/cores/windows
std::vector<std::vector<std::vector<double>>> buildWindowCache(const std::vector<WindowConfig> &windows,
                                                               const std::vector<std::vector<double>> &samples, bool use_openmp = true) {
  int n_windows = static_cast<int>(windows.size());
  std::vector<std::vector<std::vector<double>>> cache(n_windows);

  for (int i = 0; i < n_windows; i++) {
    cache[i].resize(n_windows);
    for (int j = 0; j < n_windows; j++) {
      cache[i][j].resize(samples[i].size());
    }
  }

  if (use_openmp) {
    #ifdef MULTITHREAD
    #pragma omp parallel for collapse(2) schedule(dynamic)
    #endif
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {
        for (size_t s = 0; s < samples[i].size(); s++) {
          if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kVonMises) {
            cache[i][j][s] = vonMisesWindow(samples[i][s], windows[j].center, windows[j].vonMises_kappa);
          } else if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kGaussian) {
            cache[i][j][s] = gaussianWindow(samples[i][s], windows[j].center, windows[j].width);
          } else if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kGeneralisedGaussian) {
            cache[i][j][s] = GetMulticanonicalWeightGenGaussian(samples[i][s], windows[j].center, windows[j].width);
          } else {
            MACH3LOG_ERROR("Unrecognised BiasFunction!!");
            throw MaCh3Exception(__FILE__, __LINE__);
          }
        }
      }
    }
  } else {
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {
        for (size_t s = 0; s < samples[i].size(); s++) {
          if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kVonMises) {
            cache[i][j][s] = vonMisesWindow(samples[i][s], windows[j].center, windows[j].vonMises_kappa);
          } else if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kGaussian) {
            cache[i][j][s] = gaussianWindow(samples[i][s], windows[j].center, windows[j].width);
          } else if (windows[j].umbrellaBiasFunction == M3::BiasFunction::kGeneralisedGaussian) {
            cache[i][j][s] = GetMulticanonicalWeightGenGaussian(samples[i][s], windows[j].center, windows[j].width);
          } else {
            MACH3LOG_ERROR("Unrecognised function!!!!!");
            throw MaCh3Exception(__FILE__, __LINE__);
          }
        }
      }
    }
  }

  if (debug_mode) { // this should be TH1D of cache against sample value in the
                    // variable of interest for each window, to check the cache
                    // is being built correctly and the windows look correct
                    // across the range of samples
    // TODO break this into a standalone function
    // first print the size in memory of the cache to check it is reasonable and
    // not being built incorrectly with an extra dimension or something
    size_t cache_size_bytes = 0;
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {
        cache_size_bytes += cache[i][j].size() * sizeof(double);
      }
    }
    MACH3LOG_INFO("Window cache size: {:.2f} MB", static_cast<double>(cache_size_bytes) / (1024.0 * 1024.0));
    MACH3LOG_INFO("Window cache built successfully.");
    TFile *cache_file = TFile::Open("window_cache_debug_histograms.root", "RECREATE");
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {
        std::string hist_name = "Sample" + std::to_string(i) + "_window" + std::to_string(j);
        TH1D *window_cache_hist = new TH1D(hist_name.c_str(), hist_name.c_str(), 100, -3.1415, 3.1415);
        for (size_t s = 0; s < cache[i][j].size(); s++) {
          window_cache_hist->AddBinContent(window_cache_hist->FindBin(samples[i][s]));
        }
        window_cache_hist->Write();
      }
    }

    // additionally save a TH2D of the values in the cache against the sample
    // values for a specific window to check the shape of the windows is correct
    // across the range of samples
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {
        std::string hist_name = "Sample" + std::to_string(i) + "_window" + std::to_string(j) + "_cache_2D";
        TH2D *window_cache_2D_hist = new TH2D(hist_name.c_str(), hist_name.c_str(), 100, -3.1415, 3.1415, 100, -1000, 10);
        for (size_t s = 0; s < cache[i][j].size(); s++) {
          // take the log of the likelihood values
          window_cache_2D_hist->Fill(samples[i][s], log(cache[i][j][s]));
        }
        // show overflow numbers on the hist
        window_cache_2D_hist->SetStats(1);
        window_cache_2D_hist->Write();
      }
    }

    cache_file->Close();
  }

  return cache;
}

// Main function to calculate the F matrix from a given set of samples plus weights
std::vector<std::vector<double>> calcFmatrix(std::vector<double> &z_current,
            const std::vector<WindowConfig> &windows,
            const std::vector<std::vector<double>> &samples,
            const std::vector<std::vector<std::vector<double>>> &window_cache) {
  int n_windows = static_cast<int>(windows.size());
  std::vector<std::vector<double>> F(n_windows, std::vector<double>(n_windows, 0.0));

  std::vector<double> z_inv = z_current; // make a copy to avoid modifying the original z_current
  for (size_t i = 0; i < z_current.size(); i++) {
    if (z_current[i] > 0) {
      z_inv[i] = 1.0 / z_current[i];
    } else {
      z_inv[i] = 0.0; // Handle zero or negative z values gracefully
      if (debug_mode) {
        MACH3LOG_WARN("Warning: z_current[{}] is non-positive ({}). Setting its inverse to 0 in F matrix calculation.", i, z_current[i]);
      }
    }
  }
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int i = 0; i < n_windows; i++) {
    std::vector<double> denominator_cache(samples[i].size(), 0.0);
    for (size_t s = 0; s < samples[i].size(); s++) {
      double denominator = 0.0;
      for (int k = 0; k < n_windows; k++) {
        denominator += window_cache[i][k][s] * z_inv[k];
      }
      denominator_cache[s] = 1 / denominator;
    }

    for (int j = 0; j < n_windows; j++) {
      double sum = 0.0;
      int count = 0;

      for (size_t s = 0; s < samples[i].size(); s++) {
        double sample = samples[i][s];
        double window_j = window_cache[i][j][s];
        double denominator = denominator_cache[s];

        if (denominator > 0) {
          double integrand = (window_j * z_inv[i]) * denominator;
          sum += integrand;
          count++;
        } else if (debug_mode) {
          MACH3LOG_WARN("Denominator is zero for sample {} in window {}, skipping...", sample, i);
        }
      }

      if (debug_mode) {
        MACH3LOG_INFO("F[{}][{}] sum: {}, count: {}", i, j, sum, count);
      }

      if (count > 0) {
        F[i][j] = sum / count;
      }
    }
  }

  return F;
}

// Z-solver function implementing the fixed point matrix iteration algorithm
std::vector<double> zSolver(const std::vector<double> &z_current,
        const std::vector<WindowConfig> &windows,
        const std::vector<std::vector<double>> &samples,
        const std::vector<std::vector<std::vector<double>>> &window_cache,
        bool use_openmp = true, bool verbose = false,
        [[maybe_unused]] int *total_lines = nullptr) {
  int n_windows = static_cast<int>(windows.size());
  if (verbose && !use_openmp) {
    MACH3LOG_INFO("Using single-threaded computation for F matrix...");
  }

  // F matrix and update z values
  std::vector<double> z_working = z_current;
  std::vector<std::vector<double>> F =
      calcFmatrix(z_working, windows, samples, window_cache);

  // if (verbose) {
  //     if (total_lines) *total_lines = 1; // Start counting from F matrix
  //     header std::cout << "F matrix:" << std::endl; for (int i = 0; i <
  //     n_windows; i++) {
  //         if (total_lines) (*total_lines)++; // Count each row
  //         std::cout << "[";
  //         for (int j = 0; j < n_windows; j++) {
  //             std::cout << std::setw(10) << std::fixed <<
  //             std::setprecision(5) << F[i][j]; if (j < n_windows - 1)
  //             std::cout << ", ";
  //         }
  //         std::cout << "]" << std::endl;
  //     }
  //     std::cout << std::flush;
  // }

  // Compute z_new = z_current * F
  // could do faster matrix multiplication here?
  std::vector<double> z_new(n_windows, 0.0);
  for (int i = 0; i < n_windows; i++) {
    for (int j = 0; j < n_windows; j++) {
      z_new[i] += z_current[j] * F[j][i];
    }
  }

  /// Normalize z to prevent scale drift and maintain sum(z)=1
  /// @todo: try this with magnitude = 1, the absolute scale doesn't matter but
  /// maybee this condition is preventing elements meeting their required
  /// values? also just try completely disabled
  /// double z_sum = 0.0;
  /// for (int i = 0; i < n_windows; i++) {
  ///    z_sum += z_new[i];
  ///}
  /// if (z_sum > 0) {
  ///   for (int i = 0; i < n_windows; i++) {
  ///        z_new[i] /= z_sum;
  ///    }
  //}

  // normalise the z values so that magnitude of the vector is 1
  double z_magnitude = 0.0;
  for (int i = 0; i < n_windows; i++) {
    z_magnitude += z_new[i] * z_new[i];
  }
  z_magnitude = sqrt(z_magnitude);
  if (z_magnitude > 0) {
    for (int i = 0; i < n_windows; i++) {
      z_new[i] /= z_magnitude;
    }
  }

  return z_new;
}

/// A few different convergence checks
// Check convergence
bool checkConvergence(const std::vector<double> &z_current, const std::vector<double> &z_prev, double tolerance) {
  double sum_diffs = 0.0;
  for (size_t i = 0; i < z_current.size(); i++) {
    // if (std::abs(z_current[i] - z_prev[i]) > tolerance *
    // std::max(std::abs(z_current[i]), std::abs(z_prev[i]))) { // intention
    // here: is the difference between any of the z elements more than the
    // tolerance*the magnitude of the largest element in z? in this way the
    // tolerance sets the number of decimal points below the largest element we
    // are targeting.
    //     return false;
    // }
    sum_diffs += std::abs(z_current[i] - z_prev[i]);
    // check the average difference across all elements instead
    // of requiring every element to meet the condition, this
    // should be more robust to individual elements fluctuating
    // around their target value while the overall z
    // distribution is still converging
    if (sum_diffs / static_cast<double>(z_current.size()) > tolerance) {
      return false;
    }
  }
  return true;
}

std::vector<double> getZDiffs(const std::vector<double> &z_current, const std::vector<double> &z_prev) {
  std::vector<double> diffs(z_current.size(), 0.0);
  for (size_t i = 0; i < z_current.size(); i++) {
    diffs[i] = std::abs(z_current[i] - z_prev[i]);
  }
  return diffs;
}

// This implements a check on stalling of the evolution of the window weights. Seems to be a better definition for convergence that above
// moving-average stalled-convergence check on z values
bool checkConvergenceStalled(const std::vector<double> &z_current, const std::vector<double> &z_prev, double tolerance) {
  (void)z_prev;

  static std::deque<std::vector<double>> z_history;
  static std::vector<double> previous_moving_average;
  static int stagnant_iterations = 0;

  constexpr int moving_average_window = 500;
  constexpr int stagnant_required = 500;
  const double bound = tolerance;

  if (z_current.empty()) {
    return false;
  }

  // Reset state safely if number of windows changes between solver runs.
  if (!z_history.empty() && z_history.front().size() != z_current.size()) {
    z_history.clear();
    previous_moving_average.clear();
    stagnant_iterations = 0;
  }

  z_history.push_back(z_current);
  if (static_cast<int>(z_history.size()) > moving_average_window) {
    z_history.pop_front();
  }

  // Wait until the moving-average window is fully populated.
  if (static_cast<int>(z_history.size()) < moving_average_window) {
    return false;
  }

  std::vector<double> moving_average(z_current.size(), 0.0);
  for (const auto &z_vec : z_history) {
    for (size_t i = 0; i < z_vec.size(); i++) {
      moving_average[i] += z_vec[i];
    }
  }
  for (size_t i = 0; i < moving_average.size(); i++) {
    moving_average[i] /= moving_average_window;
  }

  if (previous_moving_average.empty()) {
    previous_moving_average = moving_average;
    return false;
  }

  bool all_within_bound = true;
  for (size_t i = 0; i < moving_average.size(); i++) {
    if (std::abs(moving_average[i] - previous_moving_average[i]) > bound) {
      all_within_bound = false;
      break;
    }
  }

  if (all_within_bound) {
    stagnant_iterations++;
  } else {
    stagnant_iterations = 0;
  }

  previous_moving_average = moving_average;

  if (stagnant_iterations == stagnant_required) {
    MACH3LOG_WARN("Convergence appears stalled: moving-average change stayed within {} for {} iterations.",
                  bound, stagnant_required);
  }

  return stagnant_iterations >= stagnant_required;
}

// Main function to run the umbrella sampling solver
void UmbrellaSolver(const std::string &config_file) {
  MACH3LOG_INFO("=== Umbrella Sampling Z-Factor Solver ===");
  // Debug OpenMP status first
  MACH3LOG_INFO("Debugging OpenMP availability...");

  #ifdef MULTITHREAD
  MACH3LOG_INFO("Max threads available: {}", omp_get_max_threads());
  #else
  MACH3LOG_WARN("_OPENMP is NOT defined - OpenMP not available");
  #endif
  MACH3LOG_INFO("Loading configuration from: {}", config_file);

  // Parse configuration
  UmbrellaConfig config = parseYAMLConfig(config_file);

  if (config.windows.empty() && !config.dynamic_files) {
    MACH3LOG_ERROR("No windows defined in configuration and dynamic file loading is disabled.");
    return;
  }

  MACH3LOG_INFO("Variable of interest: {}", config.variable_of_interest);
  MACH3LOG_INFO("Output file: {}", config.output_file);

// Check OpenMP status with detailed debugging
#ifdef MULTITHREAD
  MACH3LOG_INFO("OpenMP: AVAILABLE");
  if (config.use_openmp) {
    MACH3LOG_INFO("OpenMP: ENABLED (using {} threads)", omp_get_max_threads());
  }
#else
  MACH3LOG_WARN("OpenMP: NOT AVAILABLE");
  if (config.use_openmp) {
    MACH3LOG_WARN("OpenMP: NOT AVAILABLE - falling back to single-threaded execution");
    MACH3LOG_INFO("Note: For OpenMP support, try compiling with: g++ -fopenmp ...");
    MACH3LOG_INFO("Or ensure OpenMP library is properly loaded in ROOT");
    config.use_openmp = false;
  } else {
    MACH3LOG_INFO("OpenMP: DISABLED (single-threaded execution)");
  }
#endif

  // Load data from input files
  std::vector<std::vector<double>> samples; // Declare here to ensure it exists in the scope of the entire
               // function
  if (!config.dynamic_files) {
    samples.resize(config.windows.size());
  } else {
    samples.resize(config.dynamic_n_windows);
  }
  std::vector<TFile *> input_files;
  std::vector<TTree *> input_trees;

  if (!config.dynamic_files) {
    MACH3LOG_INFO("Using static input files from configuration.");
    for (size_t i = 0; i < config.windows.size(); i++) {
      MACH3LOG_INFO("Loading file: {}", config.windows[i].input_file);
      TFile *file = M3::Open(config.windows[i].input_file.c_str(), "READ", __FILE__, __LINE__);

      TTree *tree = static_cast<TTree*>(file->Get("posteriors"));
      if (!tree) {
        MACH3LOG_ERROR("Cannot find 'posteriors' tree in {}", config.windows[i].input_file);
        file->Close();
        continue;
      }

      input_files.push_back(file);
      input_trees.push_back(tree);
    }
  } else {
    MACH3LOG_INFO("Dynamic file loading enabled. Searching for files in directory: {}", config.dynamic_pattern);
    // Use ROOT's TSystem to find all files in the directory
    TSystemDirectory dir("", config.dynamic_pattern.c_str());
    TList *files = dir.GetListOfFiles();
    int file_count = 0;
    if (files) {
      TIter next(files);
      TSystemFile *file;
      while ((file = static_cast<TSystemFile*>(next()))) {
        std::string filename = file->GetName();
        if (filename.find(".root") == std::string::npos) {
          MACH3LOG_INFO("Skipping non-root file: {}", filename);
          continue;
        }

        std::string full_path = config.dynamic_pattern + "/" + filename;
        MACH3LOG_INFO("Found file: {}", full_path);

        TFile *root_file = M3::Open(full_path, "READ", __FILE__, __LINE__);

        TTree *tree = static_cast<TTree*>(root_file->Get("posteriors"));
        if (!tree) {
          MACH3LOG_ERROR("Cannot find 'posteriors' tree in {}", full_path);
          root_file->Close();
          throw MaCh3Exception(__FILE__, __LINE__, "Missing 'posteriors' tree in file: " + full_path);
        }
        file_count++;

        input_files.push_back(root_file);
        input_trees.push_back(tree);
        MACH3LOG_INFO("Loaded tree 'posteriors' from file: {}", full_path);
      }
    } else {
      MACH3LOG_ERROR("No files found matching pattern: {}", config.dynamic_pattern);
      throw MaCh3Exception(__FILE__, __LINE__, "No files found matching pattern: " + config.dynamic_pattern);
    }

    if (file_count != config.dynamic_n_windows) {
      MACH3LOG_ERROR("Number of files found ({}) does not match expected dynamic_n_windows ({}).",
                     file_count, config.dynamic_n_windows);
      throw MaCh3Exception(__FILE__, __LINE__, "File count mismatch for dynamic loading.");
    }
  }

  for (size_t i = 0; i < input_trees.size(); i++) {
    TTree *tree = input_trees[i];
    TFile *file = input_files[i];

    MACH3LOG_INFO("Processing file {}/{}: {}", i + 1, input_trees.size(), file->GetName());
    TMacro *macro = static_cast<TMacro*>(file->Get("MaCh3_Config"));
    if (macro) {
      MACH3LOG_INFO("Found MaCh3_Config macro in file.");

      // Convert TMacro to YAML: concatenate all lines and parse
      std::stringstream yaml_text;
      TList *lines = macro->GetListOfLines();
      for (int iline = 0; iline < lines->GetEntries(); ++iline) {
        TObjString *line = static_cast<TObjString*>(lines->At(iline));
        yaml_text << line->GetString().Data() << "\n";
      }

      try {
        YAML::Node macro_yaml = YAML::Load(yaml_text.str());
        YAML::Node umbrellaConfig =
            macro_yaml["General"]["MCMC"]["Multicanonical"];
        // Extract window parameters
        config.windows[i].center = Get<double>(umbrellaConfig["Umbrella"]["UmbrellaMean"], __FILE__, __LINE__);
        MACH3LOG_INFO("Window {} center updated to {}", i, config.windows[i].center);

        // Check if using von Mises distribution
        auto biasString = Get<std::string>(umbrellaConfig["Umbrella"]["UmbrellaBiasFunction"], __FILE__, __LINE__);
        M3::BiasFunction biasMode;
        if (biasString == "gaussian") {
          biasMode = M3::BiasFunction::kGaussian;
          MACH3LOG_INFO("Window weighted with gaussian");
        } else if (biasString == "generalisedGaussian") {
          biasMode = M3::BiasFunction::kGeneralisedGaussian;
          MACH3LOG_INFO("Window weighted with generalised gaussian");
        } else if (biasString == "vonMises") {
          biasMode = M3::BiasFunction::kVonMises;
          MACH3LOG_INFO("Window weighted with vonMises");
        } else {
          MACH3LOG_ERROR("Unrecognised Bias");
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        config.windows[i].umbrellaBiasFunction = biasMode;

        if (config.windows[i].umbrellaBiasFunction == M3::BiasFunction::kVonMises) {
          // Extract von Mises sigma and compute kappa
          double vonMises_sigma = Get<double>(umbrellaConfig["Umbrella"]["UmbrellaWidth"], __FILE__, __LINE__);
          config.windows[i].vonMises_kappa = 1.0 / (vonMises_sigma * vonMises_sigma);
          config.windows[i].width = vonMises_sigma; // Store sigma in width for reference
          MACH3LOG_INFO("Window {} using von Mises: sigma = {}, kappa = {}",
                        i, vonMises_sigma, config.windows[i].vonMises_kappa);
        } else {
          // Extract Gaussian sigma
          config.windows[i].width = Get<double>(umbrellaConfig["Umbrella"]["UmbrellaWidth"], __FILE__, __LINE__);
          config.windows[i].vonMises_kappa = -1.0; // Not using von Mises
          MACH3LOG_INFO("Window {} using Gaussian: width = {}", i, config.windows[i].width);
        }
      } catch (const std::exception &e) {
        MACH3LOG_WARN("Could not parse macro as YAML: {}", e.what());
      }
    }

    double var_value;
    double logL_value;
    tree->SetBranchAddress(config.variable_of_interest.c_str(), &var_value);
    tree->SetBranchAddress("LogL", &logL_value);

    Long64_t nentries = tree->GetEntries();
    //Long64_t filtered_entries = 0;
    MACH3LOG_INFO("Window {}: {} entries", i, nentries);

    for (Long64_t entry = 0; entry < nentries; entry++) {
      tree->GetEntry(entry);
      //if (logL_value > 50.0) { // logl cut no longer needed as the posterior
      //                         // chain start has been fixed
      //  filtered_entries++;
      //  continue;
      //}
      samples[i].push_back(var_value);
    }
    //if (filtered_entries > 0) {
    //  std::cout << "Filtered " << filtered_entries
    //            << " entries with LogL > 500 from window " << i << std::endl;
    //}
  }

  // Sort by window center and keep all per-window containers aligned.
  // The final weighting stage loops over input_trees by index, so those indices
  // must track the same sorted window order used by z_current.
  for (size_t i = 0; i < config.windows.size(); i++) {
    for (size_t j = i + 1; j < config.windows.size(); j++) {
      if (config.windows[i].center > config.windows[j].center) {
        std::swap(config.windows[i], config.windows[j]);
        std::swap(samples[i], samples[j]);
        std::swap(input_trees[i], input_trees[j]);
        std::swap(input_files[i], input_files[j]);
      }
    }
  }

  // verify the order and file associations after sorting
  MACH3LOG_INFO("Final window configurations after sorting:");
  for (size_t i = 0; i < config.windows.size(); i++) {
    MACH3LOG_INFO("Window {}: center = {}, width = {}, vonMises_mode = {}, vonMises_kappa = {}, samples = {}",
      i, config.windows[i].center, config.windows[i].width,
      (config.windows[i].umbrellaBiasFunction == M3::BiasFunction::kVonMises ? "Yes" : "No"),
      config.windows[i].vonMises_kappa, samples[i].size());
  }

  // Initialize z values
  std::vector<double> z_current(config.windows.size(), 1.0);
  std::vector<double> z_prev(config.windows.size(), 1.0);

  // this should be used to pick up a solve that failed partway through due to reaching iteration max or job cancellation
  bool hacky_start = false; // Set to true to use the hardcoded starting vector, false to start with all ones
  if (hacky_start) {
    // z_current = {   0.02593,    0.02676,    0.02764,    0.03305,    0.03527,
    // 0.04275,    0.05171,    0.04808,    0.04978,    0.04560,    0.04627,
    // 0.05187,    0.04993,    0.04656,    0.04690,    0.04954,    0.04940,
    // 0.04315,    0.03616,    0.02832,    0.01998,    0.01482,    0.01167,
    // 0.00902,    0.00646,    0.00474,    0.00439,    0.00421,    0.00470,
    // 0.00553,    0.00626,    0.00822,    0.01070,    0.01477,    0.01750,
    // 0.02235};
    z_current = {
        0.023968629123202176,   0.024927005713178161,   0.026030888791054529,
        0.036203405237770721,   0.04004944137212621,    0.055993616350153479,
        0.079251266929094608,   0.065860139904686643,   0.067944181615205768,
        0.055237804689517243,   0.054778778031073304,   0.066196102148964917,
        0.059298667596959342,   0.049864361722134341,   0.048890393249559315,
        0.05284144211204099,    0.050606183191239794,   0.037329936801427918,
        0.025159187577401387,   0.01487697249782439,    0.0074395236463036573,
        0.0040535024095969992,  0.0025347709967512371,  0.0015590961074484638,
        0.00083136243905723782, 0.00047599269828616062, 0.0004454414833514952,
        0.00045205403592812914, 0.00062098732584897902, 0.00091722187629238541,
        0.0012407301891732216,  0.0023494590898020503,  0.0041640866720417773,
        0.0080599626292932776,  0.011429554834689368,   0.018117848911520615};
    z_prev = z_current; // Start with the same values for previous to avoid
                        // large initial changes
    MACH3LOG_WARN("!!!!!!!starting from hacky start vector!!!!!!");
  }

  std::vector<std::vector<double>> z_evolution;

  MACH3LOG_INFO("Starting iterative z-solver...");
  // Test OpenMP functionality once before the main loop this can probably be wrapped in debug TODO
  bool openmp_works = false;
  if (config.use_openmp) {
    MACH3LOG_INFO("Testing OpenMP parallelization...");
    #ifdef MULTITHREAD
    int max_threads = omp_get_max_threads();
    #else
    int max_threads = 1;
    #endif

    MACH3LOG_INFO("Max threads reported: {}", max_threads);

    // Test parallel region
    int actual_threads = 1;
    #ifdef MULTITHREAD
    #pragma omp parallel
    #endif
    {
      #ifdef MULTITHREAD
      #pragma omp master
      #endif
      {
        #ifdef MULTITHREAD
        actual_threads = omp_get_num_threads();
        #else
        actual_threads = 1;
        #endif
        MACH3LOG_INFO("Actual threads in parallel region: {}", actual_threads);
      }
    }

    if (actual_threads > 1) {
      MACH3LOG_INFO("OpenMP is working correctly with {} threads", actual_threads);
      openmp_works = true;
    } else {
      openmp_works = false;
    }
  }

  MACH3LOG_INFO("Precomputing window cache...");
  // TODO: this is slow as hell and causes massive memory usage, what's a smarter
  // way to do this? break out to a file? RDataFrame?
  std::vector<std::vector<std::vector<double>>> window_cache = buildWindowCache(config.windows, samples, openmp_works);

  // TFile to hold the F matrix evolution for the first 15 iterations if needed
  bool save_matrix = true; // Set to true to enable saving F matrix evolution
  // Create F_file here to avoid reopening it multiple times in the loop
  TFile *F_file = nullptr;
  if (save_matrix) {
    size_t pos = config.output_file.find(".root");
    std::string base_name = (pos != std::string::npos) ? config.output_file.substr(0, pos) : config.output_file;
    F_file = TFile::Open((base_name + "_matrix_evolution.root").c_str(), "RECREATE"); // use name from config with .root subtracted with _matrix_evolution suffix added
    if (!F_file || F_file->IsZombie()) {
      MACH3LOG_ERROR("Cannot create file {}", base_name + "_matrix_evolution.root");
      save_matrix = false; // Disable saving if file cannot be created
    }
    // add an initial FMatrix with the initial z values for reference
    std::vector<std::vector<double>> initial_F = calcFmatrix(z_current, config.windows, samples, window_cache);
    int n_windows = static_cast<int>(config.windows.size());
    TH2D initial_F_TH2D("F_matrix_initial", "Initial F matrix;Window j;Window i", n_windows, 0, n_windows, n_windows, 0, n_windows);
    for (int i = 0; i < n_windows; i++) {
      for (int j = 0; j < n_windows; j++) {initial_F_TH2D.SetBinContent(j + 1, i + 1, initial_F[i][j]); // Note the order of i and j for correct axis labeling
      }
    }
    F_file->cd();
    MACH3LOG_INFO("Saving initial F matrix to file...");
    initial_F_TH2D.Write();
  }

  // Timing variables
  auto start_time = std::chrono::high_resolution_clock::now();
  auto last_print_time = start_time;

  bool converged_robustness_check = false; // Flag to indicate if convergence check has been passed at least
             // once, used to control when to start checking for stalled
             // convergence

  if (!converged_robustness_check) {
    MACH3LOG_INFO("Starting iterative solver with convergence checks...");
  }
  TRandom3 gRandom3;
  // Iterative solver
  int total_output_lines = 0; // Track total lines printed for clearing
  for (int iteration = 0; iteration < config.max_iterations; iteration++) {
    if (iteration % config.print_frequency == 0 || iteration == 1) {
      auto current_time = std::chrono::high_resolution_clock::now();

      // Clear previous output if not the first iteration
      if (iteration > 0) {
        // Move cursor up and clear previous output
        std::cout << "\033[" << total_output_lines << "A"; // Move up
        std::cout << "\033[J"; // Clear from cursor to end of screen
      }

      // Calculate average relative change (precision metric)
      double avg_relative_change = 0.0;
      if (iteration > 0) {
        for (size_t i = 0; i < z_current.size(); i++) {
          double rel_change = std::abs(z_current[i] - z_prev[i]) / std::max(std::abs(z_current[i]), 1e-10);
          avg_relative_change += rel_change;
        }
        avg_relative_change /= static_cast<double>(z_current.size());
      }

      // Print z-values
      std::cout << "Iteration " << std::setw(6) << iteration << ", z values: [";
      for (size_t i = 0; i < z_current.size(); i++) {
        std::cout << std::setw(10) << std::fixed << std::setprecision(5) << z_current[i];
        if (i < z_current.size() - 1)
          std::cout << ", ";
      }
      if (iteration > 0) {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - last_print_time);
        double avg_time_per_iteration = static_cast<double>(duration.count()) /
                                        static_cast<double>(config.print_frequency);
        std::cout << "] (avg: " << std::setw(6) << std::setprecision(1) << avg_time_per_iteration << " ms/iter)" << std::endl;
        std::cout << "Avg relative change: " << std::scientific << std::setprecision(3) << avg_relative_change << " (target: " << config.tolerance << ")" << std::endl;
        total_output_lines = 2;
      } else {
        std::cout << "]" << std::endl;
        total_output_lines = 1;
      }
      last_print_time = current_time;
    }

    z_prev = z_current;
    z_current = zSolver(z_current, config.windows, samples, window_cache, openmp_works, iteration % config.print_frequency == 0, &total_output_lines);
    z_evolution.push_back(z_current);

    // save the first 15 iterations of the F matrix to check convergence
    // behaviour and debug if needed should be a root file with Th2D for easy
    // plotting in root, with axes of iteration number and window index, and the
    // value being the F matrix element
    if (save_matrix && (iteration < 15 || iteration % config.print_frequency == 0)) {
      std::vector<std::vector<double>> F_matrix = calcFmatrix(z_current, config.windows, samples, window_cache);
      // convert F_matrix to Th2D for saving to root file
      int n_windows = static_cast<int>(config.windows.size());
      TH2D F_TH2D(Form("F_matrix_iter_%02d", iteration),Form("F matrix at iteration %02d;Window j;Window i", iteration), n_windows, 0, n_windows, n_windows, 0, n_windows);
      for (int i = 0; i < n_windows; i++) {
        for (int j = 0; j < n_windows; j++) {
          F_TH2D.SetBinContent(j + 1, i + 1, F_matrix[i][j]); // Note the order of i and j for correct axis labeling
        }
      }
      // for the purposes of picking back up a solve after it has been
      // interrupted also save the std::vector of z_current. You can put this into hacky start to pick up
      TTree *z_tree = new TTree(Form("z_saved_iter_%02d", iteration), Form("Z vector at iteration %02d", iteration));
      z_tree->Branch("z_saved", &z_current);
      z_tree->Fill();
      F_file->cd();
      MACH3LOG_INFO("Saving F matrix for iteration {} to file...", iteration);
      F_TH2D.Write();
      z_tree->Write();
    }
    // if (save_matrix && iteration == 15) {
    //     F_file->Close();
    // }

    // after the first convergence check has been passed, randomly perturb the z
    // values to see if they return to the same values, this is a robustness
    // check to see if the solution is stable or if it is just meeting the
    // convergence criteria by chance due to small changes in z values
    bool apply_robustness_check = true;
    if (iteration % 100 == 0 &&
        (checkConvergence(z_current, z_prev, config.tolerance) || checkConvergenceStalled(z_current, z_prev, config.tolerance))) {
      if (!converged_robustness_check && apply_robustness_check) {
        MACH3LOG_INFO("Convergence check passed at iteration {}. Starting robustness check with random perturbation...", iteration);
        converged_robustness_check = true;

        // Apply random perturbation to z_current
        std::vector<double> z_perturbed = z_current;
        for (size_t i = 0; i < z_perturbed.size(); i++) {
          // Random perturbation up to 10 times the tolerance
          double perturbation = gRandom3.Uniform(-0.5, 0.5) * z_perturbed[i];
          MACH3LOG_INFO("Applying perturbation of {:.6e} to z[{}] = {:.6e}", perturbation, i, z_perturbed[i]);
          z_perturbed[i] += perturbation;
          if (z_perturbed[i] < 0)
            z_perturbed[i] = abs(z_perturbed[i]); // Ensure no negative values
        }
        z_current = z_perturbed;
          MACH3LOG_INFO("Applied random perturbation to z values for robustness check.");
      } else {
        if (checkConvergence(z_current, z_prev, config.tolerance)) {
          MACH3LOG_INFO("Convergence achieved at iteration {}", iteration);
        } else {
          MACH3LOG_WARN("Convergence appears to be stalled at iteration {}", iteration);
        }

        if (iteration == config.max_iterations - 1) {
          MACH3LOG_WARN("Reached maximum iterations without convergence.");
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double avg_time_total = static_cast<double>(total_duration.count()) /
                                static_cast<double>(iteration + 1);
        if (save_matrix) {
          F_file->Close();
        }
        MACH3LOG_INFO("Terminating at iteration {}", iteration);
        MACH3LOG_INFO("Average time per iteration: {} ms", avg_time_total);
        break;
      }
    }
  }

  std::ostringstream oss;
  oss << "Final z values: [";
  for (size_t i = 0; i < z_current.size(); i++) {
    oss << std::fixed << std::setprecision(5) << z_current[i];
    if (i < z_current.size() - 1)
      oss << ", ";
  }
  oss << "]";
  MACH3LOG_INFO("{}", oss.str());

  std::filesystem::copy_file(input_files[0]->GetName(), config.output_file,
                             std::filesystem::copy_options::overwrite_existing);
  // Create output file
  TFile *output_file = M3::Open(config.output_file.c_str(), "UPDATE", __FILE__, __LINE__);
  output_file->cd();

  TTree *input_tree = dynamic_cast<TTree*>(output_file->Get("posteriors"));
  // Create combined tree with weights
  TTree *combined_tree = input_tree->CloneTree(0);

  // Variables for the combined tree
  double umbrella_weight;
  int window_id;
  double delta_cp;

  combined_tree->Branch("umbrella_weight", &umbrella_weight, "umbrella_weight/D");
  combined_tree->Branch("window_id", &window_id, "window_id/I");

  // Fill combined tree
  for (size_t i = 0; i < input_trees.size(); i++) {
    TTree *tree = input_trees[i];

    Long64_t nentries = tree->GetEntries();

    // KS: This is to avoid warnings about missing umbrella branches...
    int oldLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    combined_tree->CopyAddresses(tree);
    gErrorIgnoreLevel = oldLevel;

    /// @todo code now assumes it is only for delta CP
    tree->SetBranchAddress("delta_cp", &delta_cp);
    // KS: SetBranchAddress above decouples the input branch address, so update the
    // copied output branch address to use the current delta_cp value.
    combined_tree->GetBranch("delta_cp")->SetAddress(&delta_cp);
    if (z_current[i] == 0) {
      MACH3LOG_WARN("Z value for window {} is zero, skipping weighting for this window to avoid division by zero.", i);
    }
    window_id = static_cast<int>(i);

    for (Long64_t entry = 0; entry < nentries; entry++) {
      tree->GetEntry(entry);

      if (z_current[i] == 0) {
        umbrella_weight = 0.0; // If z is zero, we cannot apply the umbrella weight, so we set it to 0 (completely downweigh this window's contribution)
      } else {
        // Calculate umbrella weight for this event
        // The umbrella weight corrects for the bias introduced by the window function Weight is 1 / sum of all window contributions (equation 4 from paper)
        double denominator = 1 / summedWindowsWeighted(delta_cp, config.windows, z_current);

        // umbrella_weight = z_current[i] / denominator; // with or without z_current[i] / denominator? why did I have this originally
        umbrella_weight = denominator; // This is the correct form based on the paper - the z_current[i] factor is already included in the summedWindowsWeighted function
      }

      if (combined_tree->Fill() < 0) {
        MACH3LOG_WARN("Failed writing output tree. Check disk quota/space and write permissions for: {}", config.output_file);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }
  // Write final results
  combined_tree->Write(input_tree->GetName(), TObject::kOverwrite);

  TDirectory* UmbreallaDir = output_file->mkdir("Umbrealla");
  UmbreallaDir->cd();
  // Save diagnostics
  TCanvas c1("c1", "Z Evolution", 800, 600);
  std::vector<TGraph*> z_graphs(config.windows.size());
  TLegend legend(0.7, 0.7, 0.9, 0.9);

  double ymax = 0.0;
  double ymin = std::numeric_limits<double>::max();

  for (size_t i = 0; i < config.windows.size(); i++) {
    std::vector<double> iterations, z_vals;
    for (size_t j = 0; j < z_evolution.size(); j++) {
      iterations.push_back(static_cast<int>(j));
      z_vals.push_back(z_evolution[j][i]);
    }

    for (double val : z_vals) {
      if (val > ymax)
        ymax = val;
      if (val < ymin)
        ymin = val;
    }

    z_graphs[i] = new TGraph(static_cast<int>(iterations.size()), &iterations[0], &z_vals[0]);
    z_graphs[i]->SetLineColor(static_cast<Color_t>(i + 1));
    z_graphs[i]->SetLineWidth(2);
    z_graphs[i]->SetName(Form("z_evolution_window_%lu", i));
    z_graphs[i]->SetTitle("Evolution of Z Values");

    if (i == 0) {
      z_graphs[i]->GetXaxis()->SetTitle("Iteration");
      z_graphs[i]->GetYaxis()->SetTitle("Z Value");
      z_graphs[i]->Draw("AL");
    } else {
      z_graphs[i]->Draw("L SAME");
    }

    legend.AddEntry(z_graphs[i], Form("Window %lu", i), "l");
    z_graphs[i]->Write();
  }

  legend.Draw();
  z_graphs[0]->SetMaximum(ymax);
  z_graphs[0]->SetMinimum(ymin);
  c1.Update();
  c1.SetLogy();
  c1.Write();

  // Create summary histogram of delta_cp distribution
  TH1D *h_delta_cp = new TH1D("h_delta_cp_weighted", "Weighted Delta CP Distribution", 100, -TMath::Pi(), TMath::Pi());
  TH1D *h_delta_cp_unweighted = new TH1D("h_delta_cp_unweighted", "Unweighted Delta CP Distribution", 100, -TMath::Pi(), TMath::Pi());

  combined_tree->Draw("delta_cp>>h_delta_cp_weighted", "umbrella_weight", "goff");
  combined_tree->Draw("delta_cp>>h_delta_cp_unweighted", "", "goff");

  h_delta_cp->Write();
  h_delta_cp_unweighted->Write();

  UmbreallaDir->Close();
  delete UmbreallaDir;

  // Get entry count before closing the file
  Long64_t total_entries = combined_tree->GetEntries();

  output_file->cd();

  YAML::Node yaml_config = M3OpenConfig(config_file);
  YAML::Node umbrella_config;
  umbrella_config["UmbrellaSolver"] = yaml_config["UmbrellaSolver"];

  // Convert YAML -> TMacro
  TMacro UmbrellaHeader = YAMLtoTMacro(umbrella_config, "Umbrella_Config");
  UmbrellaHeader.Write();

  output_file->Close();

  // Close input files
  for (TFile *file : input_files) {
    file->Close();
  }

  MACH3LOG_INFO("Umbrella sampling outputs created (make sure to check for convergence issues)!");
  MACH3LOG_INFO("Output written to: {}", config.output_file);
  MACH3LOG_INFO("Combined tree contains {} entries with umbrella weights.", total_entries);
}

// Main function for compiled version
int main(int argc, char *argv[]) {
  SetMaCh3LoggerFormat();
  std::string config_file = "umbrella_config.yaml";
  if (argc > 1) {
    config_file = argv[1];
  }

  MACH3LOG_INFO("Running compiled version with OpenMP support");
  try {
    UmbrellaSolver(config_file);
  } catch (const std::exception &e) {
    MACH3LOG_ERROR("Error: {}", e.what());
    return 1;
  }

  // Ensure all OpenMP threads are properly terminated
  #ifdef MULTITHREAD
  #pragma omp barrier
  #endif

  return 0;
}
