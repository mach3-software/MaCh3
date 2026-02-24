#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "Samples/SampleHandlerFD.h"

namespace py = pybind11;


// Add these includes at the top of the file
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "TH1.h"
#include "TH2.h"

namespace py = pybind11;

// Helper function to convert TH1 to numpy arrays
std::tuple<py::array_t<double>, py::array_t<double>> TH1ToNumpy(TH1* hist) {
    if (!hist) {
        throw std::runtime_error("Histogram pointer is null");
    }
    
    int nbins = hist->GetNbinsX();
    
    // Create numpy array for bin contents
    py::array_t<double> contents(nbins);
    auto contents_buf = contents.request();
    double* contents_ptr = static_cast<double*>(contents_buf.ptr);
    
    // Create numpy array for bin edges (nbins + 1 edges)
    py::array_t<double> edges(nbins + 1);
    auto edges_buf = edges.request();
    double* edges_ptr = static_cast<double*>(edges_buf.ptr);
    
    // Copy bin contents (ROOT bins start at 1, not 0)
    for (int i = 0; i < nbins; ++i) {
        std::cout << "Bin countent for bin " << i + 1 << ": " << hist->GetBinContent(i + 1) << std::endl;
        contents_ptr[i] = hist->GetBinContent(i + 1);
    }
    
    // Copy bin edges
    for (int i = 0; i <= nbins; ++i) {
        edges_ptr[i] = hist->GetBinLowEdge(i + 1);
    }
    // Add the upper edge of the last bin
    edges_ptr[nbins] = hist->GetBinLowEdge(nbins + 1) + hist->GetBinWidth(nbins + 1);
    
    return std::make_tuple(contents, edges);
}

// Helper function to convert TH2 to numpy arrays
std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> TH2ToNumpy(TH2* hist) {
    if (!hist) {
        throw std::runtime_error("Histogram pointer is null");
    }
    
    int nbinsX = hist->GetNbinsX();
    int nbinsY = hist->GetNbinsY();
    
    // Create 2D numpy array for bin contents (shape: nbinsY x nbinsX to match numpy convention)
    py::array_t<double> contents({nbinsY, nbinsX});
    auto contents_buf = contents.request();
    double* contents_ptr = static_cast<double*>(contents_buf.ptr);
    
    // Create numpy arrays for bin edges
    py::array_t<double> edgesX(nbinsX + 1);
    auto edgesX_buf = edgesX.request();
    double* edgesX_ptr = static_cast<double*>(edgesX_buf.ptr);
    
    py::array_t<double> edgesY(nbinsY + 1);
    auto edgesY_buf = edgesY.request();
    double* edgesY_ptr = static_cast<double*>(edgesY_buf.ptr);
    
    // Copy bin contents (ROOT bins start at 1, not 0)
    // Note: numpy uses row-major order (C-style), so we iterate Y then X
    for (int iy = 0; iy < nbinsY; ++iy) {
        for (int ix = 0; ix < nbinsX; ++ix) {
            contents_ptr[iy * nbinsX + ix] = hist->GetBinContent(ix + 1, iy + 1);
        }
    }
    
    // Copy X bin edges
    for (int i = 0; i <= nbinsX; ++i) {
        edgesX_ptr[i] = hist->GetXaxis()->GetBinLowEdge(i + 1);
    }
    edgesX_ptr[nbinsX] = hist->GetXaxis()->GetBinLowEdge(nbinsX + 1) + 
                         hist->GetXaxis()->GetBinWidth(nbinsX + 1);
    
    // Copy Y bin edges
    for (int i = 0; i <= nbinsY; ++i) {
        edgesY_ptr[i] = hist->GetYaxis()->GetBinLowEdge(i + 1);
    }
    edgesY_ptr[nbinsY] = hist->GetYaxis()->GetBinLowEdge(nbinsY + 1) + 
                         hist->GetYaxis()->GetBinWidth(nbinsY + 1);
    
    return std::make_tuple(contents, edgesX, edgesY);
}

// Add these bindings to the PySampleHandlerFD class definition:

/// @brief EW: As SampleHandlerBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerBase : public SampleHandlerBase {
public:
    /* Inherit the constructors */
    using SampleHandlerBase::SampleHandlerBase;

    /* Trampoline (need one for each virtual function) */
    std::string GetName() const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetName, /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::string GetSampleTitle(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetSampleTitle,       /* Name of function in C++ (must match Python name) */
            iSample               /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    int GetNOscChannels(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            int,                  /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetNOscChannels,      /* Name of function in C++ (must match Python name) */
            iSample               /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Reweight() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,              /* Return type */
            SampleHandlerBase, /* Parent class */
            "reweight",
            Reweight           /* Name of function in C++ (must match Python name) */
        );
    }


    /* Trampoline (need one for each virtual function) */
    double GetSampleLikelihood(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                /* Return type */
            SampleHandlerBase,     /* Parent class */
            "get_sample_likelihood",
            GetSampleLikelihood,   /* Name of function in C++ (must match Python name) */
            iSample                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                  /* Return type */
            SampleHandlerBase,     /* Parent class */
            "clean_memory_before_fit",
            CleanMemoryBeforeFit   /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void PrintRates(const bool DataOnly = false) override  {
        PYBIND11_OVERRIDE_PURE(
            void,                 /* Return type */
            SampleHandlerBase,    /* Parent class */
            PrintRates,           /* Name of function in C++ (must match Python name) */
            DataOnly              /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::string GetKinVarName(const int iSample, const int Dimension) const override  {
        PYBIND11_OVERRIDE_PURE(
            std::string,          /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetKinVarName,        /* Name of function in C++ (must match Python name) */
            iSample,              /* Argument(s) */
            Dimension             /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const override  {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,  /* Return type */
            SampleHandlerBase,    /* Parent class */
            GetKinVarName,        /* Name of function in C++ (must match Python name) */
            Sample,               /* Argument(s) */
            KinematicParameter    /* Argument(s) */
        );
    }

    TH1* GetDataHist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                  /* Return type */
            SampleHandlerBase,     /* Parent class */
            GetDataHist,           /* Name of function in C++ (must match Python name) */
            Sample                 /* Argument(s) */
        );
    }

    TH1* GetMCHist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                  /* Return type */
            SampleHandlerBase,     /* Parent class */
            GetMCHist,             /* Name of function in C++ (must match Python name) */
            Sample                 /* Argument(s) */
        );
    }

    TH1* GetW2Hist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                  /* Return type */
            SampleHandlerBase,     /* Parent class */
            GetW2Hist,             /* Name of function in C++ (must match Python name) */
            Sample                 /* Argument(s) */
        );
    }

    double GetLikelihood() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,            /* Return type */
            SampleHandlerBase, /* Parent class */
            "get_likelihood",  /* Python name*/
            GetLikelihood      /* Name of function in C++ (must match Python name) */
                               /* Argument(s) */
        );
    }

    TH1* Get1DVarHistByModeAndChannel(const int iSample,
                                      const std::string& ProjectionVar_Str,
                                      int kModeToFill = -1,
                                      int kChannelToFill = -1,
                                      int WeightStyle = 0,
                                      TAxis* Axis = nullptr) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                     /* Return type */
            SampleHandlerBase,        /* Parent class */
            Get1DVarHistByModeAndChannel, /* Name of function in C++ */
            iSample,
            ProjectionVar_Str,
            kModeToFill,
            kChannelToFill,
            WeightStyle,
            Axis
        );
    }

    TH2* Get2DVarHistByModeAndChannel(const int iSample,
                                      const std::string& ProjectionVar_StrX,
                                      const std::string& ProjectionVar_StrY,
                                      int kModeToFill = -1,
                                      int kChannelToFill = -1,
                                      int WeightStyle = 0,
                                      TAxis* AxisX = nullptr,
                                      TAxis* AxisY = nullptr) override {
        PYBIND11_OVERRIDE_PURE(
            TH2*,                        /* Return type */
            SampleHandlerBase,           /* Parent class */
            Get2DVarHistByModeAndChannel, /* Name of function in C++ */
            iSample,
            ProjectionVar_StrX,
            ProjectionVar_StrY,
            kModeToFill,
            kChannelToFill,
            WeightStyle,
            AxisX,
            AxisY
        );
    }

    TH1* Get1DVarHist(const int iSample,
                      const std::string &ProjectionVar,
                      const std::vector<KinematicCut> &EventSelectionVec = {},
                      int WeightStyle = 0,
                      TAxis *Axis = nullptr,
                      const std::vector<KinematicCut> &SubEventSelectionVec = {}) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                     /* Return type */
            SampleHandlerBase,        /* Parent class */
            Get1DVarHist,             /* Name of function in C++ */
            iSample,
            ProjectionVar,
            EventSelectionVec,
            WeightStyle,
            Axis,
            SubEventSelectionVec
        );
    }

    TH2* Get2DVarHist(const int iSample,
                      const std::string& ProjectionVarX,
                      const std::string& ProjectionVarY,
                      const std::vector<KinematicCut>& EventSelectionVec = {},
                      int WeightStyle = 0,
                      TAxis* AxisX = nullptr,
                      TAxis* AxisY = nullptr,
                      const std::vector<KinematicCut>& SubEventSelectionVec = {}) override {
        PYBIND11_OVERRIDE_PURE(
            TH2*,                     /* Return type */
            SampleHandlerBase,        /* Parent class */
            Get2DVarHist,             /* Name of function in C++ */
            iSample,
            ProjectionVarX,
            ProjectionVarY,
            EventSelectionVec,
            WeightStyle,
            AxisX,
            AxisY,
            SubEventSelectionVec
        );
    }

    int GetNDim(const int Sample) const override {
        PYBIND11_OVERRIDE_PURE(
            int,                      /* Return type */
            SampleHandlerBase,        /* Parent class */
            GetNDim,                  /* Name of function in C++ */
            Sample
        );
    }

    std::string GetFlavourName(const int iSample,
                               const int iChannel) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,              /* Return type */
            SampleHandlerBase,        /* Parent class */
            GetFlavourName,           /* Name of function in C++ */
            iSample,
            iChannel
        );
    }
};


/// @brief As SampleHandlerFD is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerFD : public SampleHandlerFD {
public:
    /* Inherit the constructors */
    using SampleHandlerFD::SampleHandlerFD;

    /* Trampoline (need one for each virtual function) */
    void AddAdditionalWeightPointers() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                             /* Return type */
            SampleHandlerFD,                  /* Parent class */
            "add_additional_weight_pointers", /*python name*/
            AddAdditionalWeightPointers,      /* Name of function in C++ */
                                              /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            SampleHandlerFD, /* Parent class */
            CleanMemoryBeforeFit       /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupSplines() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_splines", /*python name*/
            SetupSplines,     /* Name of function in C++ */
                             /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Init() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "init",          /*python name*/
            Init,            /* Name of function in C++ */
                             /* Argument(s) */
        );
    }
    
    /* Trampoline (need one for each virtual function) */
    int SetupExperimentMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,             /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_experiment_MC", /*python name*/
            SetupExperimentMC,     /* Name of function in C++ */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupFDMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,            /* Return type */
            SampleHandlerFD, /* Parent class */
            "setup_FD_MC",   /*python name*/
            SetupFDMC,       /* Name of function in C++ */
        );
    }

    int ReturnKinematicParameterFromString(std::string) {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_kinematic_by_name",    /* python name*/
            ReturnKinematicParameterFromString, /* Name of function in C++ (must match Python name) */
            py::arg("variable_name")
        );
    }
    
    std::string ReturnStringFromKinematicParameter(int) {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::string,                /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_kinematic_name",       /* python name*/
            ReturnStringFromKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable_id")
        );
    }

    double ReturnKinematicParameter(std::string, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }

    double ReturnKinematicParameter(int, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                     /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value",/* python name*/
            ReturnKinematicParameter,  /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }

    const double *GetPointerToKinematicParameter(std::string, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }
    const double *GetPointerToKinematicParameter(double, int) override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                   /* Return type */
            SampleHandlerFD,            /* Parent class */
            "get_event_kinematic_value_reference",/* python name*/
            GetPointerToKinematicParameter, /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")            /* Argument(s) */
        );
    }
    
    void RegisterFunctionalParameters() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,
            SampleHandlerFD,
            "register_functional_parameters",
            RegisterFunctionalParameters
        );
    }
};

void initSamples(py::module &m){

    auto m_samples = m.def_submodule("samples");

    m_samples.doc() =
        "This is a Python binding of MaCh3s C++ based samples library.";

    // Bind the systematic type enum that lets us set different types of systematics
    py::enum_<TestStatistic>(m_samples, "TestStatistic")
        .value("Poisson", TestStatistic::kPoisson)
        .value("Barlow_Beeston", TestStatistic::kBarlowBeeston)
        .value("Ice_Cube", TestStatistic::kIceCube)
        .value("Pearson", TestStatistic::kPearson)
        .value("Dembinski_Abdelmottele", TestStatistic::kDembinskiAbdelmotteleb)
        .value("N_Test_Statistics", TestStatistic::kNTestStatistics);

    py::class_<SampleHandlerBase, PySampleHandlerBase /* <--- trampoline*/>(m_samples, "SampleHandlerBase")
        .def(py::init())
        
        .def(
            "reweight", 
            &SampleHandlerBase::Reweight,
            "reweight the MC events in this sample. You will need to override this."
        )
        
        .def(
            "get_likelihood", 
            &SampleHandlerBase::GetLikelihood,
            "Get the sample likelihood at the current point in your model space. You will need to override this."
        )
        
        .def(
            "set_test_stat",
            &SampleHandlerBase::SetTestStatistic,
            "Set the test statistic that should be used when calculating likelihoods. \n\
            :param test_stat: The new test statistic to use",
            py::arg("test_stat")
        )

        .def(
            "get_bin_LLH",
            py::overload_cast<double, double, double>(&SampleHandlerBase::GetTestStatLLH, py::const_),
            "Get the LLH for a bin by comparing the data and MC. The result depends on having previously set the test statistic using :py:meth:`pyMaCh3.samples.SampleHandlerBase.set_test_stat` \n\
            :param data: The data content of the bin. \n\
            :param mc: The mc content of the bin \n\
            :param w2: The Sum(w_{i}^2) (sum of weights squared) in the bin, which is sigma^2_{MC stats}",
            py::arg("data"), 
            py::arg("mc"), 
            py::arg("w2")
        )
    ; // End of SampleHandlerBase binding

    py::class_<SampleHandlerFD, PySampleHandlerFD /* <--- trampoline*/,
               SampleHandlerBase>(m_samples, "SampleHandlerFD")
        .def(
            py::init<std::string, ParameterHandlerGeneric *>(),
            "This should never be called directly as SampleHandlerFD is an abstract base class. \n\
            However when creating a derived class, in the __init__() method, you should call the parent constructor i.e. this one by doing:: \n\
            \n\
            \tsuper(<your derived SampleHandler class>, self).__init__(*args) \n\
            \n ",
            py::arg("mc_version"), py::arg("xsec_cov"))

        .def(
            "get_mc_hist",
            [](SampleHandlerFD &self, const int Dimension) {

              //self.Reweight();
              
              // Get the histogram pointer BEFORE cloning
              TH1 *hist_original = self.GetMCHist(Dimension);
              
              // Debug: Check the original histogram
              std::cout << "=== ORIGINAL HISTOGRAM ===" << std::endl;
              if (!hist_original) {
                throw std::runtime_error("GetMCHist returned null pointer");
              }
              std::cout << "Dimension is " << Dimension << std::endl;
              std::cout << "Original histogram address: " << hist_original << std::endl;
              std::cout << "Original histogram name: " << hist_original->GetName() << std::endl;
              std::cout << "Original histogram bins: " << hist_original->GetNbinsX() << std::endl;
              std::cout << "Original histogram integral: " << hist_original->Integral() << std::endl;
              std::cout << "Original histogram entries: " << hist_original->GetEntries() << std::endl;
              
              // Print first few bin contents
              std::cout << "First 5 bin contents of original:" << std::endl;
              for (int i = 1; i <= std::min(5, hist_original->GetNbinsX()); ++i) {
                std::cout << "  Bin " << i << ": " << hist_original->GetBinContent(i) << std::endl;
              }
              
              // Now clone it
              TH1D *hist = (TH1D*)hist_original->Clone("cloned_hist");
              
              // Debug: Check the cloned histogram
              std::cout << "=== CLONED HISTOGRAM ===" << std::endl;
              std::cout << "Cloned histogram address: " << hist << std::endl;
              std::cout << "Cloned histogram name: " << hist->GetName() << std::endl;
              std::cout << "Cloned histogram integral: " << hist->Integral() << std::endl;
              std::cout << "Cloned histogram entries: " << hist->GetEntries() << std::endl;
              
              // Print first few bin contents of clone
              std::cout << "First 5 bin contents of clone:" << std::endl;
              for (int i = 1; i <= std::min(5, hist->GetNbinsX()); ++i) {
                std::cout << "  Bin " << i << ": " << hist->GetBinContent(i) << std::endl;
              }

              if (Dimension == 1) {
                // 1D histogram
                auto [contents, edgesX] = TH1ToNumpy(hist);
                auto edgesY = py::array_t<double>();
                return py::make_tuple(contents, edgesX, edgesY);
              } else if (Dimension == 2) {
                // 2D histogram - cast to TH2
                TH2 *hist2d = dynamic_cast<TH2 *>(hist);
                if (!hist2d) {
                  throw std::runtime_error("Failed to cast to TH2");
                }
                auto [contents, edgesX, edgesY] = TH2ToNumpy(hist2d);
                return py::make_tuple(contents, edgesX, edgesY);
              } else {
                throw std::invalid_argument("Dimension must be 1 or 2");
              }
            },
            py::return_value_policy::reference_internal,
            py::arg("Dimension"),
            "Get MC histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")

        .def(
            "get_data_hist",
            [](SampleHandlerFD &self, const int Dimension) {
              TH1 *hist = self.GetDataHist(Dimension);

              if (Dimension == 1) {
                // 1D histogram
                auto [contents, edgesX] = TH1ToNumpy(hist);
                auto edgesY = py::array_t<double>();
                return py::make_tuple(contents, edgesX, edgesY);
              } else if (Dimension == 2) {
                // 2D histogram - cast to TH2
                TH2 *hist2d = dynamic_cast<TH2 *>(hist);
                if (!hist2d) {
                  throw std::runtime_error("Failed to cast to TH2");
                }
                auto [contents, edgesX, edgesY] = TH2ToNumpy(hist2d);
                return py::make_tuple(contents, edgesX, edgesY);
              } else {
                throw std::invalid_argument("Dimension must be 1 or 2");
              }
            },
            py::arg("Dimension"),
            "Get Data histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")

        .def(
            "get_w2_hist",
            [](SampleHandlerFD &self, const int Dimension) {
              TH1 *hist = self.GetW2Hist(Dimension);

              if (Dimension == 1) {
                // 1D histogram
                auto [contents, edgesX] = TH1ToNumpy(hist);
                auto edgesY = py::array_t<double>();
                return py::make_tuple(contents, edgesX, edgesY);
              } else if (Dimension == 2) {
                // 2D histogram - cast to TH2
                TH2 *hist2d = dynamic_cast<TH2 *>(hist);
                if (!hist2d) {
                  throw std::runtime_error("Failed to cast to TH2");
                }
                auto [contents, edgesX, edgesY] = TH2ToNumpy(hist2d);
                return py::make_tuple(contents, edgesX, edgesY);
              } else {
                throw std::invalid_argument("Dimension must be 1 or 2");
              }
            },
            py::arg("Dimension"),
            "Get W2 histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D");

    /* Not sure if this will be needed in future versions of MaCh3 so leaving commented for now
    py::class_<fdmc_base>(m_samples, "MCstruct")
        .def(py::init())
        
        // Because a lot of the variables in fdmc_base use c style arrays,
        // we need to provide some setter functions to be able to set them using more
        // "pythony" objects, e.g. lists and numpy arrays
        .def(
            "set_event_variable_values", 
            [](fdmc_base &self, int dim, py::array_t<double, py::array::c_style> &array)
            {
                py::buffer_info bufInfo = array.request();

                if ( dim > 2 )
                    throw MaCh3Exception(__FILE__, __LINE__, "Currently only dimensions of 1 or 2 are supported sorry :(");

                if ( bufInfo.ndim != 1 )
                    throw MaCh3Exception(__FILE__, __LINE__, "Number of dimensions in parameter array must be one if setting only one of the event variable arrays!");
                
                if( dim ==1 )
                    self.x_var = array.data();
                    
                else if ( dim == 2)
                    self.y_var = array.data();
            }
        )
    ;
    */
}
