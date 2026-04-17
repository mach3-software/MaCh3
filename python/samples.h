#pragma once

/// @file samples.h
/// @author Ewan Miller

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "Samples/SampleHandlerBase.h"
#include "TH1.h"
#include "TH2.h"

namespace py = pybind11;

// Helper function to convert TH1 to numpy arrays
std::tuple<py::array_t<M3::float_t>, py::array_t<M3::float_t>> TH1ToNumpy(std::unique_ptr<TH1>& hist) {
    if (!hist) {
        throw std::runtime_error("Histogram pointer is null");
    }
    
    int nbins = hist->GetNbinsX();
    
    // Create numpy array for bin contents
    py::array_t<M3::float_t> contents(nbins);
    auto contents_buf = contents.request();
    M3::float_t* contents_ptr = static_cast<M3::float_t*>(contents_buf.ptr);
    
    // Create numpy array for bin edges (nbins + 1 edges)
    py::array_t<M3::float_t> edges(nbins + 1);
    auto edges_buf = edges.request();
    M3::float_t* edges_ptr = static_cast<M3::float_t*>(edges_buf.ptr);
    
    // Copy bin contents (ROOT bins start at 1, not 0)
    for (int i = 0; i < nbins; ++i) {
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
std::tuple<py::array_t<M3::float_t>, py::array_t<M3::float_t>, py::array_t<M3::float_t>> TH2ToNumpy(std::unique_ptr<TH2>& hist) {
    if (!hist) {
        throw std::runtime_error("Histogram pointer is null");
    }
    
    int nbinsX = hist->GetNbinsX();
    int nbinsY = hist->GetNbinsY();
    
    // Create 2D numpy array for bin contents (shape: nbinsY x nbinsX to match numpy convention)
    py::array_t<M3::float_t> contents({nbinsY, nbinsX});
    auto contents_buf = contents.request();
    M3::float_t* contents_ptr = static_cast<M3::float_t*>(contents_buf.ptr);
    
    // Create numpy arrays for bin edges
    py::array_t<M3::float_t> edgesX(nbinsX + 1);
    auto edgesX_buf = edgesX.request();
    M3::float_t* edgesX_ptr = static_cast<M3::float_t*>(edgesX_buf.ptr);
    
    py::array_t<M3::float_t> edgesY(nbinsY + 1);
    auto edgesY_buf = edgesY.request();
    M3::float_t* edgesY_ptr = static_cast<M3::float_t*>(edgesY_buf.ptr);
    
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


inline py::tuple HistToNumpy(std::unique_ptr<TH1>& hist, int dimension)
{
    if (!hist)
        throw std::runtime_error("Null histogram");

    if (dimension == 1)
    {
        auto [c, x] = TH1ToNumpy(hist);
        py::array_t<M3::float_t> y;
        return py::make_tuple(c, x, y);
    }

    if (dimension == 2)
    {
        const TH2* h2 = std::make_unique<TH2>(hist.release());
        if (!h2)
            throw std::runtime_error("Expected TH2");

        auto [c, x, y] = TH2ToNumpy(h2);
        return py::make_tuple(c, x, y);
    }

    throw std::invalid_argument("Only 1D or 2D supported");
}

// Add these bindings to the PySampleHandlerBase class definition:

/// @brief EW: As SampleHandlerBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerInterface : public SampleHandlerInterface {
public:
    /* Inherit the constructors */
    using SampleHandlerInterface::SampleHandlerInterface;

    /* Trampoline (need one for each virtual function) */
    std::string GetName() const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,               /* Return type */
            SampleHandlerInterface,    /* Parent class */
            GetName,                   /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::string GetSampleTitle(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,              /* Return type */
            SampleHandlerInterface,   /* Parent class */
            GetSampleTitle,           /* Name of function in C++ (must match Python name) */
            iSample                   /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    int GetNOscChannels(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE(
            int,                     /* Return type */
            SampleHandlerInterface,  /* Parent class */
            GetNOscChannels,         /* Name of function in C++ (must match Python name) */
            iSample                  /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Reweight() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                   /* Return type */
            SampleHandlerInterface, /* Parent class */
            "reweight",
            Reweight                /* Name of function in C++ (must match Python name) */
        );
    }


    /* Trampoline (need one for each virtual function) */
    double GetSampleLikelihood(const int iSample) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                    /* Return type */
            SampleHandlerInterface,    /* Parent class */
            "get_sample_likelihood",
            GetSampleLikelihood,       /* Name of function in C++ (must match Python name) */
            iSample                    /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                       /* Return type */
            SampleHandlerInterface,     /* Parent class */
            "clean_memory_before_fit",
            CleanMemoryBeforeFit        /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void PrintRates(const bool DataOnly = false) override  {
        PYBIND11_OVERRIDE_PURE(
            void,                   /* Return type */
            SampleHandlerInterface, /* Parent class */
            PrintRates,             /* Name of function in C++ (must match Python name) */
            DataOnly                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::string GetKinVarName(const int iSample, const int Dimension) const override  {
        PYBIND11_OVERRIDE_PURE(
            std::string,             /* Return type */
            SampleHandlerInterface,  /* Parent class */
            GetKinVarName,           /* Name of function in C++ (must match Python name) */
            iSample,                 /* Argument(s) */
            Dimension                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    std::vector<double> ReturnKinematicParameterBinning(const int Sample, const std::string &KinematicParameter) const override  {
        PYBIND11_OVERRIDE_PURE(
            std::vector<double>,     /* Return type */
            SampleHandlerInterface,  /* Parent class */
            GetKinVarName,           /* Name of function in C++ (must match Python name) */
            Sample,                  /* Argument(s) */
            KinematicParameter       /* Argument(s) */
        );
    }

    TH1* GetDataHist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                     /* Return type */
            SampleHandlerInterface,   /* Parent class */
            GetDataHist,              /* Name of function in C++ (must match Python name) */
            Sample                    /* Argument(s) */
        );
    }

    TH1* GetMCHist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                   /* Return type */
            SampleHandlerInterface, /* Parent class */
            GetMCHist,              /* Name of function in C++ (must match Python name) */
            Sample                  /* Argument(s) */
        );
    }

    TH1* GetW2Hist(const int Sample) override {
        PYBIND11_OVERRIDE_PURE(
            TH1*,                   /* Return type */
            SampleHandlerInterface, /* Parent class */
            GetW2Hist,              /* Name of function in C++ (must match Python name) */
            Sample                  /* Argument(s) */
        );
    }

    double GetLikelihood() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                 /* Return type */
            SampleHandlerInterface, /* Parent class */
            "get_likelihood",       /* Python name*/
            GetLikelihood           /* Name of function in C++ (must match Python name) */
                                    /* Argument(s) */
        );
    }

    std::unique_ptr<TH1> Get1DVarHistByModeAndChannel(const int iSample,
                                                      const std::string& ProjectionVar_Str,
                                                      const int kModeToFill = -1,
                                                      const int kChannelToFill = -1,
                                                      const int WeightStyle = 0) override {
        (void) iSample;
        (void) ProjectionVar_Str;
        (void) kModeToFill;
        (void) kChannelToFill;
        (void) WeightStyle;
        return nullptr;
    }

    std::unique_ptr<TH2> Get2DVarHistByModeAndChannel(const int iSample,
                                                      const std::string& ProjectionVar_StrX,
                                                      const std::string& ProjectionVar_StrY,
                                                      const int kModeToFill = -1,
                                                      const int kChannelToFill = -1,
                                                      const int WeightStyle = 0) override {
      (void) iSample;
      (void) ProjectionVar_StrX;
      (void) ProjectionVar_StrY;
      (void) kModeToFill;
      (void) kChannelToFill;
      (void) WeightStyle;
      return nullptr;
    }

    std::unique_ptr<TH1> Get1DVarHist(const int iSample,
                                      const std::string &ProjectionVar,
                                      const std::vector<KinematicCut> &EventSelectionVec = {},
                                      const int WeightStyle = 0,
                                      const std::vector<KinematicCut> &SubEventSelectionVec = {}) override {
      (void) iSample;
      (void) ProjectionVar;
      (void) EventSelectionVec;
      (void) WeightStyle;
      (void) SubEventSelectionVec;
      return nullptr;
    }

    std::unique_ptr<TH2> Get2DVarHist(const int iSample,
                                      const std::string& ProjectionVarX,
                                      const std::string& ProjectionVarY,
                                      const std::vector<KinematicCut>& EventSelectionVec = {},
                                      const int WeightStyle = 0,
                                      const std::vector<KinematicCut>& SubEventSelectionVec = {}) override {
      (void) iSample;
      (void) ProjectionVarX;
      (void) ProjectionVarY;
      (void) EventSelectionVec;
      (void) WeightStyle;
      (void) SubEventSelectionVec;
      return nullptr;
    }

    int GetNDim(const int Sample) const override {
        PYBIND11_OVERRIDE_PURE(
            int,                      /* Return type */
            SampleHandlerInterface,   /* Parent class */
            GetNDim,                  /* Name of function in C++ */
            Sample
        );
    }

    std::string GetFlavourName(const int iSample,
                               const int iChannel) const override {
        PYBIND11_OVERRIDE_PURE(
            std::string,              /* Return type */
            SampleHandlerInterface,   /* Parent class */
            GetFlavourName,           /* Name of function in C++ */
            iSample,
            iChannel
        );
    }
};


/// @brief As SampleHandlerBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySampleHandlerBase : public SampleHandlerBase {
public:
    /* Inherit the constructors */
    using SampleHandlerBase::SampleHandlerBase;

    /* Trampoline (need one for each virtual function) */
    void AddAdditionalWeightPointers() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                               /* Return type */
            SampleHandlerBase,                  /* Parent class */
            "add_additional_weight_pointers",   /*python name*/
            AddAdditionalWeightPointers,        /* Name of function in C++ */
                                                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void CleanMemoryBeforeFit() override {
        PYBIND11_OVERRIDE_PURE(
            void,                   /* Return type */
            SampleHandlerBase,      /* Parent class */
            CleanMemoryBeforeFit    /* Name of function in C++ (must match Python name) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupSplines() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,               /* Return type */
            SampleHandlerBase,  /* Parent class */
            "setup_splines",    /*python name*/
            SetupSplines,       /* Name of function in C++ */
                                /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void Init() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,              /* Return type */
            SampleHandlerBase, /* Parent class */
            "init",            /*python name*/
            Init,              /* Name of function in C++ */
                               /* Argument(s) */
        );
    }
    
    /* Trampoline (need one for each virtual function) */
    int SetupExperimentMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            int,                    /* Return type */
            SampleHandlerBase,      /* Parent class */
            "setup_experiment_MC",  /*python name*/
            SetupExperimentMC,      /* Name of function in C++ */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void InititialiseData() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                /* Return type */
            SampleHandlerBase,   /* Parent class */
            "inititialis_data",  /*python name*/
            InititialiseData,    /* Name of function in C++ */
        );
    }

    /* Trampoline (need one for each virtual function) */
    void SetupMC() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,               /* Return type */
            SampleHandlerBase,  /* Parent class */
            "setup_MC",         /*python name*/
            SetupMC,            /* Name of function in C++ */
        );
    }
    
    double ReturnKinematicParameter(int, int) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            double,                       /* Return type */
            SampleHandlerBase,            /* Parent class */
            "get_event_kinematic_value",  /* python name*/
            ReturnKinematicParameter,     /* Name of function in C++ (must match Python name) */
            py::arg("variable"),
            py::arg("event")              /* Argument(s) */
        );
    }
    
    const double *GetPointerToKinematicParameter(int, int) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            const double *,                           /* Return type */
            SampleHandlerBase,                        /* Parent class */
            "get_event_kinematic_value_reference",    /* python name*/
            GetPointerToKinematicParameter,           /* Name of function in C++ (must match Python name) */
            py::arg("variable"),                      /* Argument(s) */
            py::arg("event")                          /* Argument(s) */
        );
    }

    void RegisterFunctionalParameters() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                               /* Return type */
            SampleHandlerBase,                  /* Parent class */
            "register_functional_parameters",
            RegisterFunctionalParameters        /* Name of function in C++ (must match Python name) */
        );
    }
};

void initSamplesModule(py::module &m_samples){

    // Bind the systematic type enum that lets us set different types of systematics
    py::enum_<TestStatistic>(m_samples, "TestStatistic")
        .value("Poisson", TestStatistic::kPoisson)
        .value("Barlow_Beeston", TestStatistic::kBarlowBeeston)
        .value("Ice_Cube", TestStatistic::kIceCube)
        .value("Pearson", TestStatistic::kPearson)
        .value("Dembinski_Abdelmottele", TestStatistic::kDembinskiAbdelmotteleb)
        .value("N_Test_Statistics", TestStatistic::kNTestStatistics);

    py::class_<SampleHandlerInterface, PySampleHandlerInterface /* <--- trampoline*/>(m_samples, "SampleHandlerInterface")
        .def(py::init())

        .def(
            "reweight",
            &SampleHandlerInterface::Reweight,
            "reweight the MC events in this sample. You will need to override this.")

        .def(
            "get_n_samples",
            &SampleHandlerInterface::GetNSamples,
            "Get the total number of samples"
        )

        .def(
            "get_likelihood",
            &SampleHandlerInterface::GetLikelihood,
            "Get the sample likelihood at the current point in your model space. You will need to override this.")

        .def(
            "set_test_stat",
            &SampleHandlerInterface::SetTestStatistic,
            "Set the test statistic that should be used when calculating likelihoods. \n\
            :param test_stat: The new test statistic to use",
            py::arg("test_stat"))

        .def(
            "get_bin_LLH",
            py::overload_cast<double, double, double>(&SampleHandlerInterface::GetTestStatLLH, py::const_),
            "Get the LLH for a bin by comparing the data and MC. The result depends on having previously set the test statistic using :py:meth:`pyMaCh3.samples.SampleHandlerInterface.set_test_stat` \n\
            :param data: The data content of the bin. \n\
            :param mc: The mc content of the bin \n\
            :param w2: The Sum(w_{i}^2) (sum of weights squared) in the bin, which is sigma^2_{MC stats}",
            py::arg("data"),
            py::arg("mc"),
            py::arg("w2"))

        ; // End of SampleHandlerInterface binding

    py::class_<SampleHandlerBase, PySampleHandlerBase /* <--- trampoline*/,
               SampleHandlerInterface>(m_samples, "SampleHandlerBase")
        .def(
            py::init<std::string, ParameterHandlerGeneric *>(),
            "This should never be called directly as SampleHandlerBase is an abstract base class. \n\
            However when creating a derived class, in the __init__() method, you should call the parent constructor i.e. this one by doing:: \n\
            \n\
            \tsuper(<your derived SampleHandler class>, self).__init__(*args) \n\
            \n ",
            py::arg("mc_version"), py::arg("xsec_cov"))

        // ================
        // Useful getters
        // ===============
        .def(
            "get_sample_title",
            &SampleHandlerBase::GetSampleTitle,
            py::arg("sample"),
            "Get the title for a given sample"
        )

        .def(
            "get_n_dim",
            &SampleHandlerBase::GetNDim,
            py::arg("sample"),
            "Get the dimension of a given sample"
        )

        .def(
            "add_data",
            py::overload_cast<const int, const std::vector<double>&>(&SampleHandlerBase::AddData),
            py::arg("sample"),
            py::arg("data_array"),
            "Set the data for your sample handler (assumes the binning is the same as your MC!)"
        )


        // ================
        // Histogramming 
        // ================
        .def(
            "get_data_array", 
            &SampleHandlerBase::GetDataArray,
            py::arg("sample"),
            "Returns the contents of the MC histogram as a flat list"
        )

        .def(
            "get_mc_array", 
            &SampleHandlerBase::GetMCArray,
            py::arg("sample"),
            "Returns the contents of the MC histogram as a flat list"
        )

        .def(
            "get_w2_array", 
            &SampleHandlerBase::GetW2Array,
            py::arg("sample"),
            "Returns the contents of the W2 histogram as a flat list"
        )

        .def(
            "get_mc_hist",
            [](SampleHandlerBase &self, const int sample)
            {
                int Dimension = self.GetNDim(sample);
                auto hist_original = M3::Clone(std::self.GetMCHist(sample));

                auto edges = HistToNumpy(hist_original, Dimension);
                return edges;
            },
            
            py::return_value_policy::reference_internal,
            py::arg("sample"),
            "Get MC histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")

        .def(
            "get_data_hist",
            [](SampleHandlerBase &self, const int sample)
            {
                int Dimension = self.GetNDim(sample);
                auto hist_original = M3::Clone(std::self.GetDataHist(sample));

                auto edges = HistToNumpy(hist_original, Dimension);
                return edges;

            },
            py::arg("Dimension"),
            "Get Data histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")

        .def(
            "get_w2_hist",
            [](SampleHandlerBase &self, const int sample)
            {
                int Dimension = self.GetNDim(sample);
                auto hist_original = M3::Clone(std::self.GetW2Hist(sample));

                auto edges = HistToNumpy(hist_original, Dimension);
                return edges;
            },


            py::arg("sample"),
            "Get W2 histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D");
        
        .def("get_var_hist", 
            [](SampleHandlerBase &self, const int iSample,
                                    const std::string& ProjectionVarX,
                                    const std::string& ProjectionVarY="",
                                    const std::vector<KinematicCut> &EventSelectionVec = {}
                                    int WeightStyle = 0,
                                    const std::vector< KinematicCut >& SubEventSelectionVec = {})
            {

                py::array_t<M3::float_t> edgesY, edgesX, contents
                if(ProjectionVarY==""){
                    auto hist = self.Get1DVarHistByModeAndChannel(iSample, ProjectionVarX, WeightStyle, SubEventSelectionVec;
                    const auto[contents, edgesX] = TH1ToNumpy(hist.get())
                    const auto edgesY = py::array_t<M3::float_t>();
                } else{                
                    auto hist = self.Get2DVarHistByModeAndChannel(iSample, ProjectionVarX, ProjectionVarY, EventSelectionVec, WeightStyle, SubEventSelectionVec);
                    const auto [contents, edgesX, edgesY] = TH2ToNumpy(hist.get());
                }
                return py::make_tuple(contents, edgesX, edgesY);
            }
        )
        ; // End of SampleHandler Base


    py::class_<KinematicCut>(m, "KinematicCut")
        .def(py::init<>(), "Simple wrapper around Kinematic cuts")
        .def_readwrite("param_name", &KinematicCut::ParamToCutOnIt, "Parameter to cut on")
        .def_readwrite("lower_bound", &KinematicCut::LowerBound, "Lower bound")
        .def_readwrite("upper_bound", &KinematicCut::UpperBound, "Upper Bound");

        
    /* Not sure if this will be needed in future versions of MaCh3 so leaving commented for now
    py::class_<fdmc_base>(m_samples, "MCstruct")
        .def(py::init())
        
        // Because a lot of the variables in fdmc_base use c style arrays,
        // we need to provide some setter functions to be able to set them using more
        // "pythony" objects, e.g. lists and numpy arrays
        .def(
            "set_event_variable_values", 
            [](fdmc_base &self, int dim, py::array_t<M3::float_t, py::array::c_style> &array)
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
