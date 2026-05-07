#pragma once

/// @file samples.h
/// @author Ewan Miller

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "Samples/SampleHandlerBase.h"
#include "TH1.h"
#include "TH2.h"

#include "histutils.h"

namespace py = pybind11;


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
            "get_n_dim",
            &SampleHandlerInterface::GetNDim,
            py::arg("sample"),
            "Get the dimension of a given sample"
        )

        .def(
            "get_mc_hist",
            [](SampleHandlerBase &self, const int sample)
            {
                auto hist_original = M3::Clone(self.GetMCHist(sample));

                auto edges = HistToNumpy(hist_original);
                return edges;
            },

            py::return_value_policy::reference_internal,
            py::arg("sample"),
            "Get MC histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")

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

        .def(
            "add_data",
            py::overload_cast<const int, const std::vector<double>&>(&SampleHandlerBase::AddData),
            py::arg("sample"),
            py::arg("data_array"),
            "Set the data for your sample handler (assumes the binning is the same as your MC!)"
        )

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
            "get_data_array", 
            py::overload_cast<const int>(&SampleHandlerBase::GetDataArray, py::const_),
            py::arg("sample"),
            "Returns the contents of the MC histogram as a flat list"
        )

        .def(
            "get_mc_array", 
            py::overload_cast<const int>(&SampleHandlerBase::GetMCArray, py::const_),
            py::arg("sample"),
            "Returns the contents of the MC histogram as a flat list"
        )

        .def(
            "get_w2_array", 
            py::overload_cast<const int>(&SampleHandlerBase::GetW2Array, py::const_),
            py::arg("sample"),
            "Returns the contents of the W2 histogram as a flat list"
        )

        .def(
            "get_data_hist",
            [](SampleHandlerBase &self, const int sample)
            {
                auto hist_original = M3::Clone(self.GetDataHist(sample));

                auto edges = HistToNumpy(hist_original);
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
                auto hist_original = M3::Clone(self.GetW2Hist(sample));

                auto edges = HistToNumpy(hist_original);
                return edges;
            },


            py::arg("sample"),
            "Get W2 histogram as numpy arrays.\n"
            "For 1D: Returns (contents, edges)\n"
            "For 2D: Returns (contents, edgesX, edgesY)\n"
            "where contents is shape (nbinsY, nbinsX) for 2D")
        
        .def("get_var_hist", 
            [](SampleHandlerBase &self, const int iSample,
                                    const std::string& ProjectionVarX,
                                    const std::string& ProjectionVarY="",
                                    const std::vector<KinematicCut> &EventSelectionVec = {},
                                    int WeightStyle = 0,
                                    const std::vector< KinematicCut >& SubEventSelectionVec = {})
            {

                py::array_t<M3::float_t> edgesY, edgesX, contents;
                std::unique_ptr<TH1> hist;
                int dim;
                if(ProjectionVarY==""){
                    hist = self.Get1DVarHist(iSample,
                                             ProjectionVarX,
                                             EventSelectionVec,
                                             WeightStyle,
                                             SubEventSelectionVec);
                } else{
                    hist = self.Get2DVarHist(iSample, 
                                     ProjectionVarX,
                                     ProjectionVarY,
                                     EventSelectionVec,
                                     WeightStyle,
                                     SubEventSelectionVec);
                }                    
                return HistToNumpy(hist);
            }
            )
        ; // End of SampleHandler Base


    py::class_<KinematicCut>(m_samples, "KinematicCut")
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
