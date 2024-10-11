// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
// MaCh3 includes
#include "splines/SplineBase.h"
#include "splines/SplineMonolith.h"
#include "splines/SplineStructs.h"
#include "samplePDF/Structs.h" // <- The spline stuff that's in here should really be moved to splineStructs.h but I ain't doing that right now
// ROOT includes
#include "TSpline.h"

namespace py = pybind11;

// As SplineBase is an abstract base class we have to do some gymnastics to get it to get it into python
class PySplineBase : public SplineBase {
public:
    /* Inherit the constructors */
    using SplineBase::SplineBase;

    /* Trampoline (need one for each virtual function) */
    void Evaluate() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,        /* Return type */
            SplineBase,  /* Parent class */
            "evaluate",  /* Name in python*/
            Evaluate     /* Name of function in C++ (must match Python name) */
        );
    }
    
    std::string GetName() const override {
        PYBIND11_OVERRIDE_NAME(
            std::string, /* Return type */
            SplineBase,  /* Parent class */
            "get_name",  /* Name in python*/
            GetName      /* Name of function in C++ (must match Python name) */
        );
    }

    void FindSplineSegment() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,             /* Return type */
            SplineBase,       /* Parent class */
            "find_segment",   /* Name in python*/
            FindSplineSegment /* Name of function in C++ (must match Python name) */
        );
    }

    void CalcSplineWeights() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,                /* Return type */
            SplineBase,          /* Parent class */
            "calculate_weights", /* Name in python*/
            CalcSplineWeights    /* Name of function in C++ (must match Python name) */
        );
    }

    void ModifyWeights() override {
        PYBIND11_OVERRIDE_PURE_NAME(
            void,             /* Return type */
            SplineBase,       /* Parent class */
            "modify_weights", /* Name in python*/
            ModifyWeights     /* Name of function in C++ (must match Python name) */
        );
    }
};


void initSplines(py::module &m){

    auto m_splines = m.def_submodule("splines");
    m_splines.doc() = 
        "This is a Python binding of MaCh3s C++ based spline library.";

    // Bind the interpolation type enum that lets us set different interpolation types for our splines
    py::enum_<SplineInterpolation>(m_splines, "InterpolationType")
            .value("Linear", SplineInterpolation::kLinear)
            .value("Linear_Func", SplineInterpolation::kLinearFunc)
            .value("Cubic_TSpline3", SplineInterpolation::kTSpline3)
            .value("Cubic_Monotonic", SplineInterpolation::kMonotonic)
            .value("Cubic_Akima", SplineInterpolation::kAkima)
            .value("N_Interpolation_Types", SplineInterpolation::kSplineInterpolations);

    
    py::class_<SplineBase, PySplineBase /* <--- trampoline*/>(m_splines, "SplineBase");

    py::class_<TResponseFunction_red>(m_splines, "_ResponseFunctionBase")
        .doc() = "Base class of the response function, this binding only exists for consistency with the inheritance structure of the c++ code. Just pretend it doesn't exist and don't worry about it...";

    // Bind the TSpline3_red class. Decided to go with a clearer name of ResponseFunction for the python binding
    // and make the interface a bit more python-y. Additionally remove passing root stuff so we don't need to deal 
    // with root python binding and can just pass it native python objects.
    py::class_<TSpline3_red, TResponseFunction_red>(m_splines, "ResponseFunction")
        .def(
            // define a more python friendly constructor that massages the inputs and passes them
            // through to the c++ constructor
            py::init
            (
                // Just take in some vectors, then build a TSpline3 and pass this to the constructor
                [](std::vector<double> xVals, std::vector<double> yVals, SplineInterpolation interpType)
                {
                    if ( xVals.size() != yVals.size() )
                    {
                        throw MaCh3Exception(__FILE__, __LINE__, "Different number of x values and y values!");
                    }

                    int length = xVals.size();

                    TSpline3 *splineTmp = new TSpline3( "spline_tmp", xVals.data(), yVals.data(), length );

                    return std::make_unique<TSpline3_red>(splineTmp, interpType);
                }
            )
        )

        .def(
            "find_segment", 
            &TSpline3_red::FindX, 
            "Find the segment that a particular *value* lies in. \n"
            ":param value: The value to test",
            py::arg("value")
        )

        .def(
            "evaluate",
            &TSpline3_red::Eval,
            "Evaluate the response function at a particular *value*. \n"
            ":param value: The value to evaluate at.",
            py::arg("value")
        )
    ; // End of binding for ResponseFunction

    py::class_<SMonolith, SplineBase>(m_splines, "EventSplineMonolith")
        .def(
            py::init(
                [](std::vector<std::vector<TResponseFunction_red*>> &responseFns, const bool saveFlatTree)
                {
                    std::vector<RespFuncType> respFnTypes;
                    for(int i = 0; i > responseFns.size(); i++)
                    {
                        // ** WARNING **
                        // Right now I'm only pushing back TSpline3_reds as thats all thats supported right now
                        // In the future there might be more
                        // I think what would be best to do would be to store the interpolation type somehow in the ResponseFunction objects
                        // then just read them here and pass through to the constructor
                        respFnTypes.push_back(RespFuncType::kTSpline3_red);
                    }
                    return std::make_unique<SMonolith>(responseFns, respFnTypes, saveFlatTree);
                }
            ),
            "Create an EventSplineMonolith \n"
            ":param master_splines: These are the 'knot' values to make splines from. \n"
            ":param save_flat_tree: Whether we want to save monolith into speedy flat tree",
            py::arg("master_splines"),
            py::arg("save_flat_tree") = false
        )

        .def(
            py::init<std::string>(),
            "Constructor where you pass path to preprocessed root FileName which is generated by creating an EventSplineMonolith with the `save_flat_tree` flag set to True. \n"
            ":param file_name: The name of the file to read from.",
            py::arg("file_name")
        )

        .def(
            "evaluate",
            &SMonolith::Evaluate,
            "Evaluate the splines at their current values."
        )

        .def(
            "sync_mem_transfer",
            &SMonolith::SynchroniseMemTransfer,
            "This is important when running on GPU. After calculations are done on GPU we copy memory to CPU. This operation is asynchronous meaning while memory is being copied some operations are being carried. Memory must be copied before actual reweight. This function make sure all has been copied."
        )

        .def(
            "get_event_weight",
            &SMonolith::retPointer,
            py::return_value_policy::reference,
            "Get the weight of a particular event. \n"
            ":param event: The index of the event whose weight you would like.",
            py::arg("event")
        )

        .def(
            "set_param_value_array",
            // Wrap up the setSplinePointers method so that we can take in a numpy array and get 
            // pointers to it's sweet sweet data and use those pointers in the splineMonolith 
            [](SMonolith &self, py::array_t<double, py::array::c_style> &array)
            {
                py::buffer_info bufInfo = array.request();

                if ( bufInfo.ndim != 1)
                {
                    throw MaCh3Exception(__FILE__, __LINE__, "Number of dimensions in parameter array must be one!");
                }

                if ( bufInfo.shape[0] != self.GetNParams() )
                {
                    throw MaCh3Exception(__FILE__, __LINE__, "Number of entries in parameter array must equal the number of parameters!");
                }

                std::vector<const double *> paramVec;
                paramVec.resize(self.GetNParams());

                for( int idx = 0; idx < self.GetNParams(); idx++ )
                {
                    // booooo pointer arithmetic
                    paramVec[idx] = array.data() + idx;
                } 

                self.setSplinePointers(paramVec);
            },
            "Set the array that the monolith should use to read parameter values from. \n"
            "Usage of this might vary a bit from what you're used to in python. \n"
            "Rather than just setting the values here, what you're really doing is setting pointers in the underlying c++ code. \n"
            "What that means is that you pass an array to this function like:: \n"
            "\n event_spline_monolith_instance.set_param_value_array(array) \n\n"
            "Then when you set values in that array as normal, they will also be updated inside of the event_spline_monolith_instance.",
            py::arg("array")

        )

        .doc() = "This 'monolith' deals with event by event weighting using splines."

    ; // End of binding for EventSplineMonolith
}