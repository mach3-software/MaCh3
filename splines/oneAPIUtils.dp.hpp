#pragma once
#include <sycl/sycl.hpp>

struct SplineMonoUSM {

    SplineMonoUSM(sycl::queue& queue,
                  int coeff_x_size,
                  int coeff_many_size,
                  int nKnots_arr_size,
                  int paramNo_arr_size);
    ~SplineMonoUSM();

    sycl::queue& m_queue;

    // *******************
    /// CPU arrays to hold X coefficient
    float* coeff_x;

    /// CPU arrays to hold other coefficients
    float* coeff_many;

    /// CPU Number of knots per spline
    unsigned int* nKnots_arr;

    /// CPU array with the number of points per spline (not per spline point!)
    short int* paramNo_arr;

    int* param_n_knots;
};