#include "splines/oneAPIUtils.dp.hpp"
 


SplineMonoUSM::SplineMonoUSM(sycl::queue queue,
                             int coeff_x_size,
                             int coeff_many_size,
                             int nKnots_arr_size,
                             int nKnots_arr_size):m_queue(&queue){
    coeff_x_usm = sycl::malloc_shared<float>(coeff_x_size, queue);
    coeff_many_usm = sycl::malloc_shared<float>(coeff_many_size, queue);
    nKnots_arr_usm = sycl::malloc_shared<unsigned int>(nKnots_arr_size, queue);
    paramNo_arr_usm = sycl::malloc_shared<short int>(nKnots_arr_size, queue);
}

SplineMonoUSM::~SplineMonoUSM(){
    sycl::free(coeff_x_usm, m_queue);
    sycl::free(coeff_many_usm, m_queue);
    sycl::free(nKnots_arr_usm, m_queue);
    sycl::free(paramNo_arr_usm, m_queue);
}