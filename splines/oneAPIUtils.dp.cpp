#include "splines/oneAPIUtils.dp.hpp"
#include "sycl/usm.hpp"
 


SplineMonoUSM::SplineMonoUSM(sycl::queue& queue,
                             int coeff_x_size,
                             int coeff_many_size,
                             int nKnots_arr_size,
                             int paramNo_arr_size,
                             int number_of_events_size):m_queue(queue){
    coeff_x = sycl::malloc_shared<float>(coeff_x_size, queue);
    coeff_many = sycl::malloc_shared<float>(coeff_many_size, queue);
    nKnots_arr = sycl::malloc_shared<unsigned int>(nKnots_arr_size, queue);
    paramNo_arr = sycl::malloc_shared<short int>(paramNo_arr_size, queue);
    splines_per_event_arr = sycl::malloc_host<unsigned short int>(number_of_events_size, queue);
}

SplineMonoUSM::~SplineMonoUSM(){
    sycl::free(coeff_x, m_queue);
    sycl::free(coeff_many, m_queue);
    sycl::free(nKnots_arr, m_queue);
    sycl::free(paramNo_arr, m_queue);
    sycl::free(splines_per_event_arr, m_queue);
}