#pragma once

// MaCh3 Includes
#include "manager/MaCh3Logger.h"

// ROOT Includes
#include "TMatrixDSym.h"

// C++ Includes
#include <vector>

/// @brief Contains information about adaptive covariance matrix
/// @see An adaptive Metropolis algorithm, H.Haario et al., 2001 for more info!


namespace adaptive_mcmc{

///@brief struct encapsulating all adaptive MCMC information
struct AdaptiveMCMCStruct{
  // Meta variables related to adaption run time
  int start_adaptive_throw; /// When do we start throwing
  int start_adaptive_update; /// When do we stop update the adaptive matrix
  int end_adaptive_update; /// Steps between changing throw matrix
  int adaptive_update_step; /// Steps between changing throw matrix
  std::vector<int> adapt_block_matrix_indices; ///Indices for block-matrix adaption
  std::vector<int> adapt_block_sizes; ///Size of blocks for adaption

  // Variables directely linked to adaption
  std::vector<double> par_means; /// Mean values for all parameters
  TMatrixDSym*  adaptive_covariance; /// Full adaptive covariance matrix
};


///@brief Basic printer method for AMCMC struct
inline void print_adaptive_struct(AdaptiveMCMCStruct adaptive_struct){
  MACH3LOG_INFO("Adaptive MCMC Info:");
  MACH3LOG_INFO("Throwing from New Matrix from Step : {}", adaptive_struct.start_adaptive_throw);
  MACH3LOG_INFO("Adaption Matrix Start Update       : {}", adaptive_struct.start_adaptive_update);
  MACH3LOG_INFO("Adaption Matrix Ending Updates     : {}", adaptive_struct.end_adaptive_update);
  MACH3LOG_INFO("Steps Between Updates              : {}", adaptive_struct.adaptive_update_step);
}




} // adaptive_mcmc namespace
