#include "Fitters/PTMR2T2.h"

#ifdef MPIENABLED
PTMR2T2::PTMR2T2(manager *const manager) : MR2T2(manager) {
    AlgorithmName = "PTMR2T2";

    // Get the temperature parameter
    temperature = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TemperatureScale"], 3.0);

    temperature = std::pow(temperature, mpi_rank); // Each rank gets a different temperature
    n_down_swap = 0; // Number of swaps attempted
    n_up_swap = 0;


    max_temp = 100000;
    min_temp = 1.0;
    if(temperature > max_temp){
        temperature = max_temp;
    }
    if(temperature < min_temp){
        temperature = min_temp;
    }
    target_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TargetSwapRate"], 0.234);
    adaption_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TemperatureAdaptionRate"], 0.01);

    temp_vec = std::vector<double>(n_procs);
    MPI_Allgather(&temperature, 1, MPI_DOUBLE, temp_vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MACH3LOG_INFO("Using PTMR2T2 fitter with temperature {} and temperatures: ", temperature);
    for(int i = 0; i < n_procs; ++i){
        MACH3LOG_INFO("Rank {}: T = {}", i, temp_vec[i]);
    } 
}

void PTMR2T2::RunMCMC(){
    // Make a single transfer vector which will link to systematic values

    for(size_t s = 0; s < systematics.size(); ++s)
    {
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            transfer_vector_send.push_back(systematics[s]->GetParCurr(p));
        }
    }

    transfer_vector_send = std::vector<double>(transfer_vector_send.size() + 2); // Add space for header info
    transfer_vector_recv = std::vector<double>(transfer_vector_send.size());

    MR2T2::RunMCMC();
}

void PTMR2T2::DoStep()
{
    MR2T2::DoStep();
    // Now we do the swap
    Swap();
    // AdaptTemperature();
}

void PTMR2T2::ProposeStep()
{
    do_swap = false;

    MR2T2::ProposeStep();

    /// Need to scale the likelihoods by the temperature
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        syst_llh[s] /= temperature;
    }
    for(size_t s = 0; s < samples.size(); ++s){
        sample_llh[s] /= temperature;
    }

    logLProp /= temperature;
}

void PTMR2T2::Swap()
{
    // Calculate swap partner
    const int offset = (mpi_rank % 2 == static_cast<int>(step) % 2) ? 1 : -1;
    const int swap_rank = mpi_rank + offset;

    // Early exit if no valid partner
    if (swap_rank < 0 || swap_rank >= n_procs)
    {
        do_swap = false;
        return;
    }

    // ========================================================================
    // Pack metadata: [logLProp, temperature, decision (if lower rank)]
    // This reduces to a single Sendrecv for metadata exchange
    // ========================================================================

    const bool i_am_lower = (swap_rank > mpi_rank);
    double send_meta[3] = {logLProp, temperature, 0.0};
    double recv_meta[3] = {0.0, 0.0, 0.0};

    // If I'm the lower rank, compute decision before sending
    if (i_am_lower)
    {
        // We need partner's data first, so we'll do decision in next step
        // For now, just exchange basic data
    }

    MPI_Sendrecv(
        send_meta, 2, MPI_DOUBLE, swap_rank, 0,
        recv_meta, 2, MPI_DOUBLE, swap_rank, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    LogL_replica = recv_meta[0];
    temp_replica = recv_meta[1];

    // Make and communicate decision
    int decision = 0;

    if (i_am_lower)
    {
        const double swap_rand = random->Uniform(0.0, 1.0);
        const double swap_prob = std::exp((LogL_replica - logLProp) * (1.0 / temperature - 1.0 / temp_replica));

        do_swap = swap_rand < swap_prob &&
                  !out_of_bounds &&
                  LogL_replica < M3::_LARGE_LOGL_ &&
                  logLProp < M3::_LARGE_LOGL_;

        decision = do_swap ? 1 : 0;
        MPI_Send(&decision, 1, MPI_INT, swap_rank, 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&decision, 1, MPI_INT, swap_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        do_swap = (decision == 1);
    }

    // Execute swap
    if (do_swap)
    {
        i_am_lower ? n_up_swap++ : n_down_swap++;
        SwapStepInformation(swap_rank);
    }
}

/*
void PTMR2T2::Swap(){
    // We're gonna use an even/odd swapping regime    
    
    // Alternate between (0,1),(2,3)... on even steps and (1,2),(3,4)... on odd steps
    int swap_rank = mpi_rank + ((mpi_rank % 2 == static_cast<int>(step) % 2) ? 1 : -1);

    MPI_Barrier(MPI_COMM_WORLD);

    if(swap_rank > mpi_rank and swap_rank < n_procs){
        // We send info
        MPI_Send(&logLProp, 1, MPI_DOUBLE, swap_rank, 0, MPI_COMM_WORLD);
        MPI_Send(&temperature, 1, MPI_DOUBLE, swap_rank, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(swap_rank < mpi_rank and swap_rank >= 0){
        // We receive info
        MPI_Recv(&LogL_replica, 1, MPI_DOUBLE, swap_rank, 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&temp_replica, 1, MPI_DOUBLE, swap_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
    }


    if (swap_rank < mpi_rank and swap_rank >= 0){
        double swap_rand = random->Uniform(0.0, 1.0);
        double swap_prob = std::exp((LogL_replica - logLProp) * (1.0/temperature - 1.0/temp_replica));
        if (swap_rand < swap_prob and !out_of_bounds and !(LogL_replica >= M3::_LARGE_LOGL_ or logLProp >= M3::_LARGE_LOGL_)){
            do_swap = true;
        }

        MPI_Send(&do_swap, 1, MPI_C_BOOL, swap_rank, 0, MPI_COMM_WORLD);
    }
    else if(swap_rank > mpi_rank and swap_rank < n_procs){
        MPI_Recv(&do_swap, 1, MPI_C_BOOL, swap_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /// Now do the swap if we decided to
    if(do_swap){
        if(swap_rank > mpi_rank){
            n_up_swap++;
        }
        else{
            n_down_swap++;
        }
        // We swap over our step information
        SwapStepInformation(swap_rank);
    }
}
*/
void PTMR2T2::SwapStepInformation(int swap_rank)
{    
    // Pack data into single send buffer
    const size_t header_size = 2;
    transfer_vector_send[0] = logLProp;
    transfer_vector_send[1] = temperature;
    
    int idx = header_size;
    for(size_t s = 0; s < systematics.size(); ++s)
    {
        const auto& params = systematics[s]->GetParPropVec();
        std::copy(params.begin(), params.end(), transfer_vector_send.begin() + idx);
        idx += systematics[s]->GetNumParams();
    }
    
    // ========================================================================
    // MPI_Sendrecv: Simultaneous send/receive in single call
    // This is the safest and often fastest approach for paired exchanges
    // ========================================================================
    
    MPI_Sendrecv(
        transfer_vector_send.data(),                        // send buffer
        static_cast<int>(transfer_vector_send.size()),     // send count
        MPI_DOUBLE,                                         // send type
        swap_rank,                                          // destination
        0,                                                  // send tag
        transfer_vector_recv.data(),                        // receive buffer
        static_cast<int>(transfer_vector_recv.size()),     // receive count
        MPI_DOUBLE,                                         // receive type
        swap_rank,                                          // source
        0,                                                  // receive tag
        MPI_COMM_WORLD,                                     // communicator
        MPI_STATUS_IGNORE                                   // status
    );
    
    // Extract received data
    LogL_replica = transfer_vector_recv[0];
    temp_replica = transfer_vector_recv[1];
    
    // Update parameters
    idx = header_size;
    for(size_t s = 0; s < systematics.size(); ++s)
    {
        const int num_params = systematics[s]->GetNumParams();
        systematics[s]->SetParameters(
            std::vector<double>(
                transfer_vector_recv.begin() + idx,
                transfer_vector_recv.begin() + idx + num_params
            )
        );
        idx += num_params;
    }
    
    logLProp = LogL_replica * temperature / temp_replica;
}

void PTMR2T2::SynchronizeTemperatures()
{
    // Gather all temperatures to ensure the ladder is properly ordered
    MPI_Allgather(&temperature, 1, MPI_DOUBLE, temp_vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    // Rank 0 sorts and redistributes if needed
    if (mpi_rank == 0)
    {
        std::sort(temp_vec.begin(), temp_vec.end());
        temp_vec[0] = 1.0; // Ensure canonical chain
    }

    // Broadcast the sorted temperatures
    MPI_Bcast(temp_vec.data(), n_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each chain takes its temperature
    temperature = temp_vec[mpi_rank];
}

void PTMR2T2::AdaptTemperature()
{
    // Only adapt for chains that are not at the boundaries
    if (mpi_rank == 0)
    {
        // This is the canonical chain (T=1), don't adapt
        return;
    }

    // Calculate current swap acceptance rate
    int total_swaps = n_up_swap + n_down_swap;
    if (total_swaps == 0)
        return; // No swaps yet, can't adapt

    double current_swap_rate = static_cast<double>(n_up_swap) / static_cast<double>(total_swaps);

    // Adaptive step using stochastic approximation
    // If acceptance rate is too low, decrease temperature (move closer to neighbor)
    // If acceptance rate is too high, increase temperature (move further from neighbor)
    double delta = adaption_rate * (current_swap_rate - target_rate);

    // Update temperature using log scale for stability
    double log_temp = std::log(temperature);
    log_temp += delta;
    temperature = std::exp(log_temp);

    // Enforce temperature bounds
    temperature = std::max(min_temp, std::min(max_temp, temperature));

    // Optionally: Synchronize temperatures across chains periodically
    // This ensures the temperature ladder remains ordered
    if (step % 1000 == 0)
    {
        SynchronizeTemperatures();
    }
}

#endif