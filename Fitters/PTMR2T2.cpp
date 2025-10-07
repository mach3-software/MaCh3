#include "Fitters/PTMR2T2.h"

#ifdef MPIENABLED
PTMR2T2::PTMR2T2(manager *const manager) : MR2T2(manager) {
    AlgorithmName = "PTMR2T2";

    // Get the AnnealTemp parameter
    AnnealTemp = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["AnnealTempScale"], 3.0);

    AnnealTemp = std::pow(AnnealTemp, mpi_rank); // Each rank gets a different AnnealTemp

    n_down_swap = 0; // Number of swaps attempted
    n_up_swap = 0;


    max_temp = 100000000000;
    min_temp = 1.0;
    if(AnnealTemp > max_temp){
        AnnealTemp = max_temp;
    }
    if(AnnealTemp < min_temp){
        AnnealTemp = min_temp;
    }
    target_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TargetSwapRate"], 0.234);
    adaption_rate = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["AnnealTempAdaptionRate"], 0.01);

    temp_vec = std::vector<double>(n_procs);
    MPI_Allgather(&AnnealTemp, 1, MPI_DOUBLE, temp_vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MACH3LOG_INFO("Using PTMR2T2 fitter with AnnealTemp {} and AnnealTemps: ", AnnealTemp);
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
    if (step % 500 == 0 && step > 0)
        GlobalAdaptTemperatures();
}

void PTMR2T2::PrepareOutput(){
    MR2T2::PrepareOutput();
    outTree->Branch("n_down_swap", &n_down_swap);
    outTree->Branch("n_up_swap", &n_up_swap);
    outTree->Branch("AnnealTemp", &AnnealTemp);

    /// We want to also store debug information
    if (mpi_rank == 0)
    {
        outTree->Branch("Barrier_step", &Barrier_step);
        outTree->Branch("Barrier_lambda", &Barrier_lambda);
        outTree->Branch("Barrier_Lambda", &Barrier_Lambda);
        outTree->Branch("Barrier_Temps", &Barrier_Temps);
        outTree->Branch("Barrier_Betas", &Barrier_Betas);
    }
}


void PTMR2T2::ProposeStep()
{
    do_swap = false;

    MR2T2::ProposeStep();

    /// Need to scale the likelihoods by the AnnealTemp
    for (size_t s = 0; s < systematics.size(); ++s)
    {
        syst_llh[s] /= AnnealTemp;
    }
    for(size_t s = 0; s < samples.size(); ++s){
        sample_llh[s] /= AnnealTemp;
    }

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
    // Pack metadata: [logLProp, AnnealTemp, decision (if lower rank)]
    // This reduces to a single Sendrecv for metadata exchange
    // ========================================================================

    const bool i_am_lower = (swap_rank > mpi_rank);
    double send_meta[3] = {logLProp, AnnealTemp, 0.0};
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

        /// Need to do a bit of scale

        const double swap_prob = std::exp((LogL_replica - logLProp) * (1.0 / AnnealTemp - 1.0 / temp_replica));

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

void PTMR2T2::SwapStepInformation(int swap_rank)
{    
    // Pack data into single send buffer
    const size_t header_size = 2;
    transfer_vector_send[0] = logLProp;
    transfer_vector_send[1] = AnnealTemp;
    
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
    
    logLProp = LogL_replica * AnnealTemp / temp_replica;
}

void PTMR2T2::SynchronizeTemperatures()
{
    // Gather all AnnealTemps to ensure the ladder is properly ordered
    MPI_Allgather(&AnnealTemp, 1, MPI_DOUBLE, temp_vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    // Rank 0 sorts and redistributes if needed
    if (mpi_rank == 0)
    {
        std::sort(temp_vec.begin(), temp_vec.end());
        temp_vec[0] = 1.0; // Ensure canonical chain
    }

    // Broadcast the sorted AnnealTemps
    MPI_Bcast(temp_vec.data(), n_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each chain takes its AnnealTemp
    AnnealTemp = temp_vec[mpi_rank];
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
    // If acceptance rate is too low, decrease AnnealTemp (move closer to neighbor)
    // If acceptance rate is too high, increase AnnealTemp (move further from neighbor)
    double delta = adaption_rate * (current_swap_rate - target_rate);

    // Update AnnealTemp using log scale for stability
    double log_temp = std::log(AnnealTemp);
    log_temp += delta;
    AnnealTemp = std::exp(log_temp);

    // Enforce AnnealTemp bounds
    AnnealTemp = std::max(min_temp, std::min(max_temp, AnnealTemp));

    // Optionally: Synchronize AnnealTemps across chains periodically
    // This ensures the AnnealTemp ladder remains ordered
    if (step % 1000 == 0)
    {
        SynchronizeTemperatures();
    }
}

void PTMR2T2::GlobalAdaptTemperatures()
{
    // --------------------------------------------------------------------
    // 0. Gather swap acceptance rates and temperatures from all ranks
    // --------------------------------------------------------------------
    double local_accept = 0.0;

    int total_swaps = n_up_swap + n_down_swap;
    if (total_swaps > 0)
    {
        local_accept = static_cast<double>(n_up_swap) / static_cast<double>(total_swaps);
    }

    std::vector<double> all_accepts(n_procs, 0.0);
    std::vector<double> all_temps(n_procs, 0.0);

    MPI_Gather(&local_accept, 1, MPI_DOUBLE,
               all_accepts.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gather(&AnnealTemp, 1, MPI_DOUBLE,
               all_temps.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // --------------------------------------------------------------------
    // 1. Rank 0 computes new temperatures based on communication barriers
    // --------------------------------------------------------------------
    if (mpi_rank == 0)
    {
        // Sort by temperature to ensure correct ordering
        std::vector<int> indices(n_procs);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
                  [&](int a, int b)
                  { return all_temps[a] < all_temps[b]; });

        std::vector<double> temps_sorted(n_procs);
        for (int i = 0; i < n_procs; ++i)
            temps_sorted[i] = all_temps[indices[i]];

        // Convert to inverse temperatures β = 1/T
        std::vector<double> betas(n_procs);
        for (int i = 0; i < n_procs; ++i)
            betas[i] = 1.0 / temps_sorted[i];

        // Compute local communication barriers λ_i = -log(accept_i)
        std::vector<double> lambda(n_procs - 1);
        for (int i = 0; i < n_procs - 1; ++i)
        {
            double acc = std::max(1e-6, std::min(1.0 - 1e-6, all_accepts[i]));
            lambda[i] = -std::log(acc);
        }

        // Compute cumulative barrier Λ
        std::vector<double> Lambda(n_procs);
        Lambda[0] = 0.0;
        for (int i = 1; i < n_procs; ++i)
            Lambda[i] = Lambda[i - 1] + lambda[i - 1];

        // Target equal spacing in Λ-space
        std::vector<double> new_Lambda(n_procs);
        for (int i = 0; i < n_procs; ++i)
            new_Lambda[i] = i * Lambda.back() / (n_procs - 1);

        // Linear interpolation helper lambda
        auto interpolate = [&](const std::vector<double> &x,
                               const std::vector<double> &y,
                               double x_new)
        {
            if (x_new <= x.front())
                return y.front();
            if (x_new >= x.back())
                return y.back();

            auto it = std::upper_bound(x.begin(), x.end(), x_new);
            size_t j = std::distance(x.begin(), it) - 1;
            double t = (x_new - x[j]) / (x[j + 1] - x[j]);
            return y[j] + t * (y[j + 1] - y[j]);
        };

        // Interpolate new β’s and convert back to temperatures
        std::vector<double> new_betas(n_procs);
        for (int i = 0; i < n_procs; ++i)
            new_betas[i] = interpolate(Lambda, betas, new_Lambda[i]);

        for (int i = 0; i < n_procs; ++i)
            all_temps[i] = 1.0 / new_betas[i];

        // Enforce T[0] = 1.0 canonical chain
        all_temps[0] = 1.0;
        std::sort(all_temps.begin(), all_temps.end());
    }

    // --------------------------------------------------------------------
    // 2. Broadcast new temperature ladder to all ranks
    // --------------------------------------------------------------------
    MPI_Bcast(all_temps.data(), n_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    AnnealTemp = all_temps[mpi_rank];

    // --------------------------------------------------------------------
    // 3. Reset swap counters for next adaptation window
    // --------------------------------------------------------------------
    n_up_swap = 0;
    n_down_swap = 0;

    if (mpi_rank == 0)
    {
        MACH3LOG_INFO("Adapted global temperature ladder:");
        for (int i = 0; i < n_procs; ++i)
            MACH3LOG_INFO("   Rank {} -> T = {}", i, all_temps[i]);
    }
}

#endif