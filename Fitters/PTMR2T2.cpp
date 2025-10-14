#include "Fitters/PTMR2T2.h"
#include <chrono>

#ifdef MPIENABLED
PTMR2T2::PTMR2T2(manager *const manager) : MR2T2(manager) {
    AlgorithmName = "PTMR2T2";

    // Get the TempScale parameter
    TempScale = GetFromManager<double>(fitMan->raw()["General"]["MCMC"]["TempScale"], 3.0);

    TempScale = std::pow(TempScale, mpi_rank); // Each rank gets a different TempScale

    n_down_swap = 0; // Number of swaps attempted
    n_up_swap = 0;

    temp_adapt_step = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["TempAdaptRate"], 500);
    n_temp_adapts = GetFromManager<int>(fitMan->raw()["General"]["MCMC"]["NTempAdapts"], 5);

    // How many steps BEFORE we try to swap [so we can burn-in]
    steps_before_swapping = GetFromManager<unsigned int>(fitMan->raw()["General"]["MCMC"]["StepsBeforeSwapping"], 10000);

    round_trip_rate = 0.0;
    
    inv_TempScale = 1.0 / TempScale;

    temp_vec = std::vector<double>(n_procs);
    MPI_Allgather(&TempScale, 1, MPI_DOUBLE, temp_vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MACH3LOG_INFO("Using PTMR2T2 fitter with TempScale {} and AnnealTemps: ", TempScale);
    for(int i = 0; i < n_procs; ++i){
        MACH3LOG_INFO("Rank {}: T = {}", i, temp_vec[i]);
    }
    
    // Initialize timing benchmarks
}

void PTMR2T2::PrepareOutput()
{
    FitterBase::PrepareOutput();

    // Store temperature scale [used in reweighting so good to have/event]
    outTree->Branch("TempScale", &TempScale, "TempScale/D");

    // // Store timing benchmarks
    // #ifdef DEBUG
    // outTree->Branch("time_base_dostep", &time_base_dostep, "time_base_dostep/D");
    // outTree->Branch("time_swap", &time_swap, "time_swap/D");
    // outTree->Branch("time_global_adapt", &time_global_adapt, "time_global_adapt/D");
    // MACH3LOG_INFO("PTMR2T2 timing branches added to output tree");
    // #endif
}

void PTMR2T2::RunMCMC(){
    // Make a single transfer vector which will link to systematic values
    
    // First calculate the total size needed
    size_t total_params = 0;
    for(size_t s = 0; s < systematics.size(); ++s)
    {
        total_params += systematics[s]->GetNumParams();
        
        // We also want to make sure we're not starting with a tiny step scale
        // Thankfully it's not too bad to do this since temperature scales with sqrt(sigma)
        double init_scale = systematics[s]->GetGlobalStepScale();
        // Scale by sqrt(temp)
        systematics[s]->SetStepScale(init_scale * std::sqrt(TempScale));
    }

    // Allocate transfer vectors with space for header (3 elements: logLCurr, sample_logLCurr, TempScale)
    const size_t header_size = 3;
    transfer_vector_send.resize(header_size + total_params);
    transfer_vector_recv.resize(header_size + total_params);

    MR2T2::RunMCMC();
    sample_logLCurr = sample_logLProp;
}

void PTMR2T2::DoStep()
{
    // Adds in timing if debug enabled
    MR2T2::DoStep();

    if(step > steps_before_swapping){
        Swap();
    }
}


void PTMR2T2::ProposeStep()
{
    do_swap = false;
    MR2T2::ProposeStep();
    sample_logLProp = 0.0;
    for(size_t s = 0; s < samples.size(); ++s){
        logLProp -= sample_llh[s];
        // For swapping!
        sample_logLProp += sample_llh[s];
        sample_llh[s] *= inv_TempScale;
        logLProp += sample_llh[s];
    }
}

void PTMR2T2::Swap()
{
    const int offset = (mpi_rank % 2 == static_cast<int>(step) % 2) ? 1 : -1;
    const int swap_rank = mpi_rank + offset;

    if (swap_rank < 0 || swap_rank >= n_procs)
    {
        do_swap = false;
        return;
    }

    const bool i_am_lower = (swap_rank > mpi_rank);

    const int header_size = 3;

    double send_meta[header_size] = {logLCurr, sample_logLCurr, TempScale};
    double recv_meta[header_size] = {0.0, 0.0, 0.0};

    MPI_Sendrecv(
        send_meta, header_size, MPI_DOUBLE, swap_rank, 0,
        recv_meta, header_size, MPI_DOUBLE, swap_rank, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    LogLCurr_replica = recv_meta[0];
    sample_LogLCurr_replica= recv_meta[1];
    temp_replica = recv_meta[2];
    // Bit of a hack, send info as double to ensure we can do 1 message + avoid overhead

    int decision = 0;

    /// If either chain is out of bounds, we don't swap
    if (i_am_lower)
    {
        const double swap_rand = random->Uniform(0.0, 1.0);

        const double inv_temp_replica = 1.0 / temp_replica;

        /// We need to get our LogL at temperature = 0
        double logLCurr_t0 = logLCurr + (1 - inv_TempScale) * sample_logLCurr;
        double LogLCurr_replica_t0 = LogLCurr_replica + (1 - inv_temp_replica) * sample_LogLCurr_replica;

        double log_swap =-1*(LogLCurr_replica_t0 - logLCurr_t0) * (inv_TempScale - inv_temp_replica);

        // Protect against overflow in exp() for numerical stability
        const double max_log_swap = 1000.0;
        if (log_swap > max_log_swap) {
            log_swap = max_log_swap;
        }
        
        const double swap_prob = std::exp(log_swap);
        // MACH3LOG_INFO("Swap prob {}: logLCurr_t0 = {}, logL_replica = {}, T_curr = {}, T_replica = {}, swap_prob = {}, rand = {}",
        //               swap_prob, logLCurr_t0, LogLCurr_replica_t0, TempScale, temp_replica, swap_prob, swap_rand);

        do_swap = swap_rand < swap_prob &&
                  LogLCurr_replica < M3::_LARGE_LOGL_ &&
                  logLCurr < M3::_LARGE_LOGL_;  

        decision = do_swap ? 1 : 0;
        MPI_Send(&decision, 1, MPI_INT, swap_rank, 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&decision, 1, MPI_INT, swap_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        do_swap = (decision == 1);
    }

    if (do_swap)
    {
        i_am_lower ? n_up_swap++ : n_down_swap++;
        SwapStepInformation(swap_rank);
    }
}

void PTMR2T2::SwapStepInformation(int swap_rank)
{    
    const size_t header_size = 3;
    transfer_vector_send[0] = logLCurr;
    transfer_vector_send[1] = sample_logLCurr;
    transfer_vector_send[2] = TempScale;

    size_t idx = header_size;
    for(size_t s = 0; s < systematics.size(); ++s)
    {
        for (int p = 0; p < systematics[s]->GetNumParams(); ++p)
        {
            transfer_vector_send[idx++] = systematics[s]->GetParCurr(p);
        }
    }
    
    MPI_Sendrecv(
        transfer_vector_send.data(),
        static_cast<int>(transfer_vector_send.size()),
        MPI_DOUBLE,
        swap_rank,
        0,
        transfer_vector_recv.data(),
        static_cast<int>(transfer_vector_recv.size()),
        MPI_DOUBLE,
        swap_rank,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );
    
    LogLCurr_replica = transfer_vector_recv[0];
    sample_LogLCurr_replica = transfer_vector_recv[1];
    temp_replica = transfer_vector_recv[2];
    
    idx = header_size;
    for(size_t s = 0; s < systematics.size(); ++s)
    {
        const int num_params = systematics[s]->GetNumParams();
        for (int p = 0; p < num_params; ++p)
        {
            systematics[s]->SetParCurrProp(p, transfer_vector_recv[idx++]);
        }
    }

    // Now we need to swap information between chains

    #ifdef DEBUG
    double DEBUG_LLH_before = logLCurr;
    #endif
    
    // Get the penalty term first since this is not temp scaled
    double penalty = LogLCurr_replica - sample_LogLCurr_replica/temp_replica;
    if (penalty < 0)
    {
        MACH3LOG_CRITICAL("Negative Penalty term when swapping!!! logLCurr_replica: {}, sample_LogLCurr_replica: {}, temp_replica: {}, penalty {}", LogLCurr_replica, sample_LogLCurr_replica, temp_replica, penalty);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    

    // Scale sample likelihood
    sample_LogLCurr_replica *= temp_replica/TempScale;

    // Now do the swap
    sample_logLCurr = sample_LogLCurr_replica;
    
    logLCurr = sample_LogLCurr_replica + penalty;
    
    #ifdef DEBUG
    MACH3LOG_DEBUG("LLH Before swap: {} | After {}", DEBUG_LLH_before, logLCurr);
    #endif
}


void PTMR2T2::PrintProgress(){
    MR2T2::PrintProgress();

    MACH3LOG_INFO("Temperature: {:.2f}, Up Going: {:.2f}", TempScale, static_cast<double>(n_up_swap));

    // Optimized: Use MPI_Gather instead of sequential sends/receives
    std::vector<double> all_mpi_info;
    double local_mpi_info[4] = {
        static_cast<double>(mpi_rank), 
        TempScale, 
        static_cast<double>(n_up_swap), 
        static_cast<double>(n_down_swap)
    };
    
    if (mpi_rank == 0)
    {
        all_mpi_info.resize(4 * n_procs);
    }
    
    int result = MPI_Gather(local_mpi_info, 4, MPI_DOUBLE, 
                            all_mpi_info.data(), 4, MPI_DOUBLE, 
                            0, MPI_COMM_WORLD);
    
    if (result != MPI_SUCCESS)
    {
        MACH3LOG_ERROR("MPI_Gather failed in PrintProgress with error code {}", result);
        return;
    }
    
    if (mpi_rank == 0)
    {
        for (int i = 0; i < n_procs; ++i)
        {
            int rank = static_cast<int>(all_mpi_info[i * 4 + 0]);
            double temp = all_mpi_info[i * 4 + 1];
            int n_up = static_cast<int>(all_mpi_info[i * 4 + 2]);
            int n_down = static_cast<int>(all_mpi_info[i * 4 + 3]);
            int total_steps_up_down = n_up + n_down;

            double swap_rate = 0.0;
            if (total_steps_up_down > 0)
            {
                swap_rate = static_cast<double>(n_up) / static_cast<double>(total_steps_up_down);
            }
            MACH3LOG_INFO("Rank {}: Swap Rate: {:.2f} (up: {}, down: {}), Temp: {:.2f} | Swap/step = {:.2f} %", rank, swap_rate, n_up, n_down, temp, 100*total_steps_up_down/static_cast<double>(step+1-stepStart));
        }
    }
}

void PTMR2T2::AcceptStep()
{
    MR2T2::AcceptStep();
    sample_logLCurr = sample_logLProp;
}

// void PTMR2T2::GlobalAdaptTemperatures()
// {
//     // --------------------------------------------------------------------
//     // 0. Gather swap acceptance rates and temperatures from all ranks
//     // --------------------------------------------------------------------
//     double local_accept = 0.0;

//     int total_swaps = n_up_swap + n_down_swap;
//     if (total_swaps > 0)
//     {
//         local_accept = static_cast<double>(n_up_swap) / static_cast<double>(total_swaps);
//     }

//     std::vector<double> all_accepts(n_procs, 0.0);
//     std::vector<double> all_temps(n_procs, 0.0);

//     MPI_Gather(&local_accept, 1, MPI_DOUBLE,
//                all_accepts.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     MPI_Gather(&TempScale, 1, MPI_DOUBLE,
//                all_temps.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     // --------------------------------------------------------------------
//     // 1. Rank 0 computes new temperatures based on communication barriers
//     // --------------------------------------------------------------------
//     if (mpi_rank == 0)
//     {
//         // Sort by temperature to ensure correct ordering
//         std::vector<int> indices(n_procs);
//         std::iota(indices.begin(), indices.end(), 0);
//         std::sort(indices.begin(), indices.end(),
//                   [&](int a, int b)
//                   { return all_temps[a] < all_temps[b]; });

//         std::vector<double> temps_sorted(n_procs);
//         for (int i = 0; i < n_procs; ++i)
//             temps_sorted[i] = all_temps[indices[i]];

//         // Convert to inverse temperatures β = 1/T
//         std::vector<double> betas(n_procs);
//         for (int i = 0; i < n_procs; ++i)
//             betas[i] = 1.0 / temps_sorted[i];

//         // Sort acceptances to match temperature ordering
//         std::vector<double> accepts_sorted(n_procs);
//         for (int i = 0; i < n_procs; ++i)
//             accepts_sorted[i] = all_accepts[indices[i]];

//         // Compute local communication barriers λ_i = -log(accept_i)
//         std::vector<double> lambda(n_procs - 1);
//         for (int i = 0; i < n_procs - 1; ++i)
//         {
//             double acc = std::max(1e-6, std::min(1.0 - 1e-6, accepts_sorted[i]));
//             lambda[i] = -std::log(acc);
//         }

//         // Compute cumulative barrier Λ
//         std::vector<double> Lambda(n_procs);
//         Lambda[0] = 0.0;
//         for (int i = 1; i < n_procs; ++i)
//             Lambda[i] = Lambda[i - 1] + lambda[i - 1];

//         // Target equal spacing in Λ-space
//         std::vector<double> new_Lambda(n_procs);
//         for (int i = 0; i < n_procs; ++i)
//             new_Lambda[i] = i * Lambda.back() / (n_procs - 1);

//         // Linear interpolation helper lambda
//         auto interpolate = [&](const std::vector<double> &x,
//                                const std::vector<double> &y,
//                                double x_new)
//         {
//             if (x_new <= x.front())
//                 return y.front();
//             if (x_new >= x.back())
//                 return y.back();

//             auto it = std::upper_bound(x.begin(), x.end(), x_new);
//             size_t j = std::distance(x.begin(), it) - 1;
//             double t = (x_new - x[j]) / (x[j + 1] - x[j]);
//             return y[j] + t * (y[j + 1] - y[j]);
//         };

//         // Interpolate new β’s and convert back to temperatures
//         std::vector<double> new_betas(n_procs);
//         for (int i = 0; i < n_procs; ++i)
//             new_betas[i] = interpolate(Lambda, betas, new_Lambda[i]);

//         for (int i = 0; i < n_procs; ++i)
//             all_temps[i] = 1.0 / new_betas[i];

//         all_temps[0] = 1.0;
//         std::sort(all_temps.begin(), all_temps.end());
//     }

//     MPI_Bcast(all_temps.data(), n_procs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     TempScale = all_temps[mpi_rank];

//     double local_accept_for_rt = 0.0;
//     int total_swaps_rt = n_up_swap + n_down_swap;
//     if (total_swaps_rt > 0)
//         local_accept_for_rt = static_cast<double>(n_up_swap) / static_cast<double>(total_swaps_rt);

//     double global_round_trip = 0.0;
//     MPI_Allreduce(&local_accept_for_rt, &global_round_trip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     round_trip_rate = global_round_trip / static_cast<double>(n_procs);

//     if (mpi_rank == 0)
//     {
//         MACH3LOG_INFO("Adapted global temperature ladder:");
//         MACH3LOG_INFO("Global mean swap rate: {:.3f}", round_trip_rate);

//         for (int i = 0; i < n_procs; ++i)
//         {
//             double swap_rate = (i < n_procs - 1) ? all_accepts[i] : 0.0;
//             MACH3LOG_INFO("Rank {} -> T = {:.4f} | Swap Rate {:.3f}", i, all_temps[i], swap_rate);
//         }
//     }

//     n_up_swap = 0;
//     n_down_swap = 0;
// }

#endif