#include"main_header.h"

int main( const int argC, char* const argV[] ){ // arguments required for boost program options
  // ====== Initialize MPI ======
  // world_size is the number of cores, my_rank is the number of "this" core 
  int world_size, my_rank;
  MPI_Init( nullptr, nullptr );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

  // ====== Initialize Time Measurement ======
  tmm::DerivedTimeMeasure my_clock( my_rank, world_size );

  // ====== Import, Interpret and Print Parameters ======
  const ps::ParameterSpace my_pspace( argC, argV, world_size, my_rank );
  print::print_R0( my_rank, my_pspace.create_essentials_string() );

  // ====== Initialize Run Time Data and Generate Seed ======
  rtd::RunTimeData my_rtdata( my_pspace, my_rank );
  my_rtdata.generated_seed = rd::generate_seed( my_rtdata.get_seed_str(), my_rank );

  // ====== Initialize several Constant Spin Matrices on the Hilbert Space ======
  if( !my_pspace.uncoupled_spins )
  {
    init::write_general_spin_matrices( my_pspace.spin_float_list, my_pspace.num_HilbertSpaceDimension );
    init::write_cluster_Hamiltonian( my_pspace.num_Spins, my_pspace.num_HilbertSpaceDimension, my_pspace.J, my_pspace.spinspin_cmodel, my_pspace.chemical_shift, my_pspace.local_extra_interaction );
  }
  else
  {
    init::write_uncoupled_spin_matrices();
  }

  // ====== Declare the Spin Cluster Correlation Functions ======
  CluCorrTen my_spin_correlations_Re = init::generate_initial_environment_spin_correlations( my_pspace, "Re" );
  CluCorrTen my_spin_correlations_Im = init::generate_initial_environment_spin_correlations( my_pspace, "Im" );
  func::SiteFields my_spin_expvals = ( my_pspace.init_corr_mode == "import" ) ? my_pspace.imported_spin_expvals : func::SiteFields( my_pspace.num_Spins, FieldVector{0.,0.,0.} );
  my_clock.measure( "initialization", true );


  print::print_R0( my_rank, "+++++++++++++++++++++++++ Self-Consistent Iteration ++++++++++++++++++++++++++\n" );

  // Persistent pCN chain state across self-consistent iterations, kept in the basis-
  // independent full frequency representation (the per-block diagonalization changes
  // between iterations but the Matsubara frequencies do not). Empty before the first
  // iteration; carried forward thereafter to warm-start the chain.
  mvgb::EigenValuesBlocks my_full_freq_persistent;

  my_clock.enter_loop();
  do
  {
    my_clock.new_iteration();
    std::string it_str = std::to_string( my_rtdata.num_Iterations + 1 );
    print::print_R0( my_rank, "|\n|\n--------------------------- Iteration Step " + it_str + " ---------------------------\n" );

    // ====== Initialize Duration Estimator ======
    std::vector<size_t> loop_sizes{ my_rtdata.get_num_SamplesPerCore() };
    std::vector<std::string> loop_names{ "pCN production" };
    tmm::DurationEstimator my_MC_estimator( my_rank, "pCN simulation", loop_sizes, loop_names, true );

    // ====== Build the Random Generator ======
    std::mt19937 engine{ static_cast<uint>(rd::throw_seed( my_rtdata.generated_seed, my_rtdata.num_Iterations, my_rank )) };

    // ====== Initialize the Spin Cluster Correlation Functions ======
    CluCorrTen my_new_spin_correlations_Re{ my_pspace.correlation_categories, my_pspace.symmetry_type, my_pspace.num_TimePoints };
    CluCorrTen my_new_spin_correlations_Im{ my_pspace.correlation_categories, my_pspace.symmetry_type, my_pspace.num_TimePoints };
    func::SiteFields my_new_spin_expvals( my_pspace.num_Spins, FieldVector{0.,0.,0.} );

    // ====== Build the Mean-Field Correlations from the environment spin correlations ======
    auto [my_meanfield_mean, my_meanfield_correlations] = my_pspace.mf_model->self_consistency( my_spin_correlations_Re, my_spin_expvals );
    func::FrequencyCovarianceCluster my_meanfield_covariances{ my_meanfield_correlations, my_pspace.symmetry_type, my_pspace.num_Spins };

    // ====== Build the Mean-Field Distribution from the Correlations ======
    mvgb::EigenValuesBlocks my_eig;
    mvgb::OrthogonalTransformationBlocks my_ortho;
    my_meanfield_covariances.diagonalize( my_eig, my_ortho ); // diagonalization on cluster with eigen does not work yet
    my_rtdata.process_and_check_eigenvalues( my_eig );        // truncates negative eigenvalues in-place
    my_clock.measure( "mean-field distribution", true );

    // ---------------------- Preconditioned Crank-Nicolson Simulation ----------------------
    // One pCN chain per MPI core targets pi(V) propto p(V) Z(V). From the second self-
    // consistent iteration onwards the chain is warm-started from the previous iteration's
    // final state (rotated into the new diagonal basis), skipping most of the burn-in.
    const bool warm_start = !my_full_freq_persistent.empty();
    func::PCNChainCluster chain = warm_start
        ? func::PCNChainCluster( my_pspace, my_meanfield_covariances, my_eig, my_ortho, my_meanfield_mean, my_pspace.mh_step_size, engine, my_full_freq_persistent )
        : func::PCNChainCluster( my_pspace, my_meanfield_covariances, my_eig, my_ortho, my_meanfield_mean, my_pspace.mh_step_size, engine );

    // ====== Burn-In ======
    const size_t num_burn_in = warm_start
        ? static_cast<size_t>( my_pspace.mh_warm_burn_in_frac * static_cast<RealType>( my_pspace.mh_burn_in ) )
        : my_pspace.mh_burn_in;
    for( size_t s = 0; s < num_burn_in; ++s )
    {
      bool accepted = chain.step();
      my_rtdata.mh_proposed_count++;
      if( accepted ) my_rtdata.mh_accepted_count++;
    }
    my_clock.measure( "pCN burn-in", true );

    // ====== Production Sweep ======
    // Observables are accumulated from the current chain state every step after burn-in;
    // on a rejected proposal the cached state is re-accumulated (holding-time weighting).
    // The samples are split into batch-mean blocks for the error estimate; init_blocks sets
    // the block length and close_block() finalizes one block every block_length steps.
    const size_t num_SamplesPerCore = my_rtdata.get_num_SamplesPerCore();
    my_rtdata.init_blocks( my_new_spin_correlations_Re, num_SamplesPerCore );
    my_clock.enter_loop();
    for( size_t s = 0; s < num_SamplesPerCore; ++s )
    {
      bool accepted = chain.step();
      my_rtdata.mh_proposed_count++;
      if( accepted ) my_rtdata.mh_accepted_count++;

      if( my_pspace.uncoupled_spins )
      {
        func::compute_uncoupled_spin_observables_mh( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, chain.forward_uncoupled(), chain.backward_uncoupled(), chain.Z_i_list(), my_rtdata, my_pspace );
      }
      else
      {
        func::compute_spin_observables_mh( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, chain.forward(), chain.backward(), chain.Z(), my_rtdata, my_pspace );
      }
      if( ( s + 1 ) % my_rtdata.block_length == 0 ) my_rtdata.close_block();
      my_MC_estimator.obtain( my_clock.measure( "pCN production" ) );
    }
    my_clock.leave_loop();
    my_clock.measure( "pCN simulation", true );

    // Save the chain's final state (basis-independent) to seed the next iteration's warm start.
    my_full_freq_persistent = chain.full_frequency();
    // ---------------------- End of Preconditioned Crank-Nicolson Simulation ---------------

    // ====== Share the Results of Each MPI Process ======
    func::MPI_share_results_mh( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, my_rtdata );
    my_rtdata.record_mh_acceptance();
    my_clock.measure( "MPI communication", true );

    // ====== Finalize the Spin Correlations and Compute the Iteration Error ======
    const RealType N = static_cast<RealType>( my_rtdata.get_num_Samples() );
    func::normalize_mh( my_new_spin_correlations_Re, my_new_spin_correlations_Im, N );
    for( auto& spin_expval_i : my_new_spin_expvals ){ spin_expval_i /= N; }

    // ====== Compute the Standard Error from the pCN simulation ======
    // Dispatches on errmethod: 'blocking' (batch means, default; robust to autocorrelation and
    // also fills tau_int) or 'ar1' (legacy acceptance-based single-mode factor).
    my_rtdata.compute_sample_stds( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, N, my_pspace.mh_step_size );

    // Only [0, beta/2] was sampled; fill the upper half (and its std errors / tau_int) by the
    // reflection symmetry g^{ab}_{ij}(beta-tau) = g^{ab}_{ij}(tau)^* before the convergence check.
    func::symmetrize_upper_half( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_rtdata.sample_stds_Re, my_rtdata.sample_stds_Im, my_rtdata.tau_int_Re, my_rtdata.tau_int_Im );

    my_rtdata.compute_iteration_error( my_new_spin_correlations_Re, my_spin_correlations_Re );
    my_spin_correlations_Re = std::move( my_new_spin_correlations_Re );
    my_spin_correlations_Im = std::move( my_new_spin_correlations_Im );
    my_spin_expvals = std::move( my_new_spin_expvals );
    my_rtdata.finalize_iteration_step(); // increments num_Iterations and writes details to terminal
    my_clock.measure( "iteration-step finalization", true );
    print::print_R0( my_rank, "------------------------------------------------------------------\n" ) ;
  }
  while( !my_rtdata.terminate() ); 
  my_clock.leave_loop();
  my_clock.measure( "self consistency", true );
  print::print_R0( my_rank, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" );


  // ====== Store the Results and Run Time Data ======
  stc::HDF5_Storage my_data_storage( my_rank, my_pspace, my_rtdata );
  my_data_storage.store_main( my_pspace, my_rtdata, my_spin_correlations_Re, my_spin_correlations_Im, my_spin_expvals );
  
  // ====== Stop and store the Time Measurement ======
  my_clock.measure( "storing", true );
  my_clock.stop();
  my_data_storage.store_time( my_clock );

  // ====== Clean up ======
  my_data_storage.finalize();
  MPI_Finalize();
}
