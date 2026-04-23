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
  init::write_general_spin_matrices( my_pspace.spin_float_list, my_pspace.num_HilbertSpaceDimension );
  init::write_cluster_Hamiltonian( my_pspace.num_Spins, my_pspace.num_HilbertSpaceDimension, my_pspace.J, my_pspace.spinspin_cmodel, my_pspace.chemical_shift, my_pspace.local_extra_interaction );

  // ====== Declare the Spin Cluster Correlation Functions ======
  CluCorrTen my_spin_correlations_Re = init::generate_initial_environment_spin_correlations( my_pspace, "Re" );
  CluCorrTen my_spin_correlations_Im = init::generate_initial_environment_spin_correlations( my_pspace, "Im" );
  func::SiteFields my_spin_expvals = ( my_pspace.init_corr_mode == "import" ) ? my_pspace.imported_spin_expvals : func::SiteFields( my_pspace.num_Spins, FieldVector{0.,0.,0.} );
  my_clock.measure( "initialization", true );


  print::print_R0( my_rank, "+++++++++++++++++++++++++ Self-Consistent Iteration ++++++++++++++++++++++++++\n" );
  my_clock.enter_loop();
  do
  {
    my_clock.new_iteration();
    std::string it_str = std::to_string( my_rtdata.num_Iterations + 1 );
    print::print_R0( my_rank, "|\n|\n--------------------------- Iteration Step " + it_str + " ---------------------------\n" );

    // ====== Initialize Duration Estimator ======
    std::vector<size_t> loop_sizes{ my_rtdata.get_num_SetsPerCore(), my_rtdata.get_num_SamplesPerCore() };
    std::vector<std::string> loop_names{ "mean-field sampling", "time evolution" };
    tmm::DurationEstimator my_MC_estimator( my_rank, "Monte-Carlo simulation", loop_sizes, loop_names, true );

    // ====== Build the Random Generator ======
    std::mt19937 engine{ static_cast<uint>(rd::throw_seed( my_rtdata.generated_seed, my_rtdata.num_Iterations, my_rank )) };

    // ====== Initialize the Spin Cluster Correlation Functions ======
    CluCorrTen my_new_spin_correlations_Re{ my_pspace.correlation_categories, my_pspace.symmetry_type, my_pspace.num_TimePoints };
    CluCorrTen my_new_spin_correlations_Im{ my_pspace.correlation_categories, my_pspace.symmetry_type, my_pspace.num_TimePoints };
    func::SiteFields my_new_spin_expvals( my_pspace.num_Spins, FieldVector{0.,0.,0.} );
    RealType my_partition_function = RealType{0.};

    // ====== Build the Mean-Field Correlations from the environment spin correlations ====== 
    auto [my_meanfield_mean, my_meanfield_correlations] = my_pspace.mf_model->self_consistency( my_spin_correlations_Re, my_spin_expvals );
    func::FrequencyCovarianceCluster my_meanfield_covariances{ my_meanfield_correlations, my_pspace.symmetry_type, my_pspace.num_Spins };

    // ====== Build the Mean-Field Distribution from the Correlations ======
    mvgb::EigenValuesBlocks my_eig;
    mvgb::OrthogonalTransformationBlocks my_ortho;
    print::print_R0( my_rank, "Diagonalizing the mean-field covariance matrix..." );
    my_meanfield_covariances.diagonalize( my_eig, my_ortho ); // diagonalization on cluster with eigen does not work yet
    print::print_R0( my_rank, "Diagonalization finished." );
    my_rtdata.process_and_check_eigenvalues( my_eig );
    print::print_R0( my_rank, "Eigenvalues processed." );
    mvgb::DiagonalBasisNormalDistributionsBlocks my_dist{ my_eig };
    print::print_R0( my_rank, "Mean-field distribution initialized." );
    my_clock.measure( "mean-field distribution", true );

    // -------------------------- Monte-Carlo Simulation -------------------------
    size_t remaining_samples{ my_rtdata.get_num_SamplesPerCore() };
    my_clock.enter_loop();
    while( remaining_samples > 0 ) // efficiency gain through cache coherence : the noise is sampled in sets beforehand 
    {
      size_t num_SamplesInThisSet = std::min( remaining_samples, my_pspace.num_SamplesPerSet );
      remaining_samples = (size_t) std::max( (int)0, (int)remaining_samples-(int)my_pspace.num_SamplesPerSet );

      // ====== Sample the Mean-Field Noise in Matsubara/Frequency Space ======
      auto meanfield_trajectories = my_meanfield_covariances.sample_time_trajectories( my_eig, my_ortho, engine, num_SamplesInThisSet );

      my_MC_estimator.obtain( my_clock.measure( "mean-field sampling" ) );
      my_clock.enter_loop();
      for( size_t sample = 0; sample < num_SamplesInThisSet; ++sample )
      {
        // ====== Propagate in imaginary time and sample thermal correlations ======
        auto [Z, forward_propagators, backward_propagators] = func::compute_propagators( meanfield_trajectories[sample], my_meanfield_mean, my_pspace );
        my_partition_function += Z;
        func::compute_spin_observables( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, forward_propagators, backward_propagators, Z, my_rtdata, my_pspace );
        my_MC_estimator.obtain( my_clock.measure( "time evolution" ) );
      }
      my_clock.leave_loop();
      my_clock.update_time();
    }
    my_clock.leave_loop();
    my_clock.measure( "Monte-Carlo simulation", true );
    // -----------------------End of Monte-Carlo Simulation ----------------------

    // ====== Share the Results of Each MPI Process ======
    func::MPI_share_results( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expvals, my_rtdata, my_partition_function );
    my_clock.measure( "MPI communication", true );

    // ====== Finalize the Spin Correlations and Compute the Iteration Error ======
    func::normalize( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_partition_function );
    for( auto& spin_expval_i : my_new_spin_expvals ){ spin_expval_i /= my_partition_function; }

    // ====== Compute the Standard Deviation from the MC-simulation ======
    my_rtdata.compute_and_process_spin_expval_stds( my_new_spin_expvals, my_partition_function );
    my_rtdata.compute_and_process_sample_stds( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_partition_function );

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
