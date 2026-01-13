#include"main_header.h"

int main(const int argC, char* const argV[]){ // arguments required for boost program options
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

  // ====== Initialize Duration Estimator ======
  std::vector<size_t> loop_sizes{ my_pspace.num_SetsPerCore, my_pspace.num_SamplesPerCore };
  std::vector<std::string> loop_names{ "mean-field sampling", "time evolution" };
  tmm::DurationEstimator my_MC_estimator( my_rank, "Monte-Carlo simulation", loop_sizes, loop_names, true );

  // ====== Initialize several Constant Spin Matrices on the Hilbert Space ======
  init::write_general_spin_matrices( my_pspace.spin_float_list, my_pspace.num_HilbertSpaceDimension );
  init::write_cluster_Hamiltonian( my_pspace.num_Spins, my_pspace.num_HilbertSpaceDimension, my_pspace.J, my_pspace.spinspin_cmodel, my_pspace.chemical_shift, my_pspace.local_extra_interaction );

  // ====== Declare the Spin Cluster Correlation Functions ======
  std::vector<CorrTen> my_environment_spin_correlations = init::generate_environment_spin_correlations( my_pspace );
  my_clock.measure( "initialization", true );
  
  // ====== Build the Random Generator ======
  std::mt19937 engine{rd::throw_seed( my_rtdata.generated_seed, 0, my_rank )};

  // ====== Initialize the Spin Cluster Correlation Functions ======
  CluCorrTen my_new_spin_correlations{ my_pspace.compute_only_categories, my_pspace.symmetry_type, my_pspace.num_TimePoints };

  // ====== Build the Mean-Field Correlations from the environment spin correlations ====== 
  auto my_meanfield_covariances = my_pspace.mf_model->self_consistency( my_environment_spin_correlations );

  // ====== Build the Mean-Field Distribution from the Correlations ======
  mvgb::EigenValuesBlocks my_eig;
  mvgb::OrthogonalTransformationBlocks my_ortho;
  my_meanfield_covariances.diagonalize( my_eig, my_ortho );
  my_rtdata.process_and_check_eigenvalues( my_eig );
  mvgb::DiagonalBasisNormalDistributionsBlocks my_dist{ my_eig };
  my_clock.measure( "mean-field distribution", true );


  // -------------------------- Monte-Carlo Simulation -------------------------
  size_t remaining_samples{ my_rtdata.get_num_SamplesPerCore() };
  my_clock.enter_loop();
  while( remaining_samples > 0 ) // efficiency gain through cache coherence : the noise is sampled in sets beforehand 
  {
    size_t num_SamplesInThisSet = std::min( remaining_samples, my_pspace.num_SamplesPerSet );
    remaining_samples = (size_t) std::max( (int)0, (int)remaining_samples-(int)my_pspace.num_SamplesPerSet );

    // ====== Sample the Mean-Field Noise ======
    auto num_noises_per_block = mss::return_num_noises_per_block( my_pspace.symmetry_type, num_SamplesInThisSet );      
    auto noise_ptr = std::make_shared<const mvgb::GaussianNoiseVectorsBlocks>(num_noises_per_block, my_dist, engine, my_ortho);
    mss::NoiseTensorClusterReaderScheme scheme{ noise_ptr, my_pspace.symmetry_type, num_SamplesInThisSet, my_pspace.num_TimePoints, my_pspace.num_Spins }; // scheme for reading the noise

    my_MC_estimator.obtain( my_clock.measure( "mean-field sampling" ) );
    my_clock.enter_loop();
    for( size_t sample = 0; sample < num_SamplesInThisSet; ++sample )
    {
      // ====== Initialize Time Evolution Operator and Hamiltonian ======
      Operator TimeEvolutionOperator{ IDENTITY };
      SparseObservable old_Hamiltonian{}, new_Hamiltonian{};
      my_pspace.mf_model->compute_Hamiltonian( old_Hamiltonian, scheme.read(), my_pspace.symmetry_type );

      // ====== Propagate ======
      for( size_t time = 1; time < my_pspace.num_TimePoints; ++time )
      {
        my_pspace.mf_model->compute_Hamiltonian( new_Hamiltonian, scheme.read(), my_pspace.symmetry_type );
        func::propagate( TimeEvolutionOperator, old_Hamiltonian, new_Hamiltonian, my_pspace );
        old_Hamiltonian = new_Hamiltonian;
        func::compute_spin_correlations( my_new_spin_correlations, TimeEvolutionOperator, time, my_rtdata, my_pspace );
      }
      my_MC_estimator.obtain( my_clock.measure( "time evolution" ) );
    }
    my_clock.leave_loop();
    my_clock.update_time();
  }
  my_clock.leave_loop();
  my_clock.measure( "Monte-Carlo simulation", true );
  // -----------------------End of Monte-Carlo Simulation ----------------------


  // ====== Share the Results of Each MPI Process ======
  func::MPI_share_results( my_new_spin_correlations, my_rtdata );
  my_clock.measure( "MPI communication", true );

  // ====== Compute the Standard Deviation from the MC-simulation ======
  my_rtdata.compute_sample_stds( my_new_spin_correlations );

  // ====== Finalize the Spin Correlations ======
  func::normalize_and_fill_time_zero( my_new_spin_correlations, my_rtdata.get_num_Samples(), my_pspace.spin_float_list );
  my_clock.measure( "finalization", true );

  // ====== Store the Results and Run Time Data ======
  stc::HDF5_Storage my_data_storage( my_rank, my_pspace, my_rtdata );
  my_data_storage.store_main( my_pspace, my_rtdata, my_new_spin_correlations );
  
  // ====== Stop and store the Time Measurement ======
  my_clock.measure( "storing", true );
  my_clock.stop();
  my_data_storage.store_time( my_clock );

  // ====== Clean up ======
  my_data_storage.finalize();
  MPI_Finalize();
}
