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

  // ====== Initialize Hilbert Space Matrices ======
  func::initialize_matrices( my_pspace );

  // ====== Initialize Spin Correlations ======
  FieldVector my_spin_expval = my_pspace.initial_spin_expval;
  auto my_spin_correlations_Re = func::generate_initial_spin_correlations( my_pspace );
  auto my_spin_correlations_Im = func::generate_initial_spin_correlations( my_pspace );

  // ====== Initialize Run Time Data ======
  rtd::RunTimeData my_rtdata( my_pspace, my_rank );
  my_rtdata.generated_seed = rd::generate_seed( my_rtdata.get_seed_str(), my_rank );

  // ====== Initialize Duration Estimator ======
  std::vector<size_t> loop_sizes{ my_pspace.num_SetsPerCore, my_pspace.num_SamplesPerCore };
  std::vector<std::string> loop_names{ "mean-field sampling", "time evolution" };
  tmm::DurationEstimator my_MC_estimator( my_rank, "Monte-Carlo simulation", loop_sizes, loop_names, true );


  my_clock.measure( "initialization", true );
  print::print_R0( my_rank, "++++++++++++++++++++++++ Self-Consistent Iteration ++++++++++++++++++++++++\n" );
  my_clock.enter_loop();
  do
  {
    my_clock.new_iteration();
    std::string it_str = std::to_string(my_rtdata.num_Iterations + 1);
    print::print_R0( my_rank, "|\n|\n------------------------ Iteration Step " + it_str + " ------------------------\n" );
  
    // ====== Initialize new Spin Correlations ======
    RealType my_partition_function = RealType{0.};
    RealType my_partition_function_sq = RealType{0.};
    FieldVector my_new_spin_expval{0.,0.,0.}; // needs to be made time dependent when H is time dependent
    CorrTen my_new_spin_correlations_Re{ my_pspace.correlation_symmetry_type, my_pspace.num_TimePoints };
    CorrTen my_new_spin_correlations_Im{ my_pspace.correlation_symmetry_type, my_pspace.num_TimePoints };
    // ====== Determine the Mean-Field Moments Self-Consistently ======
    auto [my_meanfield_mean, my_meanfield_covariances] = func::self_consistent_equations( my_pspace, my_spin_correlations_Re, my_spin_expval );
    
    // ====== Build the Mean-Field Distribution from the Correlations ======
    fmvg::EigenValuesList my_eig;
    fmvg::OrthogonalTransformationList my_ortho;
    my_meanfield_covariances.diagonalize( my_eig, my_ortho );
    my_rtdata.process_and_check_eigenvalues( my_eig );
    fmvg::NormalDistributionsList my_dist{ my_eig };

    // ====== Build the Random Generator ======
    std::mt19937 engine{ static_cast<uint>( rd::throw_seed( my_rtdata.generated_seed, my_rtdata.num_Iterations, my_rank ) ) };
    my_clock.measure( "mean-field distribution", true );

    // ---------------------------------------------- Monte-Carlo Simulation ------------------------------------------------------- 
    size_t remaining_samples{ my_rtdata.get_num_SamplesPerCore() };
    my_clock.enter_loop();
    while( remaining_samples > 0 ) // efficiency gain through cache coherence : the noise is sampled in sets 
    {
      size_t num_SamplesInThisSet = std::min( remaining_samples, my_pspace.num_SamplesPerSet );
      remaining_samples = (size_t) std::max( (int)0, (int)remaining_samples-(int)my_pspace.num_SamplesPerSet );

      // ====== Sample the Mean-Field Noise ======
      fmvg::FrequencyNoiseVectors my_noise_vectors{ my_dist, my_ortho, engine, num_SamplesInThisSet };

      my_MC_estimator.obtain( my_clock.measure( "mean-field sampling" ) );
      my_clock.enter_loop();
      for( size_t sample = 0; sample < num_SamplesInThisSet; ++sample )
      {
        // ====== Compute Time Propagators ======
        auto sample_noise = my_noise_vectors( sample );
        auto [Z, TimePropagators, TimePropagators_inv] = func::compute_propagators( my_pspace, sample_noise, my_meanfield_mean );

        // ====== Compute Time-Evolved Spin Operators ======
        std::vector<Operator> S_x_of_t{}, S_y_of_t{}, S_z_of_t{}; 
        func::compute_S_of_t( my_pspace, TimePropagators, TimePropagators_inv, S_x_of_t, S_y_of_t, S_z_of_t ); 
        auto propagation_duration = my_clock.measure( "time propagation" );

        // ====== Compute Spin MeanValues and Correlations ======
        my_partition_function += Z;
        my_rtdata.Z_sqsum += Z*Z;
        func::compute_spin_correlations( my_pspace, my_rtdata, my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expval, Z, S_x_of_t, S_y_of_t, S_z_of_t );
        
        auto expectationvalues_duration = my_clock.measure( "expectation values" );
        my_MC_estimator.obtain( tmm::DurationQuantity{"time evolution", propagation_duration.m_duration + expectationvalues_duration.m_duration} );
      }
      my_clock.leave_loop();
      my_clock.update_time();
    }
    my_clock.leave_loop();
    my_clock.measure( "Monte-Carlo simulation", true );
    my_MC_estimator.reset();
    // ----------------------------------------------------------------------------------------------------------------------------- 

    // ====== Share the Results of Each MPI Process ======
    func::MPI_share_results( my_rtdata, my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expval, my_partition_function ); // this could be done more efficiently

    my_clock.measure( "MPI communication", true );

    // ====== Finalize the Iteration ======
    func::normalize( my_rtdata, my_new_spin_correlations_Re, my_partition_function ); // normalize with respect to the sample size
    func::normalize( my_rtdata, my_new_spin_correlations_Im, my_partition_function ); // normalize with respect to the sample size
    my_new_spin_expval /= my_partition_function;
    my_rtdata.compute_sample_stds( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_partition_function ); // statistical error of the Monte-Carlo simulation
    my_rtdata.compute_iteration_error( my_spin_correlations_Re, my_new_spin_correlations_Re );
    my_spin_correlations_Re = std::move( my_new_spin_correlations_Re );
    my_spin_correlations_Im = std::move( my_new_spin_correlations_Im );
    my_spin_expval = std::move(my_new_spin_expval);
    my_rtdata.finalize_iteration_step();
    my_clock.measure( "iteration-step finalization", true );
    print::print_R0( my_rank, "------------------------------------------------------------------\n" ) ;
  }
  while( !my_rtdata.terminate() );
  my_clock.leave_loop();
  my_clock.measure( "self consistency", true );
  print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" );


  // ====== Store the Data ======
  stoc::HDF5_Storage my_data_storage( my_rank, my_pspace, my_rtdata.termination );
  my_data_storage.store_main( my_pspace, my_rtdata, my_spin_correlations_Re, my_spin_correlations_Im, my_spin_expval ); // store mag moments as well

  // ====== Stop and Store the Time Measurement ======
  my_clock.measure( "storing", true );
  my_clock.stop();
  my_data_storage.store_time( my_clock );

  // ====== Clean up ======
  my_data_storage.finalize();
  MPI_Finalize();
}
