#include"main_header.h"

/* Three-branch finite-temperature real-time spinDMFT.

   Edge-grid V_M and correlated V_+,V_- fields are sampled from E[V V^T].
   Independent mode uses (sum N)/(sum Z_M), optionally with paired latent
   fluctuations r,-r. Ordinary pCN targets p0(r) Re Z_M(r); antithetic pCN
   targets p0(r) Re[Z_M(r)+Z_M(-r)] and measures the sign-symmetrized pair. */

int main( const int argC, char* const argV[] )
{
  int world_size, my_rank;
  MPI_Init( nullptr, nullptr );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

  tmm::DerivedTimeMeasure my_clock( my_rank, world_size );
  const ps::ParameterSpace my_pspace( argC, argV, world_size, my_rank );
  print::print_R0( my_rank, my_pspace.create_essentials_string() );
  if( my_pspace.gaussian_factorization=="svd" )
      print::print_R0( my_rank,
          "INFO: svd uses the canonical complex-SVD Takagi factorization and "
          "samples the same complete Gaussian ensemble as dense.\n" );
  func::initialize_matrices( my_pspace );
  func::ComplexMagnetizationTrajectory my_magnetization_time(
      my_pspace.num_RealTimePoints,func::ComplexFieldVector{} );
  const FieldVector initial_spin_expval=my_pspace.initial_spin_expval;
  for( auto& magnetization:my_magnetization_time )
      for( size_t c=0;c<3;++c )
          magnetization[c]=ComplexType{initial_spin_expval[c],RealType{0.}};
  CorrelationSet my_correlations = func::generate_initial_correlations(
      my_pspace, initial_spin_expval );
  CorrelationSet my_standard_errors{ my_pspace.correlation_symmetry_type,
                                     my_pspace.num_TimePoints,
                                     my_pspace.num_RealTimePoints };

  rtd::RunTimeData my_rtdata( my_pspace, my_rank );
  my_rtdata.generated_seed = rd::generate_seed( my_rtdata.get_seed_str(), my_rank );

  std::vector<size_t> loop_sizes{ my_pspace.num_SamplesPerCore };
  const std::string sampling_loop_name=my_pspace.sampling_strategy=="pcn"
      ?my_pspace.antithetic_pairs
       ?"antithetic pCN complex-field sampling"
       :"pCN complex-field sampling"
      :my_pspace.antithetic_pairs
       ?"antithetic independent complex-field sampling"
       :"independent complex-field sampling";
  std::vector<std::string> loop_names{ sampling_loop_name };
  tmm::DurationEstimator my_MC_estimator(
      my_rank, "Monte-Carlo simulation", loop_sizes, loop_names, true );

  my_clock.measure( "initialization", true );
  print::print_R0( my_rank,
      my_pspace.uses_harmonic_bath()
      ?"++++++++++++++++++++ Prescribed Harmonic-Bath Sampling +++++++++++++++++++\n"
      :"++++++++++++++++++++++++ Self-Consistent Iteration ++++++++++++++++++++++++\n" );
  my_clock.enter_loop();
  do
  {
    my_clock.new_iteration();
    const std::string it_str = std::to_string( my_rtdata.num_Iterations + 1 );
    print::print_R0( my_rank,
        std::string{"|\n|\n------------------------ "}
        +(my_pspace.uses_harmonic_bath()?"Sampling Step ":"Iteration Step ")+it_str+
        " ------------------------\n" );

    CorrelationSet new_correlations{ my_pspace.correlation_symmetry_type,
                                     my_pspace.num_TimePoints,
                                     my_pspace.num_RealTimePoints };
    CorrelationSet new_standard_errors{ my_pspace.correlation_symmetry_type,
                                        my_pspace.num_TimePoints,
                                        my_pspace.num_RealTimePoints };
    const func::SelfConsistentField field = my_pspace.uses_harmonic_bath()
        ?func::prescribed_harmonic_bath_field(my_pspace)
        :func::self_consistent_equations(
            my_pspace,my_correlations,my_magnetization_time);
    auto sampler=func::make_complex_gaussian_sampler(
        my_pspace.gaussian_factorization,field.covariance,
        my_pspace.num_TimeSteps,my_pspace.num_RealTimePoints,
        my_pspace.delta_real_t,my_pspace.fft_cross_frequency_cutoff );
    my_rtdata.record_complex_field_diagnostics(
        field.covariance_symmetry_error, field.branch_identity_error,
        sampler->reconstruction_error(),
        sampler->covariance_approximation_error(),
        sampler->latent_dimension(),sampler->largest_factorization_dimension() );

    std::mt19937 engine{ static_cast<uint>( rd::throw_seed(
        my_rtdata.generated_seed, my_rtdata.num_Iterations, my_rank ) ) };
    my_clock.measure( "mean-field distribution", true );

    my_rtdata.begin_iteration_accumulation();
    my_clock.enter_loop();
    if( my_pspace.sampling_strategy=="pcn" )
    {
      func::PCNChain chain(
          my_pspace,*sampler,field.mean_time,my_pspace.mh_step_size,engine);
      for( size_t step=0;step<my_pspace.mh_burn_in;++step ) chain.step();
      // Burn-in is not part of a production sample.  Without resetting the
      // timer here, its full duration is attributed to the first sample and
      // multiplied by num_SamplesPerCore by the duration estimator.
      my_clock.update_time();
      for( size_t sample=0;sample<my_pspace.num_SamplesPerCore;++sample )
      {
        chain.step();
        if( my_pspace.antithetic_pairs )
          func::compute_contour_pair_correlations(
              my_rtdata,chain.trajectory(),chain.antithetic_trajectory(),
              RealType{1.}/chain.real_partition(),
              my_pspace.spin_insertion_strategy);
        else
          func::compute_contour_correlations(
              my_rtdata,chain.trajectory(),RealType{1.}/chain.real_partition(),
              my_pspace.spin_insertion_strategy);
        my_MC_estimator.obtain(my_clock.measure(sampling_loop_name));
      }
      my_rtdata.record_pcn_diagnostics(
          chain.accepted(),chain.proposed(),chain.rejected_nonpositive(),
          chain.maximum_relative_imaginary_partition());
    }
    else
    {
      const auto accumulate_field=[&]( const auto& joint_field )
      {
        const auto trajectory=func::compute_contour_trajectory(
            my_pspace,joint_field,field.mean_time);
        func::compute_contour_correlations(
            my_rtdata,trajectory,RealType{1.},
            my_pspace.spin_insertion_strategy);
        my_MC_estimator.obtain(my_clock.measure(sampling_loop_name));
      };
      if( my_pspace.antithetic_pairs )
      {
        for( size_t pair=0;pair<my_pspace.num_SamplesPerCore/2;++pair )
        {
          auto latent=sampler->draw_latent(engine);
          const auto positive_field=sampler->field_from_latent(latent);
          for( auto& value:latent ) value=-value;
          const auto negative_field=sampler->field_from_latent(latent);
          accumulate_field(positive_field);
          accumulate_field(negative_field);
        }
      }
      else
      {
        for( size_t sample=0;sample<my_pspace.num_SamplesPerCore;++sample )
          accumulate_field(sampler->draw(engine));
      }
    }
    my_clock.leave_loop();
    my_clock.measure( "Monte-Carlo simulation", true );
    my_MC_estimator.reset();
    my_rtdata.mpi_reduce_and_finalize(new_correlations,new_standard_errors);
    my_clock.measure( "MPI communication", true );

    func::ComplexMagnetizationTrajectory new_magnetization_time(
        my_pspace.num_RealTimePoints,func::ComplexFieldVector{} );
    auto new_magnetization_Re_errors=my_rtdata.magnetization_time_Re_stds;
    auto new_magnetization_Im_errors=my_rtdata.magnetization_time_Im_stds;
    for( size_t t=0;t<new_magnetization_time.size();++t )
        for( size_t c=0;c<3;++c )
            new_magnetization_time[t][c]=ComplexType{
                my_rtdata.magnetization_time_Re[t][c],
                my_rtdata.magnetization_time_Im[t][c]};
    if( my_pspace.constant_magnetization_time
        &&!my_pspace.uses_harmonic_bath() )
    {
        new_magnetization_time=func::project_constant_magnetization(
            new_magnetization_time );
        new_magnetization_Re_errors.assign(
            new_magnetization_Re_errors.size(),new_magnetization_Re_errors.front());
        new_magnetization_Im_errors.assign(
            new_magnetization_Im_errors.size(),new_magnetization_Im_errors.front());
    }

    const auto residual=func::iteration_residual(
        my_correlations,new_correlations,new_standard_errors,
        my_magnetization_time,new_magnetization_time,
        new_magnetization_Re_errors,new_magnetization_Im_errors );
    my_rtdata.record_iteration_error(residual.absolute,residual.standardized);

    if( my_pspace.uses_harmonic_bath() )
    {
        // There is no fixed-point update in prescribed-bath mode.  Store the
        // measured estimator directly instead of mixing it with an arbitrary
        // spinDMFT initial guess.
        my_correlations=std::move(new_correlations);
        my_magnetization_time=std::move(new_magnetization_time);
    }
    else
    {
        my_correlations=func::mix_correlations(
            my_correlations,new_correlations,my_pspace.mixing_alpha );
        my_magnetization_time=func::mix_magnetization_trajectory(
            my_magnetization_time,new_magnetization_time,my_pspace.mixing_alpha );
    }
    my_standard_errors = std::move(new_standard_errors);
    my_rtdata.finalize_iteration_step();
    my_clock.measure( "iteration-step finalization", true );
    print::print_R0( my_rank,
        "------------------------------------------------------------------\n" );
  }
  while( !my_rtdata.terminate() );
  my_clock.leave_loop();
  my_clock.measure( "self consistency", true );
  print::print_R0( my_rank,
      "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" );

  stoc::HDF5_Storage my_data_storage(
      my_rank, my_pspace, my_rtdata.termination );
  my_data_storage.store_main(
      my_pspace, my_rtdata, my_correlations, my_standard_errors );

  my_clock.measure( "storing", true );
  my_clock.stop();
  my_data_storage.store_time( my_clock );
  my_data_storage.finalize();
  MPI_Finalize();
}
