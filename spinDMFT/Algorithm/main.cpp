#include"main_header.h"

/* spinDMFT with a preconditioned Crank-Nicolson (pCN) sampler in the diagonal
   Matsubara frequency basis.

   Instead of importance sampling V ~ p(V) and reweighting by Z(V), we run a pCN
   chain per MPI core with target density pi(V) propto p(V) Z(V):

   - State lives in the diagonal basis of the per-block covariance matrix, where p(V)
     factorizes into independent N(0, eig[block][i]) Gaussians.
   - pCN proposal: V' = sqrt(1 - beta^2) V + beta xi, with xi[block][i] drawn fresh
     from the prior N(0, sqrt(eig[block][i])). The proposal preserves p(V) exactly,
     so the Gaussian factors cancel in the M-H ratio.
   - Acceptance: min(1, Z(V')/Z(V)). Z <= 0 is treated as a forbidden region.
   - For Gaussian-prior * positive-likelihood targets pCN has dimension-independent
     autocorrelation, in contrast to random-walk M-H which scales as O(d).
   - The per-sample observable is Tr[U(beta-t) S U(t) S] / Z(V), divided sample-by-
     sample inside compute_spin_correlations. Observables are accumulated from the
     current chain state every step after burn-in. */

int main( const int argC, char* const argV[] ){ // arguments required for boost program options
  // ====== Initialize MPI ======
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
  std::vector<size_t> loop_sizes{ my_pspace.num_SamplesPerCore };
  std::vector<std::string> loop_names{ "mean-field sampling" };
  tmm::DurationEstimator my_MC_estimator( my_rank, "Monte-Carlo simulation", loop_sizes, loop_names, true );

  // Persistent chain state across self-consistent iterations, expressed in the basis-
  // independent Matsubara representation V_full = ortho * V_diag. Empty before the
  // first iteration; carried forward thereafter to skip the cold-start burn-in.
  std::vector<fmvg::Vector> my_V_full_persistent;

  my_clock.measure( "initialization", true );
  print::print_R0( my_rank, "++++++++++++++++++++++++ Self-Consistent Iteration ++++++++++++++++++++++++\n" );
  my_clock.enter_loop();
  do
  {
    my_clock.new_iteration();
    std::string it_str = std::to_string(my_rtdata.num_Iterations + 1);
    print::print_R0( my_rank, "|\n|\n------------------------ Iteration Step " + it_str + " ------------------------\n" );

    // ====== Initialize new Spin Correlations ======
    FieldVector my_new_spin_expval{0.,0.,0.};
    FieldVector my_spin_expval_sqsum{0.,0.,0.};
    CorrTen my_new_spin_correlations_Re{ my_pspace.correlation_symmetry_type, my_pspace.num_TimePoints };
    CorrTen my_new_spin_correlations_Im{ my_pspace.correlation_symmetry_type, my_pspace.num_TimePoints };

    // ====== Determine the Mean-Field Moments Self-Consistently ======
    auto [my_meanfield_mean, my_meanfield_covariances] = func::self_consistent_equations( my_pspace, my_spin_correlations_Re, my_spin_expval );

    // ====== Build the Mean-Field Distribution from the Correlations ======
    fmvg::EigenValuesList my_eig;
    fmvg::OrthogonalTransformationList my_ortho;
    my_meanfield_covariances.diagonalize( my_eig, my_ortho );
    my_rtdata.process_and_check_eigenvalues( my_eig ); // truncates negative eigenvalues to zero in-place

    // ====== Build the Random Generator ======
    std::mt19937 engine{ static_cast<uint>( rd::throw_seed( my_rtdata.generated_seed, my_rtdata.num_Iterations, my_rank ) ) };
    my_clock.measure( "mean-field distribution", true );

    // ---------------------------------------------- Preconditioned Crank-Nicolson Simulation -------------------------------------
    // The pCN chain owns its state (V in the diagonal Matsubara basis, Z, cached
    // propagators and S^a(t)) and exposes a single step() method. From the second
    // self-consistent iteration onwards we warm-start the chain from the previous
    // iteration's final V_full, rotated into the new diagonal basis.
    const bool warm_start = !my_V_full_persistent.empty();
    std::unique_ptr<func::PCNChain> chain_ptr;
    if( warm_start )
    {
      std::vector<fmvg::Vector> V_diag_init( my_ortho.size() );
      for( size_t b = 0; b < my_ortho.size(); ++b )
      {
        V_diag_init[b] = blaze::trans( my_ortho[b] ) * my_V_full_persistent[b];
      }
      chain_ptr = std::make_unique<func::PCNChain>(
          my_pspace, my_eig, my_ortho, my_meanfield_mean,
          my_pspace.mh_step_size, engine, std::move( V_diag_init ) );
    }
    else
    {
      chain_ptr = std::make_unique<func::PCNChain>(
          my_pspace, my_eig, my_ortho, my_meanfield_mean,
          my_pspace.mh_step_size, engine );
    }
    func::PCNChain& chain = *chain_ptr;

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

    // ====== Production Sweep ======
    // compute_spin_correlations accumulates the per-sample correlation
    //   Tr[U(beta-t) S U(t) S] / Z(V)
    // (no IS reweighting is needed because V is drawn from pi propto p(V) Z(V)).
    my_clock.enter_loop();
    for( size_t s = 0; s < my_pspace.num_SamplesPerCore; ++s )
    {
      bool accepted = chain.step();
      my_rtdata.mh_proposed_count++;
      if( accepted ) my_rtdata.mh_accepted_count++;

      // On a rejected proposal the cached S^a(t) is re-accumulated, which is exactly
      // the contribution of the current chain state weighted by its holding time.
      func::compute_spin_correlations( my_rtdata,
          my_new_spin_correlations_Re, my_new_spin_correlations_Im,
          my_new_spin_expval, my_spin_expval_sqsum, chain.Z(),
          chain.S_x(), chain.S_y(), chain.S_z() );
      my_MC_estimator.obtain( my_clock.measure( "mean-field sampling" ) );
    }
    my_clock.leave_loop();
    my_clock.measure( "Monte-Carlo simulation", true );
    my_MC_estimator.reset();

    // Save the chain's final state in the basis-independent Matsubara representation,
    // ready to seed the next self-consistent iteration's warm start.
    my_V_full_persistent = chain.V_full();
    // -------------------------------------------------------------------------------------------------------------------------------

    // ====== Share the Results of Each MPI Process ======
    func::MPI_share_results( my_rtdata, my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expval, my_spin_expval_sqsum );
    my_rtdata.record_mh_acceptance();
    my_clock.measure( "MPI communication", true );

    // ====== Finalize the Iteration ======
    func::normalize( my_new_spin_correlations_Re, my_pspace.num_Samples );
    func::normalize( my_new_spin_correlations_Im, my_pspace.num_Samples );
    my_new_spin_expval /= my_pspace.num_Samples;
    my_rtdata.compute_sample_stds( my_new_spin_correlations_Re, my_new_spin_correlations_Im, my_new_spin_expval, my_pspace.num_Samples, my_pspace.mh_step_size );
    my_rtdata.compute_iteration_error( my_spin_correlations_Re, my_new_spin_correlations_Re );
    my_spin_correlations_Re = std::move( my_new_spin_correlations_Re );
    my_spin_correlations_Im = std::move( my_new_spin_correlations_Im );
    my_spin_expval = std::move( my_new_spin_expval );
    my_rtdata.finalize_iteration_step();
    my_clock.measure( "iteration-step finalization", true );
    print::print_R0( my_rank, "------------------------------------------------------------------\n" );
  }
  while( !my_rtdata.terminate() );
  my_clock.leave_loop();
  my_clock.measure( "self consistency", true );
  print::print_R0( my_rank, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n" );


  // ====== Store the Data ======
  stoc::HDF5_Storage my_data_storage( my_rank, my_pspace, my_rtdata.termination );
  my_data_storage.store_main( my_pspace, my_rtdata, my_spin_correlations_Re, my_spin_correlations_Im, my_spin_expval );

  // ====== Stop and Store the Time Measurement ======
  my_clock.measure( "storing", true );
  my_clock.stop();
  my_data_storage.store_time( my_clock );

  // ====== Clean up ======
  my_data_storage.finalize();
  MPI_Finalize();
}
