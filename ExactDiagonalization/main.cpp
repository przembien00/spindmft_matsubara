#include"main_header.h"

int main( const int argC, char* const argV[] ){ // arguments required for boost program options

    // ====== Initialize MPI ======
    // world_size is the number of cores, my_rank is the number of "this" core 
    int world_size, my_rank;
    MPI_Init( nullptr, nullptr );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    // ====== Import, Interpret and Print Parameters ======
    const ps::ParameterSpace my_pspace( argC, argV, world_size, my_rank );
    print::print_R0( my_rank, my_pspace.create_essentials_string() );
    rtd::RunTimeData my_rtdata( my_pspace ); // initialize run time data

    // ====== Initialize Time Measurement ======
    tmm::SpinEDTimeMeasure my_clock( my_rank, world_size ); // initialize clock
    tmm::SimpleDurationEstimator my_estimator( my_rank, "time evolution", my_pspace.num_TimePoints, true ); // initialize duration estimator

    // ====== Construct the System ======
    my_pspace.spin_model->compute_Hamiltonian(); // Set up the Hamiltonian
    my_pspace.spin_model->compute_density_operator(); // Set up rho
    my_clock.measure( "construction", true );
  
    // ====== Diagonalize the Hamiltonian ======
    my_pspace.spin_model->diagonalize();
    my_clock.measure( "diagonalization", true );

    // ====== Transform Relevant Operators to the Diagonal Basis ======
    my_pspace.spin_model->transform_to_diagonal_basis();
    my_clock.measure( "transformation", true );

    // ====== Compute the Spin Autocorrelations ======
    my_clock.enter_loop();
    for( uint t_point = 0; t_point < my_pspace.num_TimePoints; ++t_point )
    {
        RealType time = my_pspace.delta_t * static_cast<RealType>( t_point );
        my_pspace.spin_model->compute_means_and_autocorrelations_at( time );
        my_estimator.obtain( my_clock.measure( "quantities_at" ) );
    }
    my_clock.leave_loop();
    my_clock.measure( "autocorrelations", true );

/*
    // ====== Share the Results of Each MPI Process ======
    self::MPI_share_results( my_rtdata, my_new_spin_meanvalue_holder, my_new_spin_correlation_holder ); // this could be done more efficiently
    my_clock.measure( "communication", true );
*/

    // ====== Save the Data ======
    stoc::HDF5_Storage my_data_storage( my_rank, my_pspace );
    my_data_storage.store_main( my_pspace, my_rtdata, my_pspace.spin_model->m_means, my_pspace.spin_model->m_correlations );

    // ====== Stop the Time Measurement ======
    my_clock.measure( "saving", true );
    my_clock.stop();
      
    // ====== Save the Time Measurement ======
    my_data_storage.store_time( my_clock );

    // ====== Clean up ======
    my_data_storage.finalize();
    MPI_Finalize();

}
