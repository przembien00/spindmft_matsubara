#include"Functions.h"
#include<map>
#include<random>
#include<functional>
#include<iostream>
#include<cmath>

#include<Globals/MPI_Types.h>
#include<Standard_Algorithms/Standard_Algorithms.h>
#include<Standard_Algorithms/Numerics.h>
namespace stda = Standard_Algorithms;
namespace num = Standard_Algorithms::Numerics;

#include"Error_Handling.h"
namespace error = Functions::Error_Handling;

#include<Matrices/Diagonalization.h>
namespace diag = Matrices::Diagonalization;

#include<Physics/Spin.h>
namespace sp = Physics::Spin;

#include"../Mean_Field_Models/Mean_Field_Models.h"


// return spin operator at specific site and in specific direction
inline const SparseObservable& S_at_site_in_direction( const size_t i, const size_t alpha )
{
    if( alpha == 0 ){ return S_X_LIST[i]; }
    else if( alpha == 1 ){ return S_Y_LIST[i]; }
    else{ return S_Z_LIST[i]; }
}


namespace Functions::Initialization
{

// initialization function : initialize ZERO, IDENTITY, and S_X_LIST, S_Y_LIST, S_Z_LIST for all sites
void write_general_spin_matrices( const std::vector<RealType>& spin_float_list, const size_t num_HilbertSpaceDimension )
{
  // ZERO, IDENTITY:
  std::tie(ZERO, IDENTITY) = sp::create_zero_and_identity( num_HilbertSpaceDimension );

  // S_X_LIST, S_Y_LIST, S_Z_LIST:
  size_t index = 0;
  for( auto& spin : spin_float_list )
  {
    Observable S_X{}, S_Y{}, S_Z{};
    sp::write_spin_matrices( spin, S_X, S_Y, S_Z );
    S_X_LIST.emplace_back( SparseObservable{sp::create_linear_spin_term( spin_float_list, index, S_X )} );
    S_Y_LIST.emplace_back( SparseObservable{sp::create_linear_spin_term( spin_float_list, index, S_Y )} );
    S_Z_LIST.emplace_back( SparseObservable{sp::create_linear_spin_term( spin_float_list, index, S_Z )} );
    index++;
  }
}

// initialization function : initialize H_CLUSTER
void write_cluster_Hamiltonian( const size_t num_Spins, const size_t num_HilbertSpaceDimension, const SymmMatrix& J, const ph::SpinModel& spinspin_cmodel, const ph::ChemicalShift& chemical_shift, const ph::LocalExtraInteraction& local_extra_interaction )
{
    Observable H_CLUSTER_DENSE{ZERO}; // create Hamiltonian as dense matrix

    // add local spin terms:
    for(size_t spin_i = 0; spin_i < num_Spins; spin_i++)
    { 
      // chemical shifts:
      if( chemical_shift.m_name != "none" )
      {
        H_CLUSTER_DENSE += chemical_shift.at(spin_i) * S_Z_LIST[spin_i]; // H_local += h * SIGMA^Z_i
      }

      // extra interactions:
      if( local_extra_interaction.m_name != "none" )
      {
        H_CLUSTER_DENSE += local_extra_interaction.m_term(spin_i,ZERO,S_X_LIST[spin_i],S_Y_LIST[spin_i],S_Z_LIST[spin_i]);
      }
    }

    // add spin-spin interaction terms (spin-spin interactions), H_local += sum_ij 
    auto D = spinspin_cmodel.coupling_matrix;
    for(size_t spin_i = 0; spin_i < num_Spins-1; spin_i++)
    {
      for(size_t spin_j = spin_i+1; spin_j < num_Spins; spin_j++) // spin_i is always in front of spin_j
      {
        // loop over all direction pairs to compute and add h_ij = J_ij * vS_i^T * D * vS_j
        for(size_t alpha = 0; alpha < 3; alpha++) 
        {
          for(size_t beta = 0; beta < 3; beta++) // loop over the three Pauli matrices
          {
            // add h_ij^ab = J_ij * D^ab S^alpha_i * S^beta_j
            H_CLUSTER_DENSE += J(spin_i,spin_j)*D(alpha,beta) * S_at_site_in_direction(spin_i,alpha) * S_at_site_in_direction(spin_j,beta); 
          }
        }
      }
    }

    // make it sparse
    H_CLUSTER = SparseObservable{H_CLUSTER_DENSE};
}

// initialization function : initialize the spin correlations
CluCorrTen generate_initial_environment_spin_correlations( const ps::ParameterSpace& pspace, const std::string& component )
{
  // reserve
  CluCorrTen CCT{pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints};

  if( pspace.init_corr_mode == "generate" ) // set the spin correlations to some pregiven functions
  {
    // create the zero correlation:
    corr::CorrelationVector ZC{  ph::NonDiagonalSpinCorrelation{"zero"}.create_discretization( pspace.delta_t, pspace.num_TimePoints, 0. ) };

    // insert them into the spin correlations
    CCT.iterate( [&]( CorrTen& CT, const IndexPair& ij )
    {
      if( ij[0] == ij[1] ) // for spin autocorrelations set the diagonal correlations to DC and nondiagonal correlations to NDC
      {
        // create the DC and NDC correlation:
        auto spin = pspace.spin_float_list[ij[0]]; 
        corr::CorrelationVector DC{  pspace.init_diag_corr.create_discretization( pspace.delta_t, pspace.num_TimePoints, spin ) };
        corr::CorrelationVector NDC{ pspace.init_nondiag_corr.create_discretization( pspace.delta_t, pspace.num_TimePoints, spin ) };

        // write them to the correlation tensor:
        CT.iterate( [&]( corr::CorrelationVector& C, const IndexPair& alphabeta )
        {
          C = (alphabeta[0]==alphabeta[1]) ? DC : NDC;
        } );
      }
      else // set any spin paircorrelations to ZC (zero)
      {
          std::fill( CT.begin(), CT.end(), ZC );
      }
    } );
  }
  else // interpret the imported correlations (they are still in semi-linearized form), perhaps extrapolate them if the time discretization doesn't match
  {
    const ps::CluList& imported_correlations_linearized = ( component == "Im" ) ? pspace.imported_correlations_Im_linearized : pspace.imported_correlations_Re_linearized;
    stda::for_2each( CCT.begin(), CCT.end(), imported_correlations_linearized.cbegin(), [&]( CorrTen& CT, const std::vector<RealType>& linearized_CT )
    {
      for( size_t dir=0; dir<CT.size(); ++dir ) // delinearize imported correlation tensors
      {
          // read correlation from linearized data
          const auto start_it = linearized_CT.cbegin() + dir*pspace.old_num_TimePoints;
          const auto end_it   = start_it + pspace.old_num_TimePoints;
          std::vector<RealType> correlation(start_it,end_it); // copy construction
          if( pspace.extrapolate_imported_spin_correlations ) // extrapolate if desired || otherwise the discretization matches already
          {
              correlation = num::extrapolate( correlation, pspace.num_TimePoints, pspace.old_delta_t, pspace.delta_t );
          }

          // write to correlation tensor
          CT[dir] = Corr{std::move(correlation)};
      }
    } );
  }
  return CCT;
}

};


namespace Functions
{
// compute the imaginary-time short-step propagator exp[-dt/2 (H_old + H_new)]
Operator compute_shortstep_propagator( const SparseObservable& old_Hamiltonian, const SparseObservable& new_Hamiltonian, const ps::ParameterSpace& pspace )
{
    SparseObservable Exponent = RealType{0.5} * ( old_Hamiltonian + new_Hamiltonian ); // CFET2

    if( pspace.matrix_exponential_computation == "taylor" ) // Computation via truncated taylor expansion (not recommended, large numerical errors)
    {
        return IDENTITY - pspace.delta_t * Exponent + RealType{0.5} * pspace.delta_t * pspace.delta_t * Exponent * Exponent;
    }
    else if( pspace.matrix_exponential_computation == "diagonalization" ) // Computation via diagonalization
    {
        // diagonalize:
        Operator UnitaryTransformation{};
        blaze::DynamicVector<RealType,blaze::columnVector> EigenValues{};
        Observable Exponent_DENSE{Exponent};
        diag::diagonalize_cplx( Exponent_DENSE, EigenValues, UnitaryTransformation );

        // compute exp(-i dt D) with D being the diagonalized exponent:
        DiagonalOperator DiagonalExponentMatrix{ ZERO };
        std::for_each( EigenValues.begin(), EigenValues.end(), [&,index=0]( const RealType& value ) mutable
        { 
            DiagonalExponentMatrix(index,index) = exp( ComplexType{ RealType{-1.0} * pspace.delta_t * value, RealType{0.0} } );
            index++;
        } );

        // transpose:
        return ctrans( UnitaryTransformation ) * DiagonalExponentMatrix * UnitaryTransformation;
    }
    else 
    {
        error::MATRIX_COMPUTATION( __PRETTY_FUNCTION__, pspace.matrix_exponential_computation );
        return Operator{};
    }
}

// compute all propagators U(t,0) and U(beta,t) from one sampled mean-field trajectory
std::tuple<RealType, std::vector<Operator>, std::vector<Operator>> compute_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace )
{
    std::vector<Operator> shortstep_propagators{};
    shortstep_propagators.reserve( pspace.num_TimePoints - 1 );

    SparseObservable old_Hamiltonian{}, new_Hamiltonian{};
    auto old_fields = Vs_of_t[0];
    pspace.mf_model->compute_Hamiltonian( old_Hamiltonian, std::move(old_fields), mean_fields, pspace.symmetry_type );
    for( size_t time = 1; time < pspace.num_TimePoints - 1; ++time )
    {
        auto new_fields = Vs_of_t[time];
        pspace.mf_model->compute_Hamiltonian( new_Hamiltonian, std::move(new_fields), mean_fields, pspace.symmetry_type );
        shortstep_propagators.emplace_back( compute_shortstep_propagator( old_Hamiltonian, new_Hamiltonian, pspace ) );
        old_Hamiltonian = new_Hamiltonian;
    }
    auto closing_fields = Vs_of_t[0];
    pspace.mf_model->compute_Hamiltonian( new_Hamiltonian, std::move(closing_fields), mean_fields, pspace.symmetry_type );
    shortstep_propagators.emplace_back( compute_shortstep_propagator( old_Hamiltonian, new_Hamiltonian, pspace ) );

    std::vector<Operator> forward_propagators{ IDENTITY };
    for( const auto& U_step : shortstep_propagators )
    {
        forward_propagators.emplace_back( U_step * forward_propagators.back() );
    }

    std::vector<Operator> backward_propagators( pspace.num_TimePoints, IDENTITY );
    for( size_t time = shortstep_propagators.size(); time > 0; --time )
    {
        backward_propagators[time-1] = backward_propagators[time] * shortstep_propagators[time-1];
    }

    RealType Z = std::real( blaze::trace( forward_propagators.back() ) );
    return std::make_tuple( Z, std::move( forward_propagators ), std::move( backward_propagators ) );
}

// compute an imaginary-time correlation Tr[U(beta,t) S_i^a U(t,0) S_j^b]
inline ComplexType correlation( const Operator& forward_U, const Operator& backward_U, const ten::IndexPair& ij, const ten::IndexPair& alphabeta )
{
    return blaze::trace( backward_U * S_at_site_in_direction(ij[0],alphabeta[0]) * forward_U * S_at_site_in_direction(ij[1],alphabeta[1]) );
}

// computes spin expectation values and correlations within the Monte-Carlo simulation
void compute_spin_observables( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<Operator>& forward_propagators, const std::vector<Operator>& backward_propagators, const RealType Z, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace )
{
    rtdata.Z_sqsum += Z * Z;

    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            ComplexType value = blaze::trace( forward_propagators.back() * S_at_site_in_direction(site,alpha) );
            const RealType x = std::real( value );
            spin_expvals[site][alpha] += x;
            rtdata.spin_expval_sqsum[site][alpha] += x * x;
            rtdata.spin_expval_cov[site][alpha] += x * Z;
        }
    }

    spin_CCT_Re.iterate2( rtdata.sample_sqsum_Re, [&]( auto& CT, auto& sqsum_CT, auto& ij )
    {
        auto& cov_CT = rtdata.sample_cov_Re( ij[0], ij[1] );
        CT.iterate2( sqsum_CT, [&]( auto& C, auto& sqsum_C, const auto& alphabeta )
        {
            auto cov_index = std::distance( cov_CT.m_direction_pairs.cbegin(), std::find( cov_CT.m_direction_pairs.cbegin(), cov_CT.m_direction_pairs.cend(), alphabeta ) );
            auto& cov_C = cov_CT[cov_index];
            for( size_t tau_index = 0; tau_index < pspace.num_TimePoints; ++tau_index )
            {
                ComplexType x = correlation( forward_propagators[tau_index], backward_propagators[tau_index], ij, alphabeta );
                const RealType x_re = std::real( x );
                C.at( tau_index ) += x_re;
                sqsum_C.at( tau_index ) += x_re * x_re;
                cov_C.at( tau_index ) += x_re * Z;
            }
        } );
    } );

    spin_CCT_Im.iterate2( rtdata.sample_sqsum_Im, [&]( auto& CT, auto& sqsum_CT, auto& ij )
    {
        auto& cov_CT = rtdata.sample_cov_Im( ij[0], ij[1] );
        CT.iterate2( sqsum_CT, [&]( auto& C, auto& sqsum_C, const auto& alphabeta )
        {
            auto cov_index = std::distance( cov_CT.m_direction_pairs.cbegin(), std::find( cov_CT.m_direction_pairs.cbegin(), cov_CT.m_direction_pairs.cend(), alphabeta ) );
            auto& cov_C = cov_CT[cov_index];
            for( size_t tau_index = 0; tau_index < pspace.num_TimePoints; ++tau_index )
            {
                ComplexType x = correlation( forward_propagators[tau_index], backward_propagators[tau_index], ij, alphabeta );
                const RealType x_im = std::imag( x );
                C.at( tau_index ) += x_im;
                sqsum_C.at( tau_index ) += x_im * x_im;
                cov_C.at( tau_index ) += x_im * Z;
            }
        } );
    } );
}

// shares the average and stddev results among the cores
void MPI_share_results( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, SiteFields& spin_expvals, rtd::RunTimeData& rtdata, RealType& partition_function )
{
    for( auto& spin_expval_i : spin_expvals )
    {
        std::vector<RealType> send_buf = { spin_expval_i[0], spin_expval_i[1], spin_expval_i[2] };
        std::vector<RealType> receive_buf( 3 );
        MPI_Allreduce( send_buf.data(), receive_buf.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_expval_i = FieldVector{ receive_buf[0], receive_buf[1], receive_buf[2] };
    }

    for( auto& spin_expval_i : rtdata.spin_expval_sqsum )
    {
        std::vector<RealType> send_buf = { spin_expval_i[0], spin_expval_i[1], spin_expval_i[2] };
        std::vector<RealType> receive_buf( 3 );
        MPI_Allreduce( send_buf.data(), receive_buf.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_expval_i = FieldVector{ receive_buf[0], receive_buf[1], receive_buf[2] };
    }

    for( auto& spin_expval_i : rtdata.spin_expval_cov )
    {
        std::vector<RealType> send_buf = { spin_expval_i[0], spin_expval_i[1], spin_expval_i[2] };
        std::vector<RealType> receive_buf( 3 );
        MPI_Allreduce( send_buf.data(), receive_buf.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_expval_i = FieldVector{ receive_buf[0], receive_buf[1], receive_buf[2] };
    }

    // share regular results
    std::for_each( spin_corr_Re.begin(), spin_corr_Re.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() ); 
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf; // at time zero the correlation value is set later
        } );
    } );

    std::for_each( spin_corr_Im.begin(), spin_corr_Im.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() ); 
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf; // at time zero the correlation value is set later
        } );
    } );

    // share squared average for sample average (for numerical analysis of the statistical error)
    std::for_each( rtdata.sample_sqsum_Re.begin(), rtdata.sample_sqsum_Re.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() ); 
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf; // at time zero the std is set to zero later
        } );
    } );

    std::for_each( rtdata.sample_sqsum_Im.begin(), rtdata.sample_sqsum_Im.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    std::for_each( rtdata.sample_cov_Re.begin(), rtdata.sample_cov_Re.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    std::for_each( rtdata.sample_cov_Im.begin(), rtdata.sample_cov_Im.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    RealType Z_sqsum_local = rtdata.Z_sqsum;
    MPI_Allreduce( &Z_sqsum_local, &rtdata.Z_sqsum, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    RealType Z_local = partition_function;
    MPI_Allreduce( &Z_local, &partition_function, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
}

// normalize the spin correlations with respect to the partition function
void normalize( CluCorrTen& CCT_Re, CluCorrTen& CCT_Im, const RealType partition_function )
{
  CCT_Re.iterate( [&]( CorrTen& CT, auto& )
  {
    CT.iterate( [&]( corr::CorrelationVector& corr, auto& )
    {
      std::for_each( corr.begin(), corr.end(), [partition_function]( RealType& c )
      {
        c = c / partition_function;
      } );
    } );
  } );
  CCT_Im.iterate( [&]( CorrTen& CT, auto& )
  {
    CT.iterate( [&]( corr::CorrelationVector& corr, auto& )
    {
      std::for_each( corr.begin(), corr.end(), [partition_function]( RealType& c )
      {
        c = c / partition_function;
      } );
    } );
  } );
}

};
