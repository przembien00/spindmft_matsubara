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

#include<Physics/CFET.h>
namespace cfet = Physics::CFET;

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

// initialization function : initialize 2x2 matrices for uncoupled mode
void write_uncoupled_spin_matrices()
{
    Operator Sx(2,2); Sx(0,0)=0; Sx(0,1)=0.5; Sx(1,0)=0.5; Sx(1,1)=0; S_X_2x2 = Sx;
    Operator Sy(2,2); Sy(0,0)=0; Sy(0,1)=ComplexType(0, -0.5); Sy(1,0)=ComplexType(0, 0.5); Sy(1,1)=0; S_Y_2x2 = Sy;
    Operator Sz(2,2); Sz(0,0)=0.5; Sz(0,1)=0; Sz(1,0)=0; Sz(1,1)=-0.5; S_Z_2x2 = Sz;
    Operator I(2,2); I(0,0)=1; I(0,1)=0; I(1,0)=0; I(1,1)=1; I_2x2 = I;
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
    RealType nonlocal_value{};
    if( pspace.init_pair_corr_mode == "zero" )
    {
      nonlocal_value = RealType{0.};
    }
    else if( pspace.init_pair_corr_mode == "positive" )
    {
      nonlocal_value = RealType{0.1};
    }
    else if( pspace.init_pair_corr_mode == "negative" )
    {
      nonlocal_value = RealType{-0.1};
    }
    else
      {
          error::INIT_CORRELATIONS_UNKNOWN( __PRETTY_FUNCTION__, pspace.init_pair_corr_mode );
      }
    corr::CorrelationVector NLC{ std::vector<RealType>( pspace.num_TimePoints, nonlocal_value ) };

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
      else // set any spin paircorrelations to the configured constant nonlocal seed
      {
          std::fill( CT.begin(), CT.end(), NLC );
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

namespace
{
std::vector<RealType> fourier_transform( const Corr& corr )
{
    const size_t num_frequencies = corr.size() - 1;
    std::vector<RealType> result( num_frequencies, RealType{0.} );
    for( size_t k = 0; k < num_frequencies; ++k )
    {
        for( size_t n = 0; n < num_frequencies; ++n )
        {
            result[k] += corr[n] * std::cos( ( RealType{2.} * M_PI * k * n ) / RealType(num_frequencies) );
        }
    }
    return result;
}

RealType frequency_component( const CluCorrTen& correlations, const size_t i, const size_t j, const size_t alpha, const size_t beta, const size_t frequency )
{
    return fourier_transform( correlations(i,j)(alpha,beta) )[frequency];
}

void set_frequency_covariance_entry( Multivariate_Gaussian::SymmetricMatrix& matrix, const size_t num_spins, const size_t alpha, const size_t i, const size_t beta, const size_t j, const RealType value )
{
    matrix( alpha*num_spins + i, beta*num_spins + j ) = value;
}

}

FrequencyCovarianceCluster::FrequencyCovarianceCluster( const CluCorrTen& correlations, const char symmetry_type, const size_t num_spins )
{
    fill( correlations, symmetry_type, num_spins );
}

size_t FrequencyCovarianceCluster::flat_index( const size_t frequency, const size_t block ) const
{
    return frequency*m_block_sizes.size() + block;
}

void FrequencyCovarianceCluster::fill( const CluCorrTen& correlations, const char symmetry_type, const size_t num_spins )
{
    m_symmetry_type = symmetry_type;
    m_num_spins = num_spins;
    m_num_frequencies = correlations[0][0].size() - 1;
    m_covariances.clear();
    m_fft = std::make_shared<FFT::InverseRealDFT>( m_num_frequencies ); // freq -> imaginary-time plan

    switch( m_symmetry_type )
    {
        case 'A':
        {
            m_block_sizes = std::vector<size_t>{ m_num_spins, m_num_spins, m_num_spins };
            break;
        }
        case 'B':
        {
            m_block_sizes = std::vector<size_t>{ m_num_spins, m_num_spins, m_num_spins };
            break;
        }
        case 'C':
        {
            m_block_sizes = std::vector<size_t>{ 2*m_num_spins, m_num_spins };
            break;
        }
        case 'D':
        {
            m_block_sizes = std::vector<size_t>{ 3*m_num_spins };
            break;
        }
        default:
        {
            error::MATRIX_COMPUTATION( __PRETTY_FUNCTION__, std::string(1,m_symmetry_type) );
        }
    }

    m_covariances.resize( m_num_frequencies*m_block_sizes.size() );
    for( size_t frequency = 0; frequency < m_num_frequencies; ++frequency )
    {
        for( size_t block = 0; block < m_block_sizes.size(); ++block )
        {
            m_covariances[flat_index(frequency,block)].resize( m_block_sizes[block] );
            m_covariances[flat_index(frequency,block)] = blaze::ZeroMatrix<RealType,blaze::rowMajor>( m_block_sizes[block], m_block_sizes[block] );
        }

        for( size_t i = 0; i < m_num_spins; ++i )
        {
            for( size_t j = 0; j < m_num_spins; ++j )
            {
                if( m_symmetry_type == 'A' )
                {
                    const RealType xx = frequency_component( correlations, i, j, 0, 0, frequency );
                    for( size_t block = 0; block < 3; ++block )
                    {
                        m_covariances[flat_index(frequency,block)]( i, j ) = xx;
                    }
                }
                else if( m_symmetry_type == 'B' )
                {
                    const RealType xx = frequency_component( correlations, i, j, 0, 0, frequency );
                    const RealType zz = frequency_component( correlations, i, j, 2, 2, frequency );
                    m_covariances[flat_index(frequency,0)]( i, j ) = xx;
                    m_covariances[flat_index(frequency,1)]( i, j ) = xx;
                    m_covariances[flat_index(frequency,2)]( i, j ) = zz;
                }
                else if( m_symmetry_type == 'C' )
                {
                    for( size_t alpha = 0; alpha < 2; ++alpha )
                    {
                        for( size_t beta = 0; beta < 2; ++beta )
                        {
                            set_frequency_covariance_entry( m_covariances[flat_index(frequency,0)], m_num_spins, alpha, i, beta, j, frequency_component( correlations, i, j, alpha, beta, frequency ) );
                        }
                    }
                    m_covariances[flat_index(frequency,1)]( i, j ) = frequency_component( correlations, i, j, 2, 2, frequency );
                }
                else if( m_symmetry_type == 'D' )
                {
                    for( size_t alpha = 0; alpha < 3; ++alpha )
                    {
                        for( size_t beta = 0; beta < 3; ++beta )
                        {
                            set_frequency_covariance_entry( m_covariances[flat_index(frequency,0)], m_num_spins, alpha, i, beta, j, frequency_component( correlations, i, j, alpha, beta, frequency ) );
                        }
                    }
                }
            }
        }
    }
}

void FrequencyCovarianceCluster::diagonalize( mvgb::EigenValuesBlocks& eig, mvgb::OrthogonalTransformationBlocks& ortho ) const
{
    eig.clear();
    ortho.clear();
    eig.reserve( m_covariances.size() );
    ortho.reserve( m_covariances.size() );
    std::for_each( m_covariances.cbegin(), m_covariances.cend(), [&]( const auto& matrix )
    {
        Multivariate_Gaussian::CovarianceMatrix covariance{ matrix };
        Multivariate_Gaussian::EigenValues evals{};
        Multivariate_Gaussian::OrthogonalTransformation transformation{};
        covariance.diagonalize( evals, transformation );
        eig.emplace_back( std::move(evals) );
        ortho.emplace_back( std::move(transformation) );
    } );
}

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

// compute the short-step propagators, the forward propagators U(t,0) and Z = tr U(beta,0)
// from one sampled mean-field trajectory. The backward propagators U(beta,t) are NOT built
// here: the pCN Metropolis test needs only Z, so backward (one extra N x N matmul per time
// point) is deferred to compute_backward_propagators and built only for accepted states.
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

    RealType Z = std::real( blaze::trace( forward_propagators.back() ) );
    return std::make_tuple( Z, std::move( forward_propagators ), std::move( shortstep_propagators ) );
}

// build the backward propagators U(beta,t) from the short-step propagators (coupled cluster)
std::vector<Operator> compute_backward_propagators( const std::vector<Operator>& shortstep_propagators, const size_t num_TimePoints )
{
    std::vector<Operator> backward_propagators( num_TimePoints, IDENTITY );
    for( size_t time = shortstep_propagators.size(); time > 0; --time )
    {
        backward_propagators[time-1] = backward_propagators[time] * shortstep_propagators[time-1];
    }
    return backward_propagators;
}

// ------------------------------------------------------------------------------------------------
// UNCOUPLED SPINS MODE
// ------------------------------------------------------------------------------------------------

// uncoupled counterpart of compute_propagators: returns Z_i, the forward propagators
// [time][site] and the per-site short-step propagators [site][step]. Backward propagators
// are deferred to compute_uncoupled_backward_propagators (built only for accepted states).
std::tuple<std::vector<RealType>, std::vector<std::vector<Operator>>, std::vector<std::vector<Operator>>> compute_uncoupled_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace )
{

    std::vector<std::vector<Operator>> forward_propagators( pspace.num_TimePoints, std::vector<Operator>(pspace.num_Spins) );
    std::vector<std::vector<Operator>> shortstep_per_site( pspace.num_Spins );
    std::vector<RealType> Z_i_list(pspace.num_Spins);

    Operator I2(2, 2);
    I2(0,0)=1; I2(0,1)=0; I2(1,0)=0; I2(1,1)=1;

    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        std::vector<Operator> shortstep_propagators{};
        shortstep_propagators.reserve( pspace.num_TimePoints - 1 );

        FieldVector old_fields = Vs_of_t[0][site] + mean_fields[site];
        if( pspace.chemical_shift.m_name != "none" ) old_fields[2] += pspace.chemical_shift.at(site);

        for( size_t time = 1; time < pspace.num_TimePoints - 1; ++time )
        {
            FieldVector new_fields = Vs_of_t[time][site] + mean_fields[site];
            if( pspace.chemical_shift.m_name != "none" ) new_fields[2] += pspace.chemical_shift.at(site);

            shortstep_propagators.emplace_back( cfet::CFET4opt_for_single_spin_one_half( new_fields, old_fields, pspace.delta_t, S_X_2x2, S_Y_2x2, S_Z_2x2, I_2x2 ) );
            old_fields = new_fields;
        }
        FieldVector closing_fields = Vs_of_t[0][site] + mean_fields[site];
        if( pspace.chemical_shift.m_name != "none" ) closing_fields[2] += pspace.chemical_shift.at(site);

        shortstep_propagators.emplace_back( cfet::CFET4opt_for_single_spin_one_half( closing_fields, old_fields, pspace.delta_t, S_X_2x2, S_Y_2x2, S_Z_2x2, I_2x2 ) );

        forward_propagators[0][site] = I2;
        for( size_t time = 0; time < shortstep_propagators.size(); ++time )
        {
            forward_propagators[time+1][site] = shortstep_propagators[time] * forward_propagators[time][site];
        }

        Z_i_list[site] = std::real( blaze::trace( forward_propagators.back()[site] ) );
        shortstep_per_site[site] = std::move( shortstep_propagators );
    }

    return std::make_tuple( Z_i_list, std::move( forward_propagators ), std::move( shortstep_per_site ) );
}

// build the backward propagators U(beta,t) [time][site] from the per-site short-step
// propagators (uncoupled cluster)
std::vector<std::vector<Operator>> compute_uncoupled_backward_propagators( const std::vector<std::vector<Operator>>& shortstep_per_site, const size_t num_TimePoints, const size_t num_Spins )
{
    Operator I2(2, 2);
    I2(0,0)=1; I2(0,1)=0; I2(1,0)=0; I2(1,1)=1;

    std::vector<std::vector<Operator>> backward_propagators( num_TimePoints, std::vector<Operator>(num_Spins) );
    for( size_t site = 0; site < num_Spins; ++site )
    {
        const std::vector<Operator>& shortstep_propagators = shortstep_per_site[site];
        backward_propagators[num_TimePoints - 1][site] = I2;
        for( size_t time = shortstep_propagators.size(); time > 0; --time )
        {
            backward_propagators[time-1][site] = backward_propagators[time][site] * shortstep_propagators[time-1];
        }
    }
    return backward_propagators;
}

inline const Observable& S_2x2( size_t alpha )
{
    if(alpha == 0) return S_X_2x2;
    if(alpha == 1) return S_Y_2x2;
    return S_Z_2x2;
}

inline ComplexType uncoupled_correlation( const std::vector<Operator>& forward_U_tau, const std::vector<Operator>& backward_U_tau, const std::vector<Operator>& forward_U_beta, const ten::IndexPair& ij, const ten::IndexPair& alphabeta )
{
    if( ij[0] != ij[1] )
    {
        return blaze::trace( backward_U_tau[ij[0]] * S_2x2(alphabeta[0]) * forward_U_tau[ij[0]] ) * blaze::trace( forward_U_beta[ij[1]] * S_2x2(alphabeta[1]) );
    }
    else
    {
        return blaze::trace( backward_U_tau[ij[0]] * S_2x2(alphabeta[0]) * forward_U_tau[ij[0]] * S_2x2(alphabeta[1]) );
    }
}

// ========================================================================================
// PRECONDITIONED CRANK-NICOLSON SAMPLER
// ========================================================================================

// Reconstruct one imaginary-time mean-field trajectory from a frequency-space sample given
// in the full (non-diagonal) block basis. Mirrors the per-sample logic in
// sample_time_trajectories but reads the coefficients from `full` instead of drawing them.
TimeTrajectory FrequencyCovarianceCluster::trajectory_from_full_frequency( const mvgb::EigenValuesBlocks& full ) const
{
    std::vector<std::vector<FieldVector>> frequency_values( m_num_spins, std::vector<FieldVector>( m_num_frequencies, FieldVector{0.,0.,0.} ) );
    for( size_t frequency = 0; frequency < m_num_frequencies; ++frequency )
    {
        if( m_symmetry_type == 'A' || m_symmetry_type == 'B' )
        {
            for( size_t site = 0; site < m_num_spins; ++site )
            {
                frequency_values[site][frequency][0] = full[flat_index(frequency,0)][site];
                frequency_values[site][frequency][1] = full[flat_index(frequency,1)][site];
                frequency_values[site][frequency][2] = full[flat_index(frequency,2)][site];
            }
        }
        else if( m_symmetry_type == 'C' )
        {
            for( size_t site = 0; site < m_num_spins; ++site )
            {
                frequency_values[site][frequency][0] = full[flat_index(frequency,0)][site];
                frequency_values[site][frequency][1] = full[flat_index(frequency,0)][m_num_spins + site];
                frequency_values[site][frequency][2] = full[flat_index(frequency,1)][site];
            }
        }
        else if( m_symmetry_type == 'D' )
        {
            for( size_t site = 0; site < m_num_spins; ++site )
            {
                frequency_values[site][frequency][0] = full[flat_index(frequency,0)][site];
                frequency_values[site][frequency][1] = full[flat_index(frequency,0)][m_num_spins + site];
                frequency_values[site][frequency][2] = full[flat_index(frequency,0)][2*m_num_spins + site];
            }
        }
    }

    // Inverse real DFT, performed independently for each site and each of the 3 field
    // components in O(N log N) via FFTW. FFT::InverseRealDFT reproduces the project's
    // half-complex convention exactly (DC ~ 1/sqrt(N), cosine/sine modes ~ 1/sqrt(N/2)),
    // so this is a drop-in replacement for the former O(N^2) inverse_fourier_at_time sums.
    const size_t N = m_num_frequencies;
    TimeTrajectory trajectory( N, SiteFields( m_num_spins, FieldVector{0.,0.,0.} ) );
    std::vector<RealType> comp_freq( N ), comp_time( N );
    for( size_t site = 0; site < m_num_spins; ++site )
    {
        for( size_t c = 0; c < 3; ++c )
        {
            for( size_t f = 0; f < N; ++f ) comp_freq[f] = frequency_values[site][f][c];
            m_fft->transform( comp_freq.data(), comp_time.data() );
            for( size_t t = 0; t < N; ++t ) trajectory[t][site][c] = comp_time[t];
        }
    }
    return trajectory;
}

// pCN observable accumulator (coupled cluster): fold O_raw / Z into the running sums.
void compute_spin_observables_mh( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<Operator>& forward_propagators, const std::vector<Operator>& backward_propagators, const RealType Z, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace )
{
    const RealType inv_Z = ( Z != RealType{0.} ) ? RealType{1.} / Z : RealType{0.};

    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            const RealType x = std::real( blaze::trace( forward_propagators.back() * S_at_site_in_direction(site,alpha) ) ) * inv_Z;
            spin_expvals[site][alpha] += x;
            rtdata.spin_expval_sqsum[site][alpha] += x * x;
            rtdata.cur_block_spin[site][alpha] += x;       // batch-means accumulation
        }
    }

    // Each correlation is Tr[ U(beta,t) S_i^a U(t,0) S_j^b ] = Tr[ M_{i,a} N_{j,b} ] with
    // M_{i,a} = U(beta,t) S_i^a and N_{j,b} = U(t,0) S_j^b. The left factor M depends only
    // on (i,a), the right factor N only on (j,b), so we build each once per tau and reuse it
    // across all site/direction pairs, then close with Tr[M N] = sum_{kl} M_kl N_lk =
    // sum( M % trans(N) ) -- an O(D^2) reduction instead of a third O(D^3) matrix product.
    //
    // Which (site,direction) operators are actually required is dictated by the symmetry type
    // via the stored site- and direction-pair lists: e.g. type A keeps only the (0,0) pair, so
    // only S_x enters here, while type D needs all three directions.
    const ten::IndexPairList& site_pairs = spin_CCT_Re.get_site_pairs();
    const ten::IndexPairList& dir_pairs  = spin_CCT_Re[0].get_direction_pairs();
    const size_t num_site_pairs = site_pairs.size();
    const size_t num_dir_pairs  = dir_pairs.size();

    // flag the left (M) and right (N) operators the symmetry type actually demands
    std::vector<std::array<bool,3>> need_left( pspace.num_Spins, {false,false,false} );
    std::vector<std::array<bool,3>> need_right( pspace.num_Spins, {false,false,false} );
    for( const auto& ij : site_pairs )
    {
        for( const auto& alphabeta : dir_pairs )
        {
            need_left[ ij[0] ][ alphabeta[0] ]  = true;
            need_right[ ij[1] ][ alphabeta[1] ] = true;
        }
    }

    // caches reused across tau; only the flagged (site,direction) entries are ever filled
    std::vector<std::array<Operator,3>> M( pspace.num_Spins );
    std::vector<std::array<Operator,3>> N( pspace.num_Spins );

    // Only the lower half [0, beta/2] is accumulated; the upper half is reconstructed once per
    // iteration by symmetrize_upper_half via g^{ab}_{ij}(beta-tau) = g^{ab}_{ij}(tau)^*.
    const size_t n_half = pspace.num_TimePoints - pspace.num_TimePoints / 2;
    for( size_t tau_index = 0; tau_index < n_half; ++tau_index )
    {
        const Operator& bwd = backward_propagators[tau_index];
        const Operator& fwd = forward_propagators[tau_index];
        for( size_t site = 0; site < pspace.num_Spins; ++site )
        {
            for( size_t dir = 0; dir < 3; ++dir )
            {
                if( need_left[site][dir] )  { M[site][dir] = bwd * S_at_site_in_direction(site,dir); }
                if( need_right[site][dir] ) { N[site][dir] = fwd * S_at_site_in_direction(site,dir); }
            }
        }

        for( size_t p = 0; p < num_site_pairs; ++p )
        {
            const auto& ij = site_pairs[p];
            CorrTen& CT_Re = spin_CCT_Re[p];
            CorrTen& CT_Im = spin_CCT_Im[p];
            CorrTen& sq_Re = rtdata.sample_sqsum_Re[p];
            CorrTen& sq_Im = rtdata.sample_sqsum_Im[p];
            CorrTen& blk_Re = rtdata.cur_block_Re[p];
            CorrTen& blk_Im = rtdata.cur_block_Im[p];
            for( size_t d = 0; d < num_dir_pairs; ++d )
            {
                const auto& alphabeta = dir_pairs[d];
                const ComplexType o = blaze::sum( M[ ij[0] ][ alphabeta[0] ] % blaze::trans( N[ ij[1] ][ alphabeta[1] ] ) ) * inv_Z;
                const RealType re = std::real( o );
                const RealType im = std::imag( o );
                CT_Re[d].at( tau_index ) += re;  sq_Re[d].at( tau_index ) += re * re;  blk_Re[d].at( tau_index ) += re;
                CT_Im[d].at( tau_index ) += im;  sq_Im[d].at( tau_index ) += im * im;  blk_Im[d].at( tau_index ) += im;
            }
        }
    }
}

// pCN observable accumulator (uncoupled cluster): divide each contribution site-by-site,
// i.e. by Z_i for on-site (i==j) correlations and by Z_i Z_j for inter-site correlations.
void compute_uncoupled_spin_observables_mh( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<std::vector<Operator>>& forward_propagators, const std::vector<std::vector<Operator>>& backward_propagators, const std::vector<RealType>& Z_i_list, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace )
{
    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        const RealType inv_Z_i = ( Z_i_list[site] != RealType{0.} ) ? RealType{1.} / Z_i_list[site] : RealType{0.};
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            const RealType x = std::real( blaze::trace( forward_propagators.back()[site] * S_2x2(alpha) ) ) * inv_Z_i;
            spin_expvals[site][alpha] += x;
            rtdata.spin_expval_sqsum[site][alpha] += x * x;
            rtdata.cur_block_spin[site][alpha] += x;       // batch-means accumulation
        }
    }

    // The mean, sum-of-squares and current-block clusters share the same site-pair and
    // direction-pair layout (built from the same correlation_categories / symmetry_type), so we
    // traverse them in lockstep by linear index rather than via operator()(i,j) / operator()(a,b).
    const size_t n_half = pspace.num_TimePoints - pspace.num_TimePoints / 2;
    auto accumulate = [&]( CluCorrTen& spin_CCT, CluCorrTen& sqsum_CCT, CluCorrTen& block_CCT, const bool real_part )
    {
        for( size_t p = 0; p < spin_CCT.size(); ++p )
        {
            const auto& ij = spin_CCT.get_site_pair( p );
            const RealType inv_Z_i = ( Z_i_list[ij[0]] != RealType{0.} ) ? RealType{1.} / Z_i_list[ij[0]] : RealType{0.};
            const RealType inv_Z_j = ( Z_i_list[ij[1]] != RealType{0.} ) ? RealType{1.} / Z_i_list[ij[1]] : RealType{0.};
            const RealType inv_weight = ( ij[0] == ij[1] ) ? inv_Z_i : inv_Z_i * inv_Z_j;
            CorrTen& CT       = spin_CCT[p];
            CorrTen& sqsum_CT = sqsum_CCT[p];
            CorrTen& block_CT = block_CCT[p];
            const ten::IndexPairList& dir_pairs = CT.get_direction_pairs();
            for( size_t d = 0; d < CT.size(); ++d )
            {
                const auto& alphabeta = dir_pairs[d];
                corr::CorrelationVector& C       = CT[d];
                corr::CorrelationVector& sqsum_C = sqsum_CT[d];
                corr::CorrelationVector& block_C = block_CT[d];
                // Only the lower half [0, beta/2] is accumulated; symmetrize_upper_half fills
                // the rest via g^{ab}_{ij}(beta-tau) = g^{ab}_{ij}(tau)^*.
                for( size_t tau_index = 0; tau_index < n_half; ++tau_index )
                {
                    const ComplexType val = uncoupled_correlation( forward_propagators[tau_index], backward_propagators[tau_index], forward_propagators.back(), ij, alphabeta );
                    const RealType o = ( real_part ? std::real( val ) : std::imag( val ) ) * inv_weight;
                    C.at( tau_index ) += o;
                    sqsum_C.at( tau_index ) += o * o;
                    block_C.at( tau_index ) += o;          // batch-means accumulation
                }
            }
        }
    };
    accumulate( spin_CCT_Re, rtdata.sample_sqsum_Re, rtdata.cur_block_Re, true );
    accumulate( spin_CCT_Im, rtdata.sample_sqsum_Im, rtdata.cur_block_Im, false );
}

// reduce one CluCorrTen of per-(i,j,alpha,beta,tau) sums across MPI ranks
static void mpi_sum_clucorrtensor( CluCorrTen& CCT )
{
    std::for_each( CCT.begin(), CCT.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv( corr.size() );
            MPI_Allreduce( corr.data(), rcv.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv;
        } );
    } );
}

void MPI_share_results_mh( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, SiteFields& spin_expvals, rtd::RunTimeData& rtdata )
{
    auto reduce_fieldvectors = []( std::vector<FieldVector>& v )
    {
        for( auto& f : v )
        {
            std::vector<RealType> send = { f[0], f[1], f[2] };
            std::vector<RealType> recv( 3 );
            MPI_Allreduce( send.data(), recv.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            f = FieldVector{ recv[0], recv[1], recv[2] };
        }
    };
    reduce_fieldvectors( spin_expvals );
    reduce_fieldvectors( rtdata.spin_expval_sqsum );

    mpi_sum_clucorrtensor( spin_corr_Re );
    mpi_sum_clucorrtensor( spin_corr_Im );
    mpi_sum_clucorrtensor( rtdata.sample_sqsum_Re );
    mpi_sum_clucorrtensor( rtdata.sample_sqsum_Im );
}

void normalize_mh( CluCorrTen& CCT_Re, CluCorrTen& CCT_Im, const RealType N )
{
    const RealType inv_N = RealType{1.} / N;
    auto scale = [inv_N]( CluCorrTen& CCT )
    {
        CCT.iterate( [&]( CorrTen& CT, auto& )
        {
            CT.iterate( [&]( corr::CorrelationVector& corr, auto& )
            {
                std::for_each( corr.begin(), corr.end(), [inv_N]( RealType& c ){ c *= inv_N; } );
            } );
        } );
    };
    scale( CCT_Re );
    scale( CCT_Im );
}

// Fill the upper half (tau > beta/2) of every (site-pair, direction-pair) correlation from the
// computed lower half via the reflection symmetry g^{ab}_{ij}(beta-tau) = g^{ab}_{ij}(tau)^*,
// i.e. Re even and Im odd about beta/2. The mirror partner of grid point k (tau_k = k*delta_t)
// is point N-1-k (= beta - tau_k), always in the computed lower half. Standard errors mirror
// with no sign change, since each reconstructed value equals its mirror source deterministically.
// The self-paired midpoint (present for an odd number of time points) lies in the lower half and
// keeps its raw sampled value.
void symmetrize_upper_half( CluCorrTen& Re, CluCorrTen& Im, CluCorrTen& Re_std, CluCorrTen& Im_std, CluCorrTen& Re_tau, CluCorrTen& Im_tau )
{
    const size_t n_site_pairs = Re.size();
    for( size_t p = 0; p < n_site_pairs; ++p )
    {
        const size_t n_dir_pairs = Re[p].size();
        for( size_t d = 0; d < n_dir_pairs; ++d )
        {
            corr::CorrelationVector& re   = Re[p][d];
            corr::CorrelationVector& im   = Im[p][d];
            corr::CorrelationVector& rstd = Re_std[p][d];
            corr::CorrelationVector& istd = Im_std[p][d];
            corr::CorrelationVector& rtau = Re_tau[p][d];
            corr::CorrelationVector& itau = Im_tau[p][d];
            const size_t N = re.size();          // num_TimePoints
            const size_t n_half = N - N / 2;     // computed points in [0, beta/2]
            for( size_t k = n_half; k < N; ++k )
            {
                const size_t m = N - 1 - k;      // mirror index, in the computed lower half
                re[k]   =  re[m];
                im[k]   = -im[m];
                rstd[k] =  rstd[m];
                istd[k] =  istd[m];
                rtau[k] =  rtau[m];
                itau[k] =  itau[m];
            }
        }
    }
}

// ----------------------------- PCNChainCluster implementation -----------------------------

void PCNChainCluster::build_prior()
{
    m_prior_dists.assign( m_eig.size(), {} );
    for( size_t f = 0; f < m_eig.size(); ++f )
    {
        m_prior_dists[f].resize( m_eig[f].size() );
        for( size_t i = 0; i < m_eig[f].size(); ++i )
        {
            const RealType sigma = ( m_eig[f][i] > RealType{0.} ) ? std::sqrt( m_eig[f][i] ) : RealType{0.};
            m_prior_dists[f][i] = std::normal_distribution<RealType>( RealType{0.}, sigma );
        }
    }
}

void PCNChainCluster::draw_prior_into( mvgb::EigenValuesBlocks& V_diag )
{
    V_diag.resize( m_eig.size() );
    for( size_t f = 0; f < m_eig.size(); ++f )
    {
        V_diag[f].resize( m_eig[f].size() );
        for( size_t i = 0; i < m_eig[f].size(); ++i )
        {
            V_diag[f][i] = ( m_prior_dists[f][i].stddev() > RealType{0.} ) ? m_prior_dists[f][i]( m_engine ) : RealType{0.};
        }
    }
}

mvgb::EigenValuesBlocks PCNChainCluster::rotate_to_full( const mvgb::EigenValuesBlocks& V_diag ) const
{
    mvgb::EigenValuesBlocks full( m_eig.size() );
    for( size_t f = 0; f < m_eig.size(); ++f )
    {
        full[f] = m_ortho[f] * V_diag[f];
    }
    return full;
}

void PCNChainCluster::evaluate( const mvgb::EigenValuesBlocks& V_diag, RealType& Z, RealType& logZ,
                                std::vector<Operator>& fwd, std::vector<Operator>& shortstep,
                                std::vector<RealType>& Z_i, std::vector<std::vector<Operator>>& fwd_unc,
                                std::vector<std::vector<Operator>>& shortstep_unc ) const
{
    const TimeTrajectory traj = m_cov.trajectory_from_full_frequency( rotate_to_full( V_diag ) );
    if( m_uncoupled )
    {
        std::tie( Z_i, fwd_unc, shortstep_unc ) = compute_uncoupled_propagators( traj, m_mean_fields, m_pspace );
        logZ = RealType{0.};
        bool positive = true;
        for( const RealType zi : Z_i )
        {
            if( zi > RealType{0.} ) { logZ += std::log( zi ); }
            else { positive = false; }
        }
        Z = positive ? RealType{1.} : RealType{0.}; // validity flag; observables use Z_i directly
    }
    else
    {
        std::tie( Z, fwd, shortstep ) = compute_propagators( traj, m_mean_fields, m_pspace );
        logZ = ( Z > RealType{0.} ) ? std::log( Z ) : RealType{0.};
    }
}

void PCNChainCluster::refresh_backward()
{
    if( m_uncoupled )
    {
        m_backward_unc = compute_uncoupled_backward_propagators( m_shortstep_unc, m_pspace.num_TimePoints, m_pspace.num_Spins );
    }
    else
    {
        m_backward = compute_backward_propagators( m_shortstep, m_pspace.num_TimePoints );
    }
}

PCNChainCluster::PCNChainCluster( const ps::ParameterSpace& pspace,
                                  const FrequencyCovarianceCluster& cov,
                                  const mvgb::EigenValuesBlocks& eig,
                                  const mvgb::OrthogonalTransformationBlocks& ortho,
                                  const SiteFields& mean_fields,
                                  RealType pcn_beta,
                                  std::mt19937& engine )
    : m_pspace( pspace ), m_cov( cov ), m_eig( eig ), m_ortho( ortho ), m_mean_fields( mean_fields ),
      m_engine( engine ), m_uncoupled( pspace.uncoupled_spins ),
      m_pcn_beta( pcn_beta ),
      m_pcn_root( std::sqrt( std::max( RealType{0.}, RealType{1.} - pcn_beta * pcn_beta ) ) )
{
    build_prior();

    // cold start: draw from the prior; retry on a non-positive Z (numerically degenerate)
    draw_prior_into( m_V_diag );
    evaluate( m_V_diag, m_Z, m_logZ, m_forward, m_shortstep, m_Z_i_list, m_forward_unc, m_shortstep_unc );
    for( size_t retries = 0; !( m_Z > RealType{0.} ) && retries < 32; ++retries )
    {
        draw_prior_into( m_V_diag );
        evaluate( m_V_diag, m_Z, m_logZ, m_forward, m_shortstep, m_Z_i_list, m_forward_unc, m_shortstep_unc );
    }
    refresh_backward(); // build the backward propagators for the initial accepted state
}

PCNChainCluster::PCNChainCluster( const ps::ParameterSpace& pspace,
                                  const FrequencyCovarianceCluster& cov,
                                  const mvgb::EigenValuesBlocks& eig,
                                  const mvgb::OrthogonalTransformationBlocks& ortho,
                                  const SiteFields& mean_fields,
                                  RealType pcn_beta,
                                  std::mt19937& engine,
                                  const mvgb::EigenValuesBlocks& full_freq_initial )
    : PCNChainCluster( pspace, cov, eig, ortho, mean_fields, pcn_beta, engine )
{
    // warm start: rotate the carried-over full-basis state into the new diagonal basis,
    // zeroing modes the new prior cannot support (eig <= 0), then refresh the caches.
    m_V_diag.resize( m_eig.size() );
    for( size_t f = 0; f < m_eig.size(); ++f )
    {
        const Multivariate_Gaussian::EigenValues V_diag_f = blaze::trans( m_ortho[f] ) * full_freq_initial[f];
        m_V_diag[f].resize( m_eig[f].size() );
        for( size_t i = 0; i < m_eig[f].size(); ++i )
        {
            m_V_diag[f][i] = ( m_prior_dists[f][i].stddev() > RealType{0.} ) ? V_diag_f[i] : RealType{0.};
        }
    }
    evaluate( m_V_diag, m_Z, m_logZ, m_forward, m_shortstep, m_Z_i_list, m_forward_unc, m_shortstep_unc );
    refresh_backward();
}

mvgb::EigenValuesBlocks PCNChainCluster::full_frequency() const
{
    return rotate_to_full( m_V_diag );
}

bool PCNChainCluster::step()
{
    // 1.) pCN proposal in the diagonal basis: V' = sqrt(1-beta^2) V + beta xi, xi ~ p(V)
    mvgb::EigenValuesBlocks V_diag_prop( m_eig.size() );
    for( size_t f = 0; f < m_eig.size(); ++f )
    {
        Multivariate_Gaussian::EigenValues xi( m_eig[f].size() );
        for( size_t i = 0; i < m_eig[f].size(); ++i )
        {
            xi[i] = ( m_prior_dists[f][i].stddev() > RealType{0.} ) ? m_prior_dists[f][i]( m_engine ) : RealType{0.};
        }
        V_diag_prop[f] = m_pcn_root * m_V_diag[f] + m_pcn_beta * xi;
    }

    // 2.) evaluate the target at the proposal (forward propagators + Z only; backward is
    //     skipped here since the Metropolis test below needs only Z)
    RealType Z_prop{}, logZ_prop{};
    std::vector<Operator> fwd_prop, shortstep_prop;
    std::vector<RealType> Z_i_prop;
    std::vector<std::vector<Operator>> fwd_unc_prop, shortstep_unc_prop;
    evaluate( V_diag_prop, Z_prop, logZ_prop, fwd_prop, shortstep_prop, Z_i_prop, fwd_unc_prop, shortstep_unc_prop );

    // 3.) pCN acceptance: the Gaussian prior factors cancel, leaving min(1, Z'/Z).
    //     A non-positive Z marks a forbidden region.
    if( !( Z_prop > RealType{0.} ) || !( m_Z > RealType{0.} ) ) return false;
    const RealType log_alpha = logZ_prop - m_logZ;
    if( std::log( m_uniform01( m_engine ) ) >= log_alpha ) return false;

    // 4.) accept: swap in the proposed state and its cached propagators, then build the
    //     backward propagators for the newly accepted state (needed by the observables)
    m_V_diag = std::move( V_diag_prop );
    m_Z = Z_prop;
    m_logZ = logZ_prop;
    m_forward = std::move( fwd_prop );
    m_shortstep = std::move( shortstep_prop );
    m_Z_i_list = std::move( Z_i_prop );
    m_forward_unc = std::move( fwd_unc_prop );
    m_shortstep_unc = std::move( shortstep_unc_prop );
    refresh_backward();
    return true;
}

};
