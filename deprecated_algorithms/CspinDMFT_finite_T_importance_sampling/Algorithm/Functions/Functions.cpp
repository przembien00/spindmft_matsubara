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

FieldVector inverse_fourier_at_time( const std::vector<FieldVector>& frequency_values, const size_t time )
{
    const size_t num_frequencies = frequency_values.size();
    FieldVector result{0.,0.,0.};
    result += frequency_values[0] / std::sqrt( RealType(num_frequencies) );
    for( size_t frequency = 1; frequency < num_frequencies/2 + 1; ++frequency )
    {
        result += frequency_values[frequency] * std::cos( ( RealType{2.} * M_PI * time * frequency ) / RealType(num_frequencies) ) / std::sqrt( RealType(num_frequencies) / RealType{2.} );
    }
    for( size_t frequency = num_frequencies/2 + 1; frequency < num_frequencies; ++frequency )
    {
        result += -frequency_values[frequency] * std::sin( ( RealType{2.} * M_PI * time * frequency ) / RealType(num_frequencies) ) / std::sqrt( RealType(num_frequencies) / RealType{2.} );
    }
    return result;
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

TimeTrajectories FrequencyCovarianceCluster::sample_time_trajectories( const mvgb::EigenValuesBlocks& eig, const mvgb::OrthogonalTransformationBlocks& ortho, std::mt19937& engine, const size_t num_samples ) const
{
    mvgb::DiagonalBasisNormalDistributionsBlocks distributions{ eig };
    std::vector<size_t> num_noises_per_block( eig.size(), num_samples );
    mvgb::GaussianNoiseVectorsBlocks frequency_noise{ num_noises_per_block, distributions, engine, ortho };

    TimeTrajectories trajectories( num_samples );
    for( size_t sample = 0; sample < num_samples; ++sample )
    {
        std::vector<std::vector<FieldVector>> frequency_values( m_num_spins, std::vector<FieldVector>( m_num_frequencies, FieldVector{0.,0.,0.} ) );
        for( size_t frequency = 0; frequency < m_num_frequencies; ++frequency )
        {
            if( m_symmetry_type == 'A' || m_symmetry_type == 'B' )
            {
                for( size_t site = 0; site < m_num_spins; ++site )
                {
                    frequency_values[site][frequency][0] = frequency_noise[flat_index(frequency,0)]( site, sample );
                    frequency_values[site][frequency][1] = frequency_noise[flat_index(frequency,1)]( site, sample );
                    frequency_values[site][frequency][2] = frequency_noise[flat_index(frequency,2)]( site, sample );
                }
            }
            else if( m_symmetry_type == 'C' )
            {
                for( size_t site = 0; site < m_num_spins; ++site )
                {
                    frequency_values[site][frequency][0] = frequency_noise[flat_index(frequency,0)]( site, sample );
                    frequency_values[site][frequency][1] = frequency_noise[flat_index(frequency,0)]( m_num_spins + site, sample );
                    frequency_values[site][frequency][2] = frequency_noise[flat_index(frequency,1)]( site, sample );
                }
            }
            else if( m_symmetry_type == 'D' )
            {
                for( size_t site = 0; site < m_num_spins; ++site )
                {
                    frequency_values[site][frequency][0] = frequency_noise[flat_index(frequency,0)]( site, sample );
                    frequency_values[site][frequency][1] = frequency_noise[flat_index(frequency,0)]( m_num_spins + site, sample );
                    frequency_values[site][frequency][2] = frequency_noise[flat_index(frequency,0)]( 2*m_num_spins + site, sample );
                }
            }
        }

        trajectories[sample].resize( m_num_frequencies );
        for( size_t time = 0; time < m_num_frequencies; ++time )
        {
            trajectories[sample][time].resize( m_num_spins );
            for( size_t site = 0; site < m_num_spins; ++site )
            {
                trajectories[sample][time][site] = inverse_fourier_at_time( frequency_values[site], time );
            }
        }
    }
    return trajectories;
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

// ------------------------------------------------------------------------------------------------
// UNCOUPLED SPINS MODE
// ------------------------------------------------------------------------------------------------

std::tuple<std::vector<RealType>, std::vector<std::vector<Operator>>, std::vector<std::vector<Operator>>> compute_uncoupled_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace )
{

    std::vector<std::vector<Operator>> forward_propagators( pspace.num_TimePoints, std::vector<Operator>(pspace.num_Spins) );
    std::vector<std::vector<Operator>> backward_propagators( pspace.num_TimePoints, std::vector<Operator>(pspace.num_Spins) );
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
        
        backward_propagators[pspace.num_TimePoints - 1][site] = I2;
        for( size_t time = shortstep_propagators.size(); time > 0; --time )
        {
            backward_propagators[time-1][site] = backward_propagators[time][site] * shortstep_propagators[time-1];
        }
        
        Z_i_list[site] = std::real( blaze::trace( forward_propagators.back()[site] ) );
    }
    
    return std::make_tuple( Z_i_list, std::move( forward_propagators ), std::move( backward_propagators ) );
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

void accumulate_uncoupled_Z( std::vector<RealType>& uncoupled_partition_functions, rtd::RunTimeData& rtdata, const std::vector<RealType>& Z_i_list, const ps::ParameterSpace& pspace )
{
    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        uncoupled_partition_functions[site] += Z_i_list[site];
        rtdata.uncoupled_Z_sqsum[site] += Z_i_list[site] * Z_i_list[site];
    }
    rtdata.uncoupled_Z_cov.iterate( [&]( RealType& C, const auto& ij )
    {
        C += Z_i_list[ij[0]] * Z_i_list[ij[1]];
    } );
}

void compute_uncoupled_spin_observables( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<std::vector<Operator>>& forward_propagators, const std::vector<std::vector<Operator>>& backward_propagators, const std::vector<RealType>& Z_i_list, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace )
{
    for( size_t site = 0; site < pspace.num_Spins; ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            ComplexType value = blaze::trace( forward_propagators.back()[site] * S_2x2(alpha) );
            RealType x = std::real( value );
            spin_expvals[site][alpha] += x;
            rtdata.spin_expval_sqsum[site][alpha] += x * x;
            rtdata.spin_expval_cov[site][alpha] += x * Z_i_list[site];
        }
    }

    spin_CCT_Re.iterate2( rtdata.sample_sqsum_Re, [&]( auto& CT, auto& sqsum_CT, auto& ij )
    {
        auto& cov_CT_i = rtdata.uncoupled_sample_cov_Re_i( ij[0], ij[1] );
        auto& cov_CT_j = rtdata.uncoupled_sample_cov_Re_j( ij[0], ij[1] );
        
        CT.iterate2( sqsum_CT, [&]( auto& C, auto& sqsum_C, const auto& alphabeta )
        {
            auto cov_index = std::distance( cov_CT_i.m_direction_pairs.cbegin(), std::find( cov_CT_i.m_direction_pairs.cbegin(), cov_CT_i.m_direction_pairs.cend(), alphabeta ) );
            auto& cov_C_i = cov_CT_i[cov_index];
            auto& cov_C_j = cov_CT_j[cov_index];
            
            for( size_t tau_index = 0; tau_index < pspace.num_TimePoints; ++tau_index )
            {
                ComplexType val = uncoupled_correlation( forward_propagators[tau_index], backward_propagators[tau_index], forward_propagators.back(), ij, alphabeta );
                RealType x = std::real( val );
                
                C.at( tau_index ) += x;
                sqsum_C.at( tau_index ) += x * x;
                cov_C_i.at( tau_index ) += x * Z_i_list[ij[0]];
                cov_C_j.at( tau_index ) += x * Z_i_list[ij[1]];
            }
        } );
    } );

    spin_CCT_Im.iterate2( rtdata.sample_sqsum_Im, [&]( auto& CT, auto& sqsum_CT, auto& ij )
    {
        auto& cov_CT_i = rtdata.uncoupled_sample_cov_Im_i( ij[0], ij[1] );
        auto& cov_CT_j = rtdata.uncoupled_sample_cov_Im_j( ij[0], ij[1] );
        
        CT.iterate2( sqsum_CT, [&]( auto& C, auto& sqsum_C, const auto& alphabeta )
        {
            auto cov_index = std::distance( cov_CT_i.m_direction_pairs.cbegin(), std::find( cov_CT_i.m_direction_pairs.cbegin(), cov_CT_i.m_direction_pairs.cend(), alphabeta ) );
            auto& cov_C_i = cov_CT_i[cov_index];
            auto& cov_C_j = cov_CT_j[cov_index];
            
            for( size_t tau_index = 0; tau_index < pspace.num_TimePoints; ++tau_index )
            {
                ComplexType val = uncoupled_correlation( forward_propagators[tau_index], backward_propagators[tau_index], forward_propagators.back(), ij, alphabeta );
                RealType x = std::imag( val );
                
                C.at( tau_index ) += x;
                sqsum_C.at( tau_index ) += x * x;
                cov_C_i.at( tau_index ) += x * Z_i_list[ij[0]];
                cov_C_j.at( tau_index ) += x * Z_i_list[ij[1]];
            }
        } );
    } );
}

void normalize_uncoupled( CluCorrTen& CCT_Re, CluCorrTen& CCT_Im, SiteFields& spin_expvals, const std::vector<RealType>& partition_functions, const size_t num_Samples )
{
    for( size_t site = 0; site < spin_expvals.size(); ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            spin_expvals[site][alpha] /= partition_functions[site];
        }
    }

    CCT_Re.iterate2( CCT_Im, [&partition_functions, num_Samples]( auto& CT_re, auto& CT_im, const auto& ij )
    {
        RealType Z_ij = partition_functions[ij[0]];
        if( ij[0] != ij[1] )
        {
            Z_ij = partition_functions[ij[0]] * partition_functions[ij[1]] / static_cast<RealType>(num_Samples);
        }
        
        CT_re.iterate2( CT_im, [Z_ij]( corr::CorrelationVector& c_re, corr::CorrelationVector& c_im, const auto& alphabeta )
        {
            for( auto& val : c_re ) val *= 1.0 / Z_ij;
            for( auto& val : c_im ) val *= 1.0 / Z_ij;
        } );
    } );
}

void MPI_share_uncoupled_results( std::vector<RealType>& partition_functions, rtd::RunTimeData& rtdata )
{
    std::vector<RealType> global_partition_functions( partition_functions.size(), 0.0 );
    MPI_Allreduce( partition_functions.data(), global_partition_functions.data(), partition_functions.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    partition_functions = std::move( global_partition_functions );

    std::vector<RealType> global_zsq( rtdata.uncoupled_Z_sqsum.size(), 0.0 );
    MPI_Allreduce( rtdata.uncoupled_Z_sqsum.data(), global_zsq.data(), rtdata.uncoupled_Z_sqsum.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    rtdata.uncoupled_Z_sqsum = std::move( global_zsq );

    rtdata.uncoupled_Z_cov.iterate( [&]( RealType& c, const auto& )
    {
        RealType rcv_c;
        MPI_Allreduce( &c, &rcv_c, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        c = rcv_c;
    } );

    std::for_each( rtdata.uncoupled_sample_cov_Re_i.begin(), rtdata.uncoupled_sample_cov_Re_i.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    std::for_each( rtdata.uncoupled_sample_cov_Re_j.begin(), rtdata.uncoupled_sample_cov_Re_j.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    std::for_each( rtdata.uncoupled_sample_cov_Im_i.begin(), rtdata.uncoupled_sample_cov_Im_i.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
        } );
    } );

    std::for_each( rtdata.uncoupled_sample_cov_Im_j.begin(), rtdata.uncoupled_sample_cov_Im_j.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() );
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf;
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
