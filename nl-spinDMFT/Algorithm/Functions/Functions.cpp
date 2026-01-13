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
CorrTenList generate_environment_spin_correlations( const ps::ParameterSpace& pspace )
{
  // reserve
  CorrTenList spin_corr( pspace.import_only_categories.size(), CorrTen{pspace.symmetry_type, pspace.num_TimePoints} );

  if( pspace.init_corr_mode == "generate" ) // set the spin correlations to some pregiven functions, import_only is translated to initialize_only
  {
    // create the zero correlation:
    corr::CorrelationVector ZC{  ph::NonDiagonalSpinCorrelation{"zero"}.create_discretization( pspace.delta_t, pspace.num_TimePoints, 0. ) };

    // insert them into the spin correlations
    stda::for_2each( spin_corr.begin(), spin_corr.end(), pspace.import_only_categories.cbegin(), [&]( CorrTen& CT, const IndexPair& ij )
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
    if( pspace.imported_correlations_src_directory == "spinDMFT" || pspace.imported_correlations_src_directory == "CspinDMFT" || pspace.imported_correlations_src_directory == "p-spinDMFT" )
    {
      stda::for_2each( spin_corr.begin(), spin_corr.end(), pspace.imported_correlations_linearized.cbegin(), [&]( CorrTen& CT, const std::vector<RealType>& imported_CT_linearized )
      {
        for( size_t dir=0; dir<CT.size(); ++dir ) // delinearize imported correlation tensors
        {
            // read correlation from linearized data
            const auto start_it = imported_CT_linearized.cbegin() + dir*pspace.old_num_TimePoints;
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
    else
    {
        error::IMPLEMENTATION_MISSING( __PRETTY_FUNCTION__, pspace.imported_correlations_src_directory );
    }
  }
  return spin_corr;
}

}


namespace Functions
{

// propagate the time evolution operator by a time step delta_t
void propagate( Operator& TimeEvolutionOperator, const SparseObservable& old_Hamiltonian, const SparseObservable& new_Hamiltonian, const ps::ParameterSpace& pspace )
{
    SparseObservable Exponent = RealType{0.5} * ( old_Hamiltonian + new_Hamiltonian ); // CFET2

    if( pspace.matrix_exponential_computation == "taylor" ) // Computation via truncated taylor expansion (not recommended, large numerical errors)
    {
        TimeEvolutionOperator = ( IDENTITY  +  ComplexType( 0.0, -pspace.delta_t ) * Exponent  +  ComplexType( -pspace.delta_t*pspace.delta_t /2.0 , 0.0 ) * Exponent * Exponent ) * TimeEvolutionOperator; 
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
            DiagonalExponentMatrix(index,index) = exp( ComplexType{RealType{0.0},RealType{-1.0}*pspace.delta_t} * value ); 
            index++;
        } );

        // transpose:
        TimeEvolutionOperator = ctrans( UnitaryTransformation ) * DiagonalExponentMatrix * UnitaryTransformation * TimeEvolutionOperator;
    }
    else 
    {
        error::MATRIX_COMPUTATION( __PRETTY_FUNCTION__, pspace.matrix_exponential_computation );
    }
}

// compute an autocorrelation at two specific times t1 and t2 in dependence of U, U^+, ij and alphabeta
inline RealType correlation( const Operator& U, const Operator& U_dagger, const ten::IndexPair& ij, const ten::IndexPair& alphabeta, const ps::ParameterSpace& pspace )
{
    return RealType{1.} / static_cast<RealType>(pspace.num_HilbertSpaceDimension) * num::cast_if_real_soft( blaze::trace( U_dagger * S_at_site_in_direction(ij[0],alphabeta[0]) * U * S_at_site_in_direction(ij[1],alphabeta[1]) ) ); // <S^alpha_i(time)S^beta_j(0)>
}

// computes the desired spin correlations and their stddev within the Monte-Carlo simulation
void compute_spin_correlations( CluCorrTen& spin_CCT, const Operator& TimeEvolutionOperator, const size_t time, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace )
{
    auto TimeEvolutionOperator_dagger = blaze::ctrans( TimeEvolutionOperator );

    spin_CCT.iterate2( rtdata.sample_sqsum, [&]( auto& CT, auto& sqsum_CT, auto& ij )
    {
        CT.iterate2( sqsum_CT, [&]( auto& C, auto& sqsum_C, const auto& alphabeta )
        {
            RealType x = correlation( TimeEvolutionOperator, TimeEvolutionOperator_dagger, ij, alphabeta, pspace );
            C.at( time ) += x; // add to correlation
            sqsum_C.at( time ) += x*x; // add to square sum
        } );
    } );
}

// shares the average and stddev results among the cores
void MPI_share_results( CluCorrTen& spin_corr, rtd::RunTimeData& rtdata )
{
    // share regular results
    std::for_each( spin_corr.begin(), spin_corr.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() ); 
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf; // at time zero the correlation value is set later
        } );
    } );

    // share squared average for sample average (for numerical analysis of the statistical error)
    std::for_each( rtdata.sample_sqsum.begin(), rtdata.sample_sqsum.end(), []( CorrTen& CT )
    {
        std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& corr )
        {
            std::vector<RealType> rcv_buf( corr.size() ); 
            MPI_Allreduce( corr.data(), rcv_buf.data(), corr.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
            corr = rcv_buf; // at time zero the std is set to zero later
        } );
    } );
}

// normalize the spin correlations with respect to the sample size and fill the first time point with the correct value
void normalize_and_fill_time_zero( CluCorrTen& CCT, const size_t num_Samples, const std::vector<RealType>& spin_float_list )
{
  CCT.iterate( [&]( CorrTen& CT, auto& ij )
  {
    CT.iterate( [&]( corr::CorrelationVector& corr, auto& alphabeta )
    {
      // first element of diagonal autocorrelation is initialized to S*(S+1)/3, else 0. (default value)
      if( ij[0]==ij[1] )
      {
        RealType max_val = spin_float_list[ij[0]]*(spin_float_list[ij[0]]+1.)/3.; // S(S+1)/3 is the value of av( S^alpha^2 )
        if( alphabeta[0]==alphabeta[1] )
        {
          corr.at(0) = max_val;
        }
      }

      // the other elements are normalized
      std::for_each( corr.begin()+1, corr.end(), [num_Samples]( RealType& c )
      {
        c = c / num_Samples;
      } );
    } );
  } );
}

}