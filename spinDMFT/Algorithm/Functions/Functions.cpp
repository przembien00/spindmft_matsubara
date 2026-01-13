#include"Functions.h"
#include<memory>
#include<Globals/MPI_Types.h>

#include<Physics/CFET.h>
namespace cfet = Physics::CFET;

#include<Physics/Spin.h>
namespace sp = Physics::Spin;

#include<Standard_Algorithms/Numerics.h>
namespace num = Standard_Algorithms::Numerics;

#include"Error_Handling.h"
namespace error = spinDMFT::Functions::Error_Handling;

namespace spinDMFT::Functions
{

/* initializes several important matrices on the Hilbert space */
void initialize_matrices( const ps::ParameterSpace& pspace )
{
    // initialize 0 and 1
    ZERO = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( pspace.num_HilbertSpaceDimension, pspace.num_HilbertSpaceDimension );
    IDENTITY = blaze::IdentityMatrix<ComplexType,blaze::rowMajor>( pspace.num_HilbertSpaceDimension );

    // initialize the spin matrices
    sp::write_spin_matrices( pspace.spin_float, S_X, S_Y, S_Z );

    // add extra interaction, e.g., quadrupolar interaction
    H_REST = pspace.extra_interaction.m_term( ZERO, S_X, S_Y, S_Z );
}

/* generates and returns the initial spin correlations */
CorrTen generate_initial_spin_correlations( const ps::ParameterSpace& pspace )
{
    CorrTen environment_CT{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    
    if( !pspace.load_initial_spin_correlations ) // initialize by predefined function
    {
        // build the general correlations from the given functions on the time discretizationn
        corr::CorrelationVector DC{pspace.init_diag_corr.create_discretization(pspace.delta_t, pspace.num_TimePoints, pspace.spin_float)};
        corr::CorrelationVector NDC{pspace.init_nondiag_corr.create_discretization(pspace.delta_t, pspace.num_TimePoints, pspace.spin_float)};

        // insert them into the environment correlation tensor
        environment_CT.iterate( [&](Corr &C, const auto &alphabeta)
        {
            C = (alphabeta[0]==alphabeta[1]) ? DC : NDC; // set the diagonal correlations to DC and nondiagonal correlations to NDC
        } );
    }
    else // initialize by imported linearized data and extrapolate if necessary
    {
        for( uint dir=0; dir<environment_CT.size(); ++dir ) // delinearize imported correlation tensors
        {
            size_t start = dir*pspace.old_num_TimePoints;
            size_t end   = start + pspace.old_num_TimePoints;
            std::vector<RealType> correlation{};
            std::copy( pspace.initial_correlations_linearized.cbegin() + start, pspace.initial_correlations_linearized.cbegin() + end, std::back_inserter(correlation) );
            if( pspace.extrapolate_initial_spin_correlations ) // extrapolate if necessary because the discretization doesn't match
            {
                correlation = num::extrapolate( correlation, pspace.num_TimePoints, pspace.old_delta_t, pspace.delta_t );
            }
            environment_CT[dir] = Corr{ std::move(correlation) };
        }
    }
    return environment_CT;
}

/* transform spin correlations self-consistently to meanfield correlations and write the latter into a covariance matrix */
std::tuple<FieldVector, fmvg::FrequencyCovarianceMatrix> self_consistent_equations( const ps::ParameterSpace& pspace, const CorrTen& spin_correlations, FieldVector& spin_expval )
{
    FieldVector rotated_expval = pspace.spin_model.coupling_matrix * spin_expval; // rotated spin expectation value D^ab * <S^b>
    FieldVector MF_mean = pspace.JL * rotated_expval; // mfav( V^a(t) ) = JL * sum_{b} D^ab * <S^b>
    /* this scheme applies the self-consistency equations
    mfav( V^a(t)V^b(0) ) =  JQ^2 * sum_{cd} D^ac * D^bd * <S^c(t) S^d(0)> + C^ab
    mfav( V^a(t) ) = JL * sum_{b} D^ab * <S^b>
    where JL is the linear coupling constant, JQ is the quadratic coupling constant, D is the spin-model rotation matrix and C is a potential noise 
    spin_corr represents the correlation tensor of spin correlations <S^c(t) S^d(0)> 
    ipair_new is {a,b} and ipair is {c,d} */
    auto transformation_scheme = [&]( const CorrTen& spin_corr, const IndexPair& ipair_new ) -> Corr
    {
        IndexPairList dirs = ten::get_all_direction_pairs();
        return std::pow(pspace.JQ,2) * std::accumulate( dirs.cbegin(), dirs.cend(), Corr{pspace.num_TimePoints}, [&]( Corr sum, const IndexPair& ipair )
        {
            auto D = pspace.spin_model.coupling_matrix;
            return sum + D(ipair_new[0],ipair[0]) * D(ipair_new[1],ipair[1]) * spin_corr(ipair[0],ipair[1]) + pspace.noise.m_variance_in(ipair_new[0],ipair_new[1]) + (- D(ipair_new[0],ipair[0]) * D(ipair_new[1],ipair[1]) * spin_expval[ipair[0]]*spin_expval[ipair[1]]); // D^ac * D^bd * <S^c(t) S^d(0)> - D^ac * D^bd * <S^c><S^d> + C^ab
        } ); // += JQ^2 * ...
    };
    
    // create the mean-field correlation tensor using the defined transformation scheme
    auto V_dot_V_tensor = std::make_shared<CorrTen>( spin_correlations, transformation_scheme );

    // transform CorrelationTensor to CovarianceMatrixBlocks and return
    fmvg::FrequencyCovarianceMatrix MF_Cov{ V_dot_V_tensor, pspace.correlation_symmetry_type };

    return std::make_tuple(MF_mean, MF_Cov);
}

/* helper function for spin-vector scalar product  */
Operator spin_scalar_product( const FieldVector& v )
{
    return v[0]*S_X + v[1]*S_Y + v[2]*S_Z;
}

/* computes CFET for a Hilbert space containing a single spin according to the Hamiltonian H(t) = Q(t) * S + A
T is a matrix containing the the vector Q for two adjacent time steps in its columns, T = (Q(t+) Q(t))
A is an operator which is considered only for S>1/2 (because for S=1/2 it is either constant or contains elements that can equivalently be added to the first term of H)
note that a time dependence of A is not yet considered in this implementation */
Operator CFET_for_single_spin( const std::string spin, const FieldVector& Q_new, const FieldVector& Q_old , const RealType& dt, const Operator& A )
{
    if( spin == "1/2" )
    {
        return cfet::CFET4opt_for_single_spin_one_half( Q_new, Q_old, dt, S_X, S_Y, S_Z, IDENTITY );
    }
    else
    {
        Operator H_new = spin_scalar_product(Q_new) + A;
        Operator H_old = spin_scalar_product(Q_old) + A;
        return cfet::CFET2( H_new, H_old, dt, ZERO, IDENTITY );
    }
}

/* computes the propagators from 0 to n*delta_t and returns them as a vector
according to Tm( 0, 0 ) = scheme.read_x( n, sample ); */
std::tuple<RealType, std::vector<Operator>, std::vector<Operator>> compute_propagators( const ps::ParameterSpace& pspace, std::vector<FieldVector>& noise_vectors, FieldVector& mf_mean )
{
    std::vector<Operator> shortstep_propagators{};
    Operator U_step{};
    FieldVector vV_old = noise_vectors[0]; // read the old noise vector (t_0)
    for( size_t n = 1; n < pspace.num_TimePoints - 1; ++n ) // propagator to n*delta_t
    {
        auto vV_new = noise_vectors[n]; // read the new noise vector (t_n)

        // A.) Compute the Short-Step Propagator by CFET 4 opt
        if( pspace.B.m_name == "none" ) // Q = V
        {
            U_step = CFET_for_single_spin( pspace.spin, vV_new + mf_mean, vV_old + mf_mean, pspace.delta_t, H_REST );
        }
        else // Q = V + B
        {
            U_step = CFET_for_single_spin( pspace.spin, vV_new + pspace.B.m_h + mf_mean, vV_old + pspace.B.m_h + mf_mean, pspace.delta_t, H_REST );
        }

        // B.) Compute the new Propagator and add it
        shortstep_propagators.emplace_back( U_step );

        vV_old = vV_new; // update the old noise vector (t_n-1 -> t_n)
    }

    if( pspace.B.m_name == "none" ) // Q = V
    {
        U_step = CFET_for_single_spin( pspace.spin, noise_vectors[0] + mf_mean, vV_old + mf_mean, pspace.delta_t, H_REST );
    }
    else // Q = V + B
    {
        U_step = CFET_for_single_spin( pspace.spin, noise_vectors[0] + pspace.B.m_h + mf_mean, vV_old + pspace.B.m_h + mf_mean, pspace.delta_t, H_REST );
    }
    shortstep_propagators.emplace_back( U_step );


    std::vector<Operator> propagators{ IDENTITY }; // initialization; propagator to 0 (identity)
    std::vector<Operator> propagators_inv{ IDENTITY }; // initialization; propagator to 0 (identity)
    stda::for_n_each( [&propagators, &propagators_inv]( const Operator& U_step, const Operator& U_step_inv )
    {
        propagators.emplace_back( U_step * propagators.back() );
        propagators_inv.emplace_back( propagators_inv.back() * U_step_inv );
    },
    shortstep_propagators.cbegin(), shortstep_propagators.cend(), shortstep_propagators.crbegin() );

    return std::make_tuple( std::real( blaze::trace( propagators.back() ) ), propagators, propagators_inv );
}


// compute the spin matrices over time via S(t) = U_dagger * S * U
void compute_S_of_t( const ps::ParameterSpace& pspace, const std::vector<Operator>& propagators, const std::vector<Operator>& propagators_inv, 
    std::vector<Operator>& S_x_of_t, std::vector<Operator>& S_y_of_t, std::vector<Operator>& S_z_of_t )
{
    S_x_of_t.clear(); S_x_of_t.resize( pspace.num_TimePoints );
    S_y_of_t.clear(); S_y_of_t.resize( pspace.num_TimePoints );
    S_z_of_t.clear(); S_z_of_t.resize( pspace.num_TimePoints );
    stda::for_n_each( []( const Operator& U, const Operator& U_inv, Operator& S_x, Operator& S_y, Operator& S_z )
    {
        S_x = U_inv * S_X * U;
        S_y = U_inv * S_Y * U;
        S_z = U_inv * S_Z * U;
    },
    propagators.cbegin(), propagators.cend(), propagators_inv.crbegin(), S_x_of_t.begin(), S_y_of_t.begin(), S_z_of_t.begin() );
}

// compute the autocorrelation at two specific times t1 and t2 in dependence of alphabeta
ComplexType correlation( const std::vector<Operator>& vS_t1, const std::vector<Operator>& vS_t2, const ten::IndexPair& direction )
{
    return blaze::trace( vS_t1[direction[0]] * vS_t2[direction[1]] ); // <S^alpha(t_1)S^beta(t_2)>
}

// compute the spin correlations for a single MC sample
void compute_spin_correlations( const ps::ParameterSpace& pspace, rtd::RunTimeData& rtdata,
    CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, FieldVector& spin_expval, RealType& Z,
    const std::vector<Operator>& S_x_of_t,
    const std::vector<Operator>& S_y_of_t, 
    const std::vector<Operator>& S_z_of_t )
{
    std::vector<Operator> vS_at_zero = {S_X, S_Y, S_Z}; // spin vector at time 0

    spin_expval[0] += std::real( blaze::trace(S_x_of_t[0]) ); // <Sx>
    spin_expval[1] += std::real( blaze::trace(S_y_of_t[0]) ); // <Sy>
    spin_expval[2] += std::real( blaze::trace(S_z_of_t[0]) ); // <Sz>

    // loop over all directions, then over all times
    stda::for_n_each( [&]( auto& re_g_alphabeta_of_t, auto& im_g_alphabeta_of_t, auto& re_gsq_alphabeta_of_t, auto& im_gsq_alphabeta_of_t, auto& re_cov_alphabeta_of_t, auto& im_cov_alphabeta_of_t, auto& alphabeta )
    {    
        stda::for_n_each( [&]( auto& re_g_alphabeta_at_t, auto& im_g_alphabeta_at_t, auto& re_gsq_alphabeta_at_t, auto& im_gsq_alphabeta_at_t, auto& re_cov_alphabeta_at_t, auto& im_cov_alphabeta_at_t, auto& S_x_at_t, auto& S_y_at_t, auto& S_z_at_t )
        {
            std::vector<Operator> vS_at_t = {S_x_at_t, S_y_at_t, S_z_at_t}; // spin vector at time t
            ComplexType single_sample_result = correlation( vS_at_t, vS_at_zero, alphabeta ); // compute the correlation for the single sample
            im_g_alphabeta_at_t += std::imag(single_sample_result);
            re_g_alphabeta_at_t += std::real(single_sample_result); // add the single sample results to the correlations
            re_gsq_alphabeta_at_t += std::pow(std::real(single_sample_result),2); // add the squared single sample results to the sample sqsums
            im_gsq_alphabeta_at_t += std::pow(std::imag(single_sample_result),2); // add the squared single sample results to the sample sqsums
            re_cov_alphabeta_at_t += std::real(single_sample_result) * Z; // add the covariance part for real part
            im_cov_alphabeta_at_t += std::imag(single_sample_result) * Z; // add the covariance part for imag part
        },
        re_g_alphabeta_of_t.begin(), re_g_alphabeta_of_t.end(), im_g_alphabeta_of_t.begin(), re_gsq_alphabeta_of_t.begin(), im_gsq_alphabeta_of_t.begin(), re_cov_alphabeta_of_t.begin(), im_cov_alphabeta_of_t.begin(), S_x_of_t.cbegin(), S_y_of_t.cbegin(), S_z_of_t.cbegin() ); // loop over time
    }, spin_correlations_R.begin(), spin_correlations_R.end(), spin_correlations_I.begin(), rtdata.sample_sqsum_Re.begin(), rtdata.sample_sqsum_Im.begin(), rtdata.sample_cov_Re.begin(), rtdata.sample_cov_Im.begin(), spin_correlations_R.m_direction_pairs.begin() ); 
}

// sum the results of all cores and broadcast the sum to all cores with MPI_Allreduce 
void MPI_share_results( rtd::RunTimeData& rtdata, CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, FieldVector& spin_expval, RealType& partition_function )
{
    // share mag moment results
    std::vector<RealType> rcv_buf(3);
    std::vector<RealType> snd_buf = { spin_expval[0], spin_expval[1], spin_expval[2] };
    MPI_Allreduce( snd_buf.data(), rcv_buf.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    spin_expval = { rcv_buf[0], rcv_buf[1], rcv_buf[2] };

    // share correlation results
    std::for_each( spin_correlations_R.begin(), spin_correlations_R.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share correlation results
    std::for_each( spin_correlations_I.begin(), spin_correlations_I.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share sample sqsum results
    std::for_each( rtdata.sample_sqsum_Re.begin(), rtdata.sample_sqsum_Re.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share sample sqsum results
    std::for_each( rtdata.sample_sqsum_Im.begin(), rtdata.sample_sqsum_Im.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share covariance results
    std::for_each( rtdata.sample_cov_Re.begin(), rtdata.sample_cov_Re.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );

    // share covariance results
    std::for_each( rtdata.sample_cov_Im.begin(), rtdata.sample_cov_Im.end(), []( corr::CorrelationVector& spin_c ) 
    {
        std::vector<RealType> rcv_buf( spin_c.size() ); 
        MPI_Allreduce( spin_c.data(), rcv_buf.data(), spin_c.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        spin_c = rcv_buf;
    } );    

    std::vector<RealType> send_buf = { partition_function };
    std::vector<RealType> receive_buf( 1 );
    // share partition function results
    MPI_Allreduce( send_buf.data(), receive_buf.data(), 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    partition_function = receive_buf.at(0);

    send_buf = { rtdata.Z_sqsum };
    // share partition function results
    MPI_Allreduce( send_buf.data(), receive_buf.data(), 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    rtdata.Z_sqsum = receive_buf.at(0);
}

// normalize the results of the MC-simulation by dividing them with the partition function
void normalize( rtd::RunTimeData& rtdata, CorrTen& spin_correlations, RealType& partition_function )
{
    // normalize correlation results
    std::for_each( spin_correlations.begin(), spin_correlations.end(), [&partition_function]( corr::CorrelationVector& spin_c ) 
    {
        spin_c *= 1/partition_function;
    } );
}

};