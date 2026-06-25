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
    // Only the lower half [0, beta/2] is built: the rest follows from the reflection
    // symmetry g^{ab}(beta-tau) = g^{ba}(tau)^* and is reconstructed once per iteration
    // (see symmetrize_upper_half). n_half counts the points in [0, beta/2] inclusive,
    // and the loop over propagators is truncated accordingly. The mirror partner of
    // point k is propagators_inv[N-1-k], supplied by propagators_inv.crbegin().
    const size_t n_half = pspace.num_TimePoints - pspace.num_TimePoints / 2;
    S_x_of_t.clear(); S_x_of_t.resize( n_half );
    S_y_of_t.clear(); S_y_of_t.resize( n_half );
    S_z_of_t.clear(); S_z_of_t.resize( n_half );
    stda::for_n_each( []( const Operator& U, const Operator& U_inv, Operator& S_x, Operator& S_y, Operator& S_z )
    {
        S_x = U_inv * S_X * U;
        S_y = U_inv * S_Y * U;
        S_z = U_inv * S_Z * U;
    },
    propagators.cbegin(), propagators.cbegin() + n_half, propagators_inv.crbegin(), S_x_of_t.begin(), S_y_of_t.begin(), S_z_of_t.begin() );
}

// compute the autocorrelation at two specific times t1 and t2 in dependence of alphabeta
ComplexType correlation( const std::vector<Operator>& vS_t1, const std::vector<Operator>& vS_t2, const ten::IndexPair& direction )
{
    return blaze::trace( vS_t1[direction[0]] * vS_t2[direction[1]] ); // <S^alpha(t_1)S^beta(t_2)>
}

// Accumulate per-sample contributions from one pCN chain state.
// For each direction (alpha, beta) and each time t we add
//     o = Tr[ U(beta-t; V) S^alpha U(t; V) S^beta ] / Z(V)
// to the running sum (spin_correlations), and o^2 to spin_expval_sqsum / sample_sqsum
// for a textbook standard-error-of-the-mean estimate.
void compute_spin_correlations( rtd::RunTimeData& rtdata,
    CorrTen& spin_correlations_R, CorrTen& spin_correlations_I,
    FieldVector& spin_expval, FieldVector& spin_expval_sqsum,
    const RealType Z,
    const std::vector<Operator>& S_x_of_t,
    const std::vector<Operator>& S_y_of_t,
    const std::vector<Operator>& S_z_of_t )
{
    std::vector<Operator> vS_at_zero = {S_X, S_Y, S_Z};
    const RealType inv_Z = ( Z != RealType{0.} ) ? RealType{1.} / Z : RealType{0.};

    const RealType s0 = std::real( blaze::trace(S_x_of_t[0]) ) * inv_Z;
    const RealType s1 = std::real( blaze::trace(S_y_of_t[0]) ) * inv_Z;
    const RealType s2 = std::real( blaze::trace(S_z_of_t[0]) ) * inv_Z;
    spin_expval[0] += s0;       spin_expval[1] += s1;       spin_expval[2] += s2;
    spin_expval_sqsum[0] += s0*s0; spin_expval_sqsum[1] += s1*s1; spin_expval_sqsum[2] += s2*s2;
    rtdata.cur_block_spin[0] += s0; rtdata.cur_block_spin[1] += s1; rtdata.cur_block_spin[2] += s2; // batch-means accumulation

    // Only the lower-half time points [0, beta/2] are accumulated; the inner loop is
    // driven by the (half-sized) S^a(t) caches, so re_g_of_t etc. are filled only up to
    // beta/2. The upper half is filled afterwards by symmetrize_upper_half.
    stda::for_n_each( [&]( auto& re_g_of_t, auto& im_g_of_t, auto& re_gsq_of_t, auto& im_gsq_of_t, auto& re_blk_of_t, auto& im_blk_of_t, auto& alphabeta )
    {
        stda::for_n_each( [&]( auto& S_x_at_t, auto& S_y_at_t, auto& S_z_at_t, auto& re_g_at_t, auto& im_g_at_t, auto& re_gsq_at_t, auto& im_gsq_at_t, auto& re_blk_at_t, auto& im_blk_at_t )
        {
            std::vector<Operator> vS_at_t = {S_x_at_t, S_y_at_t, S_z_at_t};
            ComplexType single_sample_result = correlation( vS_at_t, vS_at_zero, alphabeta );
            const RealType re_o = std::real(single_sample_result) * inv_Z;
            const RealType im_o = std::imag(single_sample_result) * inv_Z;
            re_g_at_t   += re_o;
            im_g_at_t   += im_o;
            re_gsq_at_t += re_o * re_o;
            im_gsq_at_t += im_o * im_o;
            re_blk_at_t += re_o;        // batch-means accumulation (current block sum)
            im_blk_at_t += im_o;
        },
        S_x_of_t.cbegin(), S_x_of_t.cend(), S_y_of_t.cbegin(), S_z_of_t.cbegin(),
        re_g_of_t.begin(), im_g_of_t.begin(), re_gsq_of_t.begin(), im_gsq_of_t.begin(),
        re_blk_of_t.begin(), im_blk_of_t.begin() );
    },
    spin_correlations_R.begin(), spin_correlations_R.end(), spin_correlations_I.begin(),
    rtdata.sample_sqsum_Re.begin(), rtdata.sample_sqsum_Im.begin(),
    rtdata.cur_block_Re.begin(), rtdata.cur_block_Im.begin(),
    spin_correlations_R.m_direction_pairs.begin() );
}

// helper: MPI_Allreduce a CorrTen of per-(alpha,beta,t) sums
static void mpi_sum_corrtensor( CorrTen& CT )
{
    std::for_each( CT.begin(), CT.end(), []( corr::CorrelationVector& v )
    {
        std::vector<RealType> rcv( v.size() );
        MPI_Allreduce( v.data(), rcv.data(), v.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
        v = rcv;
    } );
}

// in-place MPI_Allreduce of a 3-component FieldVector
static void mpi_sum_fieldvector( FieldVector& v )
{
    FieldVector rcv{ 0., 0., 0. };
    MPI_Allreduce( v.data(), rcv.data(), 3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD );
    v = rcv;
}

// sum the per-core accumulators across all MPI cores
void MPI_share_results( rtd::RunTimeData& rtdata, CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, FieldVector& spin_expval, FieldVector& spin_expval_sqsum )
{
    mpi_sum_fieldvector( spin_expval );
    mpi_sum_fieldvector( spin_expval_sqsum );
    rtdata.spin_expval_sqsum = spin_expval_sqsum;

    mpi_sum_corrtensor( spin_correlations_R );
    mpi_sum_corrtensor( spin_correlations_I );
    mpi_sum_corrtensor( rtdata.sample_sqsum_Re );
    mpi_sum_corrtensor( rtdata.sample_sqsum_Im );
}

// divide the summed-over-samples accumulator by N to obtain the Monte-Carlo mean
void normalize( CorrTen& spin_correlations, RealType N )
{
    const RealType inv_N = RealType{1.} / N;
    std::for_each( spin_correlations.begin(), spin_correlations.end(), [inv_N]( corr::CorrelationVector& spin_c )
    {
        spin_c *= inv_N;
    } );
}

// Reconstruct the upper-half (tau > beta/2) correlations from the computed lower half via
// the per-component reflection symmetry
//     g^{ab}(beta - tau) = g^{ab}(tau)^*   <=>   Re even, Im odd about beta/2,
// i.e.  Re g^{ab}(beta-tau) =  Re g^{ab}(tau),    Im g^{ab}(beta-tau) = -Im g^{ab}(tau).
// This follows from trace cyclicity g^{ab}(beta-tau) = g^{ba}(tau) together with
// Hermiticity g^{ab}(tau)^* = g^{ba}(tau). It is a SAME-component relation -- no transpose:
// using the transposed pair g^{ba} would double-apply the conjugate and (wrongly) make the
// imaginary part symmetric.  The mirror partner of grid point k (tau_k = k*delta_t) is
// point N-1-k (beta - tau_k), always in the computed lower half. The per-point standard
// errors mirror with no sign change, since each reconstructed value equals its mirror
// source deterministically.
void symmetrize_upper_half( CorrTen& Re, CorrTen& Im, CorrTen& Re_std, CorrTen& Im_std, CorrTen& tau_Re, CorrTen& tau_Im )
{
    const size_t n_pairs = Re.size();
    const size_t N = Re[0].size();
    const size_t n_half = N - N / 2; // computed points in [0, beta/2]

    for( size_t p = 0; p < n_pairs; ++p )
    {
        for( size_t k = n_half; k < N; ++k )
        {
            const size_t m = N - 1 - k; // mirror index, in the computed lower half
            Re[p][k]     =  Re[p][m];
            Im[p][k]     = -Im[p][m];
            Re_std[p][k] =  Re_std[p][m];
            Im_std[p][k] =  Im_std[p][m];
            tau_Re[p][k] =  tau_Re[p][m]; // tau_int is a chain property -> mirror with no sign change
            tau_Im[p][k] =  tau_Im[p][m];
        }
        // The self-paired midpoint tau = beta/2 (present for an odd number of time points)
        // lies in the computed lower half and keeps its raw sampled value: its imaginary
        // part is zero only on the ensemble average, so forcing it to zero per iteration
        // would bias the estimator rather than improve it.
    }
}

// --- pCN sampler helpers ---

// Inverse DFT mirroring fmvg::FrequencyNoiseVectors::fourier_back_transform but acting
// on a single sample given in the diagonal frequency basis.
std::vector<FieldVector> diag_freq_to_time( const std::vector<fmvg::Vector>& V_diag,
                                            const fmvg::OrthogonalTransformationList& ortho,
                                            FFT::InverseRealDFT& fft )
{
    size_t N = V_diag.size();
    // 1.) rotate each block back into the full Matsubara basis
    std::vector<fmvg::Vector> V_full( N );
    for( size_t b = 0; b < N; ++b )
    {
        V_full[b] = ortho[b] * V_diag[b];
    }
    // 2.) inverse real DFT, performed independently for each of the 3 field components in
    //     O(N log N) via FFTW. FFT::InverseRealDFT reproduces the project's half-complex
    //     convention exactly (DC ~ 1/sqrt(N), cosine/sine modes ~ 1/sqrt(N/2); see
    //     RealFFT.h), so this is a drop-in replacement for the former O(N^2) explicit sums.
    //     Each component's frequency coefficients are packed contiguously, transformed, then
    //     scattered back into the time-domain FieldVectors.
    std::vector<FieldVector> V_time( N );
    std::vector<RealType> comp_freq( N ), comp_time( N );
    for( size_t c = 0; c < 3; ++c )
    {
        for( size_t n = 0; n < N; ++n ) comp_freq[n] = V_full[n][c];
        fft.transform( comp_freq.data(), comp_time.data() );
        for( size_t t = 0; t < N; ++t ) V_time[t][c] = comp_time[t];
    }
    return V_time;
}

// ======================== PCNChain implementation ========================

void PCNChain::draw_prior_into( std::vector<fmvg::Vector>& V_diag )
{
    for( size_t b = 0; b < V_diag.size(); ++b )
    {
        fmvg::Vector v;
        for( size_t i = 0; i < 3; ++i ) v[i] = m_prior_dists[b][i]( m_engine );
        V_diag[b] = v;
    }
}

void PCNChain::rebuild_propagators_and_S()
{
    std::vector<FieldVector> V_time = diag_freq_to_time( m_V_diag, m_ortho, m_fft );
    std::tie( m_Z, m_propagators, m_propagators_inv ) = compute_propagators( m_pspace, V_time, m_mf_mean );
    compute_S_of_t( m_pspace, m_propagators, m_propagators_inv, m_S_x, m_S_y, m_S_z );
}

// Common construction: store references and build the per-block prior distributions
// from the (possibly truncated) eigenvalues.
PCNChain::PCNChain( const ps::ParameterSpace& pspace,
                    const fmvg::EigenValuesList& eig,
                    const fmvg::OrthogonalTransformationList& ortho,
                    const FieldVector& mf_mean,
                    RealType pcn_beta,
                    std::mt19937& engine )
    : m_pspace( pspace ),
      m_ortho( ortho ),
      m_mf_mean( mf_mean ),
      m_engine( engine ),
      m_pcn_beta( pcn_beta ),
      m_pcn_root( std::sqrt( std::max( RealType{0.}, RealType{1.} - pcn_beta * pcn_beta ) ) ),
      m_prior_dists( ortho.size() ),
      m_uniform01( RealType{0.}, RealType{1.} ),
      m_fft( ortho.size() ),
      m_V_diag( ortho.size() )
{
    // Per-block, per-component prior distribution xi ~ N(0, sqrt(eig)). Truncated
    // (non-positive) eigenvalues give a degenerate proposal at zero, freezing those
    // modes -- consistent with the prior having no support there.
    for( size_t b = 0; b < eig.size(); ++b )
    {
        for( size_t i = 0; i < 3; ++i )
        {
            const RealType sigma = ( eig[b][i] > RealType{0.} ) ? std::sqrt( eig[b][i] ) : RealType{0.};
            m_prior_dists[b][i] = std::normal_distribution<RealType>( RealType{0.}, sigma );
        }
    }

    // Cold start: draw an initial state; retry on Z <= 0 (a numerically degenerate region).
    draw_prior_into( m_V_diag );
    rebuild_propagators_and_S();
    for( size_t retries = 0; !(m_Z > RealType{0.}) && retries < 32; ++retries )
    {
        draw_prior_into( m_V_diag );
        rebuild_propagators_and_S();
    }
}

// Warm start: same construction, but seed V_diag from the supplied state instead of
// drawing from the prior. Frozen modes (eig <= 0) are zeroed for consistency.
PCNChain::PCNChain( const ps::ParameterSpace& pspace,
                    const fmvg::EigenValuesList& eig,
                    const fmvg::OrthogonalTransformationList& ortho,
                    const FieldVector& mf_mean,
                    RealType pcn_beta,
                    std::mt19937& engine,
                    std::vector<fmvg::Vector> V_diag_initial )
    : PCNChain( pspace, eig, ortho, mf_mean, pcn_beta, engine )
{
    // Replace the cold-start state with the supplied warm state, zeroing components
    // where the new prior has no support, then rebuild Z, propagators, and S^a(t).
    for( size_t b = 0; b < m_V_diag.size(); ++b )
    {
        for( size_t i = 0; i < 3; ++i )
        {
            m_V_diag[b][i] = ( eig[b][i] > RealType{0.} ) ? V_diag_initial[b][i] : RealType{0.};
        }
    }
    rebuild_propagators_and_S();
}

std::vector<fmvg::Vector> PCNChain::V_full() const
{
    std::vector<fmvg::Vector> result( m_V_diag.size() );
    for( size_t b = 0; b < m_V_diag.size(); ++b )
    {
        result[b] = m_ortho[b] * m_V_diag[b];
    }
    return result;
}

bool PCNChain::step()
{
    // 1.) pCN proposal in the diagonal basis: V' = sqrt(1-beta^2) V + beta xi, xi ~ p(V)
    std::vector<fmvg::Vector> V_diag_prop( m_V_diag.size() );
    for( size_t b = 0; b < m_V_diag.size(); ++b )
    {
        fmvg::Vector xi;
        for( size_t i = 0; i < 3; ++i ) xi[i] = m_prior_dists[b][i]( m_engine );
        V_diag_prop[b] = m_pcn_root * m_V_diag[b] + m_pcn_beta * xi;
    }

    // 2.) evaluate target at the proposal
    std::vector<FieldVector> V_time_prop = diag_freq_to_time( V_diag_prop, m_ortho, m_fft );
    auto [Z_prop, propagators_prop, propagators_inv_prop] = compute_propagators( m_pspace, V_time_prop, m_mf_mean );

    // 3.) pCN acceptance: the Gaussian prior factors cancel, leaving min(1, Z'/Z).
    //     Z <= 0 is treated as a forbidden region.
    if( !(Z_prop > RealType{0.}) || !(m_Z > RealType{0.}) ) return false;
    const RealType log_alpha = std::log( Z_prop ) - std::log( m_Z );
    if( std::log( m_uniform01( m_engine ) ) >= log_alpha ) return false;

    // 4.) on acceptance, swap in the new state and refresh the cached S^a(t)
    m_V_diag = std::move( V_diag_prop );
    m_Z = Z_prop;
    m_propagators = std::move( propagators_prop );
    m_propagators_inv = std::move( propagators_inv_prop );
    compute_S_of_t( m_pspace, m_propagators, m_propagators_inv, m_S_x, m_S_y, m_S_z );
    return true;
}

};