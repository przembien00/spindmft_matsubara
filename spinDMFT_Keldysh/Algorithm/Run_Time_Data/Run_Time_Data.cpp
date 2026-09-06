#include"Run_Time_Data.h"

#include<algorithm>
#include<array>
#include<cmath>
#include<iostream>
#include<limits>
#include<stdexcept>

#include<mpi.h>
#include<Standard_Algorithms/Print_Routines.h>

#include"Ratio_Statistics.h"

namespace spinDMFT::Run_Time_Data
{
namespace print = Print_Routines;
namespace
{
MPI_Datatype mpi_real_type()
{
#ifdef USE_FLOAT
    return MPI_FLOAT;
#else
    return MPI_DOUBLE;
#endif
}

void mpi_sum_buffer( std::vector<RealType>& values )
{
    if( values.size()>static_cast<size_t>(std::numeric_limits<int>::max()) )
        throw std::overflow_error("MPI reduction buffer exceeds the MPI int count range");
    std::vector<RealType> reduced(values.size());
    MPI_Allreduce(values.data(),reduced.data(),static_cast<int>(values.size()),
                  mpi_real_type(),MPI_SUM,MPI_COMM_WORLD);
    values=std::move(reduced);
}

bool finite( ComplexType value )
{ return std::isfinite(std::real(value))&&std::isfinite(std::imag(value)); }
}

RunTimeData::ComplexBlockSums::ComplexBlockSums(
    char symmetry, size_t num_imaginary_edge_points, size_t num_real_points )
    : correlations(symmetry,num_imaginary_edge_points,num_real_points),
      mag_Re(num_real_points,MagVec{symmetry}),
      mag_Im(num_real_points,MagVec{symmetry}),
      closure(num_real_points,ComplexType{}),
      closure_abs(num_real_points,RealType{})
{}

RunTimeData::RunTimeData( const ps::ParameterSpace& pspace, int my_rank )
    : contour_sample_stds(pspace.correlation_symmetry_type,pspace.num_TimePoints,
                          pspace.num_RealTimePoints),
      contour_tau_int(pspace.correlation_symmetry_type,pspace.num_TimePoints,
                      pspace.num_RealTimePoints),
      m_my_rank(my_rank),m_num_print_digits(pspace.num_PrintDigits),
      m_symmetry(pspace.correlation_symmetry_type),
      m_num_imaginary_edge_points(pspace.num_TimePoints),
      m_num_real_points(pspace.num_RealTimePoints),
      m_num_samples_per_core(pspace.num_SamplesPerCore),
      m_num_cores(pspace.get_num_Cores()),m_seed(pspace.seed),
      m_self_consistency(pspace.self_consistency),
      m_harmonic_bath(pspace.uses_harmonic_bath()),
      m_pcn(pspace.sampling_strategy=="pcn"),
      m_antithetic_pairs(pspace.antithetic_pairs
                         &&pspace.sampling_strategy=="independent"),
      m_closed_contour_observable_normalization(
          pspace.correlation_normalization=="closed-contour"),
      m_iteration_error_sigma_threshold(pspace.iteration_error_sigma_threshold),
      m_iteration_limit(pspace.Iteration_Limit),
      m_covariance_tolerance(pspace.covariance_tolerance),
      m_branch_identity_tolerance(pspace.branch_identity_tolerance),
      m_takagi_tolerance(pspace.takagi_tolerance),
      m_minimum_phase_magnitude(pspace.minimum_phase_magnitude),
      m_denominator_constancy_tolerance(pspace.denominator_constancy_tolerance),
      m_imaginary_magnetization_sigma(pspace.imaginary_magnetization_sigma),
      m_partition_imaginary_tolerance(pspace.partition_imaginary_tolerance),
      m_corr_directions(CorrTen{pspace.correlation_symmetry_type,
                                pspace.num_TimePoints}.get_direction_pairs()),
      m_mag_directions(mag::determine_magnetization_directions(
          pspace.correlation_symmetry_type)),
      m_total(pspace.correlation_symmetry_type,pspace.num_TimePoints,
              pspace.num_RealTimePoints),
      m_sample_squares(pspace.correlation_symmetry_type,pspace.num_TimePoints,
                       pspace.num_RealTimePoints),
      m_current_sample(pspace.correlation_symmetry_type,pspace.num_TimePoints,
                       pspace.num_RealTimePoints)
{
    if( m_num_samples_per_core==0 )
        throw std::invalid_argument("numSamplesPerCore must be positive");
    m_num_blocks=std::min(std::max<size_t>(1,pspace.num_blocks),m_num_samples_per_core);
    while( m_num_blocks>1
        &&(m_num_samples_per_core%m_num_blocks!=0
           ||(m_antithetic_pairs
              &&(m_num_samples_per_core/m_num_blocks)%2!=0)) )
        --m_num_blocks;
    m_samples_per_block=m_num_samples_per_core/m_num_blocks;
    if( m_antithetic_pairs&&m_samples_per_block%2!=0 )
        throw std::logic_error(
            "antithetic-pair jackknife blocks must contain complete pairs");
    begin_iteration_accumulation();
}

void RunTimeData::begin_iteration_accumulation()
{
    m_samples_seen=0; m_current_block=0;
    m_total=ComplexBlockSums{m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    m_sample_squares=CorrelationSet{
        m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    m_current_sample=CorrelationSet{
        m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    m_blocks.assign(m_num_blocks,ComplexBlockSums{
        m_symmetry,m_num_imaginary_edge_points,m_num_real_points});
}

void RunTimeData::begin_sample( ComplexType Z, RealType observable_normalization )
{
    if( m_samples_seen>=m_num_samples_per_core )
        throw std::runtime_error("too many independent samples accumulated in one iteration");
    if( !finite(Z)||!std::isfinite(observable_normalization) )
        throw std::runtime_error("a contour trajectory produced a non-finite sample weight");
    m_current_observable_normalization=observable_normalization;
    m_current_block=std::min(m_samples_seen/m_samples_per_block,m_num_blocks-1);
    m_total.partition+=Z; m_total.partition_abs+=std::abs(Z);
    m_total.partition_abs_sq+=std::norm(Z);
    auto& block=m_blocks[m_current_block];
    block.partition+=Z; block.partition_abs+=std::abs(Z);
    block.partition_abs_sq+=std::norm(Z);
}

void RunTimeData::accumulate_closed_contour_trace( size_t t, ComplexType value )
{
    if( t>=m_num_real_points||!finite(value) )
        throw std::runtime_error("invalid closed-contour trace");
    value*=m_current_observable_normalization;
    m_total.closure[t]+=value; m_blocks[m_current_block].closure[t]+=value;
    m_total.closure_abs[t]+=std::abs(value);
    m_blocks[m_current_block].closure_abs[t]+=std::abs(value);
}

void RunTimeData::accumulate_magnetization( size_t t, size_t c,
                                            ComplexType numerator )
{
    if( t>=m_num_real_points||c>=m_mag_directions.size()||!finite(numerator) )
        throw std::runtime_error("invalid complex magnetization numerator");
    numerator*=m_current_observable_normalization;
    m_total.mag_Re[t][c]+=std::real(numerator);
    m_total.mag_Im[t][c]+=std::imag(numerator);
    m_blocks[m_current_block].mag_Re[t][c]+=std::real(numerator);
    m_blocks[m_current_block].mag_Im[t][c]+=std::imag(numerator);
}

void RunTimeData::accumulate_edge_correlation(
    size_t t,size_t p,size_t tau,ComplexType numerator )
{
    if( !finite(numerator) )
        throw std::runtime_error("non-finite edge contour-correlation numerator");
    numerator*=m_current_observable_normalization;
    const RealType re=std::real(numerator),im=std::imag(numerator);
    m_total.correlations.Re[t][p][tau]+=re;
    m_total.correlations.Im[t][p][tau]+=im;
    m_blocks[m_current_block].correlations.Re[t][p][tau]+=re;
    m_blocks[m_current_block].correlations.Im[t][p][tau]+=im;
    if( m_pcn )
    {
        m_current_sample.Re[t][p][tau]+=re;
        m_current_sample.Im[t][p][tau]+=im;
    }
}

void RunTimeData::end_sample()
{
    if( m_pcn )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    const RealType re=m_current_sample.Re[t][p][tau];
                    const RealType im=m_current_sample.Im[t][p][tau];
                    m_sample_squares.Re[t][p][tau]+=re*re;
                    m_sample_squares.Im[t][p][tau]+=im*im;
                    m_current_sample.Re[t][p][tau]=RealType{};
                    m_current_sample.Im[t][p][tau]=RealType{};
                }
    ++m_samples_seen;
}

ComplexType RunTimeData::edge_value(
    const CorrelationSet& values,size_t t,size_t p,size_t tau )
{ return {values.Re[t][p][tau],values.Im[t][p][tau]}; }

void RunTimeData::record_complex_field_diagnostics(
    RealType symmetry_error,RealType branch_error,RealType factor_error,
    RealType approximation_error,size_t latent_dimension,
    size_t largest_factorization_dimension )
{
    covariance_symmetry_errors.push_back(symmetry_error);
    branch_identity_errors.push_back(branch_error);
    gaussian_factor_reconstruction_errors.push_back(factor_error);
    gaussian_covariance_approximation_errors.push_back(approximation_error);
    gaussian_factor_latent_dimensions.push_back(
        static_cast<unsigned int>(latent_dimension));
    gaussian_largest_factorization_dimensions.push_back(
        static_cast<unsigned int>(largest_factorization_dimension));
}

void RunTimeData::record_pcn_diagnostics(
    const size_t accepted,const size_t proposed,
    const size_t rejected_nonpositive,
    const RealType maximum_relative_imaginary_sampling_weight )
{
    std::array<RealType,3> local{
        static_cast<RealType>(accepted),static_cast<RealType>(proposed),
        static_cast<RealType>(rejected_nonpositive)};
    std::array<RealType,3> global{};
    MPI_Allreduce(local.data(),global.data(),static_cast<int>(global.size()),
                  mpi_real_type(),MPI_SUM,MPI_COMM_WORLD);
    RealType global_imaginary{};
    MPI_Allreduce(&maximum_relative_imaginary_sampling_weight,&global_imaginary,1,
                  mpi_real_type(),MPI_MAX,MPI_COMM_WORLD);
    mh_acceptance_rates.push_back(global[1]>RealType{}?global[0]/global[1]:RealType{});
    mh_nonpositive_rejection_rates.push_back(
        global[1]>RealType{}?global[2]/global[1]:RealType{});
    maximum_relative_imaginary_sampling_weights.push_back(global_imaginary);
}

void RunTimeData::mpi_reduce_and_finalize(
    CorrelationSet& correlations,CorrelationSet& standard_errors )
{
    if( m_samples_seen!=m_num_samples_per_core )
        throw std::runtime_error("the production sweep did not accumulate the requested samples");
    const size_t corr_count=m_num_real_points*m_corr_directions.size()
                           *m_num_imaginary_edge_points;
    const size_t mag_count=m_num_real_points*m_mag_directions.size();

    // Reduce every iteration total in one contiguous collective.  The previous
    // nested-tensor implementation issued O(N_t*N_components) Allreduces.
    const size_t square_count=m_pcn?2*corr_count:0;
    std::vector<RealType> totals(
        4+2*(corr_count+mag_count+m_num_real_points)
          +m_num_real_points+square_count);
    size_t total_offset{};
    totals[total_offset++]=std::real(m_total.partition);
    totals[total_offset++]=std::imag(m_total.partition);
    totals[total_offset++]=m_total.partition_abs;
    totals[total_offset++]=m_total.partition_abs_sq;
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t p=0;p<m_corr_directions.size();++p )
            for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
            {
                totals[total_offset++]=m_total.correlations.Re[t][p][tau];
                totals[total_offset++]=m_total.correlations.Im[t][p][tau];
            }
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            totals[total_offset++]=m_total.mag_Re[t][c];
            totals[total_offset++]=m_total.mag_Im[t][c];
        }
    for( const ComplexType value:m_total.closure )
    {
        totals[total_offset++]=std::real(value);
        totals[total_offset++]=std::imag(value);
    }
    for( const RealType value:m_total.closure_abs )
        totals[total_offset++]=value;
    if( m_pcn )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    totals[total_offset++]=m_sample_squares.Re[t][p][tau];
                    totals[total_offset++]=m_sample_squares.Im[t][p][tau];
                }
    if( total_offset!=totals.size() )
        throw std::logic_error("internal packed-total size mismatch");
    mpi_sum_buffer(totals);

    total_offset=0;
    m_total.partition={totals[total_offset],totals[total_offset+1]};
    total_offset+=2;
    m_total.partition_abs=totals[total_offset++];
    m_total.partition_abs_sq=totals[total_offset++];
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t p=0;p<m_corr_directions.size();++p )
            for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
            {
                m_total.correlations.Re[t][p][tau]=totals[total_offset++];
                m_total.correlations.Im[t][p][tau]=totals[total_offset++];
            }
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            m_total.mag_Re[t][c]=totals[total_offset++];
            m_total.mag_Im[t][c]=totals[total_offset++];
        }
    for( auto& value:m_total.closure )
    {
        value={totals[total_offset],totals[total_offset+1]};
        total_offset+=2;
    }
    for( auto& value:m_total.closure_abs )
        value=totals[total_offset++];
    if( m_pcn )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    m_sample_squares.Re[t][p][tau]=totals[total_offset++];
                    m_sample_squares.Im[t][p][tau]=totals[total_offset++];
                }
    if( total_offset!=totals.size() )
        throw std::logic_error("internal unpacked-total size mismatch");
    if( !m_pcn&&!denominator_resolved(m_total.partition,m_total.partition_abs) )
        throw std::runtime_error("the complex Z_M sum is statistically/numerically unresolved");
    if( m_closed_contour_observable_normalization
        &&!denominator_resolved(
            m_total.closure.back(),m_total.closure_abs.back()) )
        throw std::runtime_error(
            "the final closed-contour D(T) sum is statistically/numerically unresolved");

    partition_sum_Re.push_back(std::real(m_total.partition));
    partition_sum_Im.push_back(std::imag(m_total.partition));
    partition_abs_sums.push_back(m_total.partition_abs);
    average_phase_magnitudes.push_back(std::abs(m_total.partition)/m_total.partition_abs);
    effective_sample_sizes.push_back(m_total.partition_abs_sq>RealType{0.}
        ?std::norm(m_total.partition)/m_total.partition_abs_sq:RealType{0.});

    const RealType global_sample_count=static_cast<RealType>(
        m_num_samples_per_core*m_num_cores);
    auto estimate=[&]( const ComplexType numerator )
    {
        return m_pcn?numerator/global_sample_count
                    :exact_complex_ratio(
                        numerator,m_total.partition,m_total.partition_abs);
    };
    auto observable_estimate=[&]( const ComplexType numerator )
    {
        return m_closed_contour_observable_normalization
            ?exact_complex_ratio(
                numerator,m_total.closure.back(),m_total.closure_abs.back())
            :estimate(numerator);
    };

    correlations=CorrelationSet{m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t p=0;p<m_corr_directions.size();++p )
        {
            for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
            {
                const auto ratio=observable_estimate(
                    edge_value(m_total.correlations,t,p,tau));
                correlations.Re[t][p][tau]=std::real(ratio);
                correlations.Im[t][p][tau]=std::imag(ratio);
            }
        }

    magnetization_time_Re.assign(m_num_real_points,FieldVector{});
    magnetization_time_Im.assign(m_num_real_points,FieldVector{});
    magnetization_time_Re_stds.assign(m_num_real_points,FieldVector{});
    magnetization_time_Im_stds.assign(m_num_real_points,FieldVector{});
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            const ComplexType numerator{m_total.mag_Re[t][c],m_total.mag_Im[t][c]};
            const auto ratio=observable_estimate(numerator);
            const size_t direction=m_mag_directions[c];
            magnetization_time_Re[t][direction]=std::real(ratio);
            magnetization_time_Im[t][direction]=std::imag(ratio);
        }
    closed_contour_ratio_Re.resize(m_num_real_points);
    closed_contour_ratio_Im.resize(m_num_real_points);
    RealType closure_residual{};
    for( size_t t=0;t<m_num_real_points;++t )
    {
        const auto ratio=estimate(m_total.closure[t]);
        closed_contour_ratio_Re[t]=std::real(ratio);
        closed_contour_ratio_Im[t]=std::imag(ratio);
        closure_residual=std::max(closure_residual,std::abs(ratio-ComplexType{1.,0.}));
    }
    denominator_constancy_residuals.push_back(closure_residual);

    standard_errors=CorrelationSet{m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    closed_contour_residual_Re_sample_stds.assign(m_num_real_points,RealType{});
    closed_contour_residual_Im_sample_stds.assign(m_num_real_points,RealType{});
    closed_contour_residual_abs_sample_stds.assign(m_num_real_points,RealType{});

    const size_t pooled=m_num_blocks*m_num_cores;
    if( pooled<2 )
    {
        if( m_pcn )
        {
            blocking_curve_block_lengths.assign(
                1,static_cast<RealType>(m_samples_per_block));
            blocking_curve_mean_errors.assign(1,RealType{});
            blocking_curve_max_errors.assign(1,RealType{});
            maximum_tau_int.push_back(RealType{});
            block_length_to_tau_ratios.push_back(
                std::numeric_limits<RealType>::infinity());
            effective_sample_sizes.back()=global_sample_count;
        }
        contour_sample_stds=standard_errors;
        return;
    }

    if( m_pcn )
    {
        // Contiguous batch means are the pCN error estimator.  Starting from
        // the configured blocks, merge adjacent blocks by powers of two and
        // retain the largest standard error among scales with at least eight
        // pooled batches.  A still-rising stored blocking curve signals that
        // the production chain is too short for a resolved error estimate.
        constexpr size_t moment_stride=4;
        constexpr size_t minimum_pooled_batches=8;
        const size_t corr_base=0;
        const size_t mag_base=corr_count;
        const size_t closure_base=mag_base+mag_count;
        const size_t quantity_count=closure_base+m_num_real_points;
        blocking_curve_block_lengths.clear();
        blocking_curve_mean_errors.clear();
        blocking_curve_max_errors.clear();

        auto add_value=[]( std::vector<RealType>& moments,
                           const size_t quantity,const ComplexType value )
        {
            const size_t q=moment_stride*quantity;
            moments[q]+=std::real(value);
            moments[q+2]+=std::imag(value);
        };
        auto component_error=[]( const RealType centered_square_sum,
                                 const RealType count )
        {
            if( count<RealType{2.} ) return RealType{};
            return std::sqrt(std::max(RealType{},centered_square_sum)
                             /((count-RealType{1.})*count));
        };

        for( size_t merge=1;m_num_blocks/merge>=2;merge*=2 )
        {
            const size_t local_groups=m_num_blocks/merge;
            const size_t global_groups=local_groups*m_num_cores;
            std::vector<RealType> moments(moment_stride*quantity_count,RealType{});
            const RealType inverse_group_samples=RealType{1.}/static_cast<RealType>(
                merge*m_samples_per_block);
            auto group_observable_estimate=[&](
                const ComplexType numerator,const size_t first )
            {
                if( !m_closed_contour_observable_normalization )
                    return inverse_group_samples*numerator;
                ComplexType denominator{};
                RealType denominator_abs{};
                for( size_t b=0;b<merge;++b )
                {
                    denominator+=m_blocks[first+b].closure.back();
                    denominator_abs+=m_blocks[first+b].closure_abs.back();
                }
                return exact_complex_ratio(
                    numerator,denominator,denominator_abs);
            };
            for( size_t group=0;group<local_groups;++group )
            {
                const size_t first=group*merge;
                size_t flat{};
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t p=0;p<m_corr_directions.size();++p )
                        for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau,++flat )
                        {
                            ComplexType sum{};
                            for( size_t b=0;b<merge;++b )
                                sum+=edge_value(
                                    m_blocks[first+b].correlations,t,p,tau);
                            add_value(moments,corr_base+flat,
                                      group_observable_estimate(sum,first));
                        }
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t c=0;c<m_mag_directions.size();++c )
                    {
                        RealType re{},im{};
                        for( size_t b=0;b<merge;++b )
                        {
                            re+=m_blocks[first+b].mag_Re[t][c];
                            im+=m_blocks[first+b].mag_Im[t][c];
                        }
                        const size_t flat_mag=t*m_mag_directions.size()+c;
                        add_value(moments,mag_base+flat_mag,
                            group_observable_estimate(
                                ComplexType{re,im},first));
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    ComplexType sum{};
                    for( size_t b=0;b<merge;++b )
                        sum+=m_blocks[first+b].closure[t];
                    add_value(moments,closure_base+t,
                              inverse_group_samples*sum);
                }
            }
            mpi_sum_buffer(moments);
            const RealType count=static_cast<RealType>(global_groups);

            // Compute the variance in a second pass around the global batch
            // mean.  The former E[x^2]-E[x]^2 expression loses all precision
            // for nearly deterministic contour entries and can round a small
            // positive variance to zero.  That spurious zero subsequently
            // turns an otherwise finite standardized residual into infinity.
            std::vector<RealType> centered_moments(
                moment_stride*quantity_count,RealType{});
            auto add_centered_value=[&]( const size_t quantity,
                                         const ComplexType value )
            {
                const size_t q=moment_stride*quantity;
                const RealType delta_re=std::real(value)-moments[q]/count;
                const RealType delta_im=std::imag(value)-moments[q+2]/count;
                centered_moments[q+1]+=delta_re*delta_re;
                centered_moments[q+3]+=delta_im*delta_im;
            };
            for( size_t group=0;group<local_groups;++group )
            {
                const size_t first=group*merge;
                size_t flat{};
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t p=0;p<m_corr_directions.size();++p )
                        for( size_t tau=0;tau<m_num_imaginary_edge_points;
                             ++tau,++flat )
                        {
                            ComplexType sum{};
                            for( size_t b=0;b<merge;++b )
                                sum+=edge_value(
                                    m_blocks[first+b].correlations,t,p,tau);
                            add_centered_value(
                                corr_base+flat,
                                group_observable_estimate(sum,first));
                        }
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t c=0;c<m_mag_directions.size();++c )
                    {
                        RealType re{},im{};
                        for( size_t b=0;b<merge;++b )
                        {
                            re+=m_blocks[first+b].mag_Re[t][c];
                            im+=m_blocks[first+b].mag_Im[t][c];
                        }
                        add_centered_value(
                            mag_base+t*m_mag_directions.size()+c,
                            group_observable_estimate(
                                ComplexType{re,im},first));
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    ComplexType sum{};
                    for( size_t b=0;b<merge;++b )
                        sum+=m_blocks[first+b].closure[t];
                    add_centered_value(
                        closure_base+t,inverse_group_samples*sum);
                }
            }
            mpi_sum_buffer(centered_moments);
            auto errors=[&]( const size_t quantity )
            {
                const size_t q=moment_stride*quantity;
                return std::array<RealType,2>{
                    component_error(centered_moments[q+1],count),
                    component_error(centered_moments[q+3],count)};
            };
            const bool use_for_report=merge==1||global_groups>=minimum_pooled_batches;
            RealType curve_sum{},curve_max{}; size_t curve_count{};
            size_t flat{};
            for( size_t t=0;t<m_num_real_points;++t )
                for( size_t p=0;p<m_corr_directions.size();++p )
                    for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau,++flat )
                    {
                        const auto value=errors(corr_base+flat);
                        if( use_for_report )
                        {
                            standard_errors.Re[t][p][tau]=std::max(
                                standard_errors.Re[t][p][tau],value[0]);
                            standard_errors.Im[t][p][tau]=std::max(
                                standard_errors.Im[t][p][tau],value[1]);
                        }
                        for( const RealType component:value )
                            if( component>RealType{} )
                            {
                                curve_sum+=component;
                                curve_max=std::max(curve_max,component);
                                ++curve_count;
                            }
                    }
            if( use_for_report )
            {
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t c=0;c<m_mag_directions.size();++c )
                    {
                        const auto value=errors(
                            mag_base+t*m_mag_directions.size()+c);
                        const size_t direction=m_mag_directions[c];
                        magnetization_time_Re_stds[t][direction]=std::max(
                            magnetization_time_Re_stds[t][direction],value[0]);
                        magnetization_time_Im_stds[t][direction]=std::max(
                            magnetization_time_Im_stds[t][direction],value[1]);
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    const auto value=errors(closure_base+t);
                    closed_contour_residual_Re_sample_stds[t]=std::max(
                        closed_contour_residual_Re_sample_stds[t],value[0]);
                    closed_contour_residual_Im_sample_stds[t]=std::max(
                        closed_contour_residual_Im_sample_stds[t],value[1]);
                    closed_contour_residual_abs_sample_stds[t]=std::hypot(
                        closed_contour_residual_Re_sample_stds[t],
                        closed_contour_residual_Im_sample_stds[t]);
                }
            }
            blocking_curve_block_lengths.push_back(static_cast<RealType>(
                merge*m_samples_per_block));
            blocking_curve_mean_errors.push_back(
                curve_count?curve_sum/static_cast<RealType>(curve_count):RealType{});
            blocking_curve_max_errors.push_back(curve_max);
        }

        contour_tau_int=CorrelationSet{
            m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
        RealType largest_tau{};
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    const RealType mean_re=correlations.Re[t][p][tau];
                    const RealType mean_im=correlations.Im[t][p][tau];
                    const RealType var_re=std::max(RealType{},
                        m_sample_squares.Re[t][p][tau]/global_sample_count
                        -mean_re*mean_re);
                    const RealType var_im=std::max(RealType{},
                        m_sample_squares.Im[t][p][tau]/global_sample_count
                        -mean_im*mean_im);
                    const RealType tau_re=var_re>RealType{}
                        ?RealType{0.5}*standard_errors.Re[t][p][tau]
                         *standard_errors.Re[t][p][tau]*global_sample_count/var_re
                        :RealType{};
                    const RealType tau_im=var_im>RealType{}
                        ?RealType{0.5}*standard_errors.Im[t][p][tau]
                         *standard_errors.Im[t][p][tau]*global_sample_count/var_im
                        :RealType{};
                    contour_tau_int.Re[t][p][tau]=tau_re;
                    contour_tau_int.Im[t][p][tau]=tau_im;
                    largest_tau=std::max({largest_tau,tau_re,tau_im});
                }
        maximum_tau_int.push_back(largest_tau);
        effective_sample_sizes.back()=largest_tau>RealType{}
            ?global_sample_count/(RealType{2.}*largest_tau)
            :global_sample_count;
        block_length_to_tau_ratios.push_back(largest_tau>RealType{}
            ?static_cast<RealType>(m_samples_per_block)/largest_tau
            :std::numeric_limits<RealType>::infinity());
        contour_sample_stds=standard_errors;
        return;
    }

    // Each rank evaluates its own leave-one-block replicates and contributes
    // only six moments per reported quantity.  This is algebraically identical
    // to pooling all rank-local blocks, without replicating the block tensors.
    constexpr size_t moment_stride=6;
    const size_t corr_base=0;
    const size_t mag_base=corr_base+corr_count;
    const size_t closure_base=mag_base+mag_count;
    const size_t quantity_count=closure_base+m_num_real_points;
    std::vector<RealType> moments(moment_stride*quantity_count,RealType{});
    auto add_moment=[&]( const size_t quantity, const ComplexType value,
                         const ComplexType component_center )
    {
        const size_t q=moment_stride*quantity;
        const RealType re=std::real(value-component_center);
        const RealType im=std::imag(value-component_center);
        // Center magnitudes as well; their variance is shift-invariant and the
        // centered moments avoid cancellation for nearly deterministic data.
        const RealType magnitude=std::abs(value)-std::abs(component_center);
        moments[q]+=re; moments[q+1]+=re*re;
        moments[q+2]+=im; moments[q+3]+=im*im;
        moments[q+4]+=magnitude; moments[q+5]+=magnitude*magnitude;
    };
    auto block_ratio=[&]( const ComplexType total_numerator,
                          const ComplexType block_numerator,
                          const ComplexBlockSums& block )
    {
        return exact_complex_ratio(
            total_numerator-block_numerator,
            m_total.partition-block.partition,
            m_total.partition_abs-block.partition_abs);
    };
    auto observable_block_ratio=[&](
        const ComplexType total_numerator,
        const ComplexType block_numerator,const ComplexBlockSums& block )
    {
        return m_closed_contour_observable_normalization
            ?exact_complex_ratio(
                total_numerator-block_numerator,
                m_total.closure.back()-block.closure.back(),
                m_total.closure_abs.back()-block.closure_abs.back())
            :block_ratio(total_numerator,block_numerator,block);
    };

    for( const auto& block:m_blocks )
    {
        size_t flat{};
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau,++flat )
                {
                    const ComplexType total=edge_value(
                        m_total.correlations,t,p,tau);
                    const ComplexType replicate=observable_block_ratio(
                        total,edge_value(block.correlations,t,p,tau),block);
                    add_moment(corr_base+flat,replicate,
                        observable_estimate(total));
                }

        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            for( size_t t=0;t<m_num_real_points;++t )
            {
                const size_t flat_mag=t*m_mag_directions.size()+c;
                const ComplexType total{
                    m_total.mag_Re[t][c],m_total.mag_Im[t][c]};
                const ComplexType block_value{
                    block.mag_Re[t][c],block.mag_Im[t][c]};
                const ComplexType current=observable_block_ratio(
                    total,block_value,block);
                const ComplexType center=observable_estimate(total);
                add_moment(mag_base+flat_mag,current,center);
            }
        }

        for( size_t t=0;t<m_num_real_points;++t )
        {
            const ComplexType ratio=block_ratio(
                m_total.closure[t],block.closure[t],block);
            const ComplexType value=ratio-ComplexType{1.,0.};
            const ComplexType center=exact_complex_ratio(
                m_total.closure[t],m_total.partition,m_total.partition_abs)
                -ComplexType{1.,0.};
            add_moment(closure_base+t,value,center);
        }
    }
    mpi_sum_buffer(moments);

    const RealType block_count=static_cast<RealType>(pooled);
    auto moment_error=[&]( const RealType sum, const RealType square_sum )
    {
        const RealType centered=std::max(
            RealType{},square_sum-sum*sum/block_count);
        return std::sqrt((block_count-RealType{1.})/block_count*centered);
    };
    auto errors=[&]( const size_t quantity )
    {
        const size_t q=moment_stride*quantity;
        return std::array<RealType,3>{
            moment_error(moments[q],moments[q+1]),
            moment_error(moments[q+2],moments[q+3]),
            moment_error(moments[q+4],moments[q+5])};
    };

    size_t flat{};
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t p=0;p<m_corr_directions.size();++p )
            for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau,++flat )
            {
                const auto value=errors(corr_base+flat);
                standard_errors.Re[t][p][tau]=value[0];
                standard_errors.Im[t][p][tau]=value[1];
            }
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            const size_t direction=m_mag_directions[c];
            const size_t flat_mag=t*m_mag_directions.size()+c;
            const auto mag_errors=errors(mag_base+flat_mag);
            magnetization_time_Re_stds[t][direction]=mag_errors[0];
            magnetization_time_Im_stds[t][direction]=mag_errors[1];
        }
    for( size_t t=0;t<m_num_real_points;++t )
    {
        const auto value=errors(closure_base+t);
        closed_contour_residual_Re_sample_stds[t]=value[0];
        closed_contour_residual_Im_sample_stds[t]=value[1];
        closed_contour_residual_abs_sample_stds[t]=value[2];
    }
    contour_sample_stds=standard_errors;
}

void RunTimeData::record_iteration_error(
    RealType absolute_error,RealType standardized_error )
{
    absolute_iteration_errors.push_back(absolute_error);
    standardized_iteration_errors.push_back(standardized_error);
}

bool RunTimeData::diagnostics_pass() const
{
    if( covariance_symmetry_errors.empty()||branch_identity_errors.empty()
        ||gaussian_factor_reconstruction_errors.empty()||average_phase_magnitudes.empty()
        ||denominator_constancy_residuals.empty() ) return false;
    if( covariance_symmetry_errors.back()>m_covariance_tolerance
        ||branch_identity_errors.back()>m_branch_identity_tolerance
        ||gaussian_factor_reconstruction_errors.back()>m_takagi_tolerance
        ||average_phase_magnitudes.back()<m_minimum_phase_magnitude
        ||denominator_constancy_residuals.back()>m_denominator_constancy_tolerance ) return false;
    // if( m_pcn&&(maximum_relative_imaginary_sampling_weights.empty()
    //     ||maximum_relative_imaginary_sampling_weights.back()
    //       >m_partition_imaginary_tolerance) ) return false;
    for( const size_t direction:m_mag_directions )
        if( std::abs(magnetization_time_Im.front()[direction])>
            m_imaginary_magnetization_sigma*
            magnetization_time_Im_stds.front()[direction] ) return false;
    return true;
}

void RunTimeData::finalize_iteration_step()
{
    if( m_my_rank==0 )
    {
        std::cout<<"Iteration step finished:\n"
          <<print::quantity_to_output_line(36,m_harmonic_bath
             ?"initial-guess difference (unused)":"raw complex fixed-point residual",
             print::round_value_to_string(absolute_iteration_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"largest standardized residual",
             print::round_value_to_string(standardized_iteration_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"standardized residual threshold",
             print::round_value_to_string(m_iteration_error_sigma_threshold,m_num_print_digits))
          <<print::quantity_to_output_line(36,m_pcn
             ?"sampled-Z phase magnitude":"average phase magnitude",
             print::round_value_to_string(average_phase_magnitudes.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,m_pcn
             ?"minimum autocorrelation N_eff":"complex-weight N_eff",
             print::round_value_to_string(effective_sample_sizes.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"branch transpose residual",
             print::round_value_to_string(branch_identity_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"Gaussian covariance truncation",
             print::round_value_to_string(
                 gaussian_covariance_approximation_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"largest Gaussian factor block",
             std::to_string(gaussian_largest_factorization_dimensions.back()))
          <<print::quantity_to_output_line(36,"closed-contour trace residual",
             print::round_value_to_string(denominator_constancy_residuals.back(),m_num_print_digits));
        if( m_pcn )
        {
            std::cout
              <<print::quantity_to_output_line(36,"MH acceptance rate",
                 print::round_value_to_string(mh_acceptance_rates.back(),m_num_print_digits))
              <<print::quantity_to_output_line(36,"nonpositive-weight rejection rate",
                 print::round_value_to_string(
                     mh_nonpositive_rejection_rates.back(),m_num_print_digits));
        }
    }
    ++num_Iterations;
}

bool RunTimeData::terminate()
{
    if( !m_self_consistency )
    {
        termination=m_harmonic_bath
            ?"prescribed harmonic bath":"no self consistency";
        return true;
    }
    if( standardized_iteration_errors.back()<m_iteration_error_sigma_threshold
        &&diagnostics_pass() )
    {
        print::print_R0(m_my_rank,"\033[1;32mTerminating: fixed point and contour diagnostics converged.\033[0m\n");
        termination="by convergence"; return true;
    }
    if( num_Iterations>=m_iteration_limit )
    {
        print::print_R0(m_my_rank,"\033[1;31mTerminating because the iteration limit was reached.\033[0m\n");
        termination="by iteration limit"; return true;
    }
    return false;
}

}
