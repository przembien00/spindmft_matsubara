#include"Run_Time_Data.h"
#include"../Functions/Functions.h"

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
      magnetization(symmetry,num_real_points),
      closure(num_real_points,ComplexType{})
{}

RunTimeData::RunTimeData(const ps::ParameterSpace& pspace,int my_rank)
    : m_real_time_points(pspace.num_RealTimePoints),
      m_real_samples_per_core(pspace.num_SamplesPerCore),
      m_real_iteration_limit(pspace.Iteration_Limit),
      m_preparation_samples_per_core(pspace.num_PreparationSamplesPerCore),
      m_preparation_iteration_limit(pspace.Preparation_Iteration_Limit),
      m_requested_num_blocks(pspace.num_blocks),
      m_my_rank(my_rank),m_num_print_digits(pspace.num_PrintDigits),
      m_symmetry(pspace.correlation_symmetry_type),
      m_num_imaginary_edge_points(pspace.num_TimePoints),
      m_num_cores(pspace.get_num_Cores()),m_seed(pspace.seed),
      m_self_consistency(pspace.self_consistency),
      m_harmonic_bath(pspace.uses_harmonic_bath()),
      m_pcn(pspace.sampling_strategy=="pcn"),
      m_matsubara_iteration_error_threshold(pspace.num_PreparationSamples>0
          ?pspace.iteration_error_sigma_threshold
            /std::sqrt(RealType{3.}*static_cast<RealType>(pspace.num_PreparationSamples))
            *pspace.spin_float*(pspace.spin_float+RealType{1.})/RealType{3.}
          :RealType{}),
      m_iteration_error_sigma_threshold(pspace.iteration_error_sigma_threshold),
      m_real_time_residual_scale(std::max(pspace.spin_float,pspace.spin_float*pspace.spin_float)),
      m_covariance_tolerance(pspace.covariance_tolerance),
      m_takagi_tolerance(pspace.takagi_tolerance),
      m_minimum_phase_magnitude(pspace.minimum_phase_magnitude),
      m_corr_directions(CorrTen{pspace.correlation_symmetry_type,
                                pspace.num_TimePoints}.get_direction_pairs()),
      m_mag_directions(mag::determine_magnetization_directions(
          pspace.correlation_symmetry_type))
{
}

void RunTimeData::configure_iteration(bool matsubara_preparation)
{
    m_matsubara_preparation=matsubara_preparation;
    m_num_real_points=matsubara_preparation?size_t{1}:m_real_time_points;
    m_num_samples_per_core=matsubara_preparation
        ?m_preparation_samples_per_core:m_real_samples_per_core;
    m_iteration_limit=matsubara_preparation?m_preparation_iteration_limit:m_real_iteration_limit;
    if(m_num_samples_per_core==0)
        throw std::invalid_argument("samples per core must be positive for each loop");
    m_num_blocks=std::min(std::max<size_t>(1,m_requested_num_blocks),m_num_samples_per_core);
    while(m_num_blocks>1&&m_num_samples_per_core%m_num_blocks!=0)--m_num_blocks;
    m_samples_per_block=m_num_samples_per_core/m_num_blocks;
    auto& statistics=loop();
    statistics.num_blocks=m_num_blocks;
    statistics.samples_per_block=m_samples_per_block;
    if(statistics.contour_sample_stds.Re.empty())
    {
        statistics.contour_sample_stds=CorrelationSet{
            m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
        statistics.contour_tau_int=CorrelationSet{
            m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    }
    m_samples_seen=0; m_current_block=0;
    m_total=TotalSampleSums{};
    m_total.closure.assign(m_num_real_points,ComplexType{});
    m_correlations=nullptr;m_magnetization=nullptr;
    if(m_pcn)
        m_sample_squares=CorrelationSet{
            m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    m_current_sample=CorrelationSet{};
    m_partition_second_moments={};
    if(m_pcn)m_ratio_cross_moments.assign(
        m_num_real_points*m_corr_directions.size()*m_num_imaginary_edge_points,{});
    m_blocks.assign(m_num_blocks,ComplexBlockSums{
        m_symmetry,m_num_imaginary_edge_points,m_num_real_points});
}

void RunTimeData::begin_iteration_accumulation(CorrelationSet& correlations,
    Magnetization& magnetization,bool matsubara_preparation)
{
    configure_iteration(matsubara_preparation);
    correlations=CorrelationSet{m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    magnetization=Magnetization{m_symmetry,m_num_real_points};
    m_correlations=&correlations;
    m_magnetization=&magnetization;
}

void RunTimeData::begin_independent_sample(ComplexType Z)
{
    if(m_pcn)throw std::invalid_argument("fused independent measurement cannot accumulate pCN states");
    begin_sample(Z);
}

void RunTimeData::begin_sample( ComplexType Z, RealType observable_normalization )
{
    if(!m_correlations||!m_magnetization)
        throw std::logic_error("begin accumulation into caller-owned observables before sampling");
    if( m_samples_seen>=m_num_samples_per_core )
        throw std::runtime_error("too many independent samples accumulated in one iteration");
    if( !finite(Z)||!std::isfinite(observable_normalization) )
        throw std::runtime_error("a contour trajectory produced a non-finite sample weight");
    m_current_observable_normalization=observable_normalization;
    m_current_block=std::min(m_samples_seen/m_samples_per_block,m_num_blocks-1);
    // In pCN mode observable_normalization is the inverse positive likelihood
    // (|Z_M| for the staged Matsubara chain).
    // Its unknown normalization must cancel between the reweighted numerator
    // and denominator. Accumulating raw Z here loses that complex denominator
    // at finite sample count and breaks identities such as <(S^a)^2>=1/4.
    // Independent sampling passes a normalization of one and is unchanged.
    Z*=observable_normalization;
    if(!finite(Z))throw std::runtime_error("reweighted partition function is nonfinite");
    m_current_partition=Z;
    if(m_pcn)
    {
        const RealType re=std::real(Z),im=std::imag(Z);
        m_partition_second_moments[0]+=re*re;
        m_partition_second_moments[1]+=im*im;
        m_partition_second_moments[2]+=re*im;
    }
    m_total.partition+=Z; m_total.partition_abs+=std::abs(Z);
    m_total.partition_abs_sq+=std::norm(Z);
    auto& block=m_blocks[m_current_block];
    block.partition+=Z; block.partition_abs+=std::abs(Z);
}

void RunTimeData::accumulate_closed_contour_trace( size_t t, ComplexType value )
{
    if( t>=m_num_real_points||!finite(value) )
        throw std::runtime_error("invalid closed-contour trace");
    value*=m_current_observable_normalization;
    m_total.closure[t]+=value; m_blocks[m_current_block].closure[t]+=value;
}

void RunTimeData::accumulate_magnetization( size_t t, size_t c,
                                            ComplexType numerator )
{
    if( t>=m_num_real_points||c>=m_mag_directions.size()||!finite(numerator) )
        throw std::runtime_error("invalid complex magnetization numerator");
    numerator*=m_current_observable_normalization;
    m_magnetization->Re[t][c]+=std::real(numerator);
    m_magnetization->Im[t][c]+=std::imag(numerator);
    m_blocks[m_current_block].magnetization.Re[t][c]+=std::real(numerator);
    m_blocks[m_current_block].magnetization.Im[t][c]+=std::imag(numerator);
}

void RunTimeData::accumulate_edge_correlation(
    size_t t,size_t p,size_t tau,ComplexType numerator )
{
    if( !finite(numerator) )
        throw std::runtime_error("non-finite edge contour-correlation numerator");
    numerator*=m_current_observable_normalization;
    const RealType re=std::real(numerator),im=std::imag(numerator);
    m_correlations->Re[t][p][tau]+=re;
    m_correlations->Im[t][p][tau]+=im;
    m_blocks[m_current_block].correlations.Re[t][p][tau]+=re;
    m_blocks[m_current_block].correlations.Im[t][p][tau]+=im;
    if( m_pcn )
    {
        if(m_current_sample.Re.empty())
            m_current_sample=CorrelationSet{
                m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
        m_current_sample.Re[t][p][tau]+=re;
        m_current_sample.Im[t][p][tau]+=im;
    }
}

void RunTimeData::end_sample()
{
    if( m_pcn && !m_current_sample.Re.empty() )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    const RealType re=m_current_sample.Re[t][p][tau];
                    const RealType im=m_current_sample.Im[t][p][tau];
                    m_sample_squares.Re[t][p][tau]+=re*re;
                    m_sample_squares.Im[t][p][tau]+=im*im;
                    accumulate_ratio_cross_moments(
                        (t*m_corr_directions.size()+p)*m_num_imaginary_edge_points+tau,re,im);
                    m_current_sample.Re[t][p][tau]=RealType{};
                    m_current_sample.Im[t][p][tau]=RealType{};
                }
    ++m_samples_seen;
}

void RunTimeData::accumulate_ratio_cross_moments(size_t flat,RealType re,RealType im)
{
    auto& moments=m_ratio_cross_moments[flat];
    const RealType zr=std::real(m_current_partition),zi=std::imag(m_current_partition);
    moments[0]+=re*im;moments[1]+=re*zr;moments[2]+=re*zi;
    moments[3]+=im*zr;moments[4]+=im*zi;
}

void RunTimeData::accumulate_sample( const spinDMFT::Functions::MeasuredSample& sample,
                                      const RealType normalization )
{
    const size_t nt=m_num_real_points,nc=m_corr_directions.size();
    const size_t nm=m_mag_directions.size(),ni=m_num_imaginary_edge_points;
    if( sample.correlations.size()!=nt*nc*ni
        ||sample.magnetization.size()!=nt*nm||sample.closure.size()!=nt )
        throw std::invalid_argument("cached observable dimensions differ from accumulator");
    begin_sample(sample.partition,normalization);
    auto& block=m_blocks[m_current_block];
    size_t index{};
    for( size_t t=0;t<nt;++t )
    {
        accumulate_closed_contour_trace(t,sample.closure[t]);
        for( size_t c=0;c<nm;++c )
            accumulate_magnetization(t,c,sample.magnetization[t*nm+c]);
        for( size_t c=0;c<nc;++c )
        {
            RealType* total_re=m_correlations->Re[t][c].data();
            RealType* total_im=m_correlations->Im[t][c].data();
            RealType* block_re=block.correlations.Re[t][c].data();
            RealType* block_im=block.correlations.Im[t][c].data();
            RealType* square_re=m_pcn?m_sample_squares.Re[t][c].data():nullptr;
            RealType* square_im=m_pcn?m_sample_squares.Im[t][c].data():nullptr;
            for( size_t tau=0;tau<ni;++tau,++index )
            {
                const ComplexType raw=sample.correlations[index];
                if( !finite(raw) )
                    throw std::runtime_error("non-finite edge contour-correlation numerator");
                const ComplexType value=raw*normalization;
                if(!finite(value))throw std::runtime_error("reweighted correlation numerator is nonfinite");
                const RealType re=std::real(value),im=std::imag(value);
                total_re[tau]+=re; total_im[tau]+=im;
                block_re[tau]+=re; block_im[tau]+=im;
                if( m_pcn )
                {
                    square_re[tau]+=re*re;
                    square_im[tau]+=im*im;
                    accumulate_ratio_cross_moments(index,re,im);
                }
            }
        }
    }
    // The complete state was accumulated above, including its sample square.
    // Rejected pCN steps take this same path and keep their block position.
    ++m_samples_seen;
}

ComplexType RunTimeData::edge_value(
    const CorrelationSet& values,size_t t,size_t p,size_t tau )
{ return {values.Re[t][p][tau],values.Im[t][p][tau]}; }

void RunTimeData::record_complex_field_diagnostics(
    RealType symmetry_error,RealType factor_error,
    RealType approximation_error,size_t latent_dimension,
    size_t largest_factorization_dimension )
{
    loop().covariance_symmetry_errors.push_back(symmetry_error);
    loop().gaussian_factor_reconstruction_errors.push_back(factor_error);
    loop().gaussian_covariance_approximation_errors.push_back(approximation_error);
    loop().gaussian_factor_latent_dimensions.push_back(
        static_cast<unsigned int>(latent_dimension));
    loop().gaussian_largest_factorization_dimensions.push_back(
        static_cast<unsigned int>(largest_factorization_dimension));
}

void RunTimeData::record_pcn_diagnostics(
    const size_t accepted,const size_t proposed,
    const size_t rejected_zero_weight )
{
    std::array<RealType,3> local{
        static_cast<RealType>(accepted),static_cast<RealType>(proposed),
        static_cast<RealType>(rejected_zero_weight)};
    std::array<RealType,3> global{};
    MPI_Allreduce(local.data(),global.data(),static_cast<int>(global.size()),
                  mpi_real_type(),MPI_SUM,MPI_COMM_WORLD);
    loop().mh_acceptance_rates.push_back(global[1]>RealType{}?global[0]/global[1]:RealType{});
    loop().mh_nonpositive_rejection_rates.push_back(
        global[1]>RealType{}?global[2]/global[1]:RealType{});
}

void RunTimeData::mpi_reduce_and_normalize()
{
    if(!m_correlations||!m_magnetization)
        throw std::logic_error("no observable sums are bound for MPI reduction");
    if( m_samples_seen!=m_num_samples_per_core )
        throw std::runtime_error("the production sweep did not accumulate the requested samples");
    const size_t corr_count=m_num_real_points*m_corr_directions.size()
                           *m_num_imaginary_edge_points;
    const size_t mag_count=m_num_real_points*m_mag_directions.size();

    // Reduce every iteration total in one contiguous collective.  The previous
    // nested-tensor implementation issued O(N_t*N_components) Allreduces.
    const size_t square_count=m_pcn?2*corr_count:0;
    const size_t cross_count=m_pcn?3+5*corr_count:0;
    std::vector<RealType> totals(
        4+2*(corr_count+mag_count+m_num_real_points)
          +square_count+cross_count);
    size_t total_offset{};
    totals[total_offset++]=std::real(m_total.partition);
    totals[total_offset++]=std::imag(m_total.partition);
    totals[total_offset++]=m_total.partition_abs;
    totals[total_offset++]=m_total.partition_abs_sq;
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t p=0;p<m_corr_directions.size();++p )
            for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
            {
                totals[total_offset++]=m_correlations->Re[t][p][tau];
                totals[total_offset++]=m_correlations->Im[t][p][tau];
            }
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            totals[total_offset++]=m_magnetization->Re[t][c];
            totals[total_offset++]=m_magnetization->Im[t][c];
        }
    for( const ComplexType value:m_total.closure )
    {
        totals[total_offset++]=std::real(value);
        totals[total_offset++]=std::imag(value);
    }
    if( m_pcn )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    totals[total_offset++]=m_sample_squares.Re[t][p][tau];
                    totals[total_offset++]=m_sample_squares.Im[t][p][tau];
                }
    if(m_pcn)
    {
        for(const auto value:m_partition_second_moments)totals[total_offset++]=value;
        for(const auto& moments:m_ratio_cross_moments)
            for(const auto value:moments)totals[total_offset++]=value;
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
                m_correlations->Re[t][p][tau]=totals[total_offset++];
                m_correlations->Im[t][p][tau]=totals[total_offset++];
            }
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            m_magnetization->Re[t][c]=totals[total_offset++];
            m_magnetization->Im[t][c]=totals[total_offset++];
        }
    for( auto& value:m_total.closure )
    {
        value={totals[total_offset],totals[total_offset+1]};
        total_offset+=2;
    }
    if( m_pcn )
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    m_sample_squares.Re[t][p][tau]=totals[total_offset++];
                    m_sample_squares.Im[t][p][tau]=totals[total_offset++];
                }
    if(m_pcn)
    {
        for(auto& value:m_partition_second_moments)value=totals[total_offset++];
        for(auto& moments:m_ratio_cross_moments)
            for(auto& value:moments)value=totals[total_offset++];
    }
    if( total_offset!=totals.size() )
        throw std::logic_error("internal unpacked-total size mismatch");
    if( !denominator_resolved(m_total.partition,m_total.partition_abs) )
        throw std::runtime_error("the complex Z_M sum is statistically/numerically unresolved");
    loop().partition_sum_Re.push_back(std::real(m_total.partition));
    loop().partition_sum_Im.push_back(std::imag(m_total.partition));
    loop().partition_abs_sums.push_back(m_total.partition_abs);
    loop().average_phase_magnitudes.push_back(std::abs(m_total.partition)/m_total.partition_abs);
    loop().effective_sample_sizes.push_back(m_total.partition_abs_sq>RealType{0.}
        ?std::norm(m_total.partition)/m_total.partition_abs_sq:RealType{0.});

    const RealType global_sample_count=static_cast<RealType>(
        m_num_samples_per_core*m_num_cores);
    auto estimate=[&]( const ComplexType numerator )
    {
        return exact_complex_ratio(
            numerator,m_total.partition,m_total.partition_abs);
    };
    loop().magnetization_time_stds=Magnetization{m_symmetry,m_num_real_points};
    loop().closed_contour_ratio_Re.resize(m_num_real_points);
    loop().closed_contour_ratio_Im.resize(m_num_real_points);
    RealType closure_residual{};
    for( size_t t=0;t<m_num_real_points;++t )
    {
        const auto ratio=estimate(m_total.closure[t]);
        loop().closed_contour_ratio_Re[t]=std::real(ratio);
        loop().closed_contour_ratio_Im[t]=std::imag(ratio);
        closure_residual=std::max(closure_residual,std::abs(ratio-ComplexType{1.,0.}));
    }
    loop().denominator_constancy_residuals.push_back(closure_residual);

    loop().contour_sample_stds=CorrelationSet{
        m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
    auto& standard_errors=loop().contour_sample_stds;
    loop().closed_contour_residual_Re_sample_stds.assign(m_num_real_points,RealType{});
    loop().closed_contour_residual_Im_sample_stds.assign(m_num_real_points,RealType{});
    loop().closed_contour_residual_abs_sample_stds.assign(m_num_real_points,RealType{});

    const size_t pooled=m_num_blocks*m_num_cores;
    if( pooled<2 )
    {
        if( m_pcn )
        {
            loop().blocking_curve_block_lengths.assign(
                1,static_cast<RealType>(m_samples_per_block));
            loop().blocking_curve_mean_errors.assign(1,RealType{});
            loop().blocking_curve_max_errors.assign(1,RealType{});
            loop().maximum_tau_int.push_back(RealType{});
            loop().block_length_to_tau_ratios.push_back(
                std::numeric_limits<RealType>::infinity());
            loop().effective_sample_sizes.back()=global_sample_count;
        }
        normalize_observables();
        return;
    }

    if( m_pcn )
    {
        // Contiguous batch means are the pCN error estimator.  Starting from
        // the configured blocks, merge adjacent blocks by powers of two and
        // retain the largest standard error among scales with at least eight
        // pooled batches.  A still-rising stored blocking curve signals that
        // the production chain is too short for a resolved error estimate.
        constexpr size_t moment_stride=2;
        constexpr size_t minimum_pooled_batches=8;
        const size_t corr_base=0;
        const size_t mag_base=corr_count;
        const size_t closure_base=mag_base+mag_count;
        const size_t quantity_count=closure_base+m_num_real_points;
        loop().blocking_curve_block_lengths.clear();
        loop().blocking_curve_mean_errors.clear();
        loop().blocking_curve_max_errors.clear();

        auto add_value=[]( std::vector<RealType>& moments,
                           const size_t quantity,const ComplexType value )
        {
            const size_t q=moment_stride*quantity;
            moments[q]+=std::real(value);
            moments[q+1]+=std::imag(value);
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
            std::vector<ComplexType> group_partitions(local_groups);
            std::vector<RealType> group_absolute_weights(local_groups,RealType{});
            for(size_t group=0;group<local_groups;++group)
                for(size_t b=0;b<merge;++b)
                {
                    const auto& block=m_blocks[group*merge+b];
                    group_partitions[group]+=block.partition;
                    group_absolute_weights[group]+=block.partition_abs;
                }
            const auto group_partition_estimate=[&](ComplexType numerator,size_t group)
            {
                return exact_complex_ratio(numerator,group_partitions[group],
                    group_absolute_weights[group]);
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
                                      group_partition_estimate(sum,group));
                        }
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t c=0;c<m_mag_directions.size();++c )
                    {
                        RealType re{},im{};
                        for( size_t b=0;b<merge;++b )
                        {
                            re+=m_blocks[first+b].magnetization.Re[t][c];
                            im+=m_blocks[first+b].magnetization.Im[t][c];
                        }
                        const size_t flat_mag=t*m_mag_directions.size()+c;
                        add_value(moments,mag_base+flat_mag,
                            group_partition_estimate(
                                ComplexType{re,im},group));
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    ComplexType sum{};
                    for( size_t b=0;b<merge;++b )
                        sum+=m_blocks[first+b].closure[t];
                    add_value(moments,closure_base+t,
                              group_partition_estimate(sum,group));
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
                const RealType delta_im=std::imag(value)-moments[q+1]/count;
                centered_moments[q]+=delta_re*delta_re;
                centered_moments[q+1]+=delta_im*delta_im;
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
                                group_partition_estimate(sum,group));
                        }
                for( size_t t=0;t<m_num_real_points;++t )
                    for( size_t c=0;c<m_mag_directions.size();++c )
                    {
                        RealType re{},im{};
                        for( size_t b=0;b<merge;++b )
                        {
                            re+=m_blocks[first+b].magnetization.Re[t][c];
                            im+=m_blocks[first+b].magnetization.Im[t][c];
                        }
                        add_centered_value(
                            mag_base+t*m_mag_directions.size()+c,
                            group_partition_estimate(
                                ComplexType{re,im},group));
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    ComplexType sum{};
                    for( size_t b=0;b<merge;++b )
                        sum+=m_blocks[first+b].closure[t];
                    add_centered_value(
                        closure_base+t,group_partition_estimate(sum,group));
                }
            }
            mpi_sum_buffer(centered_moments);
            auto errors=[&]( const size_t quantity )
            {
                const size_t q=moment_stride*quantity;
                return std::array<RealType,2>{
                    component_error(centered_moments[q],count),
                    component_error(centered_moments[q+1],count)};
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
                        loop().magnetization_time_stds.Re[t][c]=std::max(
                            loop().magnetization_time_stds.Re[t][c],value[0]);
                        loop().magnetization_time_stds.Im[t][c]=std::max(
                            loop().magnetization_time_stds.Im[t][c],value[1]);
                    }
                for( size_t t=0;t<m_num_real_points;++t )
                {
                    const auto value=errors(closure_base+t);
                    loop().closed_contour_residual_Re_sample_stds[t]=std::max(
                        loop().closed_contour_residual_Re_sample_stds[t],value[0]);
                    loop().closed_contour_residual_Im_sample_stds[t]=std::max(
                        loop().closed_contour_residual_Im_sample_stds[t],value[1]);
                    loop().closed_contour_residual_abs_sample_stds[t]=std::hypot(
                        loop().closed_contour_residual_Re_sample_stds[t],
                        loop().closed_contour_residual_Im_sample_stds[t]);
                }
            }
            loop().blocking_curve_block_lengths.push_back(static_cast<RealType>(
                merge*m_samples_per_block));
            loop().blocking_curve_mean_errors.push_back(
                curve_count?curve_sum/static_cast<RealType>(curve_count):RealType{});
            loop().blocking_curve_max_errors.push_back(curve_max);
        }

        loop().contour_tau_int=CorrelationSet{
            m_symmetry,m_num_imaginary_edge_points,m_num_real_points};
        RealType largest_tau{};
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau )
                {
                    // The ratio influence includes fluctuations of the
                    // complex phase denominator. Numerator variance alone
                    // is not invariant to likelihood normalization and
                    // gives incorrect tau/ESS for phase reweighting.
                    const size_t flat=(t*m_corr_directions.size()+p)*m_num_imaginary_edge_points+tau;
                    const auto& cross=m_ratio_cross_moments[flat];
                    const ComplexType inv=global_sample_count/m_total.partition;
                    const ComplexType g=estimate(edge_value(*m_correlations,t,p,tau));
                    const ComplexType gi=g*inv;
                    const std::array<std::array<RealType,4>,4> second{{
                        {{m_sample_squares.Re[t][p][tau],cross[0],cross[1],cross[2]}},
                        {{cross[0],m_sample_squares.Im[t][p][tau],cross[3],cross[4]}},
                        {{cross[1],cross[3],m_partition_second_moments[0],m_partition_second_moments[2]}},
                        {{cross[2],cross[4],m_partition_second_moments[2],m_partition_second_moments[1]}}
                    }};
                    const auto influence_variance=[&](const std::array<RealType,4>& coefficients)
                    {
                        long double variance{},scale{};
                        for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)
                        {
                            const long double term=static_cast<long double>(coefficients[i])
                                *coefficients[j]*second[i][j]/global_sample_count;
                            variance+=term;scale+=std::abs(term);
                        }
                        if(variance<=100*std::numeric_limits<RealType>::epsilon()*scale)return RealType{};
                        return static_cast<RealType>(variance);
                    };
                    const RealType var_re=influence_variance({std::real(inv),-std::imag(inv),-std::real(gi),std::imag(gi)});
                    const RealType var_im=influence_variance({std::imag(inv),std::real(inv),-std::imag(gi),-std::real(gi)});
                    const RealType tau_re=var_re>RealType{}
                        ?RealType{0.5}*standard_errors.Re[t][p][tau]
                         *standard_errors.Re[t][p][tau]*global_sample_count/var_re
                        :RealType{};
                    const RealType tau_im=var_im>RealType{}
                        ?RealType{0.5}*standard_errors.Im[t][p][tau]
                         *standard_errors.Im[t][p][tau]*global_sample_count/var_im
                        :RealType{};
                    loop().contour_tau_int.Re[t][p][tau]=tau_re;
                    loop().contour_tau_int.Im[t][p][tau]=tau_im;
                    largest_tau=std::max({largest_tau,tau_re,tau_im});
                }
        loop().maximum_tau_int.push_back(largest_tau);
        loop().effective_sample_sizes.back()=largest_tau>RealType{}
            ?global_sample_count/(RealType{2.}*largest_tau)
            :global_sample_count;
        loop().block_length_to_tau_ratios.push_back(largest_tau>RealType{}
            ?static_cast<RealType>(m_samples_per_block)/largest_tau
            :std::numeric_limits<RealType>::infinity());
        normalize_observables();
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
    for( const auto& block:m_blocks )
    {
        size_t flat{};
        for( size_t t=0;t<m_num_real_points;++t )
            for( size_t p=0;p<m_corr_directions.size();++p )
                for( size_t tau=0;tau<m_num_imaginary_edge_points;++tau,++flat )
                {
                    const ComplexType total=edge_value(
                        (*m_correlations),t,p,tau);
                    const ComplexType replicate=block_ratio(
                        total,edge_value(block.correlations,t,p,tau),block);
                    add_moment(corr_base+flat,replicate,
                        estimate(total));
                }

        for( size_t c=0;c<m_mag_directions.size();++c )
        {
            for( size_t t=0;t<m_num_real_points;++t )
            {
                const size_t flat_mag=t*m_mag_directions.size()+c;
                const ComplexType total{
                    m_magnetization->Re[t][c],m_magnetization->Im[t][c]};
                const ComplexType block_value{
                    block.magnetization.Re[t][c],block.magnetization.Im[t][c]};
                const ComplexType current=block_ratio(
                    total,block_value,block);
                const ComplexType center=estimate(total);
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
            const size_t flat_mag=t*m_mag_directions.size()+c;
            const auto mag_errors=errors(mag_base+flat_mag);
            loop().magnetization_time_stds.Re[t][c]=mag_errors[0];
            loop().magnetization_time_stds.Im[t][c]=mag_errors[1];
        }
    for( size_t t=0;t<m_num_real_points;++t )
    {
        const auto value=errors(closure_base+t);
        loop().closed_contour_residual_Re_sample_stds[t]=value[0];
        loop().closed_contour_residual_Im_sample_stds[t]=value[1];
        loop().closed_contour_residual_abs_sample_stds[t]=value[2];
    }
    normalize_observables();
}

void RunTimeData::normalize_observables()
{
    for(size_t t=0;t<m_num_real_points;++t)
    {
        for(size_t p=0;p<m_corr_directions.size();++p)
            for(size_t tau=0;tau<m_num_imaginary_edge_points;++tau)
            {
                const ComplexType value=exact_complex_ratio(
                    edge_value(*m_correlations,t,p,tau),m_total.partition,m_total.partition_abs);
                m_correlations->Re[t][p][tau]=std::real(value);
                m_correlations->Im[t][p][tau]=std::imag(value);
            }
        for(size_t c=0;c<m_mag_directions.size();++c)
        {
            const ComplexType value=exact_complex_ratio(
                ComplexType{m_magnetization->Re[t][c],m_magnetization->Im[t][c]},
                m_total.partition,m_total.partition_abs);
            m_magnetization->Re[t][c]=std::real(value);
            m_magnetization->Im[t][c]=std::imag(value);
        }
    }
    m_correlations=nullptr;m_magnetization=nullptr;
}

void RunTimeData::record_iteration_error(
    RealType absolute_error,RealType standardized_error )
{
    loop().absolute_iteration_errors.push_back(absolute_error);
    loop().standardized_iteration_errors.push_back(standardized_error);
}

void RunTimeData::record_matsubara_iteration_error(
    const CorrelationSet& old_correlations, const CorrelationSet& raw_correlations,
    const MagTen& old_magnetization_Re, const MagTen& raw_magnetization_Re,
    RealType standardized_error )
{
    // Match spinDMFT_cplx_field: maximum over component-wise mean absolute
    // correlation changes (Re/Im separately) and real magnetization changes.
    // Sampling errors do not set the gate.
    RealType largest{};
    auto scan=[&]( const CorrTen& old_values, const CorrTen& raw_values )
    {
        for( size_t c=0;c<old_values.size();++c )
        {
            RealType mean_abs{};
            for( size_t tau=0;tau<old_values[c].size();++tau )
                mean_abs+=std::abs(raw_values[c][tau]-old_values[c][tau]);
            mean_abs/=static_cast<RealType>(old_values[c].size());
            largest=std::isfinite(mean_abs)?std::max(largest,mean_abs)
                :std::numeric_limits<RealType>::infinity();
        }
    };
    scan(old_correlations.Re.front(),raw_correlations.Re.front());
    scan(old_correlations.Im.front(),raw_correlations.Im.front());
    for( size_t c=0;c<old_magnetization_Re.num_components();++c )
    {
        const RealType change=std::abs(
            raw_magnetization_Re.front()[c]-old_magnetization_Re.front()[c]);
        largest=std::isfinite(change)?std::max(largest,change)
            :std::numeric_limits<RealType>::infinity();
    }
    record_iteration_error(largest,standardized_error);
}

bool RunTimeData::diagnostics_pass( const MagTen& magnetization_Im ) const
{
    if( loop().covariance_symmetry_errors.empty()
        ||loop().gaussian_factor_reconstruction_errors.empty()||loop().average_phase_magnitudes.empty() ) return false;
    if( loop().covariance_symmetry_errors.back()>m_covariance_tolerance
        ||loop().gaussian_factor_reconstruction_errors.back()>m_takagi_tolerance
        ||loop().average_phase_magnitudes.back()<m_minimum_phase_magnitude ) return false;
    for( size_t c=0;c<m_mag_directions.size();++c )
        if( std::abs(magnetization_Im.front()[c])>
            m_iteration_error_sigma_threshold*
            loop().magnetization_time_stds.Im.front()[c] ) return false;
    return true;
}

void RunTimeData::finalize_iteration_step()
{
    if( m_my_rank==0 )
    {
        std::cout<<"Iteration step finished:\n"
          <<print::quantity_to_output_line(36,m_harmonic_bath
             ?"initial-guess difference (unused)":m_matsubara_preparation
                 ?"raw Matsubara absolute residual":"raw complex fixed-point residual",
             print::round_value_to_string(loop().absolute_iteration_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,"largest standardized residual",
             print::round_value_to_string(loop().standardized_iteration_errors.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,m_matsubara_preparation
             ?"absolute Matsubara threshold":"standardized residual threshold",
             print::round_value_to_string(m_matsubara_preparation
                 ?m_matsubara_iteration_error_threshold:m_iteration_error_sigma_threshold,
                 m_num_print_digits))
          <<print::quantity_to_output_line(36,m_pcn
             ?"sampled-Z phase magnitude":"average phase magnitude",
             print::round_value_to_string(loop().average_phase_magnitudes.back(),m_num_print_digits))
          <<print::quantity_to_output_line(36,m_pcn
             ?"minimum autocorrelation N_eff":"complex-weight N_eff",
             print::round_value_to_string(loop().effective_sample_sizes.back(),m_num_print_digits));
        std::cout
          <<print::quantity_to_output_line(36,"largest Gaussian factor block",
             std::to_string(loop().gaussian_largest_factorization_dimensions.back()))
          <<print::quantity_to_output_line(36,"closed-contour trace residual",
             print::round_value_to_string(loop().denominator_constancy_residuals.back(),m_num_print_digits));
        if( m_pcn )
        {
            std::cout
              <<print::quantity_to_output_line(36,"MH acceptance rate",
                 print::round_value_to_string(loop().mh_acceptance_rates.back(),m_num_print_digits))
              <<print::quantity_to_output_line(36,"zero-modulus rejection rate",
                 print::round_value_to_string(
                     loop().mh_nonpositive_rejection_rates.back(),m_num_print_digits));
        }
    }
    ++loop().num_Iterations;
}

bool RunTimeData::terminate( const MagTen& magnetization_Im )
{
    if( !m_self_consistency )
    {
        loop().termination=m_harmonic_bath
            ?"prescribed harmonic bath":"no self consistency";
        return true;
    }
    const bool fixed_point_converged=m_matsubara_preparation
        ?loop().absolute_iteration_errors.back()<m_matsubara_iteration_error_threshold
        // An exploding Monte Carlo variance can make any update look small
        // in sigma units. Also require a sub-spin-scale absolute update.
        :loop().standardized_iteration_errors.back()<m_iteration_error_sigma_threshold
            &&loop().absolute_iteration_errors.back()<m_real_time_residual_scale;
    if( fixed_point_converged&&diagnostics_pass(magnetization_Im) )
    {
        print::print_R0(m_my_rank,"\033[1;32mTerminating: fixed point and contour diagnostics converged.\033[0m\n");
        loop().termination="by convergence"; return true;
    }
    if( loop().num_Iterations>=m_iteration_limit )
    {
        print::print_R0(m_my_rank,"\033[1;31mTerminating because the iteration limit was reached.\033[0m\n");
        loop().termination="by iteration limit"; return true;
    }
    return false;
}

}
