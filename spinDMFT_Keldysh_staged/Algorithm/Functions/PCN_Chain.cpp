#include"PCN_Chain.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<stdexcept>

namespace spinDMFT::Functions
{

PCNChain::PCNChain( const ps::ParameterSpace& pspace,
                    JointComplexGaussianSampler& sampler,
                    const MeanFieldTrajectory& mean_field_time,
                    const RealType step_size,
                    std::mt19937& engine )
    : m_pspace(pspace),m_sampler(sampler),m_mean_field_time(mean_field_time),
      m_engine(engine),m_step_size(step_size),
      m_retention(std::sqrt(std::max(
          RealType{},RealType{1.}-step_size*step_size)))
{
    if(dynamic_cast<WeightedDenseComplexGaussianSampler*>(&sampler))
        throw std::invalid_argument(
            "weighted-dense requires independent sampling; pCN weight positivity is not established for this ensemble");
    if( step_size<=RealType{}||step_size>RealType{1.} )
        throw std::invalid_argument("pCN step size must lie in (0,1]");
    constexpr size_t maximum_initialization_attempts=128;
    for( size_t attempt=0;attempt<maximum_initialization_attempts;++attempt )
    {
        auto latent=m_sampler.draw_latent(m_engine);
        ContourTrajectory trajectory;
        RealType sampling_weight_real{};
        if( evaluate(latent,trajectory,
                     sampling_weight_real) )
        {
            m_latent=std::move(latent);
            m_trajectory=std::move(trajectory);
            if( m_pspace.correlation_normalization!="closed-contour" )
                complete_contour_trajectory(m_pspace,m_proposed_field,m_mean_field_time,
                                           m_trajectory,m_workspace);
            m_sampling_weight_real=sampling_weight_real;
            return;
        }
    }
    const std::string denominator=m_pspace.correlation_normalization
        =="closed-contour"?"D(T)":"Z_M";
    throw std::runtime_error("pCN could not initialize a finite state with positive Re "+denominator);
}

bool PCNChain::evaluate( const LatentVector& latent,
                         ContourTrajectory& trajectory,
                         RealType& sampling_weight_real )
{
    m_proposed_field=m_sampler.contour_field_from_latent(latent,m_pspace.uses_cf4());
    build_contour_trajectory(m_pspace,m_proposed_field,m_mean_field_time,
        trajectory,m_workspace,m_pspace.correlation_normalization!="closed-contour");
    return finite_positive_sampling_weight(sampling_weight(trajectory),sampling_weight_real);
}

ComplexType PCNChain::sampling_weight( const ContourTrajectory& trajectory ) const
{
    return m_pspace.correlation_normalization=="closed-contour"
        ?trajectory.final_closed_contour_trace:trajectory.partition_function;
}

bool PCNChain::finite_positive_sampling_weight(
    const ComplexType weight, RealType& weight_real )
{
    weight_real=std::real(weight);
    if( std::isfinite(weight_real)&&std::isfinite(std::imag(weight))
        &&weight_real>RealType{} )
    {
        const RealType scale=std::max(
            std::abs(weight_real),std::numeric_limits<RealType>::min());
        m_maximum_relative_imaginary_sampling_weight=std::max(
            m_maximum_relative_imaginary_sampling_weight,
            std::abs(std::imag(weight))/scale);
    }
    return std::isfinite(weight_real)&&std::isfinite(std::imag(weight))
           &&weight_real>RealType{};
}

bool PCNChain::step()
{
    ++m_proposed;
    const auto innovation=m_sampler.draw_latent(m_engine);
    m_proposal.resize(m_latent.size());
    auto& proposal=m_proposal;
    for( size_t i=0;i<proposal.size();++i )
        proposal[i]=m_retention*m_latent[i]+m_step_size*innovation[i];

    auto& proposed_trajectory=m_proposed_trajectory;
    RealType proposed_partition{};
    if( !evaluate(proposal,proposed_trajectory,proposed_partition) )
    {
        ++m_rejected_nonpositive;
        return false;
    }

    const RealType log_alpha=std::log(proposed_partition)
                            -std::log(m_sampling_weight_real);
    const RealType log_uniform=std::log(m_uniform01(m_engine));
    if( log_uniform>=std::min(RealType{},log_alpha) ) return false;

    if( m_pspace.correlation_normalization!="closed-contour" )
        complete_contour_trajectory(m_pspace,m_proposed_field,m_mean_field_time,
                                   proposed_trajectory,m_workspace);
    std::swap(m_latent,proposal);
    std::swap(m_trajectory,proposed_trajectory);
    m_sampling_weight_real=proposed_partition;
    ++m_accepted;
    return true;
}

}
