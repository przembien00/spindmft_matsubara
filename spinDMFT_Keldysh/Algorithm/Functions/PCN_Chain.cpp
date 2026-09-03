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
          RealType{},RealType{1.}-step_size*step_size))),
      m_antithetic_pairs(pspace.antithetic_pairs)
{
    if( step_size<=RealType{}||step_size>RealType{1.} )
        throw std::invalid_argument("pCN step size must lie in (0,1]");
    constexpr size_t maximum_initialization_attempts=128;
    for( size_t attempt=0;attempt<maximum_initialization_attempts;++attempt )
    {
        auto latent=m_sampler.draw_latent(m_engine);
        ContourTrajectory trajectory;
        ContourTrajectory antithetic_trajectory;
        RealType partition_real{};
        if( evaluate(latent,trajectory,antithetic_trajectory,partition_real) )
        {
            m_latent=std::move(latent);
            m_trajectory=std::move(trajectory);
            m_antithetic_trajectory=std::move(antithetic_trajectory);
            m_partition_real=partition_real;
            return;
        }
    }
    throw std::runtime_error(
        m_antithetic_pairs
        ?"antithetic pCN could not initialize a finite state with positive Re Z_M for both signs"
        :"pCN could not initialize a finite state with positive Re Z_M");
}

bool PCNChain::evaluate( const LatentVector& latent,
                         ContourTrajectory& trajectory,
                         ContourTrajectory& antithetic_trajectory,
                         RealType& partition_real )
{
    const auto field=m_sampler.field_from_latent(latent);
    trajectory=compute_contour_trajectory(m_pspace,field,m_mean_field_time);
    RealType positive_partition{};
    if( !finite_positive_partition(trajectory,positive_partition) ) return false;
    if( !m_antithetic_pairs )
    {
        partition_real=positive_partition;
        return true;
    }

    auto antithetic_latent=latent;
    for( auto& value:antithetic_latent ) value=-value;
    const auto antithetic_field=m_sampler.field_from_latent(antithetic_latent);
    antithetic_trajectory=compute_contour_trajectory(
        m_pspace,antithetic_field,m_mean_field_time);
    RealType negative_partition{};
    if( !finite_positive_partition(antithetic_trajectory,negative_partition) )
        return false;
    partition_real=positive_partition+negative_partition;
    return std::isfinite(partition_real)&&partition_real>RealType{};
}

bool PCNChain::finite_positive_partition(
    const ContourTrajectory& trajectory, RealType& partition_real )
{
    const ComplexType Z=trajectory.partition_function;
    partition_real=std::real(Z);
    if( std::isfinite(partition_real)&&std::isfinite(std::imag(Z))
        &&partition_real>RealType{} )
    {
        const RealType scale=std::max(
            std::abs(partition_real),std::numeric_limits<RealType>::min());
        m_maximum_relative_imaginary_partition=std::max(
            m_maximum_relative_imaginary_partition,std::abs(std::imag(Z))/scale);
    }
    return std::isfinite(partition_real)&&std::isfinite(std::imag(Z))
           &&partition_real>RealType{};
}

const ContourTrajectory& PCNChain::antithetic_trajectory() const
{
    if( !m_antithetic_pairs )
        throw std::logic_error(
            "antithetic trajectory requested from an ordinary pCN chain");
    return m_antithetic_trajectory;
}

bool PCNChain::step()
{
    ++m_proposed;
    const auto innovation=m_sampler.draw_latent(m_engine);
    LatentVector proposal(m_latent.size());
    for( size_t i=0;i<proposal.size();++i )
        proposal[i]=m_retention*m_latent[i]+m_step_size*innovation[i];

    ContourTrajectory proposed_trajectory;
    ContourTrajectory proposed_antithetic_trajectory;
    RealType proposed_partition{};
    if( !evaluate(proposal,proposed_trajectory,
                  proposed_antithetic_trajectory,proposed_partition) )
    {
        ++m_rejected_nonpositive;
        return false;
    }

    const RealType log_alpha=std::log(proposed_partition)-std::log(m_partition_real);
    const RealType log_uniform=std::log(m_uniform01(m_engine));
    if( log_uniform>=std::min(RealType{},log_alpha) ) return false;

    m_latent=std::move(proposal);
    m_trajectory=std::move(proposed_trajectory);
    m_antithetic_trajectory=std::move(proposed_antithetic_trajectory);
    m_partition_real=proposed_partition;
    ++m_accepted;
    return true;
}

}
