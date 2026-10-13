#pragma once

#include<cstddef>
#include<random>

#include"Functions.h"

namespace spinDMFT::Functions
{

// Metropolized pCN chain whose likelihood is the real part of the selected
// observable denominator: Z_M for partition-function normalization and the
// fixed final D(T) for closed-contour normalization.
class PCNChain
{
 public:
    using LatentVector = JointComplexGaussianSampler::LatentVector;

    PCNChain( const ps::ParameterSpace& pspace,
              JointComplexGaussianSampler& sampler,
              const MeanFieldTrajectory& mean_field_time,
              RealType step_size,
              std::mt19937& engine );

    bool step();
    const ContourTrajectory& trajectory() const { return m_trajectory; }
    RealType real_sampling_weight() const { return m_sampling_weight_real; }
    size_t proposed() const { return m_proposed; }
    size_t accepted() const { return m_accepted; }
    size_t rejected_nonpositive() const { return m_rejected_nonpositive; }
    RealType maximum_relative_imaginary_sampling_weight() const
    { return m_maximum_relative_imaginary_sampling_weight; }

 private:
    bool evaluate( const LatentVector& latent, ContourTrajectory& trajectory,
                   RealType& sampling_weight_real );
    ComplexType sampling_weight( const ContourTrajectory& trajectory ) const;
    bool finite_positive_sampling_weight( ComplexType weight,
                                          RealType& weight_real );

    const ps::ParameterSpace& m_pspace;
    JointComplexGaussianSampler& m_sampler;
    const MeanFieldTrajectory& m_mean_field_time;
    std::mt19937& m_engine;
    RealType m_step_size{};
    RealType m_retention{};
    std::uniform_real_distribution<RealType> m_uniform01{RealType{0.},RealType{1.}};
    LatentVector m_latent{};
    ContourTrajectory m_trajectory{};
    ContourTrajectory m_proposed_trajectory{};
    JointComplexGaussianSampler::ContourFieldSample m_proposed_field{};
    TrajectoryWorkspace m_workspace{};
    LatentVector m_proposal{};
    RealType m_sampling_weight_real{};
    size_t m_proposed{};
    size_t m_accepted{};
    size_t m_rejected_nonpositive{};
    RealType m_maximum_relative_imaginary_sampling_weight{};
};

}
