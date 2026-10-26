#pragma once
#include"Functions.h"
#include"Complex_Gaussian.h"

namespace spinDMFT::Functions
{
// pCN acts only on the real M latents, targeting p0(x)|Z_M(Fx)|.
// Production observables and their complex denominator must both be divided
// by |Z_M|. Real residual latents are not members of this chain.
class MatsubaraPCNChain
{
 public:
    using LatentVector=JointComplexGaussianSampler::LatentVector;
    MatsubaraPCNChain(const ps::ParameterSpace& p,JointComplexGaussianSampler& sampler,
        const MeanFieldTrajectory& mean,std::mt19937& engine);
    bool step();
    const LatentVector& latent() const { return m_latent; }
    ComplexType partition() const { return m_partition; }
    RealType weight() const { return m_weight; }
    size_t accepted() const { return m_accepted; }
    size_t proposed() const { return m_proposed; }
    size_t rejected_zero_weight() const { return m_rejected_zero; }
 private:
    RealType evaluate(const LatentVector& x,ComplexType& partition);
    const ps::ParameterSpace& m_parameters;
    JointComplexGaussianSampler& m_sampler;
    MeanFieldTrajectory m_mean;
    std::mt19937& m_engine;
    RealType m_step{},m_retention{},m_weight{};
    ComplexType m_partition{};
    LatentVector m_latent;
    ContourTrajectory m_trajectory;
    TrajectoryWorkspace m_workspace;
    std::uniform_real_distribution<RealType> m_uniform{RealType{},RealType{1.}};
    size_t m_accepted{},m_proposed{},m_rejected_zero{};
};
// Branch-equivalent norm: |V+|^2+|V-|^2=2|eta|^2+|nu|^2/2.
std::vector<RealType> staggered_sampling_weights(const ps::ParameterSpace& p);
std::unique_ptr<JointComplexGaussianSampler> make_staged_real_sampler(
    const ps::ParameterSpace& p,const ComplexDynamicMatrix& covariance);
ComplexDynamicMatrix matsubara_block(const SelfConsistentField& field,size_t points);
// Construct a genuine Matsubara trajectory: only M is drawn; real fields
// are zero placeholders at the single origin vertex and are never propagated.
template<typename Matrix>
void measure_matsubara_sample(const ps::ParameterSpace& preparation,
    JointComplexGaussianSampler& sampler,std::mt19937& engine,
    const MeanFieldTrajectory& mean,rtd::RunTimeData& data,
    BasicContourTrajectory<Matrix>& trajectory,BasicTrajectoryWorkspace<Matrix>& workspace,
    MeasuredSample& sample,MeasurementWorkspace& measurement);
// Replace feedback at the origin only; raw measured observables stay intact.
void freeze_origin(CorrelationSet& correlations,Magnetization& magnetization,
    const CorrelationSet& prepared,const Magnetization& prepared_magnetization);
void seed_from_preparation(const ps::ParameterSpace& real_parameters,
    const CorrelationSet& prepared,const Magnetization& prepared_magnetization,
    CorrelationSet& correlations,Magnetization& magnetization);
}
