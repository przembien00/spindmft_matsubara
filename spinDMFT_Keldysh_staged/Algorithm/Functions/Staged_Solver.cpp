#include"Staged_Solver.h"
#include<cmath>
#include<limits>

namespace spinDMFT::Functions
{
MatsubaraPCNChain::MatsubaraPCNChain(const ps::ParameterSpace& p,
    JointComplexGaussianSampler& sampler,const MeanFieldTrajectory& mean,std::mt19937& engine)
    :m_preparation(p),m_sampler(sampler),m_engine(engine),m_step(p.mh_step_size)
{
    if(!std::isfinite(m_step)||m_step<=RealType{}||m_step>RealType{1.})
        throw std::invalid_argument("Matsubara pCN step size must lie in (0,1]");
    if(mean.empty()||sampler.size()!=3*p.num_TimePoints)
        throw std::invalid_argument("Matsubara pCN requires an imaginary-only sampler and mean");
    m_preparation.num_RealTimeSteps=0;m_preparation.num_RealTimePoints=1;
    m_preparation.Tmax=RealType{};
    m_mean={mean.front()};
    m_retention=std::sqrt(std::max(RealType{},RealType{1.}-m_step*m_step));
    for(size_t attempt=0;attempt<128;++attempt)
    {
        auto x=m_sampler.draw_latent(m_engine);
        const RealType w=evaluate(x,m_partition);
        if(w>RealType{}){m_latent=std::move(x);m_weight=w;return;}
    }
    throw std::runtime_error("Matsubara pCN could not initialize a nonzero |Z_M| state");
}

RealType MatsubaraPCNChain::evaluate(const LatentVector& x,ComplexType& partition)
{
    const auto M=m_sampler.field_from_latent(x);
    JointComplexGaussianSampler::ContourFieldSample field;
    field.edge_field.resize(M.size()+6);field.edge_field=ComplexType{};
    for(size_t i=0;i<M.size();++i)field.edge_field[i]=M[i];
    build_contour_trajectory(m_preparation,field,m_mean,m_trajectory,m_workspace,true);
    partition=m_trajectory.partition_function;
    const RealType w=std::abs(partition);
    if(!std::isfinite(std::real(partition))||!std::isfinite(std::imag(partition))||!std::isfinite(w))
        throw std::runtime_error("Matsubara pCN encountered nonfinite Z_M");
    if(w>RealType{}&&!std::isfinite(RealType{1.}/w))
        throw std::runtime_error("Matsubara pCN inverse |Z_M| overflowed");
    return w;
}

bool MatsubaraPCNChain::step()
{
    ++m_proposed;
    const auto xi=m_sampler.draw_latent(m_engine);
    LatentVector proposal(m_latent.size());
    for(size_t i=0;i<proposal.size();++i)
        proposal[i]=m_retention*m_latent[i]+m_step*xi[i];
    ComplexType partition;
    const RealType w=evaluate(proposal,partition);
    if(w==RealType{}){++m_rejected_zero;return false;}
    const RealType log_alpha=std::log(w)-std::log(m_weight);
    if(std::log(m_uniform(m_engine))>=std::min(RealType{},log_alpha))return false;
    m_latent=std::move(proposal);m_weight=w;m_partition=partition;++m_accepted;
    return true;
}

std::vector<RealType> staggered_sampling_weights(const ps::ParameterSpace& p)
{
    if(!p.uses_staggered_strang())return {};
    const Contour::StaggeredLayout layout{p.num_TimePoints,
        p.num_RealTimeSteps*p.real_time_steps_per_interval()};
    std::vector<RealType> weights(layout.dimension());
    for(size_t i=0;i<weights.size();++i)
    {
        const auto sector=layout.decode(i).sector;
        weights[i]=sector==Contour::StaggeredSector::Matsubara?RealType{1.}
            :sector==Contour::StaggeredSector::Eta?RealType{2.}:RealType{0.5};
    }
    return weights;
}

std::unique_ptr<JointComplexGaussianSampler> make_staged_real_sampler(
    const ps::ParameterSpace& p,const ComplexDynamicMatrix& covariance)
{
    const auto weights=staggered_sampling_weights(p);
    if(!weights.empty()&&covariance.rows()!=weights.size())
        throw std::invalid_argument("staggered covariance does not match its sampling grid");
    if(!weights.empty())return std::make_unique<DenseComplexGaussianSampler>(covariance,weights);
    return std::make_unique<DenseComplexGaussianSampler>(covariance);
}

ComplexDynamicMatrix matsubara_block(const SelfConsistentField& field,size_t points)
{
    ComplexDynamicMatrix A(3*points,3*points);
    for(size_t i=0;i<A.rows();++i)for(size_t j=0;j<A.columns();++j)
        A(i,j)=field.covariance_source(i,j);
    return A;
}

template<typename Matrix>
void measure_matsubara_sample(const ps::ParameterSpace& p,
    JointComplexGaussianSampler& sampler,std::mt19937& engine,
    const MeanFieldTrajectory& mean,rtd::RunTimeData& data,
    BasicContourTrajectory<Matrix>& trajectory,BasicTrajectoryWorkspace<Matrix>& workspace,
    rtd::MeasuredSample& sample,MeasurementWorkspace& measurement)
{
    if(p.num_RealTimePoints!=1||p.num_RealTimeSteps!=0)
        throw std::invalid_argument("Matsubara preparation must have only the real-time origin");
    const auto M=sampler.draw(engine);
    JointComplexGaussianSampler::ContourFieldSample field;
    field.edge_field.resize(M.size()+6);field.edge_field=ComplexType{};
    for(size_t i=0;i<M.size();++i)field.edge_field[i]=M[i];
    build_contour_trajectory(p,field,mean,trajectory,workspace);
    measure_contour_observables(data,trajectory,sample,measurement,"prefix");
}

template void measure_matsubara_sample<Operator>(const ps::ParameterSpace&,
    JointComplexGaussianSampler&,std::mt19937&,const MeanFieldTrajectory&,rtd::RunTimeData&,
    BasicContourTrajectory<Operator>&,BasicTrajectoryWorkspace<Operator>&,
    rtd::MeasuredSample&,MeasurementWorkspace&);
template void measure_matsubara_sample<SpinHalfOperator>(const ps::ParameterSpace&,
    JointComplexGaussianSampler&,std::mt19937&,const MeanFieldTrajectory&,rtd::RunTimeData&,
    BasicContourTrajectory<SpinHalfOperator>&,BasicTrajectoryWorkspace<SpinHalfOperator>&,
    rtd::MeasuredSample&,MeasurementWorkspace&);

void freeze_origin(CorrelationSet& q,MagTen& re,MagTen& im,
    const CorrelationSet& prepared,const MagTen& prepared_re,const MagTen& prepared_im)
{
    q.Re.front()=prepared.Re.front();q.Im.front()=prepared.Im.front();
    re.front()=prepared_re.front();im.front()=prepared_im.front();
}

void seed_from_preparation(const ps::ParameterSpace& p,
    const CorrelationSet& prepared,const MagTen& prepared_re,const MagTen& prepared_im,
    CorrelationSet& q,MagTen& re,MagTen& im)
{
    re=MagTen{p.correlation_symmetry_type,p.num_RealTimePoints};
    im=MagTen{p.correlation_symmetry_type,p.num_RealTimePoints};
    for(size_t t=0;t<p.num_RealTimePoints;++t)
    {re[t]=prepared_re.front();im[t]=prepared_im.front();}
    // Reuse the source solver's selected real-time seed envelope/transpose.
    // Supplying the prepared t=0 primitive through the loader fields avoids
    // discarding its nonzero imaginary/transverse components.
    auto seed_parameters=p;
    seed_parameters.load_initial_spin_correlations=true;
    seed_parameters.old_num_TimePoints=p.num_TimePoints;
    seed_parameters.old_delta_t=p.delta_t;
    seed_parameters.initial_correlations_linearized.clear();
    seed_parameters.initial_correlations_imag_linearized.clear();
    for(size_t c=0;c<prepared.Re.front().size();++c)
        for(size_t tau=0;tau<p.num_TimePoints;++tau)
        {
            seed_parameters.initial_correlations_linearized.push_back(prepared.Re.front()[c][tau]);
            seed_parameters.initial_correlations_imag_linearized.push_back(prepared.Im.front()[c][tau]);
        }
    q=generate_initial_correlations(seed_parameters,re);
    freeze_origin(q,re,im,prepared,prepared_re,prepared_im);
}
}
