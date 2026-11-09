#include "../main_header.h"
namespace CspinDMFT::Functions {
void freeze_origin(ContourCorrelations& correlations,TimeMagnetization& moments,
                   const MatsubaraCorrelations& prepared,const ClusterMagnetization& prepared_moments) {
    for(size_t i=0;i<correlations.sites;++i) for(size_t a:moments.directions)
        moments.set(0,3*i+a,prepared_moments.m(3*i+a));
    for(auto ij:correlations.categories) for(auto ab:correlations.directions)
        for(size_t k=0;k<correlations.imaginary_points;++k)
            correlations.set(0,3*ij[0]+ab[0],3*ij[1]+ab[1],k,prepared.x(3*ij[0]+ab[0],3*ij[1]+ab[1],k));
}
void seed_from_preparation(const MatsubaraCorrelations& prepared,const ClusterMagnetization& prepared_moments,
                           ContourCorrelations& correlations,TimeMagnetization& moments) {
    for(size_t t=0;t<correlations.times;++t) {
        for(size_t i=0;i<correlations.sites;++i) for(size_t a:moments.directions)
            moments.set(t,3*i+a,prepared_moments.m(3*i+a));
        for(auto ij:correlations.categories) for(auto ab:correlations.directions)
            for(size_t k=0;k<correlations.imaginary_points;++k)
                correlations.set(t,3*ij[0]+ab[0],3*ij[1]+ab[1],k,prepared.x(3*ij[0]+ab[0],3*ij[1]+ab[1],k));
    }
}
void constrain_feedback(ContourCorrelations& correlations,TimeMagnetization& moments,
    const MatsubaraCorrelations& prepared,const ClusterMagnetization& prepared_moments,const ParameterSpace& p) {
    freeze_origin(correlations,moments,prepared,prepared_moments);
    if(p.constant_magnetization_time) for(size_t t=0;t<moments.times;++t) for(size_t i=0;i<moments.sites;++i)
        for(size_t a:moments.directions) moments.set(t,3*i+a,prepared_moments.m(3*i+a));
}
}

namespace CspinDMFT::Functions {
MatsubaraPCNChain::MatsubaraPCNChain(const Model& model,const ParameterSpace& p,Sampler& sampler,
    const std::vector<C>& mean,std::mt19937& engine)
    :m_model(model),m_parameters(p),m_sampler(sampler),m_engine(engine),
     m_marginal(sampler,p.preparation_layout.dimension()),
     m_mean(mean.begin(),mean.begin()+3*model.sites) {
    for(size_t retry=0;retry<1000 && !(m_weight>0.);++retry) {
        m_active=m_marginal.draw_latent(m_engine); m_weight=evaluate(m_active);
    }
    if(!(m_weight>0.)) throw std::runtime_error("could not initialize a nonzero pCN partition");
}
double MatsubaraPCNChain::evaluate(const Sampler::LatentVector& latent) {
    const auto field=m_marginal.field_from_latent(latent);
    const C z=measure<MatsubaraCorrelations>(m_model,m_parameters,m_parameters.preparation_layout,field,m_mean,true).partition;
    if(!std::isfinite(z.real()) || !std::isfinite(z.imag())) throw std::runtime_error("nonfinite pCN Matsubara partition");
    return std::abs(z);
}
void MatsubaraPCNChain::step() {
    const auto innovation=m_marginal.draw_latent(m_engine);
    Sampler::LatentVector proposal=std::sqrt(1.-m_parameters.mh_step_size*m_parameters.mh_step_size)*m_active
        +m_parameters.mh_step_size*innovation;
    const double next=evaluate(proposal); ++m_proposed;
    if(next==0.) ++m_rejected_zero;
    if(next>0. && (next>=m_weight || m_uniform(m_engine)<next/m_weight)) {
        m_active=std::move(proposal); m_weight=next; ++m_accepted;
    }
}
Field MatsubaraPCNChain::draw_field() {
    return m_sampler.field_from_latent(m_marginal.draw_joint_latent(m_active,m_engine));
}
}
