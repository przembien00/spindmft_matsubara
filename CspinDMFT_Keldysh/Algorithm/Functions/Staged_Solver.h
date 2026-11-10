#pragma once
#include "../Run_Time_Data/Run_Time_Data.h"
#include <Complex_Gaussian/Joint_Matsubara_PCN.h>
namespace CspinDMFT::Functions {
// One pCN chain on the current joint Gaussian law's Matsubara row space.
// Main controls burn-in and the production loop. Each production draw also
// refreshes the orthogonal real-field latents.
class MatsubaraPCNChain {
public:
    MatsubaraPCNChain(const Model&,const ParameterSpace&,Sampler&,const std::vector<C>& mean,std::mt19937&);
    void step();
    Field draw_field();
    double weight() const { return m_weight; }
    double acceptance_rate() const { return m_proposed?double(m_accepted)/double(m_proposed):0.; }
    double zero_modulus_rejection_rate() const { return m_proposed?double(m_rejected_zero)/double(m_proposed):0.; }
    void reset_statistics() { m_accepted=m_proposed=m_rejected_zero=0; }
private:
    double evaluate(const Sampler::LatentVector&);
    const Model& m_model;
    const ParameterSpace& m_parameters;
    Sampler& m_sampler;
    std::mt19937& m_engine;
    spinDMFT::Functions::JointMatsubaraSampler m_marginal;
    std::vector<C> m_mean;
    Sampler::LatentVector m_active;
    double m_weight=0.;
    size_t m_accepted=0,m_proposed=0,m_rejected_zero=0;
    std::uniform_real_distribution<double> m_uniform{0.,1.};
};
}
