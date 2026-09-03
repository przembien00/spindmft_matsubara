#include"../Functions/Complex_Gaussian.h"

#include<algorithm>
#include<cmath>
#include<iostream>
#include<memory>
#include<stdexcept>

namespace func = spinDMFT::Functions;

namespace
{

int require( bool condition, const char* message )
{
    if( condition ) return 0;
    std::cerr << "FAILED: " << message << '\n';
    return 1;
}

RealType relative_residual( const func::ComplexDynamicMatrix& actual,
                            const func::ComplexDynamicMatrix& expected )
{
    RealType residual{},scale{};
    for( size_t i=0;i<actual.rows();++i )
        for( size_t j=0;j<actual.columns();++j )
        {
            residual+=std::norm(actual(i,j)-expected(i,j));
            scale+=std::norm(expected(i,j));
        }
    return scale>RealType{0.}?std::sqrt(residual/scale):std::sqrt(residual);
}

func::ComplexDynamicMatrix empirical_pseudo_covariance(
    func::JointComplexGaussianSampler& sampler,
    std::mt19937& engine,const size_t draws )
{
    func::ComplexDynamicMatrix result(sampler.size(),sampler.size(),ComplexType{});
    for( size_t sample=0;sample<draws;++sample )
    {
        const auto z=sampler.draw(engine);
        for( size_t i=0;i<z.size();++i )
            for( size_t j=0;j<z.size();++j ) result(i,j)+=z[i]*z[j];
    }
    result/=static_cast<RealType>(draws);
    return result;
}

func::ComplexDynamicMatrix hermitian_covariance( const func::TakagiFactor& factor )
{
    return factor.L*blaze::ctrans(factor.L);
}

func::ComplexDynamicMatrix make_structured_keldysh_covariance()
{
    constexpr size_t M=3,R=3,plus=M,minus=M+R;
    func::ComplexDynamicMatrix A(M,M,ComplexType{}),S(R,R,ComplexType{});
    func::ComplexDynamicMatrix B(M,R,ComplexType{}),response(R,R,ComplexType{});
    A(0,0)={1.2,0.2}; A(1,1)={-0.7,0.1}; A(2,2)={0.8,-0.3};
    A(0,1)=A(1,0)={0.2,-0.1}; A(1,2)=A(2,1)={-0.3,0.2};
    S(0,0)={0.9,-0.1}; S(1,1)={0.6,0.3}; S(2,2)={-0.5,0.2};
    S(0,2)=S(2,0)={0.15,0.1};
    B(0,0)={0.3,0.2}; B(0,2)={-0.1,0.4}; B(1,1)={0.2,-0.3};
    B(2,0)={-0.25,0.1};
    response(0,0)={0.1,0.2}; response(0,1)={-0.3,0.1};
    response(1,0)={0.2,-0.1}; response(1,2)={0.15,0.25};
    response(2,1)={-0.2,0.05};

    func::ComplexDynamicMatrix covariance(M+2*R,M+2*R,ComplexType{});
    for( size_t i=0;i<M;++i )
        for( size_t j=0;j<M;++j ) covariance(i,j)=A(i,j);
    for( size_t i=0;i<M;++i )
        for( size_t j=0;j<R;++j )
        {
            covariance(i,plus+j)=covariance(i,minus+j)=B(i,j);
            covariance(plus+j,i)=covariance(minus+j,i)=B(i,j);
        }
    for( size_t i=0;i<R;++i )
        for( size_t j=0;j<R;++j )
        {
            const auto Rij=response(i,j),Rji=response(j,i);
            covariance(plus+i,plus+j)=S(i,j)+RealType{0.5}*(Rij+Rji);
            covariance(plus+i,minus+j)=S(i,j)+RealType{0.5}*(-Rij+Rji);
            covariance(minus+i,plus+j)=S(i,j)+RealType{0.5}*(Rij-Rji);
            covariance(minus+i,minus+j)=S(i,j)-RealType{0.5}*(Rij+Rji);
        }
    return covariance;
}

func::ComplexDynamicMatrix make_periodic_doubled_frequency_covariance()
{
    constexpr size_t NM=1,NR=2,L=2*NR,M=3,R=3*NR,plus=M,minus=M+R;
    func::ComplexDynamicMatrix covariance(M+2*R,M+2*R,ComplexType{});
    const RealType two_pi=RealType{2.}*std::acos(RealType{-1.});
    func::ComplexDynamicMatrix F(L,L);
    for( size_t frequency=0;frequency<L;++frequency )
        for( size_t time=0;time<L;++time )
        {
            const RealType phase=-two_pi*static_cast<RealType>(frequency*time)
                                 /static_cast<RealType>(L);
            F(frequency,time)=std::exp(ComplexType{0.,phase})
                             /std::sqrt(static_cast<RealType>(L));
        }

    for( size_t component=0;component<3;++component )
    {
        func::ComplexDynamicMatrix frequency_covariance(1+2*L,1+2*L,ComplexType{});
        frequency_covariance(0,0)=ComplexType{1.1,0.2};
        const auto add_group=[&](const std::vector<size_t>& frequencies,
                                 const func::ComplexDynamicMatrix& S,
                                 const func::ComplexDynamicMatrix& response,
                                 const std::vector<ComplexType>& B)
        {
            for( size_t i=0;i<frequencies.size();++i )
            {
                const size_t pi=1+frequencies[i],mi=1+L+frequencies[i];
                frequency_covariance(0,pi)=frequency_covariance(pi,0)=B[i];
                frequency_covariance(0,mi)=frequency_covariance(mi,0)=B[i];
                for( size_t j=0;j<frequencies.size();++j )
                {
                    const size_t pj=1+frequencies[j],mj=1+L+frequencies[j];
                    const ComplexType Rij=response(i,j),Rji=response(j,i),Sij=S(i,j);
                    frequency_covariance(pi,pj)=Sij+RealType{0.5}*(Rij+Rji);
                    frequency_covariance(pi,mj)=Sij+RealType{0.5}*(-Rij+Rji);
                    frequency_covariance(mi,pj)=Sij+RealType{0.5}*(Rij-Rji);
                    frequency_covariance(mi,mj)=Sij-RealType{0.5}*(Rij+Rji);
                }
            }
        };
        func::ComplexDynamicMatrix S0(1,1),response0(1,1),S2(1,1),response2(1,1);
        S0(0,0)={0.4,0.1}; response0(0,0)={0.2,-0.1};
        S2(0,0)={-0.25,0.08}; response2(0,0)={0.14,0.03};
        add_group({0},S0,response0,{{0.1,0.05}});
        add_group({2},S2,response2,{{-0.08,0.04}});
        func::ComplexDynamicMatrix S13(2,2,ComplexType{}),response13(2,2,ComplexType{});
        S13(0,1)=S13(1,0)={0.12,0.02};
        response13(0,1)={0.18,-0.04}; response13(1,0)={-0.09,0.07};
        add_group({1,3},S13,response13,{{0.07,-0.03},{-0.02,0.06}});

        covariance(component,component)=frequency_covariance(0,0);
        for( size_t time=0;time<NR;++time )
        {
            ComplexType mixed_plus{},mixed_minus{};
            for( size_t frequency=0;frequency<L;++frequency )
            {
                mixed_plus+=frequency_covariance(0,1+frequency)*std::conj(F(frequency,time));
                mixed_minus+=frequency_covariance(0,1+L+frequency)*std::conj(F(frequency,time));
            }
            covariance(component,plus+3*time+component)=mixed_plus;
            covariance(plus+3*time+component,component)=mixed_plus;
            covariance(component,minus+3*time+component)=mixed_minus;
            covariance(minus+3*time+component,component)=mixed_minus;
            for( size_t other_time=0;other_time<NR;++other_time )
                for( size_t branch_i=0;branch_i<2;++branch_i )
                    for( size_t branch_j=0;branch_j<2;++branch_j )
                    {
                        ComplexType value{};
                        const size_t fi_offset=1+branch_i*L;
                        const size_t fj_offset=1+branch_j*L;
                        for( size_t fi=0;fi<L;++fi )
                            for( size_t fj=0;fj<L;++fj )
                                value+=std::conj(F(fi,time))*frequency_covariance(
                                    fi_offset+fi,fj_offset+fj)*std::conj(F(fj,other_time));
                        const size_t row=(branch_i?minus:plus)+3*time+component;
                        const size_t column=(branch_j?minus:plus)+3*other_time+component;
                        covariance(row,column)=value;
                    }
        }
    }
    return covariance;
}

func::ComplexDynamicMatrix make_doubled_frequency_covariance()
{
    // Start from one transformed Matsubara point and insert an independent,
    // untransformed beta- endpoint before the two real branches.
    const auto periodic=make_periodic_doubled_frequency_covariance();
    constexpr size_t old_matsubara_size=3,endpoint_size=3;
    func::ComplexDynamicMatrix covariance(
        periodic.rows()+endpoint_size,periodic.columns()+endpoint_size,
        ComplexType{});
    const auto expanded_index=[=](const size_t old_index)
    {
        return old_index<old_matsubara_size
            ?old_index:old_index+endpoint_size;
    };
    for( size_t i=0;i<periodic.rows();++i )
        for( size_t j=0;j<periodic.columns();++j )
            covariance(expanded_index(i),expanded_index(j))=periodic(i,j);
    for( size_t component=0;component<3;++component )
    {
        const size_t endpoint=old_matsubara_size+component;
        covariance(endpoint,endpoint)=ComplexType{
            RealType{0.8}+RealType{0.1}*component,
            RealType{-0.12}+RealType{0.03}*component};
        covariance(endpoint,component)=covariance(component,endpoint)=
            ComplexType{RealType{0.13}*(component+1),RealType{-0.04}};
        for( size_t old_real=old_matsubara_size;
             old_real<periodic.rows();++old_real )
        {
            const size_t real=expanded_index(old_real);
            const ComplexType value=RealType{0.7}*periodic(component,old_real)
                +(old_real%3==component
                    ?ComplexType{RealType{0.025}*(component+1),RealType{0.015}}
                    :ComplexType{});
            covariance(endpoint,real)=covariance(real,endpoint)=value;
        }
    }
    return covariance;
}

}

int main()
{
    int failures{};
    {
        func::ComplexDynamicMatrix Gamma( 1, 1 );
        Gamma(0,0) = ComplexType{-2.,0.};
        const auto factor = func::autonne_takagi( Gamma );
        failures += require( factor.numerical_rank == 1, "negative scalar rank" );
        failures += require( factor.reconstruction_error < RealType{1e-12},
                             "negative scalar reconstruction" );
        failures += require( std::abs(std::real(factor.L(0,0))) < RealType{1e-12},
                             "negative scalar becomes an imaginary field direction" );
    }
    {
        func::ComplexDynamicMatrix Gamma( 3, 3, ComplexType{} );
        Gamma(0,0) = ComplexType{1.2,0.3};
        Gamma(0,1) = Gamma(1,0) = ComplexType{-0.4,0.2};
        Gamma(0,2) = Gamma(2,0) = ComplexType{0.1,-0.5};
        Gamma(1,1) = ComplexType{-0.7,0.4};
        Gamma(1,2) = Gamma(2,1) = ComplexType{0.6,0.1};
        Gamma(2,2) = ComplexType{0.9,-0.2};
        const auto factor = func::autonne_takagi( Gamma );
        failures += require( factor.reconstruction_error < RealType{1e-12},
                             "generic complex symmetric reconstruction" );

        func::DenseComplexGaussianSampler sampler( Gamma );
        failures += require( sampler.size() == 3, "dense sampler size" );
        failures += require( sampler.numerical_rank() == factor.numerical_rank,
                             "dense sampler rank" );
        std::mt19937 engine{1234};
        const auto field = sampler.draw(engine);
        failures += require( field.size() == 3, "dense sampler draw size" );

        func::ComplexDynamicMatrix empirical( 3, 3, ComplexType{} );
        constexpr size_t draws = 50000;
        for( size_t sample = 0; sample < draws; ++sample )
        {
            const auto z = sampler.draw(engine);
            for( size_t i = 0; i < 3; ++i )
                for( size_t j = 0; j < 3; ++j )
                    empirical(i,j) += z[i]*z[j]; // unconjugated pseudo-covariance
        }
        empirical /= static_cast<RealType>(draws);
        RealType residual{}, scale{};
        for( size_t i = 0; i < 3; ++i )
            for( size_t j = 0; j < 3; ++j )
            {
                residual += std::norm(empirical(i,j)-Gamma(i,j));
                scale += std::norm(Gamma(i,j));
            }
        failures += require( std::sqrt(residual/scale) < RealType{0.02},
                             "empirical unconjugated covariance" );

        auto latent=sampler.draw_latent(engine);
        const auto positive=sampler.field_from_latent(latent);
        for( auto& value:latent ) value=-value;
        const auto negative=sampler.field_from_latent(latent);
        RealType antithetic_residual{};
        for( size_t i=0;i<positive.size();++i )
            antithetic_residual=std::max(
                antithetic_residual,std::abs(positive[i]+negative[i]));
        failures += require( antithetic_residual<RealType{1e-13},
                             "negated latent vector produces the antithetic field" );
    }
    {
        func::ComplexDynamicMatrix nonsymmetric( 2, 2, ComplexType{} );
        nonsymmetric(0,1) = ComplexType{1.,0.};
        bool rejected{};
        try { (void)func::autonne_takagi(nonsymmetric); }
        catch( const std::invalid_argument& ) { rejected = true; }
        failures += require( rejected, "nonsymmetric input is rejected" );
    }
    {
        // A response field may have E[nu^2]=0 while remaining nonzero sample by
        // sample through its correlation with eta.
        func::ComplexDynamicMatrix keldysh( 2, 2, ComplexType{} );
        keldysh(0,0)=ComplexType{1.,0.};
        keldysh(0,1)=keldysh(1,0)=ComplexType{0.,0.4};
        func::DenseComplexGaussianSampler sampler(keldysh);
        std::mt19937 engine{9981};
        ComplexType nu_square{}; RealType nu_abs_square{};
        constexpr size_t draws=100000;
        for( size_t sample=0;sample<draws;++sample )
        {
            const auto eta_nu=sampler.draw(engine);
            nu_square+=eta_nu[1]*eta_nu[1];
            nu_abs_square+=std::norm(eta_nu[1]);
        }
        nu_square/=static_cast<RealType>(draws);
        nu_abs_square/=static_cast<RealType>(draws);
        failures += require( std::abs(nu_square)<RealType{0.01},
                             "zero nu-nu pseudo-covariance" );
        failures += require( nu_abs_square>RealType{0.1},
                             "zero pseudo-variance does not set sampled nu to zero" );
    }
    {
        const auto covariance=make_structured_keldysh_covariance();
        func::SVDComplexGaussianSampler sampler(covariance);
        failures+=require(sampler.size()==9,"SVD sampler size");
        failures+=require(sampler.reconstruction_error()<RealType{1e-12},
                          "SVD exact structured reconstruction");
        failures+=require(sampler.latent_dimension()==9,
                          "SVD sampler latent dimension");
        std::mt19937 engine{7712};
        const auto empirical=empirical_pseudo_covariance(sampler,engine,100000);
        failures+=require(relative_residual(empirical,covariance)<RealType{0.025},
                          "SVD empirical pseudo-covariance");

        auto dense=func::make_complex_gaussian_sampler("dense",covariance,1,1);
        auto svd=func::make_complex_gaussian_sampler("svd",covariance,1,1);
        failures+=require(dense->size()==svd->size(),"sampler factory alternatives");

        const auto dense_factor=func::autonne_takagi(covariance);
        const auto svd_factor=func::svd_takagi(covariance);
        failures+=require(relative_residual(
            hermitian_covariance(svd_factor),hermitian_covariance(dense_factor))
            <RealType{1e-11},"SVD and dense Hermitian covariance equivalence");
    }
    {
        // Rank-zero input must not create null-space noise.
        func::ComplexDynamicMatrix zero(9,9,ComplexType{});
        func::SVDComplexGaussianSampler sampler(zero);
        std::mt19937 engine{9012};
        const auto field=sampler.draw(engine);
        RealType norm{};
        for( const auto value:field ) norm+=std::norm(value);
        failures+=require(norm==RealType{0.},
                          "zero SVD covariance produces zero field");
    }
    {
        // Exercises exact singular-value degeneracies and a larger covariance.
        const auto covariance=make_doubled_frequency_covariance();
        const auto dense_factor=func::autonne_takagi(covariance);
        const auto svd_factor=func::svd_takagi(covariance);
        failures+=require(svd_factor.reconstruction_error<RealType{1e-11},
                          "degenerate SVD Takagi reconstruction");
        failures+=require(relative_residual(
            hermitian_covariance(svd_factor),hermitian_covariance(dense_factor))
            <RealType{1e-11},"degenerate SVD Hermitian covariance equivalence");
    }
    {
        const std::array<size_t,5> expected{4,3,2,1,0};
        for( size_t point=0;point<expected.size();++point )
            failures+=require(func::one_sided_edge_reflection_index(
                                  point,expected.size()-1)==expected[point],
                              "one-sided reflection exchanges 0+ and beta-");
    }
    {
        // The FFT path changes only representation: factorize and sample the
        // doubled covariance blockwise in frequency space, inverse transform,
        // and retain the original physical interval.
        const auto covariance=make_doubled_frequency_covariance();
        failures+=require(std::abs(covariance(0,6)-covariance(3,6))
                                  >RealType{1e-8},
                          "FFT test keeps beta- covariance distinct from 0+");
        func::FFTDenseComplexGaussianSampler sampler(covariance,1,2,1.,-1.);
        failures+=require(sampler.size()==covariance.rows(),"FFT sampler physical size");
        failures+=require(sampler.reconstruction_error()<RealType{1e-11},
                          "FFT doubled-grid reconstruction");
        std::mt19937 engine{4189};
        const auto empirical=empirical_pseudo_covariance(sampler,engine,150000);
        failures+=require(relative_residual(empirical,covariance)<RealType{0.03},
                          "FFT empirical physical pseudo-covariance");

        auto factory=func::make_complex_gaussian_sampler("fft",covariance,1,2);
        failures+=require(factory->size()==covariance.rows(),"FFT sampler factory");
        failures+=require(factory->latent_dimension()>factory->size(),
                          "FFT sampler keeps doubled-grid latent field");

        func::FFTDenseComplexGaussianSampler truncated(covariance,1,2,1.,1.);
        failures+=require(truncated.reconstruction_error()<RealType{1e-11},
                          "FFT frequency-block reconstruction");
        failures+=require(truncated.covariance_approximation_error()>RealType{0.},
                          "FFT high-frequency cross covariance is discarded");
        failures+=require(truncated.largest_factorization_dimension()<27,
                          "FFT cutoff reduces the largest dense factorization");
        failures+=require(truncated.draw(engine).size()==covariance.rows(),
                          "truncated FFT sampler physical draw size");
    }
    {
        // Every zero-rank frequency block must explicitly clear its persistent
        // output buffer rather than retain data from an earlier draw.
        func::ComplexDynamicMatrix zero(18,18,ComplexType{});
        func::FFTDenseComplexGaussianSampler sampler(zero,1,2,1.,1.);
        std::mt19937 engine{6619};
        const auto field=sampler.draw(engine);
        RealType norm{};
        for( const auto value:field ) norm+=std::norm(value);
        failures+=require(sampler.latent_dimension()==0,
                          "zero FFT covariance has zero latent rank");
        failures+=require(norm==RealType{0.},
                          "zero FFT covariance produces zero field");
    }
    {
        bool rejected{};
        try
        {
            func::ComplexDynamicMatrix covariance(9,9,ComplexType{});
            (void)func::make_complex_gaussian_sampler("unknown",covariance,1,1);
        }
        catch( const std::invalid_argument& ) { rejected=true; }
        failures+=require(rejected,"unknown sampler algorithm is rejected");
    }
    return failures == 0 ? 0 : 1;
}
