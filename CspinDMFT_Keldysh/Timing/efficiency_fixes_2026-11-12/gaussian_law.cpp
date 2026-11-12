#include <Complex_Gaussian/Complex_Gaussian.h>
#include <iomanip>
#include <iostream>

namespace f=spinDMFT::Functions;
void emit(const f::ComplexDynamicMatrix& matrix) {
    std::cout<<matrix.rows()<<' '<<matrix.columns()<<'\n';
    for(size_t i=0;i<matrix.rows();++i)for(size_t j=0;j<matrix.columns();++j)
        std::cout<<matrix(i,j).real()<<' '<<matrix(i,j).imag()<<'\n';
}
int main() {
    std::cout<<std::setprecision(17);
    std::mt19937 rng(7193); std::normal_distribution<double> normal;
    for(size_t n:{size_t{1},size_t{9},size_t{12}})for(bool split:{false,true}) {
        f::ComplexDynamicMatrix covariance(n,n,ComplexType{});
        for(size_t i=0;i<n;++i)for(size_t j=i;j<n;++j) {
            const bool same_block=!split || i/3==j/3;
            if(same_block) covariance(i,j)=covariance(j,i)=ComplexType(normal(rng),normal(rng));
        }
        // Repeated blocks plus a locally retained mode removed globally.
        if(split && n>=9) {
            for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j) covariance(3+i,3+j)=covariance(i,j);
            for(size_t i=6;i<n;++i)for(size_t j=6;j<n;++j) covariance(i,j)=0.;
            covariance(n-1,n-1)=1e-13;
        }
        if(!split && n>1) covariance(0,1)+=ComplexType(1e-14,-1e-14);
        for(int objective=0;objective<3;++objective) {
            std::vector<double> weights(n);
            for(size_t i=0;i<n;++i) weights[i]=i<n/3?1.:i<2*n/3?2.:.5;
            std::unique_ptr<f::DenseComplexGaussianSampler> sampler;
            if(objective==0) sampler=std::make_unique<f::DenseComplexGaussianSampler>(covariance);
            if(objective==1) sampler=std::make_unique<f::DenseComplexGaussianSampler>(covariance,weights);
            if(objective==2) sampler=std::make_unique<f::DenseComplexGaussianSampler>(covariance,std::max(size_t{1},n/3),n/3,1e-8);
            const auto L=sampler->leading_factor_rows(n);
            emit(L); // Exact real-linear Gaussian map, including latent order.
            emit(f::ComplexDynamicMatrix(L*blaze::trans(L))); // E[z z^T]
            emit(f::ComplexDynamicMatrix(L*blaze::ctrans(L))); // E[z z^dagger]
            std::mt19937 engine(42);
            f::ComplexDynamicMatrix draws(n,19);
            for(size_t s=0;s<draws.columns();++s) {
                const auto field=sampler->draw(engine);
                for(size_t i=0;i<n;++i) draws(i,s)=field[i];
            }
            emit(draws);
        }
    }
    f::ComplexDynamicMatrix underflow(3,3,ComplexType{});
    underflow(0,1)=underflow(1,0)=.01; underflow(2,2)=1.;
    bool rejected=false;
    try {
        const double tiny=std::numeric_limits<double>::denorm_min();
        f::DenseComplexGaussianSampler sampler(underflow,std::vector<double>{tiny,tiny,1.});
    } catch(const std::runtime_error&) { rejected=true; }
    if(!rejected) return 1;
    std::cout<<1<<' '<<0<<'\n'; // Both versions must reject the lost physical coupling.

}
