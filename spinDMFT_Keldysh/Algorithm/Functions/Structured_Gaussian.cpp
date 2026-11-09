#include"Structured_Gaussian.h"
#include"Gaussian_FFT.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<numeric>
#include<stdexcept>

namespace spinDMFT::Functions
{
namespace
{
using Block=StructuredContourCovariance::Block;
using Field=JointComplexGaussianSampler::FieldVector;
Block transpose(const Block& b)
{ Block t{};for(size_t a=0;a<3;++a)for(size_t c=0;c<3;++c)t[3*a+c]=b[3*c+a];return t; }

class Toeplitz3
{
 public:
    Toeplitz3(std::vector<Block> positive,std::vector<Block> negative)
        : m_positive(std::move(positive)),m_negative(std::move(negative))
    {
        m_n=m_positive.size();
        if(!m_n||m_negative.size()!=m_n)throw std::invalid_argument("invalid lag tables");
        m_nonzero=false;
        for(const auto* table:{&m_positive,&m_negative})for(const auto& b:*table)
            for(auto z:b)m_nonzero=m_nonzero||z!=ComplexType{};
        if(m_n<=16||!m_nonzero)return;
        m_length=gaussian_detail::fft_length(2*m_n-1);
        for(size_t ab=0;ab<9;++ab)
        {
            auto& f=m_fourier[ab];f.resize(m_length);
            f[0]=m_positive[0][ab];
            for(size_t lag=1;lag<m_n;++lag)
            { f[lag]=m_positive[lag][ab];f[m_length-lag]=m_negative[lag][ab]; }
            gaussian_detail::fft(f);
        }
    }
    void add(const Field& x,size_t input_offset,Field& y,size_t output_offset) const
    {
        if(!m_nonzero)return;
        if(!m_length)
        {
            for(size_t i=0;i<m_n;++i)for(size_t j=0;j<m_n;++j)
            {
                const auto& b=i>=j?m_positive[i-j]:m_negative[j-i];
                for(size_t a=0;a<3;++a)for(size_t c=0;c<3;++c)
                    y[output_offset+3*i+a]+=b[3*a+c]*x[input_offset+3*j+c];
            }
            return;
        }
        std::array<std::vector<ComplexType>,3> input;
        for(size_t c=0;c<3;++c)
        {
            input[c].resize(m_length);
            for(size_t j=0;j<m_n;++j)input[c][j]=x[input_offset+3*j+c];
            gaussian_detail::fft(input[c]);
        }
        std::vector<ComplexType> output(m_length);
        for(size_t a=0;a<3;++a)
        {
            std::fill(output.begin(),output.end(),ComplexType{});
            for(size_t c=0;c<3;++c)for(size_t f=0;f<m_length;++f)
                output[f]+=m_fourier[3*a+c][f]*input[c][f];
            gaussian_detail::fft(output,true);
            for(size_t i=0;i<m_n;++i)y[output_offset+3*i+a]+=output[i];
        }
    }
    // Sum of absolute values in each finite Toeplitz row, using prefix sums.
    std::vector<RealType> row_bounds() const
    {
        std::vector<RealType> rows(3*m_n);
        for(size_t a=0;a<3;++a)
        {
            std::vector<RealType> pos(m_n),neg(m_n);
            for(size_t lag=0;lag<m_n;++lag)
            {
                for(size_t b=0;b<3;++b)
                { pos[lag]+=std::abs(m_positive[lag][3*a+b]);neg[lag]+=std::abs(m_negative[lag][3*a+b]); }
                if(lag){pos[lag]+=pos[lag-1];neg[lag]+=neg[lag-1];}
            }
            for(size_t i=0;i<m_n;++i)rows[3*i+a]=pos[i]+neg[m_n-1-i]-neg[0];
        }
        return rows;
    }
    size_t stored() const
    { return 18*m_n+9*m_length; }
 private:
    size_t m_n{},m_length{};
    bool m_nonzero{};
    std::vector<Block> m_positive,m_negative;
    std::array<std::vector<ComplexType>,9> m_fourier;
};

std::vector<Block> transposed_lags(const std::vector<Block>& values)
{ auto result=values;for(auto& b:result)b=transpose(b);return result; }

// Exact Chebyshev series coefficients of sqrt(max(t,0)), t in [-1,1].
// Integral I_j=int_0^(pi/2) sqrt(cos(theta))*cos(j theta)dtheta obeys
// I_(j+2)/I_j=(1/2-j)/(5/2+j); the constant term has half normalization.
std::vector<RealType> coefficients(size_t degree)
{
    const long double pi=std::acos(-1.L);
    std::vector<long double> c(degree+1);
    c[0]=std::tgamma(.75L)/(2*std::sqrt(pi)*std::tgamma(1.25L));
    c[1]=std::tgamma(1.25L)/(std::sqrt(pi)*std::tgamma(1.75L));
    c[2]=.4L*c[0];
    for(size_t j=1;j+2<=degree;++j)c[j+2]=c[j]*(.5L-j)/(2.5L+j);
    return {c.begin(),c.end()};
}

// A uniform bound on |p(t)^2-max(t,0)|: compute the Chebyshev
// coefficients of p^2 by an alias-free Laurent convolution, subtract the
// known series of max(t,0), then sum absolute coefficients plus its EXACT
// remaining tail 1/[pi(2k+1)]. No sampled spectral grid is used as a bound.
RealType squared_function_bound(const std::vector<RealType>& c)
{
    const size_t k=c.size()-1,length=gaussian_detail::fft_length(4*k+1);
    std::vector<std::complex<long double>> a(length);
    a[0]=c[0];
    for(size_t j=1;j<=k;++j)a[j]=a[length-j]=static_cast<long double>(c[j])/2;
    gaussian_detail::fft(a);for(auto& x:a)x*=x;gaussian_detail::fft(a,true);
    const long double pi=std::acos(-1.L);
    long double bound=1/(pi*(2*k+1));
    for(size_t j=0;j<=2*k;++j)
    {
        long double target{};
        if(j==0)target=1/pi;
        else if(j==1)target=.5L;
        else if(j%2==0)
        {
            const long double l=j/2;
            target=2/pi*((j/2)%2?1:-1)/(4*l*l-1);
        }
        bound+=std::abs((j?2:1)*a[j].real()-target);
    }
    // Conservative roundoff allowance, separate from the analytic tail.
    bound+=1000*std::numeric_limits<long double>::epsilon()*length;
    bound+=100*std::numeric_limits<RealType>::epsilon()*k*k;
    return static_cast<RealType>(bound);
}
}

struct StructuredContourCovariance::Impl
{
    struct Stencil{size_t first,count;std::array<RealType,4> weights{};};
    size_t m,N,n,eta_offset,nu_offset;
    Toeplitz3 M,S,R,RT;
    std::vector<Block> mixed;
    std::vector<Stencil> interpolation;
    RealType bound{};
    std::vector<bool> active;
    Impl(size_t imaginary,size_t measured,size_t q,std::vector<Block> ml,
        std::vector<Block> sl,std::vector<Block> rl,std::vector<Block> b)
        : m(imaginary),N(measured),n(measured*q),eta_offset(3*m),nu_offset(3*(m+n)),
          M(ml,transposed_lags(ml)),S(sl,transposed_lags(sl)),
          R(rl,std::vector<Block>(n)),
          RT([&]{std::vector<Block> t(n);t[0]=transpose(rl[0]);return t;}(),transposed_lags(rl)),
          mixed(std::move(b)),active(3*(m+2*n+1))
    {
        if(!m||!N||!q||ml.size()!=m||sl.size()!=n||rl.size()!=n||mixed.size()!=(N+1)*m)
            throw std::invalid_argument("invalid compressed contour covariance dimensions");
        std::vector<RealType> counts(N+1);
        for(size_t j=0;j<n;++j)
        {
            const RealType x=(RealType(j)+RealType{.5})/RealType(q);
            const size_t interval=static_cast<size_t>(std::floor(x)),count=std::min(size_t{4},N+1);
            Stencil s{std::min(interval==0?0:interval-1,N+1-count),count,{}};
            for(size_t l=0;l<count;++l)
            {
                s.weights[l]=1;
                for(size_t k=0;k<count;++k)if(k!=l)
                    s.weights[l]*=(x-RealType(s.first+k))/(RealType(l)-RealType(k));
                counts[s.first+l]+=std::abs(s.weights[l]);
            }
            interpolation.push_back(s);
        }
        std::vector<RealType> rows(active.size()),coarse_rows(3*(N+1));
        const auto add_rows=[&](const Toeplitz3& t,size_t offset)
        {const auto v=t.row_bounds();for(size_t i=0;i<v.size();++i)rows[offset+i]+=v[i];};
        add_rows(M,0);add_rows(S,eta_offset);add_rows(R,eta_offset);add_rows(RT,nu_offset);
        for(size_t t=0;t<=N;++t)for(size_t k=0;k<m;++k)
            for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)
            {
                const RealType z=std::abs(mixed[t*m+k][3*a+b]);
                coarse_rows[3*t+a]+=z;rows[3*k+b]+=counts[t]*z;
            }
        for(size_t j=0;j<n;++j)
        {
            const auto& s=interpolation[j];
            for(size_t l=0;l<s.count;++l)for(size_t a=0;a<3;++a)
                rows[eta_offset+3*j+a]+=std::abs(s.weights[l])*coarse_rows[3*(s.first+l)+a];
        }
        for(size_t i=0;i<rows.size();++i)
        {
            if(!std::isfinite(rows[i]))throw std::invalid_argument("non-finite compressed covariance");
            active[i]=rows[i]>RealType{};bound=std::max(bound,rows[i]);
        }
        bound*=RealType{1}+RealType{100}*std::numeric_limits<RealType>::epsilon();
    }
};

StructuredContourCovariance::StructuredContourCovariance(size_t m,size_t N,size_t q,
    std::vector<Block> ml,std::vector<Block> sl,std::vector<Block> rl,std::vector<Block> b)
{
    const size_t limit=std::numeric_limits<size_t>::max()/3;
    if(!m||!N||!q||m>=limit||N>(limit-m-1)/2/q
        ||ml.size()!=m||sl.size()!=N*q||rl.size()!=N*q
        ||N==std::numeric_limits<size_t>::max()||m>std::numeric_limits<size_t>::max()/(N+1)
        ||b.size()!=(N+1)*m)
        throw std::invalid_argument("invalid compressed contour covariance dimensions");
    m_impl=std::make_unique<Impl>(m,N,q,std::move(ml),std::move(sl),std::move(rl),std::move(b));
}
StructuredContourCovariance::~StructuredContourCovariance()=default;
size_t StructuredContourCovariance::size() const {return m_impl->active.size();}
RealType StructuredContourCovariance::spectral_bound() const {return m_impl->bound;}
const std::vector<bool>& StructuredContourCovariance::active_rows() const {return m_impl->active;}
size_t StructuredContourCovariance::stored_complex_values() const
{const auto& p=*m_impl;return p.M.stored()+p.S.stored()+p.R.stored()+p.RT.stored()+9*p.mixed.size();}

StructuredContourCovariance::FieldVector StructuredContourCovariance::apply(const FieldVector& x) const
{
    const auto& p=*m_impl;
    if(x.size()!=size())throw std::invalid_argument("compressed covariance vector size mismatch");
    FieldVector y(size(),ComplexType{});
    p.M.add(x,0,y,0);p.S.add(x,p.eta_offset,y,p.eta_offset);
    p.R.add(x,p.nu_offset,y,p.eta_offset);p.RT.add(x,p.eta_offset,y,p.nu_offset);
    // B=E B0, B^T=B0^T E^T. Only the coarse mixed table is stored.
    FieldVector bx(3*(p.N+1),ComplexType{}),etx(bx.size(),ComplexType{});
    for(size_t j=0;j<p.n;++j)
    {
        const auto& s=p.interpolation[j];
        for(size_t l=0;l<s.count;++l)for(size_t a=0;a<3;++a)
            etx[3*(s.first+l)+a]+=s.weights[l]*x[p.eta_offset+3*j+a];
    }
    for(size_t t=0;t<=p.N;++t)for(size_t k=0;k<p.m;++k)
        for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)
        {
            const auto value=p.mixed[t*p.m+k][3*a+b];
            bx[3*t+a]+=value*x[3*k+b];y[3*k+b]+=value*etx[3*t+a];
        }
    for(size_t j=0;j<p.n;++j)
    {
        const auto& s=p.interpolation[j];
        for(size_t l=0;l<s.count;++l)for(size_t a=0;a<3;++a)
            y[p.eta_offset+3*j+a]+=s.weights[l]*bx[3*(s.first+l)+a];
    }
    // FFT roundoff must never populate an identically zero physical row.
    for(size_t i=0;i<size();++i)if(!p.active[i])y[i]=ComplexType{};
    return y;
}

ChebyshevComplexGaussianSampler::ChebyshevComplexGaussianSampler(
    std::shared_ptr<StructuredContourCovariance> covariance,RealType tolerance,size_t max_degree)
    : m_covariance(std::move(covariance))
{
    if(!m_covariance||!std::isfinite(tolerance)||tolerance<=RealType{}||tolerance>=RealType{1}
        ||max_degree<16||max_degree>1048576)
        throw std::invalid_argument("invalid Chebyshev sampling tolerance or maximum degree");
    if(spectral_bound()==RealType{})return;
    for(size_t degree=16;;degree=std::min(2*degree,max_degree))
    {
        auto c=coefficients(degree);
        const RealType bound=RealType{2}*squared_function_bound(c);
        if(bound<=tolerance)
        {
            m_coefficients=std::move(c);m_moment_bound=spectral_bound()*bound;break;
        }
        if(degree==max_degree)
            throw std::runtime_error("Chebyshev covariance error bound exceeds gaussianChebyshevTolerance at gaussianChebyshevMaxDegree; increase the maximum degree or use dense sampling");
    }
}

ChebyshevComplexGaussianSampler::LatentVector
ChebyshevComplexGaussianSampler::apply_scaled_lift(const LatentVector& x) const
{
    FieldVector u(size());
    for(size_t i=0;i<size();++i)u[i]={x[i],-x[size()+i]};
    const auto v=m_covariance->apply(u);
    LatentVector result(2*size());
    for(size_t i=0;i<size();++i)
    {result[i]=std::real(v[i])/spectral_bound();result[size()+i]=std::imag(v[i])/spectral_bound();}
    return result;
}

ChebyshevComplexGaussianSampler::LatentVector
ChebyshevComplexGaussianSampler::draw_latent(std::mt19937& engine)
{ LatentVector r(latent_dimension());for(auto& x:r)x=m_normal(engine);return r; }

ChebyshevComplexGaussianSampler::FieldVector
ChebyshevComplexGaussianSampler::field_from_latent(const LatentVector& latent)
{
    if(latent.size()!=latent_dimension())throw std::invalid_argument("Chebyshev latent dimension mismatch");
    FieldVector field(size(),ComplexType{});
    if(m_coefficients.empty())return field;
    // Forward Chebyshev recurrence; coefficients and degree never depend on r.
    LatentVector previous=latent,current=apply_scaled_lift(latent);
    LatentVector sum=m_coefficients[0]*previous+m_coefficients[1]*current;
    for(size_t j=2;j<m_coefficients.size();++j)
    {
        LatentVector next=RealType{2}*apply_scaled_lift(current)-previous;
        sum+=m_coefficients[j]*next;previous=std::move(current);current=std::move(next);
    }
    const RealType scale=std::sqrt(spectral_bound());
    for(size_t i=0;i<size();++i)if(m_covariance->active_rows()[i])
        field[i]=scale*ComplexType{sum[i],sum[size()+i]};
    return field;
}

ChebyshevComplexGaussianSampler::FieldVector
ChebyshevComplexGaussianSampler::draw(std::mt19937& engine)
{return field_from_latent(draw_latent(engine));}
}
