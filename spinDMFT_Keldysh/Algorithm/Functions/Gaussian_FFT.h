#pragma once

#include<algorithm>
#include<cmath>
#include<complex>
#include<vector>

namespace spinDMFT::Functions::gaussian_detail
{
// Radix-two FFT for linear convolutions. No periodic physical boundary is
// introduced: callers pad to at least the full linear-convolution length.
template<class Real>
void fft(std::vector<std::complex<Real>>& a,bool inverse=false)
{
    const size_t n=a.size();
    for(size_t i=1,j=0;i<n;++i)
    {
        size_t bit=n>>1;
        for(;j&bit;bit>>=1)j^=bit;
        j^=bit;
        if(i<j)std::swap(a[i],a[j]);
    }
    const Real pi=std::acos(Real{-1});
    for(size_t length=2;length<=n;length*=2)
    {
        const Real angle=(inverse?Real{2}:Real{-2})*pi/Real(length);
        const std::complex<Real> step{std::cos(angle),std::sin(angle)};
        for(size_t start=0;start<n;start+=length)
        {
            std::complex<Real> w{1};
            for(size_t j=0;j<length/2;++j)
            {
                const auto u=a[start+j],v=w*a[start+j+length/2];
                a[start+j]=u+v;a[start+j+length/2]=u-v;w*=step;
            }
        }
    }
    if(inverse)for(auto& x:a)x/=Real(n);
}

inline size_t fft_length(size_t minimum)
{ size_t n=1;while(n<minimum)n*=2;return n; }
}
