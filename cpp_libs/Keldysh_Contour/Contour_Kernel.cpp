#include"Contour_Kernel.h"

#include<algorithm>
#include<cmath>
#include<stdexcept>

namespace spinDMFT::Contour
{
namespace
{
void validate_primitive( const CorrelationSet& correlations )
{
    if( correlations.Re.empty() || correlations.Im.size()!=correlations.Re.size() )
        throw std::invalid_argument("inconsistent contour primitive");
}
}

ComplexType edge_value( const CorrelationSet& correlations, size_t real_time,
                        size_t first_spin, size_t second_spin,
                        size_t imaginary_edge )
{
    validate_primitive(correlations);
    if( real_time >= correlations.Re.size() )
        throw std::out_of_range("real-time primitive index out of range");
    const auto& re=correlations.Re[real_time](
        static_cast<uint>(first_spin),static_cast<uint>(second_spin));
    const auto& im=correlations.Im[real_time](
        static_cast<uint>(first_spin),static_cast<uint>(second_spin));
    if( imaginary_edge >= re.size() )
        throw std::out_of_range("imaginary edge index out of range");
    return {re[imaginary_edge],im[imaginary_edge]};
}

ComplexType greater_value( const CorrelationSet& correlations,
                           size_t a, size_t b, std::ptrdiff_t difference )
{
    const size_t beta=correlations.imaginary_points-1;
    if( difference >= 0 )
        return edge_value(correlations,static_cast<size_t>(difference),a,b,beta);
    return lesser_value(correlations,b,a,-difference);
}

ComplexType lesser_value( const CorrelationSet& correlations,
                          size_t a, size_t b, std::ptrdiff_t difference )
{
    if( difference >= 0 )
        return edge_value(correlations,static_cast<size_t>(difference),a,b,0);
    return greater_value(correlations,b,a,-difference);
}

ComplexType branch_correlation( const CorrelationSet& correlations,
                                const ContourLayout& layout,
                                const ContourIndex& first,
                                const ContourIndex& second )
{
    if( first.point>=layout.branch_size(first.branch)
        || second.point>=layout.branch_size(second.branch)
        || first.component>=3 || second.component>=3 )
        throw std::out_of_range("branch correlation index out of range");

    const bool first_M=first.branch==Branch::Matsubara;
    const bool second_M=second.branch==Branch::Matsubara;
    const auto equal_time_value=[&]()
    {
        // theta(0)=1/2: a same-branch contact is the symmetric half-sum of
        // G greater and G lesser.  Resolve it in canonical composite-index
        // order and reuse it for the transpose so E[V V^T] remains complex
        // symmetric even before Monte-Carlo noise has converged perfectly.
        const bool canonical=layout.flat(first.branch,first.point,first.component)
                           <=layout.flat(second.branch,second.point,second.component);
        const ContourIndex& left=canonical?first:second;
        const ContourIndex& right=canonical?second:first;
        if(first_M)
            return RealType{0.5}*(greater_value(correlations,left.component,right.component,0)
                                 +lesser_value(correlations,left.component,right.component,0));
        const RealType t=RealType(left.point);
        const ComplexType greater=real_time_value(
            correlations,true,t,t,left.component,right.component);
        const ComplexType lesser=real_time_value(
            correlations,false,t,t,left.component,right.component);
        return RealType{0.5}*(greater+lesser);
    };
    if( first_M && second_M )
    {
        if( first.point>second.point )
            return edge_value(correlations,0,second.component,first.component,
                              first.point-second.point);
        if( first.point<second.point )
            return edge_value(correlations,0,first.component,second.component,
                              second.point-first.point);
        return equal_time_value();
    }
    if( !first_M && second_M )
        return edge_value(correlations,first.point,first.component,
                          second.component,second.point);
    if( first_M && !second_M )
        return edge_value(correlations,second.point,second.component,
                          first.component,first.point);

    const std::ptrdiff_t difference=static_cast<std::ptrdiff_t>(first.point)
                                  -static_cast<std::ptrdiff_t>(second.point);
    const ComplexType greater=real_time_value(correlations,true,
        RealType(first.point),RealType(second.point),first.component,second.component);
    const ComplexType lesser=real_time_value(correlations,false,
        RealType(first.point),RealType(second.point),first.component,second.component);
    if( first.branch==Branch::Forward && second.branch==Branch::Forward )
    {
        if( difference>0 ) return greater;
        if( difference<0 ) return lesser;
        return equal_time_value();
    }
    if( first.branch==Branch::Backward && second.branch==Branch::Backward )
    {
        if( difference>0 ) return lesser;
        if( difference<0 ) return greater;
        return equal_time_value();
    }
    if( first.branch==Branch::Forward && second.branch==Branch::Backward )
        return lesser;
    return greater; // -+
}

ComplexType interpolated_edge_value(const CorrelationSet& correlations,
    RealType x,size_t a,size_t b,size_t tau)
{
    validate_primitive(correlations);
    const size_t n=correlations.Re.size();
    if(!std::isfinite(x)||x<RealType{}||x>RealType(n-1))
        throw std::out_of_range("staggered kernel time is outside the primitive grid");
    const size_t interval=static_cast<size_t>(std::floor(x));
    if(x==RealType(interval))return edge_value(correlations,interval,a,b,tau);
    const size_t count=std::min(size_t{4},n);
    const size_t first=std::min(interval==0?0:interval-1,n-count);
    ComplexType result{};
    for(size_t j=0;j<count;++j)
    {
        RealType weight{1.};
        for(size_t k=0;k<count;++k)if(k!=j)
            weight*=(x-RealType(first+k))/(RealType(j)-RealType(k));
        result+=weight*edge_value(correlations,first+j,a,b,tau);
    }
    return result;
}

ComplexType real_time_value(const CorrelationSet& correlations,bool greater,
    RealType t,RealType s,size_t a,size_t b)
{
    if(!correlations.full_two_time)
    {
        if(t<s)return real_time_value(correlations,!greater,s,t,b,a);
        return interpolated_edge_value(correlations,t-s,a,b,
            greater?correlations.imaginary_points-1:0);
    }
    const size_t n=correlations.Re.size();
    const auto interpolate=[&](RealType x,const auto& evaluate)
    {
        if(!std::isfinite(x)||x<RealType{}||x>RealType(n-1))
            throw std::out_of_range("two-time coordinate outside measurement grid");
        const size_t interval=static_cast<size_t>(std::floor(x));
        if(x==RealType(interval))return evaluate(interval);
        const size_t count=std::min(size_t{4},n);
        const size_t first=std::min(interval==0?0:interval-1,n-count);
        ComplexType result{};
        for(size_t j=0;j<count;++j)
        {
            RealType weight{1.};
            for(size_t k=0;k<count;++k)if(k!=j)
                weight*=(x-RealType(first+k))/(RealType(j)-RealType(k));
            result+=weight*evaluate(first+j);
        }
        return result;
    };
    return interpolate(t,[&](size_t i)
    {
        return interpolate(s,[&](size_t j)
            {return edge_value(correlations,i,a,b,correlations.real_offset(greater)+j);});
    });
}

ComplexType staggered_correlation(const CorrelationSet& correlations,
    const StaggeredIndex& first,const StaggeredIndex& second,size_t q)
{
    if(q==0)throw std::invalid_argument("staggered covariance needs positive substeps");
    using Sector=StaggeredSector;
    const auto f=first.sector,s=second.sector;
    if(f==Sector::Nu)
        return s==Sector::Eta?staggered_correlation(correlations,second,first,q):ComplexType{};
    if(s==Sector::Nu)
    {
        if(f!=Sector::Eta || second.point>first.point)return ComplexType{};
        const RealType t=correlations.full_two_time
            ?(RealType(first.point)+RealType{0.5})/RealType(q)
            :(RealType(first.point-second.point)+RealType{0.5})/RealType(q);
        const RealType second_time=correlations.full_two_time?RealType(second.point)/RealType(q):RealType{};
        return real_time_value(correlations,true,t,second_time,first.component,second.component)
              -real_time_value(correlations,false,t,second_time,first.component,second.component);
    }
    if(f==Sector::Matsubara&&s==Sector::Matsubara)
    {
        const ContourLayout layout{correlations.imaginary_points,correlations.Re.size()};
        return branch_correlation(correlations,layout,
            {Branch::Matsubara,first.point,first.component},
            {Branch::Matsubara,second.point,second.component});
    }
    if(f==Sector::Matsubara)return staggered_correlation(correlations,second,first,q);
    if(s==Sector::Matsubara)
        return interpolated_edge_value(correlations,(RealType(first.point)+RealType{0.5})/RealType(q),
            first.component,second.component,second.point);
    // N_ab(t,s)=1/2(C^>+C^<); negative lags exchange the spin indices.
    if(first.point<second.point || (first.point==second.point&&first.component>second.component))
        return staggered_correlation(correlations,second,first,q);
    const RealType t=correlations.full_two_time
        ?(RealType(first.point)+RealType{0.5})/RealType(q)
        :RealType(first.point-second.point)/RealType(q);
    const RealType second_time=correlations.full_two_time
        ?(RealType(second.point)+RealType{0.5})/RealType(q):RealType{};
    return RealType{0.5}*(real_time_value(correlations,true,t,second_time,first.component,second.component)
                        +real_time_value(correlations,false,t,second_time,first.component,second.component));
}

}
