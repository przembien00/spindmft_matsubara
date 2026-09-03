#include"Contour_Kernel.h"

#include<algorithm>
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
    const size_t beta=correlations.Re.front()[0].size()-1;
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
        const ComplexType greater=greater_value(
            correlations,left.component,right.component,0);
        const ComplexType lesser=lesser_value(
            correlations,left.component,right.component,0);
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
    const ComplexType greater=greater_value(
        correlations,first.component,second.component,difference);
    const ComplexType lesser=lesser_value(
        correlations,first.component,second.component,difference);
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

}
