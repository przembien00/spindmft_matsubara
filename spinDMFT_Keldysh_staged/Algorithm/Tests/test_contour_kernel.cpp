#include"../Contour/Contour_Kernel.h"

#include<algorithm>
#include<cmath>
#include<iostream>

namespace contour = spinDMFT::Contour;

namespace
{
int require( bool condition,const char* message )
{
    if( condition ) return 0;
    std::cerr<<"FAILED: "<<message<<'\n'; return 1;
}
bool close( ComplexType a,ComplexType b,RealType tolerance=RealType{1e-12} )
{ return std::abs(a-b)<tolerance; }
}

int main()
{
    int failures{};
    const contour::ContourLayout layout{4,3};
    failures+=require(layout.dimension()==30,"M,+,- layout dimension");
    for( size_t i=0;i<layout.dimension();++i )
    {
        const auto index=layout.decode(i);
        failures+=require(layout.flat(index.branch,index.point,index.component)==i,
                          "layout encode/decode round trip");
    }

    contour::CorrelationSet primitive{'D',4,3};
    for( size_t t=0;t<3;++t )
        for( size_t p=0;p<primitive.Re[t].size();++p )
        {
            const auto ab=primitive.Re[t].get_direction_pair(p);
            for( size_t tau=0;tau<4;++tau )
            {
                const ComplexType value=ComplexType{
                    RealType{1.}+RealType{0.2}*t+RealType{0.03}*tau
                      +RealType{0.01}*(ab[0]+ab[1]),
                    RealType{0.04}*(static_cast<int>(ab[0])-static_cast<int>(ab[1]))
                      +RealType{0.02}*t};
                primitive.Re[t][p][tau]=std::real(value);
                primitive.Im[t][p][tau]=std::imag(value);
            }
        }
    // Equal-time greater/lesser contact is symmetric, as required for a
    // classical pseudo-covariance at theta(0)=1/2.
    for( size_t a=0;a<3;++a ) for( size_t b=0;b<3;++b )
    {
        const ComplexType contact{RealType{0.5}+RealType{0.1}*(a+b),RealType{0.}};
        for( size_t p=0;p<primitive.Re[0].size();++p )
        {
            const auto ab=primitive.Re[0].get_direction_pair(p);
            if( ab[0]==a&&ab[1]==b )
            {
                primitive.Re[0][p][0]=primitive.Re[0][p][3]=std::real(contact);
                primitive.Im[0][p][0]=primitive.Im[0][p][3]=std::imag(contact);
            }
        }
    }

    using B=contour::Branch;
    const auto branch=[&](const contour::ContourIndex& first,
                          const contour::ContourIndex& second)
    { return contour::branch_correlation(primitive,layout,first,second); };
    const contour::ContourIndex plus_late{B::Forward,2,0};
    const contour::ContourIndex plus_early{B::Forward,0,1};
    const contour::ContourIndex minus_late{B::Backward,2,0};
    const contour::ContourIndex minus_early{B::Backward,0,1};
    const ComplexType greater=contour::greater_value(primitive,0,1,2);
    const ComplexType lesser=contour::lesser_value(primitive,0,1,2);
    failures+=require(close(branch(plus_late,plus_early),greater),
                      "++ is time ordered");
    failures+=require(close(branch(minus_late,minus_early),lesser),
                      "-- is anti-time ordered");
    failures+=require(close(branch(plus_late,minus_early),lesser),
                      "+- is lesser");
    failures+=require(close(branch(minus_late,plus_early),greater),
                      "-+ is greater");
    const contour::ContourIndex matsubara{B::Matsubara,1,1};
    failures+=require(close(branch(plus_late,matsubara),
                            contour::edge_value(primitive,2,0,1,1)),
                      "mixed branch uses the imaginary edge grid");

    const ComplexType cpp=branch(plus_late,plus_early);
    const ComplexType cpm=branch(plus_late,minus_early);
    const ComplexType cmp=branch(minus_late,plus_early);
    const ComplexType cmm=branch(minus_late,minus_early);
    failures+=require(close(RealType{0.25}*(cpp+cpm+cmp+cmm),
                            RealType{0.5}*(greater+lesser)),"eta-eta covariance");
    failures+=require(close(RealType{0.5}*(cpp-cpm+cmp-cmm),greater-lesser),
                      "eta-nu is the causal response block");
    failures+=require(close(cpp-cpm-cmp+cmm,ComplexType{}),"nu-nu pseudo-covariance vanishes");

    for( size_t i=0;i<layout.dimension();++i )
        for( size_t j=0;j<layout.dimension();++j )
            failures+=require(close(branch(layout.decode(i),layout.decode(j)),
                                    branch(layout.decode(j),layout.decode(i))),
                "composite contour transpose identity");

    contour::CorrelationSet contact{'D',4,3};
    const auto set_contact=[&](size_t tau,ComplexType value)
    {
        for( size_t p=0;p<contact.Re[0].size();++p )
            if( contact.Re[0].get_direction_pair(p)==std::array<size_t,2>{0,1} )
            {
                contact.Re[0][p][tau]=std::real(value);
                contact.Im[0][p][tau]=std::imag(value);
            }
    };
    const ComplexType contact_greater{RealType{2.},RealType{0.5}};
    const ComplexType contact_lesser{RealType{-1.},RealType{-0.25}};
    set_contact(3,contact_greater);
    set_contact(0,contact_lesser);
    const contour::ContourIndex plus_equal_a{B::Forward,1,0};
    const contour::ContourIndex plus_equal_b{B::Forward,1,1};
    const contour::ContourIndex minus_equal_a{B::Backward,1,0};
    const contour::ContourIndex minus_equal_b{B::Backward,1,1};
    const ComplexType contact_half=RealType{0.5}*(contact_greater+contact_lesser);
    const ComplexType contact_plus_plus=contour::branch_correlation(
        contact,layout,plus_equal_a,plus_equal_b);
    const ComplexType contact_plus_minus=contour::branch_correlation(
        contact,layout,plus_equal_a,minus_equal_b);
    const ComplexType contact_minus_plus=contour::branch_correlation(
        contact,layout,minus_equal_a,plus_equal_b);
    const ComplexType contact_minus_minus=contour::branch_correlation(
        contact,layout,minus_equal_a,minus_equal_b);
    failures+=require(close(contact_plus_plus,contact_half),
        "++ coincident time uses theta(0)=1/2");
    failures+=require(close(contour::branch_correlation(
        contact,layout,plus_equal_b,plus_equal_a),contact_half),
        "coincident transpose reuses the symmetric half-contact");
    failures+=require(close(contact_minus_minus,contact_half),
        "-- coincident time uses theta(0)=1/2");
    failures+=require(close(RealType{0.25}*(contact_plus_plus+contact_plus_minus
                                           +contact_minus_plus+contact_minus_minus),
                            contact_half),"coincident eta-eta is symmetrized");
    failures+=require(close(RealType{0.5}*(contact_plus_plus-contact_plus_minus
                                           +contact_minus_plus-contact_minus_minus),
                            RealType{0.5}*(contact_greater-contact_lesser)),
                      "coincident eta-nu keeps half the contact response");
    failures+=require(close(contact_plus_plus-contact_plus_minus
                            -contact_minus_plus+contact_minus_minus,ComplexType{}),
                      "coincident nu-nu pseudo-covariance vanishes");
    return failures==0?0:1;
}
