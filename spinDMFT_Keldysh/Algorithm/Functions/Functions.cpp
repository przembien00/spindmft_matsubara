#include"Functions.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<stdexcept>

#include<Physics/CFET.h>
#include<Physics/Spin.h>
#include<Standard_Algorithms/Numerics.h>

namespace cfet = Physics::CFET;
namespace num = Standard_Algorithms::Numerics;
namespace sp = Physics::Spin;

namespace spinDMFT::Functions
{
namespace
{

Operator spin_scalar_product( const ComplexFieldVector& field )
{
    return field[0]*S_X+field[1]*S_Y+field[2]*S_Z;
}

ComplexFieldVector total_field( const ps::ParameterSpace& pspace,
                                const ComplexFieldVector& fluctuating,
                                const FieldVector& mean )
{
    ComplexFieldVector result=fluctuating;
    for( size_t c=0;c<3;++c )
        result[c]+=ComplexType{mean[c]+pspace.B.m_h[c],RealType{0.}};
    return result;
}

Operator cfet4_step( const ps::ParameterSpace& pspace,
                     const ComplexFieldVector& fluctuating_new,
                     const FieldVector& mean_new,
                     const ComplexFieldVector& fluctuating_old,
                     const FieldVector& mean_old,
                     const ComplexType& contour_step )
{
    const ComplexFieldVector field_new=total_field(
        pspace,fluctuating_new,mean_new);
    const ComplexFieldVector field_old=total_field(
        pspace,fluctuating_old,mean_old);
    std::array<Operator,3> exponentials{};
    std::transform(cfet::BETA_LIST.cbegin(),cfet::BETA_LIST.cend(),
                   exponentials.begin(),[&]( const cfet::BetaComponent& beta )
    {
        ComplexFieldVector weighted_field{};
        for( size_t c=0;c<3;++c )
            weighted_field[c]=beta[0]*field_new[c]+beta[1]*field_old[c];
        if( pspace.num_HilbertSpaceDimension==2 && H_REST_IS_SCALAR )
            return spin_half_field_exponential(
                weighted_field,(beta[0]+beta[1])*H_REST_SCALAR,contour_step);
        const Operator weighted_hamiltonian=
            spin_scalar_product(weighted_field)+(beta[0]+beta[1])*H_REST;
        return general_matrix_exponential(
            contour_step*weighted_hamiltonian);
    });
    return exponentials[0]*exponentials[1]*exponentials[2];
}

Operator gauss_cfet4_step( const ps::ParameterSpace& pspace,
                           const ComplexFieldVector& fluctuating_early,
                           const FieldVector& mean_early,
                           const ComplexFieldVector& fluctuating_late,
                           const FieldVector& mean_late,
                           const ComplexType& contour_step )
{
    // Two-node, two-exponential fourth-order commutator-free Magnus step.
    // The left exponential is weighted toward the later Gauss node, while
    // the right exponential is weighted toward the earlier node.
    const RealType root_three=std::sqrt(RealType{3.});
    const RealType a_early=(RealType{3.}-RealType{2.}*root_three)/RealType{12.};
    const RealType a_late =(RealType{3.}+RealType{2.}*root_three)/RealType{12.};
    const ComplexFieldVector field_early=total_field(
        pspace,fluctuating_early,mean_early);
    const ComplexFieldVector field_late=total_field(
        pspace,fluctuating_late,mean_late);

    const auto exponential=[&]( const RealType early_weight,
                                const RealType late_weight )
    {
        const ComplexFieldVector weighted_field=
            early_weight*field_early+late_weight*field_late;
        const RealType rest_weight=early_weight+late_weight;
        if( pspace.num_HilbertSpaceDimension==2 && H_REST_IS_SCALAR )
            return spin_half_field_exponential(
                weighted_field,rest_weight*H_REST_SCALAR,contour_step);
        const Operator weighted_hamiltonian=
            spin_scalar_product(weighted_field)+rest_weight*H_REST;
        return general_matrix_exponential(contour_step*weighted_hamiltonian);
    };
    return exponential(a_early,a_late)*exponential(a_late,a_early);
}

template<typename Vector>
Vector cubic_interval_value( const std::vector<Vector>& values,
                             const size_t interval,
                             const RealType fraction )
{
    if( values.size()<4 )
        throw std::invalid_argument(
            "Gauss-node CF4 interpolation needs at least four edge points");
    if( interval+1>=values.size() )
        throw std::out_of_range("CF4 interpolation interval is out of range");

    const size_t first=interval==0?0
        :interval+2>=values.size()?values.size()-4:interval-1;
    const RealType x=static_cast<RealType>(interval)+fraction;
    Vector result{};
    for( size_t j=0;j<4;++j )
    {
        const RealType xj=static_cast<RealType>(first+j);
        RealType weight{1.};
        for( size_t k=0;k<4;++k )
            if( k!=j )
            {
                const RealType xk=static_cast<RealType>(first+k);
                weight*=(x-xk)/(xj-xk);
            }
        result+=weight*values[first+j];
    }
    return result;
}

RealType complex_matrix_symmetry_error( const ComplexDynamicMatrix& matrix )
{
    RealType numerator{},denominator{};
    for( size_t i=0;i<matrix.rows();++i )
        for( size_t j=0;j<matrix.columns();++j )
        {
            numerator+=std::norm(matrix(i,j)-matrix(j,i));
            denominator+=std::norm(matrix(i,j));
        }
    return denominator>RealType{0.}?std::sqrt(numerator/denominator)
                                   :std::sqrt(numerator);
}

ComplexType trace_product( const Operator& left, const Operator& right )
{
    if( left.rows()!=left.columns() || right.rows()!=right.columns()
        ||left.rows()!=right.rows() )
        throw std::invalid_argument("trace product needs equal square matrices");
    ComplexType result{};
    for( size_t i=0;i<left.rows();++i )
        for( size_t j=0;j<left.columns();++j )
            result+=left(i,j)*right(j,i);
    return result;
}

size_t tensor_direction_index( const CorrTen& tensor, size_t a, size_t b )
{
    for( size_t p=0;p<tensor.size();++p )
    {
        const auto direction=tensor.get_direction_pair(p);
        if( direction[0]==a && direction[1]==b ) return p;
    }
    throw std::logic_error("requested spin direction is absent from the correlation symmetry");
}

RealType tensor_max_difference( const CorrTen& lhs, const CorrTen& rhs )
{
    if( lhs.size()!=rhs.size() )
        throw std::invalid_argument("correlation tensors have different symmetries");
    RealType largest{};
    for( size_t p=0;p<lhs.size();++p )
    {
        if( lhs[p].size()!=rhs[p].size() )
            throw std::invalid_argument("correlation tensors have different grids");
        for( size_t i=0;i<lhs[p].size();++i )
            largest=std::max(largest,std::abs(lhs[p][i]-rhs[p][i]));
    }
    return largest;
}

void mix_tensor( CorrTen& output, const CorrTen& old_values,
                 const CorrTen& raw_values, RealType alpha )
{
    for( size_t p=0;p<output.size();++p )
        for( size_t i=0;i<output[p].size();++i )
            output[p][i]=(RealType{1.}-alpha)*old_values[p][i]
                         +alpha*raw_values[p][i];
}


}

Operator general_matrix_exponential( const Operator& A )
{
    if( A.rows()!=A.columns() )
        throw std::invalid_argument("matrix exponential requires a square matrix");
    const size_t n=A.rows();
    if( n==0 ) return Operator{};

    RealType norm_one{};
    for( size_t j=0;j<n;++j )
    {
        RealType column_sum{};
        for( size_t i=0;i<n;++i ) column_sum+=std::abs(A(i,j));
        norm_one=std::max(norm_one,column_sum);
    }
    constexpr RealType theta13=RealType{5.371920351148152};
    const int squarings=norm_one>theta13
        ?std::max(0,static_cast<int>(std::ceil(std::log2(norm_one/theta13)))):0;
    const RealType scale=std::ldexp(RealType{1.},squarings);
    const Operator B=A/scale;
    const Operator B2=B*B,B4=B2*B2,B6=B4*B2;
    const Operator I=blaze::IdentityMatrix<ComplexType,blaze::rowMajor>(n);
    const std::array<RealType,14> c{
        RealType{64764752532480000.},RealType{32382376266240000.},
        RealType{7771770303897600.}, RealType{1187353796428800.},
        RealType{129060195264000.},  RealType{10559470521600.},
        RealType{670442572800.},     RealType{33522128640.},
        RealType{1323241920.},      RealType{40840800.},
        RealType{960960.},          RealType{16380.},
        RealType{182.},             RealType{1.}
    };
    const Operator U=B*(B6*(c[13]*B6+c[11]*B4+c[9]*B2)
                       +c[7]*B6+c[5]*B4+c[3]*B2+c[1]*I);
    const Operator V=B6*(c[12]*B6+c[10]*B4+c[8]*B2)
                    +c[6]*B6+c[4]*B4+c[2]*B2+c[0]*I;
    Operator result=blaze::inv(V-U)*(V+U);
    for( int s=0;s<squarings;++s ) result=result*result;
    return result;
}

Operator spin_half_field_exponential(
    const ComplexFieldVector& field, const ComplexType scalar,
    const ComplexType contour_step )
{
    // For S=sigma/2, write contour_step*H=a*I+w dot sigma.  Since
    // (w dot sigma)^2=(w dot w)I, with no complex conjugation, the exponential
    // is exact for arbitrary complex fields.
    const ComplexType half_step=RealType{0.5}*contour_step;
    const ComplexType wx=half_step*field[0];
    const ComplexType wy=half_step*field[1];
    const ComplexType wz=half_step*field[2];
    const ComplexType q_squared=wx*wx+wy*wy+wz*wz;
    const ComplexType q=std::sqrt(q_squared);

    ComplexType sinhc{};
    if( std::abs(q)<RealType{1e-4} )
    {
        // sinh(q)/q = 1 + q^2/3! + q^4/5! + q^6/7! + ... .
        const ComplexType q_fourth=q_squared*q_squared;
        sinhc=ComplexType{1.,0.}+q_squared/RealType{6.}
             +q_fourth/RealType{120.}
             +q_fourth*q_squared/RealType{5040.};
    }
    else sinhc=std::sinh(q)/q;

    const ComplexType prefactor=std::exp(contour_step*scalar);
    const ComplexType diagonal=std::cosh(q);
    const ComplexType imaginary_unit{0.,1.};
    Operator result(2,2,ComplexType{});
    result(0,0)=prefactor*(diagonal+sinhc*wz);
    result(0,1)=prefactor*sinhc*(wx-imaginary_unit*wy);
    result(1,0)=prefactor*sinhc*(wx+imaginary_unit*wy);
    result(1,1)=prefactor*(diagonal-sinhc*wz);
    return result;
}

void initialize_matrices( const ps::ParameterSpace& pspace )
{
    ZERO=blaze::ZeroMatrix<ComplexType,blaze::rowMajor>(
        pspace.num_HilbertSpaceDimension,pspace.num_HilbertSpaceDimension);
    IDENTITY=blaze::IdentityMatrix<ComplexType,blaze::rowMajor>(
        pspace.num_HilbertSpaceDimension);
    sp::write_spin_matrices(pspace.spin_float,S_X,S_Y,S_Z);
    H_REST=pspace.extra_interaction.m_term(ZERO,S_X,S_Y,S_Z);
    H_REST_SCALAR=ComplexType{};
    H_REST_IS_SCALAR=false;
    if( pspace.num_HilbertSpaceDimension==2 )
    {
        const ComplexType h00=H_REST(0,0);
        const ComplexType h01=H_REST(0,1);
        const ComplexType h10=H_REST(1,0);
        const ComplexType h11=H_REST(1,1);
        H_REST_SCALAR=RealType{0.5}*(h00+h11);
        const RealType scale=std::max(
            RealType{1.},std::max(std::abs(h00),std::abs(h11)));
        const RealType tolerance=RealType{100.}
            *std::numeric_limits<RealType>::epsilon()*scale;
        H_REST_IS_SCALAR=std::abs(h01)<=tolerance
            &&std::abs(h10)<=tolerance
            &&std::abs(h00-h11)<=tolerance;
    }
}

std::pair<MagTen,MagTen> generate_initial_magnetization(
    const ps::ParameterSpace& pspace )
{
    MagTen magnetization_Re{
        pspace.correlation_symmetry_type,pspace.num_RealTimePoints};
    MagTen magnetization_Im{
        pspace.correlation_symmetry_type,pspace.num_RealTimePoints};
    if( pspace.load_initial_spin_correlations )
    {
        if( pspace.initial_magnetization_linearized.size()<3 )
            throw std::invalid_argument(
                "imported initial magnetization needs x, y, and z components");
        FieldVector imported_full{};
        std::copy_n(pspace.initial_magnetization_linearized.cbegin(),
                    imported_full.size(),imported_full.begin());
        const MagVec imported{
            pspace.correlation_symmetry_type,imported_full};
        for( auto& magnetization : magnetization_Re ) magnetization=imported;
    }
    return {std::move(magnetization_Re),std::move(magnetization_Im)};
}

CorrelationSet generate_initial_correlations( const ps::ParameterSpace& pspace,
                                              const MagTen& magnetization_Re )
{
    if( magnetization_Re.empty()
        ||magnetization_Re.size()!=pspace.num_RealTimePoints
        ||magnetization_Re.get_symmetry()!=pspace.correlation_symmetry_type )
        throw std::invalid_argument(
            "initial correlations need a matching nonempty magnetization tensor");
    const FieldVector spin_expectation=magnetization_Re.front().expand();
    CorrelationSet result{pspace.correlation_symmetry_type,pspace.num_TimePoints,
                          pspace.num_RealTimePoints};
    CorrTen edge_Re{pspace.correlation_symmetry_type,pspace.num_TimePoints};
    CorrTen edge_Im{pspace.correlation_symmetry_type,pspace.num_TimePoints};

    if( !pspace.load_initial_spin_correlations )
    {
        const Corr diagonal{pspace.init_diag_corr.create_discretization(
            pspace.delta_t,pspace.num_TimePoints,pspace.spin_float)};
        const Corr nondiagonal{pspace.init_nondiag_corr.create_discretization(
            pspace.delta_t,pspace.num_TimePoints,pspace.spin_float)};
        edge_Re.iterate([&](Corr& values,const auto& direction)
        { values=direction[0]==direction[1]?diagonal:nondiagonal; });
    }
    else
    {
        for( size_t p=0;p<edge_Re.size();++p )
        {
            const size_t start=p*pspace.old_num_TimePoints;
            const size_t end=start+pspace.old_num_TimePoints;
            std::vector<RealType> re(pspace.initial_correlations_linearized.cbegin()+start,
                                     pspace.initial_correlations_linearized.cbegin()+end);
            std::vector<RealType> im(pspace.initial_correlations_imag_linearized.cbegin()+start,
                                     pspace.initial_correlations_imag_linearized.cbegin()+end);
            if( pspace.extrapolate_initial_spin_correlations )
            {
                re=num::extrapolate(re,pspace.num_TimePoints,
                                    pspace.old_delta_t,pspace.delta_t);
                im=num::extrapolate(im,pspace.num_TimePoints,
                                    pspace.old_delta_t,pspace.delta_t);
            }
            edge_Re[p]=Corr{std::move(re)};
            edge_Im[p]=Corr{std::move(im)};
        }
    }

    const RealType rate=std::max(std::abs(pspace.JQ),RealType{1e-6});
    for( size_t t_index=0;t_index<result.Re.size();++t_index )
    {
        const RealType t=static_cast<RealType>(t_index)*pspace.delta_real_t;
        const RealType envelope=std::exp(-RealType{0.5}*std::pow(rate*t,2));
        result.Re[t_index].iterate([&](Corr& values,const auto& direction)
        {
            const RealType product=spin_expectation[direction[0]]
                                  *spin_expectation[direction[1]];
            const size_t source=tensor_direction_index(
                edge_Re,direction[1],direction[0]);
            for( size_t tau=0;tau<values.size();++tau )
                values[tau]=product+(edge_Re[source][tau]-product)*envelope;
        });
        result.Im[t_index].iterate([&](Corr& values,const auto& direction)
        {
            const size_t source=tensor_direction_index(
                edge_Im,direction[1],direction[0]);
            for( size_t tau=0;tau<values.size();++tau )
                values[tau]=edge_Im[source][tau]*envelope;
        });
    }
    return result;
}

CorrelationSet connected_contour_primitive(
    const CorrelationSet& correlations,
    const MagTen& magnetization_Re,
    const MagTen& magnetization_Im )
{
    if( correlations.Re.empty() || correlations.Im.size()!=correlations.Re.size()
        ||magnetization_Re.size()!=correlations.Re.size()
        ||magnetization_Im.size()!=correlations.Re.size()
        ||magnetization_Re.get_symmetry()!=magnetization_Im.get_symmetry()
        ||magnetization_Re.get_symmetry()!=correlations.Re.front().get_symmetry()
        ||magnetization_Re.get_directions()!=magnetization_Im.get_directions() )
        throw std::invalid_argument(
            "connected contour primitive needs matching correlation and magnetization grids" );

    CorrelationSet connected=correlations;
    const auto full_Re=magnetization_Re.expand();
    const auto full_Im=magnetization_Im.expand();
    for( size_t t=0;t<connected.Re.size();++t )
    {
        if( connected.Im[t].size()!=connected.Re[t].size() )
            throw std::invalid_argument("inconsistent complex contour primitive");
        for( size_t p=0;p<connected.Re[t].size();++p )
        {
            const auto direction=connected.Re[t].get_direction_pair(p);
            if( connected.Im[t][p].size()!=connected.Re[t][p].size() )
                throw std::invalid_argument("inconsistent complex contour primitive grid");
            const ComplexType disconnected=
                ComplexType{full_Re[t][direction[0]],full_Im[t][direction[0]]}
               *ComplexType{full_Re.front()[direction[1]],full_Im.front()[direction[1]]};
            for( size_t tau=0;tau<connected.Re[t][p].size();++tau )
            {
                const ComplexType full{correlations.Re[t][p][tau],
                                       correlations.Im[t][p][tau]};
                const ComplexType value=full-disconnected;
                connected.Re[t][p][tau]=std::real(value);
                connected.Im[t][p][tau]=std::imag(value);
            }
        }
    }
    return connected;
}

namespace
{
SelfConsistentField covariance_from_primitive( const ps::ParameterSpace& pspace,
    CorrelationSet primitive,const bool prescribed,const bool materialize )
{
    const contour::ContourLayout layout{pspace.num_TimePoints,pspace.num_RealTimePoints};
    struct Term {size_t c,d;RealType coefficient;};
    std::array<std::vector<Term>,9> rotations;
    std::array<RealType,9> noise{};
    for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)
    {
        if(prescribed)rotations[3*a+b].push_back({a,b,RealType{1.}});
        else
        {
            const auto& D=pspace.spin_model.coupling_matrix;
            for(size_t c=0;c<3;++c)for(size_t d=0;d<3;++d)
            {
                const RealType coefficient=D(a,c)*D(b,d);
                if(coefficient!=RealType{})rotations[3*a+b].push_back({c,d,coefficient});
            }
            noise[3*a+b]=pspace.noise.m_variance_in(a,b);
        }
    }
    const RealType scale=prescribed?RealType{1.}:pspace.JQ*pspace.JQ;
    auto values=std::make_shared<CorrelationSet>(std::move(primitive));
    // Capture owned values and rotation coefficients, never a ParameterSpace
    // reference: the source remains valid independently of its caller.
    const auto raw=[values,layout,rotations,noise,scale](size_t row,size_t col)
    {
        const auto first=layout.decode(row),second=layout.decode(col);
        const size_t pair=3*first.component+second.component;
        ComplexType value{};
        for(const auto& term:rotations[pair])
            value+=term.coefficient*contour::branch_correlation(*values,layout,
                {first.branch,first.point,term.c},{second.branch,second.point,term.d});
        return scale*value+ComplexType{noise[pair],RealType{}};
    };
    SelfConsistentField result;
    const size_t n=layout.dimension();
    if(materialize)result.covariance.resize(n,n,false);
    RealType difference{},norm{};
    // Diagnose the original unsymmetrized kernel while filling only the
    // canonical triangle. No second dense raw covariance is allocated.
    for(size_t i=0;i<n;++i)for(size_t j=i;j<n;++j)
    {
        const ComplexType upper=raw(i,j);
        if(i==j)norm+=std::norm(upper);
        else
        {
            const ComplexType lower=raw(j,i);
            norm+=std::norm(upper)+std::norm(lower);
            difference+=RealType{2.}*std::norm(upper-lower);
        }
        if(materialize){result.covariance(i,j)=upper;result.covariance(j,i)=upper;}
    }
    result.branch_identity_error=norm>RealType{}?std::sqrt(difference/norm):std::sqrt(difference);
    result.covariance_symmetry_error=RealType{};
    result.covariance_source=CovarianceSource(n,[raw](size_t i,size_t j)
    {return i<=j?raw(i,j):raw(j,i);});
    return result;
}
}

SelfConsistentField self_consistent_equations(
    const ps::ParameterSpace& pspace, const CorrelationSet& correlations,
    const MagTen& magnetization_Re, const MagTen& magnetization_Im,
    const bool materialize_covariance )
{
    SelfConsistentField result=covariance_from_primitive(pspace,
        connected_contour_primitive(correlations,magnetization_Re,magnetization_Im),
        false,materialize_covariance);
    const auto physical_magnetization=magnetization_Re.expand();
    result.mean_time.resize(magnetization_Re.size());
    for(size_t t=0;t<magnetization_Re.size();++t)
    {
        result.mean_time[t]=pspace.JL*(
            pspace.spin_model.coupling_matrix*physical_magnetization[t]);
    }
    return result;
}

CorrelationSet harmonic_bath_primitive( const ps::ParameterSpace& pspace )
{
    if( !pspace.uses_harmonic_bath() )
        throw std::invalid_argument(
            "harmonic_bath_primitive requires bath=harmonic");
    if( pspace.beta<=RealType{0.} || pspace.bath_frequency<=RealType{0.}
        || pspace.bath_coupling<RealType{0.} )
        throw std::invalid_argument("invalid prescribed harmonic-bath parameters");
    const std::string axes="xyz";
    const size_t component=axes.find(pspace.bath_component);
    if( component==std::string::npos )
        throw std::invalid_argument("invalid harmonic-bath spin component");

    CorrelationSet result{pspace.correlation_symmetry_type,
                          pspace.num_TimePoints,pspace.num_RealTimePoints};
    const RealType omega=pspace.bath_frequency;
    const RealType n_plus_one=-RealType{1.}/std::expm1(-pspace.beta*omega);
    const RealType amplitude=pspace.bath_coupling*pspace.bath_coupling*n_plus_one;
    for( size_t t_index=0;t_index<pspace.num_RealTimePoints;++t_index )
    {
        const RealType t=static_cast<RealType>(t_index)*pspace.delta_real_t;
        const ComplexType positive_phase=std::exp(ComplexType{0.,omega*t});
        for( size_t p=0;p<result.Re[t_index].size();++p )
        {
            const auto direction=result.Re[t_index].get_direction_pair(p);
            if( direction[0]!=component || direction[1]!=component ) continue;
            for( size_t tau_index=0;tau_index<pspace.num_TimePoints;++tau_index )
            {
                const RealType tau=static_cast<RealType>(tau_index)*pspace.delta_t;
                const ComplexType value=amplitude*(
                    std::exp(-omega*tau)*positive_phase
                   +std::exp(-omega*(pspace.beta-tau))*std::conj(positive_phase));
                result.Re[t_index][p][tau_index]=std::real(value);
                result.Im[t_index][p][tau_index]=std::imag(value);
            }
        }
    }
    return result;
}

SelfConsistentField prescribed_harmonic_bath_field(
    const ps::ParameterSpace& pspace,const bool materialize_covariance )
{
    auto result=covariance_from_primitive(pspace,harmonic_bath_primitive(pspace),true,materialize_covariance);
    result.mean_time.assign(pspace.num_RealTimePoints,FieldVector{});
    return result;
}

ContourTrajectory compute_contour_trajectory(
    const ps::ParameterSpace& pspace,
    const DenseComplexGaussianSampler::FieldVector& joint_field,
    const MeanFieldTrajectory& mean_field_time )
{
    JointComplexGaussianSampler::ContourFieldSample field_sample{};
    field_sample.edge_field=joint_field;
    return compute_contour_trajectory(pspace,field_sample,mean_field_time);
}

ContourTrajectory compute_contour_trajectory(
    const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field_sample,
    const MeanFieldTrajectory& mean_field_time )
{
    ContourTrajectory result;
    TrajectoryWorkspace workspace;
    build_contour_trajectory(pspace,field_sample,mean_field_time,result,workspace);
    return result;
}

void build_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field_sample,
    const MeanFieldTrajectory& mean_field_time, ContourTrajectory& result,
    TrajectoryWorkspace& workspace, const bool imaginary_only )
{
    const contour::ContourLayout layout{
        pspace.num_TimePoints,pspace.num_RealTimePoints};
    const auto& joint_field=field_sample.edge_field;
    if( joint_field.size()!=layout.dimension() )
        throw std::invalid_argument("joint Keldysh field has the wrong dimension");
    if( mean_field_time.size()!=pspace.num_RealTimePoints )
        throw std::invalid_argument("mean-field trajectory has the wrong real-time grid");
    const size_t q=pspace.real_time_substeps;
    if(pspace.num_RealTimeSteps==0||pspace.num_RealTimeSteps>=std::numeric_limits<size_t>::max()/6
        ||q>(std::numeric_limits<size_t>::max()/6-1)/pspace.num_RealTimeSteps)
        throw std::invalid_argument("invalid real-time substep grid");
    if(q>1)
    {
        if(pspace.gaussian_factorization!="fft"&&pspace.gaussian_factorization!="dense"
            &&pspace.gaussian_factorization!="weighted-dense")
            throw std::invalid_argument("real-time substeps require dense, weighted-dense, or FFT sampling");
        if(pspace.gaussian_factorization=="fft"&&!field_sample.has_real_gauss_fields())
            throw std::invalid_argument("FFT substepping requires sampled Gauss-node fields");
    }
    if( pspace.uses_cf4()
        &&(field_sample.real_gauss_fields[0].size()!=0
           ||field_sample.real_gauss_fields[1].size()!=0) )
    {
        const size_t expected=6*pspace.num_RealTimeSteps*q;
        if( field_sample.real_gauss_fields[0].size()!=expected
            ||field_sample.real_gauss_fields[1].size()!=expected )
            throw std::invalid_argument(
                "CF4 trajectory received malformed Gauss-node real-field grids");
    }

    auto& imaginary_fields=workspace.imaginary_fields;
    auto& forward_fields=workspace.forward_fields;
    auto& backward_fields=workspace.backward_fields;
    imaginary_fields.resize(pspace.num_TimePoints);
    forward_fields.resize(pspace.num_RealTimePoints);
    backward_fields.resize(pspace.num_RealTimePoints);
    for( size_t k=0;k<pspace.num_TimePoints;++k )
        for( size_t c=0;c<3;++c )
            imaginary_fields[k][c]=joint_field[layout.flat(contour::Branch::Matsubara,k,c)];
    for( size_t t=0;t<forward_fields.size();++t )
        for( size_t c=0;c<3;++c )
        {
            forward_fields[t][c]=joint_field[layout.flat(contour::Branch::Forward,t,c)];
            backward_fields[t][c]=joint_field[layout.flat(contour::Branch::Backward,t,c)];
        }

    auto& imaginary_steps=workspace.imaginary_steps;
    imaginary_steps.resize(pspace.num_TimeSteps);
    for( size_t k=0;k<pspace.num_TimeSteps;++k )
    {
        if( pspace.uses_cf4() )
        {
            const RealType offset=std::sqrt(RealType{3.})/RealType{6.};
            const auto early=cubic_interval_value(
                imaginary_fields,k,RealType{0.5}-offset);
            const auto late=cubic_interval_value(
                imaginary_fields,k,RealType{0.5}+offset);
            imaginary_steps[k]=gauss_cfet4_step(
                pspace,early,mean_field_time.front(),
                late,mean_field_time.front(),
                ComplexType{-pspace.delta_t,RealType{0.}});
        }
        else
            imaginary_steps[k]=cfet4_step(
                pspace,imaginary_fields[k+1],mean_field_time.front(),
                imaginary_fields[k],mean_field_time.front(),
                ComplexType{-pspace.delta_t,RealType{0.}});
    }
    auto& prefix=workspace.prefix;
    prefix.resize(pspace.num_TimeSteps+1); prefix.front()=IDENTITY;
    for( size_t k=0;k<pspace.num_TimeSteps;++k )
        prefix[k+1]=imaginary_steps[k]*prefix[k];
    result.imaginary_density_operator=prefix.back();
    result.partition_function=blaze::trace(result.imaginary_density_operator);

    if( !imaginary_only )
        complete_contour_trajectory(pspace,field_sample,mean_field_time,result,workspace);
}

void complete_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field_sample,
    const MeanFieldTrajectory& mean_field_time, ContourTrajectory& result,
    TrajectoryWorkspace& workspace )
{
    const auto& prefix=workspace.prefix;
    const auto& imaginary_steps=workspace.imaginary_steps;
    const auto& forward_fields=workspace.forward_fields;
    const auto& backward_fields=workspace.backward_fields;
    auto& suffix=workspace.suffix;
    suffix.resize(pspace.num_TimeSteps+1); suffix.back()=IDENTITY;
    for( size_t k=pspace.num_TimeSteps;k-- >0; )
        suffix[k]=suffix[k+1]*imaginary_steps[k];
    const std::array<const Observable*,3> spins{&S_X,&S_Y,&S_Z};
    for( size_t c=0;c<3;++c )
    {
        result.imaginary_edge_insertions[c].resize(pspace.num_TimePoints);
        for( size_t tau=0;tau<pspace.num_TimePoints;++tau )
            result.imaginary_edge_insertions[c][tau]=
                suffix[tau]*(*spins[c])*prefix[tau];
    }

    result.forward_steps.resize(pspace.num_RealTimePoints);
    result.backward_steps.resize(pspace.num_RealTimePoints);
    result.forward_steps.front()=IDENTITY; result.backward_steps.front()=IDENTITY;
    const size_t q=pspace.real_time_steps_per_interval();
    const RealType h=pspace.delta_real_t/static_cast<RealType>(q);
    for( size_t t=1;t<pspace.num_RealTimePoints;++t )
    {
        const size_t interval=t-1;
        for(size_t j=0;j<q;++j)
        {
            const size_t micro=interval*q+j;
            Operator forward_step,backward_step;
            if( pspace.uses_cf4() )
            {
                const RealType offset=std::sqrt(RealType{3.})/RealType{6.};
                const RealType early=(static_cast<RealType>(j)+RealType{0.5}-offset)
                                    /static_cast<RealType>(q);
                const RealType late=(static_cast<RealType>(j)+RealType{0.5}+offset)
                                   /static_cast<RealType>(q);
                const FieldVector mean_early=cubic_interval_value(mean_field_time,interval,early);
                const FieldVector mean_late=cubic_interval_value(mean_field_time,interval,late);
                ComplexFieldVector forward_early{},forward_late{};
                ComplexFieldVector backward_early{},backward_late{};
                if( field_sample.has_real_gauss_fields() )
                {
                    for( size_t c=0;c<3;++c )
                    {
                        forward_early[c]=field_sample.real_gauss_fields[0][6*micro+c];
                        forward_late[c]=field_sample.real_gauss_fields[1][6*micro+c];
                        backward_early[c]=field_sample.real_gauss_fields[0][6*micro+3+c];
                        backward_late[c]=field_sample.real_gauss_fields[1][6*micro+3+c];
                    }
                }
                else
                {
                    // Evaluate the same native-edge cubic at each substep's nodes;
                    // refining propagation does not introduce new random fields.
                    forward_early=cubic_interval_value(forward_fields,interval,early);
                    forward_late=cubic_interval_value(forward_fields,interval,late);
                    backward_early=cubic_interval_value(backward_fields,interval,early);
                    backward_late=cubic_interval_value(backward_fields,interval,late);
                }
                forward_step=gauss_cfet4_step(
                    pspace,forward_early,mean_early,forward_late,mean_late,
                    ComplexType{RealType{0.},-h});
                // Reverse both nodes and contour sign on the backward branch.
                backward_step=gauss_cfet4_step(
                    pspace,backward_late,mean_late,backward_early,mean_early,
                    ComplexType{RealType{0.},+h});
            }
            else
            {
                forward_step=cfet4_step(
                    pspace,forward_fields[t],mean_field_time[t],
                    forward_fields[t-1],mean_field_time[t-1],
                    ComplexType{RealType{0.},-h});
                backward_step=cfet4_step(
                    pspace,backward_fields[t-1],mean_field_time[t-1],
                    backward_fields[t],mean_field_time[t],
                    ComplexType{RealType{0.},+h});
            }
            // Store only the composed native interval; measurements stay on its edges.
            if(j==0)
            {
                result.forward_steps[t]=forward_step;
                result.backward_steps[t]=backward_step;
            }
            else
            {
                result.forward_steps[t]=forward_step*result.forward_steps[t];
                result.backward_steps[t]=result.backward_steps[t]*backward_step;
            }
        }
    }
    Operator final_density=result.imaginary_density_operator;
    for( size_t t=1;t<pspace.num_RealTimePoints;++t )
        final_density=result.forward_steps[t]*final_density
                    *result.backward_steps[t];
    result.final_closed_contour_trace=blaze::trace(final_density);
}

namespace
{
template<typename Matrix>
void measure_observables_impl( const rtd::RunTimeData& rtdata,
    const ContourTrajectory& trajectory, rtd::MeasuredSample& sample,
    MeasurementWorkspace& workspace, std::vector<Matrix>& left,
    const bool prefix_insertion )
{
    const size_t nt=rtdata.num_real_time_points();
    const size_t ni=rtdata.num_imaginary_edge_points();
    const size_t nc=rtdata.num_correlation_components();
    const size_t nm=rtdata.num_magnetization_components();
    sample.partition=trajectory.partition_function;
    sample.correlations.resize(nt*nc*ni);
    sample.magnetization.resize(nt*nm);
    sample.closure.resize(nt);
    std::array<bool,3> measured{},inserted{};
    for( size_t p=0;p<nc;++p )
    {
        const auto direction=rtdata.correlation_direction(p);
        measured[direction[0]]=true; inserted[direction[1]]=true;
    }
    for( size_t c=0;c<nm;++c ) measured[rtdata.magnetization_direction(c)]=true;
    if constexpr( std::is_same_v<Matrix,SpinHalfOperator> )
        for( size_t c=0;c<3;++c ) if( inserted[c] )
        {
            auto& packed=workspace.spin_half_insertions[c]; packed.resize(ni);
            for( size_t tau=0;tau<ni;++tau )
            {
                const auto& op=trajectory.imaginary_edge_insertions[c][tau];
                packed[tau]={op(0,0),op(0,1),op(1,0),op(1,1)};
            }
        }
    const Matrix identity(IDENTITY);
    const std::array<Matrix,3> spins{Matrix(S_X),Matrix(S_Y),Matrix(S_Z)};
    if( !prefix_insertion )
    {
        Matrix backward_total=identity;
        for( size_t t=1;t<nt;++t )
            backward_total=backward_total*Matrix(trajectory.backward_steps[t]);
        left.resize(nt); left.back()=backward_total;
        for( size_t t=nt-1;t>0;--t )
            left[t-1]=left[t]*Matrix(trajectory.forward_steps[t]);
    }
    const Matrix rho(trajectory.imaginary_density_operator);
    Matrix density=rho,forward=identity,backward_prefix=identity;
    std::array<Matrix,3> measured_spins{};
    const auto contract=[](const auto& a,const auto& b)
    {
        if constexpr( std::is_same_v<Matrix,SpinHalfOperator> )
            return a(0,0)*b(0,0)+a(0,1)*b(1,0)
                  +a(1,0)*b(0,1)+a(1,1)*b(1,1);
        else return trace_product(a,b);
    };
    for( size_t t=0;t<nt;++t )
    {
        if( t>0 )
        {
            const Matrix u(trajectory.forward_steps[t]),b(trajectory.backward_steps[t]);
            forward=u*forward;
            if( prefix_insertion ) backward_prefix=backward_prefix*b;
            density=u*density*b;
        }
        sample.closure[t]=blaze::trace(density);
        for( size_t c=0;c<3;++c ) if( measured[c] )
            measured_spins[c]=(prefix_insertion?backward_prefix:left[t])*spins[c]*forward;
        for( size_t c=0;c<nm;++c )
            sample.magnetization[t*nm+c]=contract(rho,measured_spins[rtdata.magnetization_direction(c)]);
        for( size_t p=0;p<nc;++p )
        {
            const auto direction=rtdata.correlation_direction(p);
            const Matrix& b=measured_spins[direction[0]];
            ComplexType* output=sample.correlations.data()+(t*nc+p)*ni;
            if constexpr( std::is_same_v<Matrix,SpinHalfOperator> )
            {
                const auto& packed=workspace.spin_half_insertions[direction[1]];
                const ComplexType b00=b(0,0),b10=b(1,0),b01=b(0,1),b11=b(1,1);
                for( size_t tau=0;tau<ni;++tau )
                {
                    const auto& a=packed[tau];
                    output[tau]=a[0]*b00+a[1]*b10+a[2]*b01+a[3]*b11;
                }
            }
            else
                for( size_t tau=0;tau<ni;++tau )
                    output[tau]=contract(trajectory.imaginary_edge_insertions[direction[1]][tau],b);
        }
    }
}
}

void measure_contour_observables( const rtd::RunTimeData& layout,
    const ContourTrajectory& trajectory, rtd::MeasuredSample& sample,
    MeasurementWorkspace& workspace, const std::string& insertion_strategy )
{
    const size_t nt=layout.num_real_time_points();
    if( nt==0||trajectory.forward_steps.size()!=nt||trajectory.backward_steps.size()!=nt )
        throw std::invalid_argument("trajectory and measurement real-time grids differ");
    if( insertion_strategy!="closed-contour"&&insertion_strategy!="prefix" )
        throw std::invalid_argument("unknown spin insertion strategy");
    for( const auto& insertions:trajectory.imaginary_edge_insertions )
        if( insertions.size()!=layout.num_imaginary_edge_points() )
            throw std::invalid_argument("trajectory and measurement imaginary-time grids differ");
    const bool prefix=insertion_strategy=="prefix";
    if( trajectory.imaginary_density_operator.rows()==2 )
        measure_observables_impl(layout,trajectory,sample,workspace,workspace.spin_half_left,prefix);
    else
        measure_observables_impl(layout,trajectory,sample,workspace,workspace.general_left,prefix);
}

void compute_contour_correlations( rtd::RunTimeData& rtdata,
    const ContourTrajectory& trajectory, const RealType observable_normalization,
    const std::string& insertion_strategy )
{
    rtd::MeasuredSample sample;
    MeasurementWorkspace workspace;
    measure_contour_observables(rtdata,trajectory,sample,workspace,insertion_strategy);
    rtdata.accumulate_sample(sample,observable_normalization);
}

CorrTen imaginary_time_slice( const contour::ContourCorrelation& correlations )
{
    const char symmetry=correlations.front().get_symmetry();
    CorrTen result{symmetry,correlations.front()[0].size()};
    result.iterate([&](Corr& values,const auto& direction)
    {
        values=correlations.front()(
            static_cast<uint>(direction[1]),static_cast<uint>(direction[0]));
    });
    return result;
}

CorrTen real_time_slice( const contour::ContourCorrelation& correlations )
{
    const char symmetry=correlations.front().get_symmetry();
    const size_t beta_index=correlations.front()[0].size()-1;
    CorrTen result{symmetry,correlations.size()};
    result.iterate([&](Corr& values,const auto& direction)
    {
        for( size_t t=0;t<correlations.size();++t )
            values[t]=correlations[t](static_cast<uint>(direction[0]),
                                      static_cast<uint>(direction[1]))[beta_index];
    });
    return result;
}

RealType max_contour_difference( const CorrelationSet& old_values,
                                 const CorrelationSet& new_values )
{
    RealType largest{};
    auto compare=[&](const contour::ContourCorrelation& lhs,
                     const contour::ContourCorrelation& rhs)
    {
        if( lhs.size()!=rhs.size() )
            throw std::invalid_argument("contour primitives have different real-time grids");
        for( size_t t=0;t<lhs.size();++t )
            largest=std::max(largest,tensor_max_difference(lhs[t],rhs[t]));
    };
    compare(old_values.Re,new_values.Re);
    compare(old_values.Im,new_values.Im);
    return largest;
}

CorrelationSet mix_correlations( const CorrelationSet& old_values,
                                 const CorrelationSet& raw_values,
                                 RealType alpha )
{
    if( alpha<=RealType{0.} || alpha>RealType{1.} )
        throw std::invalid_argument("fixed-point mixing alpha must lie in (0,1]");
    CorrelationSet result=raw_values;
    for( size_t t=0;t<result.Re.size();++t )
    {
        mix_tensor(result.Re[t],old_values.Re[t],raw_values.Re[t],alpha);
        mix_tensor(result.Im[t],old_values.Im[t],raw_values.Im[t],alpha);
    }
    return result;
}

MagTen mix_magnetization_tensor( const MagTen& old_values,
                                 const MagTen& raw_values,
                                 RealType alpha )
{
    if( alpha<=RealType{0.} || alpha>RealType{1.} )
        throw std::invalid_argument("fixed-point mixing alpha must lie in (0,1]");
    if( old_values.size()!=raw_values.size()
        ||old_values.get_symmetry()!=raw_values.get_symmetry()
        ||old_values.get_directions()!=raw_values.get_directions() )
        throw std::invalid_argument("magnetization trajectories have different grids");
    MagTen result=(RealType{1.}-alpha)*old_values;
    result+=alpha*raw_values;
    return result;
}

MagTen project_constant_magnetization( const MagTen& values )
{
    if( values.empty() )
        throw std::invalid_argument("cannot project an empty magnetization trajectory");
    MagTen result=values;
    for( size_t t=1;t<result.size();++t ) result[t]=result.front();
    return result;
}

IterationResidual iteration_residual(
    const CorrelationSet& old_correlations,
    const CorrelationSet& raw_correlations,
    const CorrelationSet& standard_errors,
    const MagTen& old_magnetization_Re,
    const MagTen& old_magnetization_Im,
    const MagTen& raw_magnetization_Re,
    const MagTen& raw_magnetization_Im,
    const MagTen& magnetization_Re_errors,
    const MagTen& magnetization_Im_errors )
{
    if( old_correlations.Re.size()!=raw_correlations.Re.size()
        ||old_correlations.Im.size()!=raw_correlations.Im.size()
        ||standard_errors.Re.size()!=raw_correlations.Re.size()
        ||standard_errors.Im.size()!=raw_correlations.Im.size()
        ||old_magnetization_Re.size()!=raw_magnetization_Re.size()
        ||old_magnetization_Im.size()!=raw_magnetization_Im.size()
        ||raw_magnetization_Im.size()!=raw_magnetization_Re.size()
        ||magnetization_Re_errors.size()!=raw_magnetization_Re.size()
        ||magnetization_Im_errors.size()!=raw_magnetization_Re.size() )
        throw std::invalid_argument("iteration residual inputs have different grids");

    IterationResidual result{};
    auto update=[&]( const RealType difference, const RealType error,
                     const RealType zero_error_roundoff_scale=RealType{} )
    {
        if( !std::isfinite(difference)||!std::isfinite(error) )
        {
            result.absolute=std::numeric_limits<RealType>::infinity();
            result.standardized=std::numeric_limits<RealType>::infinity();
            return;
        }
        result.absolute=std::max(result.absolute,std::abs(difference));
        if( error<RealType{0.} )
            throw std::invalid_argument("iteration residual standard errors must be non-negative");
        if( error==RealType{0.} )
        {
            const RealType roundoff_limit=RealType{64.}
                *std::numeric_limits<RealType>::epsilon()
                *zero_error_roundoff_scale;
            if( difference!=RealType{0.}
                &&!(zero_error_roundoff_scale>RealType{}
                    &&std::abs(difference)<=roundoff_limit) )
                result.standardized=std::numeric_limits<RealType>::infinity();
            return;
        }
        result.standardized=std::max(result.standardized,std::abs(difference)/error);
    };
    auto update_contour=[&]( const contour::ContourCorrelation& old_values,
                             const contour::ContourCorrelation& raw_values,
                             const contour::ContourCorrelation& errors )
    {
        for( size_t t=0;t<raw_values.size();++t )
        {
            if( old_values[t].size()!=raw_values[t].size()
                ||errors[t].size()!=raw_values[t].size() )
                throw std::invalid_argument("iteration residual inputs have different symmetries");
            for( size_t p=0;p<raw_values[t].size();++p )
            {
                if( old_values[t][p].size()!=raw_values[t][p].size()
                    ||errors[t][p].size()!=raw_values[t][p].size() )
                    throw std::invalid_argument("iteration residual inputs have different grids");
                for( size_t tau=0;tau<raw_values[t][p].size();++tau )
                {
                    // At t=0 the tau=0 and tau=beta values can be exact spin
                    // identities with genuinely zero sampling variance.  Do
                    // not turn last-bit ratio arithmetic at those two
                    // endpoints into an infinite standardized residual.
                    const bool exact_endpoint=t==0
                        &&(tau==0||tau+1==raw_values[t][p].size());
                    const RealType roundoff_scale=exact_endpoint
                        ?std::max({RealType{1.},
                                   std::abs(old_values[t][p][tau]),
                                   std::abs(raw_values[t][p][tau])})
                        :RealType{};
                    update(raw_values[t][p][tau]-old_values[t][p][tau],
                           errors[t][p][tau],roundoff_scale);
                }
            }
        }
    };
    update_contour(old_correlations.Re,raw_correlations.Re,standard_errors.Re);
    update_contour(old_correlations.Im,raw_correlations.Im,standard_errors.Im);

    const auto& directions=raw_magnetization_Re.get_directions();
    const char symmetry=raw_magnetization_Re.get_symmetry();
    if( old_magnetization_Re.get_symmetry()!=symmetry
        ||old_magnetization_Im.get_symmetry()!=symmetry
        ||raw_magnetization_Im.get_symmetry()!=symmetry
        ||magnetization_Re_errors.get_symmetry()!=symmetry
        ||magnetization_Im_errors.get_symmetry()!=symmetry
        ||old_magnetization_Re.get_directions()!=directions
        ||old_magnetization_Im.get_directions()!=directions
        ||raw_magnetization_Im.get_directions()!=directions
        ||magnetization_Re_errors.get_directions()!=directions
        ||magnetization_Im_errors.get_directions()!=directions )
        throw std::invalid_argument("iteration residual magnetization inputs have different symmetries");
    for( size_t t=0;t<raw_magnetization_Re.size();++t )
        for( size_t c=0;c<magnetization_Re_errors.num_components();++c )
        {
            const RealType difference_Re=
                raw_magnetization_Re[t][c]-old_magnetization_Re[t][c];
            const RealType difference_Im=
                raw_magnetization_Im[t][c]-old_magnetization_Im[t][c];
            result.absolute=std::max(
                result.absolute,std::hypot(difference_Re,difference_Im));
            update(difference_Re,magnetization_Re_errors[t][c]);
            update(difference_Im,magnetization_Im_errors[t][c]);
        }
    return result;
}

}
