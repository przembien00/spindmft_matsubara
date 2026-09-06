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

ComplexType rotated_connected_entry(
    const ps::ParameterSpace& pspace,
    const contour::CorrelationSet& connected_correlations,
    const contour::ContourLayout& layout,
    const contour::ContourIndex& first, const contour::ContourIndex& second )
{
    const auto& D=pspace.spin_model.coupling_matrix;
    ComplexType rotated{};
    for( size_t c=0;c<3;++c )
        for( size_t d=0;d<3;++d )
        {
            const ComplexType connected=contour::branch_correlation(
                connected_correlations,layout,{first.branch,first.point,c},
                                              {second.branch,second.point,d});
            rotated+=D(first.component,c)*D(second.component,d)*connected;
        }
    return pspace.JQ*pspace.JQ*rotated
        +ComplexType{pspace.noise.m_variance_in(first.component,second.component),
                     RealType{0.}};
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

CorrelationSet generate_initial_correlations( const ps::ParameterSpace& pspace,
                                              const FieldVector& spin_expectation )
{
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
    const ComplexMagnetizationTrajectory& magnetization_time )
{
    if( correlations.Re.empty() || correlations.Im.size()!=correlations.Re.size()
        || magnetization_time.size()!=correlations.Re.size() )
        throw std::invalid_argument(
            "connected contour primitive needs matching correlation and magnetization grids" );

    CorrelationSet connected=correlations;
    const ComplexFieldVector& magnetization_zero=magnetization_time.front();
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
                magnetization_time[t][direction[0]]
               *magnetization_zero[direction[1]];
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

SelfConsistentField self_consistent_equations(
    const ps::ParameterSpace& pspace, const CorrelationSet& correlations,
    const ComplexMagnetizationTrajectory& magnetization_time )
{
    const contour::ContourLayout layout{
        pspace.num_TimePoints,pspace.num_RealTimePoints};
    const CorrelationSet connected_correlations=
        connected_contour_primitive(correlations,magnetization_time);
    ComplexDynamicMatrix raw(layout.dimension(),layout.dimension(),ComplexType{});
    for( size_t row=0;row<layout.dimension();++row )
    {
        const auto first=layout.decode(row);
        for( size_t col=0;col<layout.dimension();++col )
            raw(row,col)=rotated_connected_entry(
                pspace,connected_correlations,layout,first,layout.decode(col));
    }

    ComplexDynamicMatrix covariance(layout.dimension(),layout.dimension(),ComplexType{});
    for( size_t row=0;row<layout.dimension();++row )
        for( size_t col=row;col<layout.dimension();++col )
        {
            covariance(row,col)=raw(row,col);
            covariance(col,row)=raw(row,col);
        }

    SelfConsistentField result{};
    result.mean_time.resize(magnetization_time.size());
    for( size_t t=0;t<magnetization_time.size();++t )
    {
        FieldVector physical_magnetization{};
        for( size_t c=0;c<3;++c )
            physical_magnetization[c]=std::real(magnetization_time[t][c]);
        result.mean_time[t]=pspace.JL
            *(pspace.spin_model.coupling_matrix*physical_magnetization);
    }
    result.branch_identity_error=complex_matrix_symmetry_error(raw);
    result.covariance_symmetry_error=complex_matrix_symmetry_error(covariance);
    result.covariance=std::move(covariance);
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
    const ps::ParameterSpace& pspace )
{
    const contour::ContourLayout layout{
        pspace.num_TimePoints,pspace.num_RealTimePoints};
    const CorrelationSet bath=harmonic_bath_primitive(pspace);
    ComplexDynamicMatrix raw(layout.dimension(),layout.dimension(),ComplexType{});
    for( size_t row=0;row<layout.dimension();++row )
        for( size_t col=0;col<layout.dimension();++col )
            raw(row,col)=contour::branch_correlation(
                bath,layout,layout.decode(row),layout.decode(col));

    ComplexDynamicMatrix covariance(layout.dimension(),layout.dimension(),ComplexType{});
    for( size_t row=0;row<layout.dimension();++row )
        for( size_t col=row;col<layout.dimension();++col )
        {
            covariance(row,col)=raw(row,col);
            covariance(col,row)=raw(row,col);
        }

    SelfConsistentField result{};
    result.mean_time.assign(pspace.num_RealTimePoints,FieldVector{});
    result.branch_identity_error=complex_matrix_symmetry_error(raw);
    result.covariance_symmetry_error=complex_matrix_symmetry_error(covariance);
    result.covariance=std::move(covariance);
    return result;
}

ContourTrajectory compute_contour_trajectory(
    const ps::ParameterSpace& pspace,
    const DenseComplexGaussianSampler::FieldVector& joint_field,
    const MeanFieldTrajectory& mean_field_time )
{
    const contour::ContourLayout layout{
        pspace.num_TimePoints,pspace.num_RealTimePoints};
    if( joint_field.size()!=layout.dimension() )
        throw std::invalid_argument("joint Keldysh field has the wrong dimension");
    if( mean_field_time.size()!=pspace.num_RealTimePoints )
        throw std::invalid_argument("mean-field trajectory has the wrong real-time grid");

    std::vector<ComplexFieldVector> imaginary_fields(pspace.num_TimePoints);
    std::vector<ComplexFieldVector> forward_fields(pspace.num_RealTimePoints);
    std::vector<ComplexFieldVector> backward_fields(pspace.num_RealTimePoints);
    for( size_t k=0;k<pspace.num_TimePoints;++k )
        for( size_t c=0;c<3;++c )
            imaginary_fields[k][c]=joint_field[layout.flat(contour::Branch::Matsubara,k,c)];
    for( size_t t=0;t<forward_fields.size();++t )
        for( size_t c=0;c<3;++c )
        {
            forward_fields[t][c]=joint_field[layout.flat(contour::Branch::Forward,t,c)];
            backward_fields[t][c]=joint_field[layout.flat(contour::Branch::Backward,t,c)];
        }

    ContourTrajectory result{};
    std::vector<Operator> imaginary_steps(pspace.num_TimeSteps);
    for( size_t k=0;k<pspace.num_TimeSteps;++k )
    {
        imaginary_steps[k]=cfet4_step(
            pspace,imaginary_fields[k+1],mean_field_time.front(),
            imaginary_fields[k],mean_field_time.front(),
            ComplexType{-pspace.delta_t,RealType{0.}});
    }
    std::vector<Operator> prefix(pspace.num_TimeSteps+1,IDENTITY);
    std::vector<Operator> suffix(pspace.num_TimeSteps+1,IDENTITY);
    for( size_t k=0;k<pspace.num_TimeSteps;++k )
        prefix[k+1]=imaginary_steps[k]*prefix[k];
    for( size_t k=pspace.num_TimeSteps;k-- >0; )
        suffix[k]=suffix[k+1]*imaginary_steps[k];
    result.imaginary_density_operator=prefix.back();
    result.partition_function=blaze::trace(result.imaginary_density_operator);

    const std::array<const Observable*,3> spins{&S_X,&S_Y,&S_Z};
    for( size_t c=0;c<3;++c )
    {
        result.imaginary_edge_insertions[c].resize(pspace.num_TimePoints);
        for( size_t tau=0;tau<pspace.num_TimePoints;++tau )
            result.imaginary_edge_insertions[c][tau]=
                suffix[tau]*(*spins[c])*prefix[tau];
    }

    result.forward_steps.assign(pspace.num_RealTimePoints,IDENTITY);
    result.backward_steps.assign(pspace.num_RealTimePoints,IDENTITY);
    for( size_t t=1;t<pspace.num_RealTimePoints;++t )
    {
        result.forward_steps[t]=cfet4_step(
            pspace,forward_fields[t],mean_field_time[t],
            forward_fields[t-1],mean_field_time[t-1],
            ComplexType{RealType{0.},-pspace.delta_real_t});
        result.backward_steps[t]=cfet4_step(
            pspace,backward_fields[t-1],mean_field_time[t-1],
            backward_fields[t],mean_field_time[t],
            ComplexType{RealType{0.},+pspace.delta_real_t});
    }
    Operator final_density=result.imaginary_density_operator;
    for( size_t t=1;t<pspace.num_RealTimePoints;++t )
        final_density=result.forward_steps[t]*final_density
                    *result.backward_steps[t];
    result.final_closed_contour_trace=blaze::trace(final_density);
    return result;
}

namespace
{
void accumulate_contour_observables(
    rtd::RunTimeData& rtdata,
    const ContourTrajectory& trajectory,
    const std::string& insertion_strategy )
{
    const std::array<const Observable*,3> spins{&S_X,&S_Y,&S_Z};
    const size_t num_real_points=rtdata.num_real_time_points();
    if( num_real_points==0 || trajectory.forward_steps.size()!=num_real_points
        || trajectory.backward_steps.size()!=num_real_points )
        throw std::invalid_argument("trajectory and measurement real-time grids differ");
    if( insertion_strategy!="closed-contour" && insertion_strategy!="prefix" )
        throw std::invalid_argument("unknown spin insertion strategy");

    // left[t] = B_-(T,0) U_+(t,T), where U_+(t,T) denotes the forward
    // continuation t -> T: U_N ... U_{t+1}, not an inverse propagator.
    // Build it right-to-left so every forward/backward step occurs once in
    // left[t] S U_+(t,0), including for t=0. No branch cancellation is assumed.
    Operator backward_total=IDENTITY;
    for( size_t t=1;t<num_real_points;++t )
        backward_total=backward_total*trajectory.backward_steps[t];
    std::vector<Operator> left(num_real_points,IDENTITY);
    left.back()=backward_total;
    for( size_t t=num_real_points-1;t>0;--t )
        left[t-1]=left[t]*trajectory.forward_steps[t];

    // Preserve the prefix closure trajectory D(t). Its final value is also the
    // fixed denominator and pCN likelihood when closed-contour normalization
    // is selected; intermediate values remain diagnostics.
    Operator density=trajectory.imaginary_density_operator;
    Operator forward=IDENTITY;
    Operator backward_prefix=IDENTITY;
    std::array<bool,3> measured_directions{};
    for( size_t p=0;p<rtdata.num_correlation_components();++p )
        measured_directions[rtdata.correlation_direction(p)[0]]=true;
    for( size_t c=0;c<rtdata.num_magnetization_components();++c )
        measured_directions[rtdata.magnetization_direction(c)]=true;
    std::array<Operator,3> measured_spins{};
    for( size_t t=0;t<num_real_points;++t )
    {
        if( t>0 )
        {
            forward=trajectory.forward_steps[t]*forward;
            backward_prefix=backward_prefix*trajectory.backward_steps[t];
            density=trajectory.forward_steps[t]*density*trajectory.backward_steps[t];
        }
        rtdata.accumulate_closed_contour_trace(t,blaze::trace(density));
        for( size_t direction=0;direction<3;++direction )
            if( measured_directions[direction] )
                measured_spins[direction]=insertion_strategy=="prefix"
                    ?backward_prefix*(*spins[direction])*forward
                    :left[t]*(*spins[direction])*forward;
        for( size_t c=0;c<rtdata.num_magnetization_components();++c )
        {
            const size_t direction=rtdata.magnetization_direction(c);
            rtdata.accumulate_magnetization(
                t,c,trace_product(trajectory.imaginary_density_operator,
                                  measured_spins[direction]));
        }
        for( size_t p=0;p<rtdata.num_correlation_components();++p )
        {
            const auto direction=rtdata.correlation_direction(p);
            for( size_t tau=0;tau<rtdata.num_imaginary_edge_points();++tau )
                rtdata.accumulate_edge_correlation(t,p,tau,trace_product(
                    trajectory.imaginary_edge_insertions[direction[1]][tau],
                    measured_spins[direction[0]]));
        }
    }
}
}

void compute_contour_correlations( rtd::RunTimeData& rtdata,
                                   const ContourTrajectory& trajectory,
                                   const RealType observable_normalization,
                                   const std::string& insertion_strategy )
{
    rtdata.begin_sample(
        trajectory.partition_function,observable_normalization);
    accumulate_contour_observables(rtdata,trajectory,insertion_strategy);
    rtdata.end_sample();
}

void compute_contour_pair_correlations(
    rtd::RunTimeData& rtdata,
    const ContourTrajectory& positive_trajectory,
    const ContourTrajectory& negative_trajectory,
    const RealType observable_normalization,
    const std::string& insertion_strategy )
{
    rtdata.begin_sample(
        positive_trajectory.partition_function
            +negative_trajectory.partition_function,
        observable_normalization);
    accumulate_contour_observables(
        rtdata,positive_trajectory,insertion_strategy);
    accumulate_contour_observables(
        rtdata,negative_trajectory,insertion_strategy);
    rtdata.end_sample();
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

ComplexMagnetizationTrajectory mix_magnetization_trajectory(
    const ComplexMagnetizationTrajectory& old_values,
    const ComplexMagnetizationTrajectory& raw_values,
    RealType alpha )
{
    if( alpha<=RealType{0.} || alpha>RealType{1.} )
        throw std::invalid_argument("fixed-point mixing alpha must lie in (0,1]");
    if( old_values.size()!=raw_values.size() )
        throw std::invalid_argument("magnetization trajectories have different grids");
    ComplexMagnetizationTrajectory result=raw_values;
    for( size_t t=0;t<result.size();++t )
        for( size_t c=0;c<3;++c )
            result[t][c]=(RealType{1.}-alpha)*old_values[t][c]
                        +alpha*raw_values[t][c];
    return result;
}

ComplexMagnetizationTrajectory project_constant_magnetization(
    const ComplexMagnetizationTrajectory& values )
{
    if( values.empty() )
        throw std::invalid_argument("cannot project an empty magnetization trajectory");
    return ComplexMagnetizationTrajectory(values.size(),values.front());
}

RealType max_magnetization_difference(
    const ComplexMagnetizationTrajectory& old_values,
    const ComplexMagnetizationTrajectory& raw_values )
{
    if( old_values.size()!=raw_values.size() )
        throw std::invalid_argument("magnetization trajectories have different grids");
    RealType largest{};
    for( size_t t=0;t<old_values.size();++t )
        for( size_t c=0;c<3;++c )
            largest=std::max(largest,std::abs(raw_values[t][c]-old_values[t][c]));
    return largest;
}

IterationResidual iteration_residual(
    const CorrelationSet& old_correlations,
    const CorrelationSet& raw_correlations,
    const CorrelationSet& standard_errors,
    const ComplexMagnetizationTrajectory& old_magnetization,
    const ComplexMagnetizationTrajectory& raw_magnetization,
    const std::vector<FieldVector>& magnetization_Re_errors,
    const std::vector<FieldVector>& magnetization_Im_errors )
{
    if( old_correlations.Re.size()!=raw_correlations.Re.size()
        ||old_correlations.Im.size()!=raw_correlations.Im.size()
        ||standard_errors.Re.size()!=raw_correlations.Re.size()
        ||standard_errors.Im.size()!=raw_correlations.Im.size()
        ||old_magnetization.size()!=raw_magnetization.size()
        ||magnetization_Re_errors.size()!=raw_magnetization.size()
        ||magnetization_Im_errors.size()!=raw_magnetization.size() )
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

    for( size_t t=0;t<raw_magnetization.size();++t )
        for( size_t c=0;c<3;++c )
        {
            const ComplexType difference=
                raw_magnetization[t][c]-old_magnetization[t][c];
            result.absolute=std::max(result.absolute,std::abs(difference));
            update(std::real(difference),magnetization_Re_errors[t][c]);
            update(std::imag(difference),magnetization_Im_errors[t][c]);
        }
    return result;
}

}
