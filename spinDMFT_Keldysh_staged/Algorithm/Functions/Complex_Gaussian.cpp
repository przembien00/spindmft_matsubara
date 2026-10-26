#include"Complex_Gaussian.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<numeric>
#include<stdexcept>
#include<vector>
#include<unordered_map>

// Use the same Fortran BLAS ABI and integer/character types as Blaze's
// existing LAPACK wrappers. This needs no additional CBLAS headers/library.
#if !defined(INTEL_MKL_VERSION)
extern "C" {
void sgemm_(char*,char*,blaze::blas_int_t*,blaze::blas_int_t*,blaze::blas_int_t*,
    float*,float*,blaze::blas_int_t*,float*,blaze::blas_int_t*,float*,float*,
    blaze::blas_int_t*,blaze::fortran_charlen_t,blaze::fortran_charlen_t);
void dgemm_(char*,char*,blaze::blas_int_t*,blaze::blas_int_t*,blaze::blas_int_t*,
    double*,double*,blaze::blas_int_t*,double*,blaze::blas_int_t*,double*,double*,
    blaze::blas_int_t*,blaze::fortran_charlen_t,blaze::fortran_charlen_t);
}
#endif

namespace spinDMFT::Functions
{

struct GaussianBlockFactors
{
    using BatchMatrix=blaze::DynamicMatrix<RealType,blaze::columnMajor>;
    struct Block
    {
        std::vector<size_t> rows;
        std::shared_ptr<TakagiFactor> factor;
        std::shared_ptr<BatchMatrix> batch_factor;
    };
    std::vector<Block> blocks;
    size_t size{},rank{},largest{};
    RealType error{};
};

namespace
{

// Setup-only reconstructions, shared by identical blocks and released after
// the physical-basis diagnostic. Never retain these in the sampling state.
using FactorReconstructions=std::unordered_map<const TakagiFactor*,ComplexDynamicMatrix>;
std::shared_ptr<GaussianBlockFactors> symmetry_factors(
    const ComplexDynamicMatrix& covariance,
    FactorReconstructions* reconstructions=nullptr);

RealType frobenius_norm( const ComplexDynamicMatrix& matrix )
{
    RealType square_sum{};
    for( size_t i=0;i<matrix.rows();++i )
        for( size_t j=0;j<matrix.columns();++j ) square_sum+=std::norm(matrix(i,j));
    return std::sqrt(square_sum);
}

RealType transpose_symmetry_error( const ComplexDynamicMatrix& matrix )
{
    if( matrix.rows()!=matrix.columns() )
        throw std::invalid_argument("Takagi factorization requires a square matrix");
    const RealType denominator=frobenius_norm(matrix);
    const RealType numerator=frobenius_norm(matrix-blaze::trans(matrix));
    return denominator>RealType{0.}?numerator/denominator:numerator;
}

RealType symmetry_tolerance( const size_t size )
{
    return std::numeric_limits<RealType>::epsilon()*RealType{1000.}
          *static_cast<RealType>(std::max(size_t{1},size));
}


}

static TakagiFactor autonne_takagi_impl( const ComplexDynamicMatrix& Gamma_input,
    ComplexDynamicMatrix* reconstruction )
{
    if( Gamma_input.rows()!=Gamma_input.columns() )
        throw std::invalid_argument("Autonne--Takagi factorization requires a square matrix");
    const size_t n=Gamma_input.rows();
    if( transpose_symmetry_error(Gamma_input)>symmetry_tolerance(n) )
        throw std::invalid_argument("Autonne--Takagi input is not complex symmetric");

    TakagiFactor result{};
    if( n==0 ) return result;
    using RealMatrix=blaze::DynamicMatrix<RealType,blaze::rowMajor>;
    using SymmetricRealMatrix=blaze::SymmetricMatrix<RealMatrix>;
    SymmetricRealMatrix lift(2*n);
    for( size_t i=0;i<n;++i )
        for( size_t j=0;j<n;++j )
        {
            const RealType re=std::real(Gamma_input(i,j));
            const RealType im=std::imag(Gamma_input(i,j));
            lift(i,j)=re;
            lift(i,n+j)=im;
            lift(n+i,j)=im;
            lift(n+i,n+j)=-re;
        }

    blaze::DynamicVector<RealType,blaze::columnVector> eigenvalues(2*n);
    RealMatrix eigenvectors(2*n,2*n);
    // Blaze stores one normalized eigenvector in each row for row-major output.
    blaze::eigen(lift,eigenvalues,eigenvectors);
    const RealType sigma_max=std::max(RealType{0.},eigenvalues[2*n-1]);
    const RealType tolerance=std::numeric_limits<RealType>::epsilon()*RealType{100.}
        *static_cast<RealType>(2*n)*std::max(RealType{1.},sigma_max);
    size_t rank{};
    for( const RealType value:eigenvalues ) if( value>tolerance ) ++rank;
    result.L.resize(n,rank,false);
    result.numerical_rank=rank;
    size_t column{};
    for( size_t ev=0;ev<2*n;++ev )
    {
        const RealType sigma=eigenvalues[ev];
        if( sigma<=tolerance ) continue;
        result.singular_values.push_back(sigma);
        const RealType root=std::sqrt(sigma);
        for( size_t i=0;i<n;++i )
            result.L(i,column)=root*ComplexType{eigenvectors(ev,i),eigenvectors(ev,n+i)};
        ++column;
    }
    ComplexDynamicMatrix reconstructed=result.L*blaze::trans(result.L);
    const RealType denominator=frobenius_norm(Gamma_input);
    result.reconstruction_error=denominator>RealType{0.}
        ?frobenius_norm(Gamma_input-reconstructed)/denominator
        :frobenius_norm(reconstructed);
    if(reconstruction)*reconstruction=std::move(reconstructed);
    return result;
}

TakagiFactor autonne_takagi(const ComplexDynamicMatrix& Gamma_input,
    RealType takagi_tolerance)
{
    auto factor=autonne_takagi_impl(Gamma_input,nullptr);
    if(!std::isfinite(factor.reconstruction_error)
        ||factor.reconstruction_error>takagi_tolerance)
        throw std::runtime_error("Autonne--Takagi factor failed pseudo-covariance reconstruction");
    return factor;
}

JointComplexGaussianSampler::ContourFieldSample
JointComplexGaussianSampler::contour_field_from_latent(
    const LatentVector& latent )
{
    // Return the directly sampled field without additional propagation nodes.
    ContourFieldSample result{};
    result.edge_field=field_from_latent(latent);
    return result;
}

JointComplexGaussianSampler::ContourFieldSample
JointComplexGaussianSampler::draw_contour_field(
    std::mt19937& engine )
{
    return contour_field_from_latent(
        draw_latent(engine));
}

namespace
{
std::shared_ptr<GaussianBlockFactors> symmetry_factors(
    const ComplexDynamicMatrix& covariance,
    FactorReconstructions* reconstructions )
{
    const size_t n=covariance.rows();
    if( n!=covariance.columns()||transpose_symmetry_error(covariance)>symmetry_tolerance(n) )
        throw std::invalid_argument("Gaussian covariance is not complex symmetric");
    auto result=std::make_shared<GaussianBlockFactors>(); result->size=n;
    std::vector<size_t> parent(n); std::iota(parent.begin(),parent.end(),size_t{});
    const auto root=[&](size_t i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
    for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j)
        // Never discard small but nonzero couplings to manufacture symmetry.
        if(covariance(i,j)!=ComplexType{}||covariance(j,i)!=ComplexType{})
            parent[root(j)]=root(i);
    std::vector<std::vector<size_t>> rows(n);
    for(size_t i=0;i<n;++i)rows[root(i)].push_back(i);
    std::vector<ComplexDynamicMatrix> unique_covariances;
    std::vector<std::shared_ptr<TakagiFactor>> unique_factors;
    RealType sigma_max{};
    for(auto& indices:rows)if(!indices.empty())
    {
        ComplexDynamicMatrix block(indices.size(),indices.size());
        for(size_t i=0;i<indices.size();++i)for(size_t j=0;j<indices.size();++j)
            block(i,j)=covariance(indices[i],indices[j]);
        size_t identical=0;
        for(;identical<unique_covariances.size();++identical)
        {
            const auto& previous=unique_covariances[identical];
            if(previous.rows()!=block.rows())continue;
            bool equal=true;
            for(size_t i=0;i<block.rows()&&equal;++i)for(size_t j=0;j<block.rows();++j)
                if(block(i,j)!=previous(i,j)){equal=false;break;}
            if(equal)break;
        }
        std::shared_ptr<TakagiFactor> factor;
        if(identical<unique_factors.size())factor=unique_factors[identical];
        else
        {
            factor=std::make_shared<TakagiFactor>();
            auto* reconstructed=reconstructions?&(*reconstructions)[factor.get()]:nullptr;
            *factor=autonne_takagi_impl(block,reconstructed);
            for(const auto sigma:factor->singular_values)sigma_max=std::max(sigma_max,sigma);
            unique_covariances.push_back(std::move(block));unique_factors.push_back(factor);
        }
        result->largest=std::max(result->largest,indices.size());
        result->blocks.push_back({std::move(indices),factor});
    }
    // Local factorizations retain a superset of the whole-matrix numerical
    // rank. Apply its original global threshold to all blocks consistently.
    const RealType tolerance=std::numeric_limits<RealType>::epsilon()*RealType{100.}
        *static_cast<RealType>(2*n)*std::max(RealType{1.},sigma_max);
    for(size_t k=0;k<unique_factors.size();++k)
    {
        auto& factor=*unique_factors[k];
        size_t rank{};for(auto sigma:factor.singular_values)if(sigma>tolerance)++rank;
        if(rank!=factor.numerical_rank)
        {
            ComplexDynamicMatrix L(factor.L.rows(),rank);size_t column{};
            std::vector<RealType> retained;
            for(size_t j=0;j<factor.numerical_rank;++j)if(factor.singular_values[j]>tolerance)
            {
                for(size_t i=0;i<L.rows();++i)L(i,column)=factor.L(i,j);
                retained.push_back(factor.singular_values[j]);++column;
            }
            factor.L=std::move(L);factor.singular_values=std::move(retained);factor.numerical_rank=rank;
            // Only a changed factor invalidates the reconstruction and its
            // already-computed local residual. Preserve the global cutoff.
            ComplexDynamicMatrix reconstructed=factor.L*blaze::trans(factor.L);
            const auto& target=unique_covariances[k];
            const RealType scale=frobenius_norm(target);
            factor.reconstruction_error=scale>RealType{}?frobenius_norm(target-reconstructed)/scale:frobenius_norm(reconstructed);
            if(reconstructions)(*reconstructions)[&factor]=std::move(reconstructed);
        }
    }
    RealType residual{},scale{};
    for(const auto& block:result->blocks)
    {
        result->rank+=block.factor->numerical_rank;
        RealType norm{};
        for(auto i:block.rows)for(auto j:block.rows)norm+=std::norm(covariance(i,j));
        scale+=norm;residual+=norm*block.factor->reconstruction_error*block.factor->reconstruction_error;
    }
    result->error=scale>RealType{}?std::sqrt(residual/scale):std::sqrt(residual);
    return result;
}

JointComplexGaussianSampler::FieldVector block_field(
    const GaussianBlockFactors& factors,const JointComplexGaussianSampler::LatentVector& latent )
{
    if(latent.size()!=factors.rank)throw std::invalid_argument("Gaussian latent-state rank mismatch");
    JointComplexGaussianSampler::FieldVector field(factors.size,ComplexType{});
    size_t offset{};
    for(const auto& block:factors.blocks)
    {
        const auto& L=block.factor->L;
        for(size_t i=0;i<L.rows();++i)
        {
            ComplexType value{};
            for(size_t j=0;j<L.columns();++j)value+=L(i,j)*latent[offset+j];
            field[block.rows[i]]=value;
        }
        offset+=L.columns();
    }
    return field;
}

std::vector<JointComplexGaussianSampler::ContourFieldSample> block_batch(
    JointComplexGaussianSampler& sampler,GaussianBlockFactors& factors,
    std::mt19937& engine,const size_t count )
{
    using Matrix=GaussianBlockFactors::BatchMatrix;
    Matrix latent(factors.rank,count);
    for(size_t sample=0;sample<count;++sample)
    {
        const auto r=sampler.draw_latent(engine);
        for(size_t i=0;i<r.size();++i)latent(i,sample)=r[i];
    }
    std::vector<JointComplexGaussianSampler::ContourFieldSample> fields(count);
    for(auto& field:fields)field.edge_field.resize(factors.size,false);
    size_t offset{};
    for(auto& block:factors.blocks)
    {
        const auto& L=block.factor->L;
        if(!block.batch_factor)
        {
            for(const auto& previous:factors.blocks)
                if(previous.factor==block.factor&&previous.batch_factor)
                {block.batch_factor=previous.batch_factor;break;}
            if(!block.batch_factor)
            {
                block.batch_factor=std::make_shared<Matrix>(2*L.rows(),L.columns());
                for(size_t j=0;j<L.columns();++j)for(size_t i=0;i<L.rows();++i)
                {
                    (*block.batch_factor)(2*i,j)=std::real(L(i,j));
                    (*block.batch_factor)(2*i+1,j)=std::imag(L(i,j));
                }
            }
        }
        auto& A=*block.batch_factor;
        Matrix values(2*L.rows(),count,RealType{});
        if(L.columns()!=0&&count!=0)
        {
            const auto integer=[](size_t value)
            {
                if(value>static_cast<size_t>(std::numeric_limits<blaze::blas_int_t>::max()))
                    throw std::overflow_error("Gaussian batch exceeds BLAS integer range");
                return static_cast<blaze::blas_int_t>(value);
            };
            auto m=integer(A.rows()),n=integer(count),k=integer(A.columns());
            auto lda=integer(A.spacing()),ldb=integer(latent.spacing()),ldc=integer(values.spacing());
            RealType one{1.},zero{};char no_transpose='N';
            // Stack Re L and Im L to use one real GEMM with the same real
            // latent coordinates. No artificial complex latent variables.
#ifdef USE_FLOAT
            sgemm_(&no_transpose,&no_transpose,&m,&n,&k,&one,A.data(),&lda,
                latent.data()+offset,&ldb,&zero,values.data(),&ldc
#else
            dgemm_(&no_transpose,&no_transpose,&m,&n,&k,&one,A.data(),&lda,
                latent.data()+offset,&ldb,&zero,values.data(),&ldc
#endif
#if !defined(INTEL_MKL_VERSION)
                ,blaze::fortran_charlen_t{1},blaze::fortran_charlen_t{1}
#endif
            );
        }
        for(size_t i=0;i<L.rows();++i)for(size_t sample=0;sample<count;++sample)
            fields[sample].edge_field[block.rows[i]]={values(2*i,sample),values(2*i+1,sample)};
        offset+=L.columns();
    }
    return fields;
}
}

std::vector<JointComplexGaussianSampler::ContourFieldSample>
JointComplexGaussianSampler::draw_contour_batch(
    std::mt19937& engine,const size_t count )
{
    std::vector<ContourFieldSample> fields;fields.reserve(count);
    for(size_t i=0;i<count;++i)fields.push_back(draw_contour_field(engine));
    return fields;
}

DenseComplexGaussianSampler::DenseComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance,RealType takagi_tolerance)
    : m_factors(symmetry_factors(covariance))
{
    validate_reconstruction(takagi_tolerance);
}
DenseComplexGaussianSampler::DenseComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance,const std::vector<RealType>& weights,
    RealType takagi_tolerance)
{
    if(covariance.rows()!=covariance.columns() || weights.size()!=covariance.rows() || weights.empty())
        throw std::invalid_argument("diagonal Gaussian weights must match the covariance");
    const RealType maximum=*std::max_element(weights.begin(),weights.end());
    for(const auto weight:weights)
    {
        if(!std::isfinite(weight)||weight<=RealType{})
            throw std::invalid_argument("Gaussian weights must be finite and positive");
        const RealType root=std::sqrt(weight/maximum)*std::sqrt(RealType{2.});
        if(root==RealType{})throw std::invalid_argument("Gaussian weight ratio underflows");
        m_weight_roots.push_back(root);
    }
    ComplexDynamicMatrix scaled=covariance;
    for(size_t i=0;i<scaled.rows();++i)for(size_t j=0;j<scaled.columns();++j)
        scaled(i,j)*=m_weight_roots[i]*m_weight_roots[j];
    FactorReconstructions reconstructions;
    m_factors=symmetry_factors(scaled,&reconstructions);
    reset(scaled);
    for(const auto& block:m_factors->blocks)
    {
        const auto& reconstructed=reconstructions.at(block.factor.get());
        for(size_t i=0;i<block.rows.size();++i)for(size_t j=0;j<block.rows.size();++j)
            scaled(block.rows[i],block.rows[j])=reconstructed(i,j)
                /(m_weight_roots[block.rows[i]]*m_weight_roots[block.rows[j]]);
    }
    const RealType norm=frobenius_norm(covariance);
    m_unweighted_reconstruction_error=norm>RealType{}?frobenius_norm(scaled-covariance)/norm
        :frobenius_norm(scaled);
    validate_reconstruction(takagi_tolerance);
}
void DenseComplexGaussianSampler::validate_reconstruction(RealType tolerance) const
{
    const RealType error=reconstruction_error();
    if(!std::isfinite(error)||error>tolerance)
        throw std::runtime_error("joint Gaussian factor failed pseudo-covariance reconstruction");
}
RealType DenseComplexGaussianSampler::reconstruction_error() const
{ return m_weight_roots.empty()?m_factors->error:m_unweighted_reconstruction_error; }
size_t DenseComplexGaussianSampler::latent_dimension() const { return m_factors->rank; }
size_t DenseComplexGaussianSampler::size() const { return m_factors->size; }
size_t DenseComplexGaussianSampler::largest_factorization_dimension() const { return m_factors->largest; }
DenseComplexGaussianSampler::LatentVector DenseComplexGaussianSampler::draw_latent(std::mt19937& engine)
{
    LatentVector latent(latent_dimension());
    for(auto& value:latent)value=m_standard_normal(engine);
    return latent;
}
ComplexDynamicMatrix DenseComplexGaussianSampler::leading_factor_rows(size_t rows) const
{
    if(rows>size())throw std::invalid_argument("factor row count exceeds field dimension");
    ComplexDynamicMatrix result(rows,latent_dimension(),ComplexType{});
    size_t offset{};
    for(const auto& block:m_factors->blocks)
    {
        const auto& L=block.factor->L;
        for(size_t i=0;i<L.rows();++i)if(block.rows[i]<rows)
            for(size_t j=0;j<L.columns();++j)
                result(block.rows[i],offset+j)=L(i,j)
                    /(m_weight_roots.empty()?RealType{1.}:m_weight_roots[block.rows[i]]);
        offset+=L.columns();
    }
    return result;
}

DenseComplexGaussianSampler::FieldVector DenseComplexGaussianSampler::field_from_latent(const LatentVector& latent)
{
    auto field=block_field(*m_factors,latent);
    if(!m_weight_roots.empty())for(size_t i=0;i<field.size();++i)field[i]/=m_weight_roots[i];
    return field;
}
DenseComplexGaussianSampler::FieldVector DenseComplexGaussianSampler::draw(std::mt19937& engine)
{ return field_from_latent(draw_latent(engine)); }
std::vector<JointComplexGaussianSampler::ContourFieldSample>
DenseComplexGaussianSampler::draw_contour_batch(std::mt19937& engine,size_t count)
{
    auto fields=block_batch(*this,*m_factors,engine,count);
    if(!m_weight_roots.empty())for(auto& sample:fields)
        for(size_t i=0;i<sample.edge_field.size();++i)sample.edge_field[i]/=m_weight_roots[i];
    return fields;
}

}
