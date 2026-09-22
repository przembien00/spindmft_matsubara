#include"Complex_Gaussian.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<numeric>
#include<stdexcept>
#include<vector>
#include<unordered_map>

#include<fftw3.h>

// Use the same Fortran BLAS ABI and integer/character types as Blaze's
// existing LAPACK wrappers. This needs no additional CBLAS headers/library.
#if !defined(INTEL_MKL_VERSION)
extern "C" {
void sgemv_(char*,blaze::blas_int_t*,blaze::blas_int_t*,float*,float*,
    blaze::blas_int_t*,float*,blaze::blas_int_t*,float*,float*,blaze::blas_int_t*,
    blaze::fortran_charlen_t);
void dgemv_(char*,blaze::blas_int_t*,blaze::blas_int_t*,double*,double*,
    blaze::blas_int_t*,double*,blaze::blas_int_t*,double*,double*,blaze::blas_int_t*,
    blaze::fortran_charlen_t);
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

std::shared_ptr<GaussianBlockFactors> symmetry_factors(
    const ComplexDynamicMatrix& covariance,bool use_svd);

// Pack a complex-by-real product as one real BLAS product. Adjacent rows
// hold Re L and Im L; columns are the original independent real latent modes.
using RealFactorMatrix=GaussianBlockFactors::BatchMatrix;

blaze::blas_int_t blas_integer(const size_t value)
{
    if(value>static_cast<size_t>(std::numeric_limits<blaze::blas_int_t>::max()))
        throw std::overflow_error("FFT Gaussian product exceeds BLAS integer range");
    return static_cast<blaze::blas_int_t>(value);
}

void multiply_real_factor(RealFactorMatrix& factor,const RealType* latent,
    const size_t count,const size_t latent_stride,RealFactorMatrix& values)
{
    values.resize(factor.rows(),count,false);
    // GEMV/GEMM quick returns need not overwrite the output for an empty
    // inner dimension. Clear it explicitly, including reused batch buffers.
    if(factor.columns()==0||count==0)
    {
        reset(values);
        return;
    }
    auto m=blas_integer(factor.rows()),k=blas_integer(factor.columns());
    auto lda=blas_integer(factor.spacing());
    RealType one{1.},zero{};
    char no_transpose='N';
    // The legacy Fortran signatures are mutable; BLAS does not modify x/B.
    auto* input=const_cast<RealType*>(latent);
    if(count==1)
    {
        blaze::blas_int_t increment=1;
#ifdef USE_FLOAT
        sgemv_(&no_transpose,&m,&k,&one,factor.data(),&lda,input,&increment,
               &zero,values.data(),&increment
#else
        dgemv_(&no_transpose,&m,&k,&one,factor.data(),&lda,input,&increment,
               &zero,values.data(),&increment
#endif
#if !defined(INTEL_MKL_VERSION)
               ,blaze::fortran_charlen_t{1}
#endif
        );
    }
    else
    {
        auto n=blas_integer(count),ldb=blas_integer(latent_stride),ldc=blas_integer(values.spacing());
#ifdef USE_FLOAT
        sgemm_(&no_transpose,&no_transpose,&m,&n,&k,&one,factor.data(),&lda,
               input,&ldb,&zero,values.data(),&ldc
#else
        dgemm_(&no_transpose,&no_transpose,&m,&n,&k,&one,factor.data(),&lda,
               input,&ldb,&zero,values.data(),&ldc
#endif
#if !defined(INTEL_MKL_VERSION)
               ,blaze::fortran_charlen_t{1},blaze::fortran_charlen_t{1}
#endif
        );
    }
}

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

#ifdef USE_FLOAT
using FFTPlan=fftwf_plan;
using FFTComplex=fftwf_complex;
FFTPlan make_fft_plan( const size_t count,const size_t internal,
                       ComplexType* data,const int direction )
{
    int length=static_cast<int>(count);
    return fftwf_plan_many_dft(1,&length,static_cast<int>(internal),
        reinterpret_cast<FFTComplex*>(data),nullptr,static_cast<int>(internal),1,
        reinterpret_cast<FFTComplex*>(data),nullptr,static_cast<int>(internal),1,
        direction,FFTW_ESTIMATE);
}
void execute_fft( const FFTPlan plan ) { fftwf_execute(plan); }
void destroy_fft_plan( const FFTPlan plan ) { if(plan) fftwf_destroy_plan(plan); }
#else
using FFTPlan=fftw_plan;
using FFTComplex=fftw_complex;
FFTPlan make_fft_plan( const size_t count,const size_t internal,
                       ComplexType* data,const int direction )
{
    int length=static_cast<int>(count);
    return fftw_plan_many_dft(1,&length,static_cast<int>(internal),
        reinterpret_cast<FFTComplex*>(data),nullptr,static_cast<int>(internal),1,
        reinterpret_cast<FFTComplex*>(data),nullptr,static_cast<int>(internal),1,
        direction,FFTW_ESTIMATE);
}
void execute_fft( const FFTPlan plan ) { fftw_execute(plan); }
void destroy_fft_plan( const FFTPlan plan ) { if(plan) fftw_destroy_plan(plan); }
#endif

static_assert(sizeof(ComplexType)==sizeof(FFTComplex),
              "ComplexType and FFTW complex storage must have identical sizes");

void unitary_fft_matrix_group( ComplexDynamicMatrix& matrix,const size_t offset,
                               const size_t count,const size_t internal,
                               const bool transform_rows )
{
    std::vector<ComplexType> buffer(count*internal);
    const FFTPlan plan=make_fft_plan(count,internal,buffer.data(),FFTW_FORWARD);
    if( !plan ) throw std::runtime_error("failed to create forward FFTW plan");
    const RealType normalization=RealType{1.}/std::sqrt(static_cast<RealType>(count));
    const size_t outer=transform_rows?matrix.columns():matrix.rows();
    for( size_t fixed=0;fixed<outer;++fixed )
    {
        for( size_t t=0;t<count;++t )
            for( size_t component=0;component<internal;++component )
            {
                const size_t index=offset+t*internal+component;
                buffer[t*internal+component]=transform_rows
                    ?matrix(index,fixed):matrix(fixed,index);
            }
        execute_fft(plan);
        for( size_t t=0;t<count;++t )
            for( size_t component=0;component<internal;++component )
            {
                const size_t index=offset+t*internal+component;
                const ComplexType value=normalization*buffer[t*internal+component];
                if( transform_rows ) matrix(index,fixed)=value;
                else matrix(fixed,index)=value;
            }
    }
    destroy_fft_plan(plan);
}

void validate_fft_covariance( const CovarianceSource& physical,
                              const size_t num_matsubara_intervals,
                              const size_t num_real )
{
    const size_t expected_size=3*(num_matsubara_intervals+1)+6*num_real;
    if( physical.rows()!=expected_size||physical.columns()!=expected_size )
        throw std::invalid_argument("FFT Gaussian covariance has incompatible contour dimensions");
}

ComplexType doubled_contour_value_unchecked(
    const CovarianceSource& physical,const size_t num_matsubara_intervals,
    const size_t num_real,const size_t first,const size_t second )
{
    const size_t num_matsubara_points=num_matsubara_intervals+1;
    const size_t matsubara_size=3*num_matsubara_points;
    const size_t embedded_real=2*num_real;
    const auto physical_real_index=[&](const size_t branch,const size_t time,
                                       const size_t component)
    {
        return matsubara_size+3*(branch*num_real+time)+component;
    };
    if( first<matsubara_size&&second<matsubara_size )
        return physical(first,second);

    // The embedding is complex symmetric, so only the Matsubara-real
    // orientation needs to be evaluated explicitly.
    if( first>=matsubara_size&&second<matsubara_size )
        return doubled_contour_value_unchecked(physical,num_matsubara_intervals,
                                               num_real,second,first);
    if( first<matsubara_size )
    {
        const size_t tau=first/3;
        const size_t a=first%3;
        const size_t relative=second-matsubara_size;
        const size_t t=relative/6;
        const size_t branch=(relative%6)/3;
        const size_t b=relative%3;
        if( t<num_real )
            return physical(first,physical_real_index(branch,t,b));
        if( t==num_real ) return ComplexType{};
        const size_t positive_time=embedded_real-t;
        const size_t reflected_matsubara=3*one_sided_edge_reflection_index(
            tau,num_matsubara_intervals)+b;
        return physical(reflected_matsubara,
                        physical_real_index(branch,positive_time,a));
    }

    const size_t first_relative=first-matsubara_size;
    const size_t second_relative=second-matsubara_size;
    const size_t first_time=first_relative/6;
    const size_t second_time=second_relative/6;
    const size_t first_branch=(first_relative%6)/3;
    const size_t first_component=first_relative%3;
    const size_t second_branch=(second_relative%6)/3;
    const size_t second_component=second_relative%3;
    const size_t lag=(first_time+embedded_real-second_time)%embedded_real;
    if( lag<num_real )
        return physical(
            physical_real_index(first_branch,lag,first_component),
            physical_real_index(second_branch,0,second_component));
    if( lag>num_real )
        return physical(
            physical_real_index(first_branch,0,first_component),
            physical_real_index(second_branch,embedded_real-lag,second_component));
    return ComplexType{};
}

RealType physical_marginal_error( const CovarianceSource& physical,
                                  const size_t num_matsubara_intervals,
                                  const size_t num_real )
{
    validate_fft_covariance(physical,num_matsubara_intervals,num_real);
    const size_t num_matsubara_points=num_matsubara_intervals+1;
    const size_t matsubara_size=3*num_matsubara_points;
    RealType residual{},scale{};
    const auto embedded_index=[&](const size_t physical_index)
    {
        if( physical_index<matsubara_size ) return physical_index;
        const size_t relative=physical_index-matsubara_size;
        const size_t branch=relative/(3*num_real);
        const size_t within=relative%(3*num_real);
        return matsubara_size+6*(within/3)+3*branch+within%3;
    };
    for( size_t i=0;i<physical.rows();++i )
        for( size_t j=0;j<physical.columns();++j )
        {
            residual+=std::norm(doubled_contour_value_unchecked(
                physical,num_matsubara_intervals,num_real,
                embedded_index(i),embedded_index(j))-physical(i,j));
            scale+=std::norm(physical(i,j));
        }
    return scale>RealType{0.}?std::sqrt(residual/scale):std::sqrt(residual);
}

struct FrequencyBlockFactorization
{
    struct Block
    {
        std::vector<size_t> rows{};
        std::shared_ptr<TakagiFactor> factor{};
    };

    std::vector<Block> blocks{};
    size_t total_rank{};
    RealType approximation_error{};
    RealType reconstruction_error{};
    size_t largest_dimension{};
};

FrequencyBlockFactorization factor_frequency_blocks(
    const CovarianceSource& physical,const size_t num_matsubara_intervals,
    const size_t num_real_points,const RealType delta_real_time,const RealType cutoff )
{
    if( delta_real_time<=RealType{0.} )
        throw std::invalid_argument("FFT Gaussian sampler needs positive real-time spacing");
    const size_t matsubara_size=3*(num_matsubara_intervals+1);
    const size_t embedded_real_points=2*num_real_points;
    const size_t real_size=6*embedded_real_points;
    validate_fft_covariance(physical,num_matsubara_intervals,num_real_points);

    std::vector<size_t> low_rows(matsubara_size);
    std::iota(low_rows.begin(),low_rows.end(),size_t{});
    std::vector<std::vector<size_t>> high_groups;
    std::vector<bool> visited(embedded_real_points,false);
    const RealType two_pi=RealType{2.}*std::acos(RealType{-1.});
    for( size_t mode=0;mode<embedded_real_points;++mode )
    {
        if( visited[mode] ) continue;
        const size_t partner=(embedded_real_points-mode)%embedded_real_points;
        visited[mode]=true;
        visited[partner]=true;
        std::vector<size_t> group;
        for( size_t component=0;component<6;++component )
            group.push_back(matsubara_size+6*mode+component);
        if( partner!=mode )
            for( size_t component=0;component<6;++component )
                group.push_back(matsubara_size+6*partner+component);
        const size_t absolute_mode=std::min(mode,embedded_real_points-mode);
        const RealType omega=two_pi*static_cast<RealType>(absolute_mode)
            /(static_cast<RealType>(embedded_real_points)*delta_real_time);
        if( cutoff<RealType{0.}||omega<=cutoff )
            low_rows.insert(low_rows.end(),group.begin(),group.end());
        else high_groups.push_back(std::move(group));
    }

    std::vector<std::vector<size_t>> blocks;
    blocks.push_back(std::move(low_rows));
    blocks.insert(blocks.end(),high_groups.begin(),high_groups.end());
    const size_t frequency_size=matsubara_size+real_size;
    std::vector<size_t> block_id(frequency_size),block_position(frequency_size);
    for( size_t block=0;block<blocks.size();++block )
        for( size_t position=0;position<blocks[block].size();++position )
        {
            const size_t row=blocks[block][position];
            block_id[row]=block;
            block_position[row]=position;
        }

    // Allocate only the retained diagonal frequency blocks.  In particular,
    // no full doubled-grid or full transformed covariance is materialized.
    std::vector<ComplexDynamicMatrix> block_covariances;
    block_covariances.reserve(blocks.size());
    for( const auto& rows:blocks )
        block_covariances.emplace_back(rows.size(),rows.size(),ComplexType{});
    RealType full_square{},discarded_square{};

    ComplexDynamicMatrix matsubara(matsubara_size,matsubara_size);
    for( size_t i=0;i<matsubara_size;++i )
        for( size_t j=0;j<matsubara_size;++j )
        {
            matsubara(i,j)=physical(i,j);
            full_square+=std::norm(matsubara(i,j));
        }
    for( const bool rows:{true,false} )
        unitary_fft_matrix_group(
            matsubara,0,num_matsubara_intervals,3,rows);
    for( size_t i=0;i<matsubara_size;++i )
        for( size_t j=0;j<matsubara_size;++j )
            block_covariances[0](block_position[i],block_position[j])=matsubara(i,j);

    ComplexDynamicMatrix mixed(matsubara_size,real_size);
    for( size_t i=0;i<matsubara_size;++i )
        for( size_t j=0;j<real_size;++j )
        {
            mixed(i,j)=doubled_contour_value_unchecked(
                physical,num_matsubara_intervals,num_real_points,
                i,matsubara_size+j);
            full_square+=RealType{2.}*std::norm(mixed(i,j));
        }
    unitary_fft_matrix_group(
        mixed,0,num_matsubara_intervals,3,true);
    unitary_fft_matrix_group(
        mixed,0,embedded_real_points,6,false);
    for( size_t i=0;i<matsubara_size;++i )
        for( size_t j=0;j<real_size;++j )
        {
            const size_t global=matsubara_size+j;
            if( block_id[global]!=0 )
            {
                // Both transpose-related mixed sectors are discarded.
                discarded_square+=RealType{2.}*std::norm(mixed(i,j));
                continue;
            }
            const size_t row=block_position[i];
            const size_t column=block_position[global];
            block_covariances[0](row,column)=mixed(i,j);
            block_covariances[0](column,row)=mixed(i,j);
        }

    // A block-circulant real covariance obeys
    // (F Gamma F^T)_{k,l}=delta_{k,-l} FFT[K](k).  One batched FFT of
    // the 6x6 lag kernel therefore populates every retained {k,-k} block.
    std::vector<ComplexType> real_spectra(embedded_real_points*36);
    for( size_t lag=0;lag<embedded_real_points;++lag )
        for( size_t first_component=0;first_component<6;++first_component )
            for( size_t second_component=0;second_component<6;++second_component )
            {
                const size_t first=matsubara_size+6*lag+first_component;
                const size_t second=matsubara_size+second_component;
                const ComplexType value=doubled_contour_value_unchecked(
                    physical,num_matsubara_intervals,num_real_points,first,second);
                real_spectra[36*lag+6*first_component+second_component]=value;
                full_square+=static_cast<RealType>(embedded_real_points)*std::norm(value);
            }
    const FFTPlan real_forward=make_fft_plan(
        embedded_real_points,36,real_spectra.data(),FFTW_FORWARD);
    if( !real_forward ) throw std::runtime_error("failed to create real covariance FFTW plan");
    execute_fft(real_forward);
    destroy_fft_plan(real_forward);
    for( size_t mode=0;mode<embedded_real_points;++mode )
    {
        const size_t partner=(embedded_real_points-mode)%embedded_real_points;
        const size_t block=block_id[matsubara_size+6*mode];
        for( size_t first_component=0;first_component<6;++first_component )
            for( size_t second_component=0;second_component<6;++second_component )
            {
                const size_t row=matsubara_size+6*mode+first_component;
                const size_t column=matsubara_size+6*partner+second_component;
                block_covariances[block](block_position[row],block_position[column])
                    =real_spectra[36*mode+6*first_component+second_component];
            }
    }

    std::vector<FrequencyBlockFactorization::Block> factored;
    factored.reserve(blocks.size());
    size_t total_rank{};
    RealType factor_target_square{},factor_residual_square{};
    size_t largest_dimension{};
    for( size_t block_index=0;block_index<blocks.size();++block_index )
    {
        const auto& rows=blocks[block_index];
        ComplexDynamicMatrix block(rows.size(),rows.size());
        // Copy one triangle into both halves so FFT roundoff is not
        // magnified by the relative check on very small high-frequency blocks.
        for( size_t i=0;i<rows.size();++i )
            for( size_t j=i;j<rows.size();++j )
            {
                const ComplexType value=RealType{0.5}*(
                    block_covariances[block_index](i,j)
                    +block_covariances[block_index](j,i));
                block(i,j)=value;
                block(j,i)=value;
            }
        const RealType block_norm=frobenius_norm(block);
        const auto split=symmetry_factors(block,true);
        total_rank+=split->rank;
        largest_dimension=std::max(largest_dimension,split->largest);
        factor_target_square+=block_norm*block_norm;
        factor_residual_square+=std::pow(split->error*block_norm,2);
        for(const auto& component:split->blocks)
        {
            std::vector<size_t> component_rows;
            for(auto row:component.rows)component_rows.push_back(rows[row]);
            factored.push_back({std::move(component_rows),component.factor});
        }
    }

    FrequencyBlockFactorization result{};
    result.blocks=std::move(factored);
    result.total_rank=total_rank;
    result.reconstruction_error=factor_target_square>RealType{0.}
        ?std::sqrt(factor_residual_square/factor_target_square)
        :std::sqrt(factor_residual_square);
    result.approximation_error=full_square>RealType{0.}
        ?std::sqrt(discarded_square/full_square):std::sqrt(discarded_square);
    result.largest_dimension=largest_dimension;
    return result;
}

}

size_t one_sided_edge_reflection_index( const size_t point,
                                        const size_t num_intervals )
{
    if( num_intervals==0 )
        throw std::invalid_argument("one-sided edge grid needs an interval");
    if( point>num_intervals )
        throw std::out_of_range("one-sided edge-grid point out of range");
    return num_intervals-point;
}

TakagiFactor autonne_takagi( const ComplexDynamicMatrix& Gamma_input )
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
    const ComplexDynamicMatrix reconstructed=result.L*blaze::trans(result.L);
    const RealType denominator=frobenius_norm(Gamma_input);
    result.reconstruction_error=denominator>RealType{0.}
        ?frobenius_norm(Gamma_input-reconstructed)/denominator
        :frobenius_norm(reconstructed);
    return result;
}

TakagiFactor svd_takagi( const ComplexDynamicMatrix& Gamma_input )
{
    if( Gamma_input.rows()!=Gamma_input.columns() )
        throw std::invalid_argument("SVD Takagi factorization requires a square matrix");
    const size_t n=Gamma_input.rows();
    if( transpose_symmetry_error(Gamma_input)>symmetry_tolerance(n) )
        throw std::invalid_argument("SVD Takagi input is not complex symmetric");

    TakagiFactor result{};
    if( n==0 ) return result;
    ComplexDynamicMatrix U,V;
    blaze::DynamicVector<RealType,blaze::columnVector> singular_values;
    blaze::svd(Gamma_input,U,singular_values,V); // Gamma=U diag(s) V; V stores V^H.
    const RealType sigma_max=singular_values.size()?singular_values[0]:RealType{};
    const RealType rank_tolerance=std::numeric_limits<RealType>::epsilon()*RealType{100.}
        *static_cast<RealType>(n)*std::max(RealType{1.},sigma_max);
    size_t rank{};
    while( rank<singular_values.size() && singular_values[rank]>rank_tolerance ) ++rank;
    result.L.resize(n,rank,false);
    result.numerical_rank=rank;
    result.singular_values.assign(singular_values.begin(),singular_values.begin()+rank);

    // D=V^H conj(U) contains the unitary phase rotation within every
    // degenerate singular-value subspace. If D=R R^T, then T=U R is the
    // corresponding Takagi basis and Gamma=T Sigma T^T.
    const RealType degeneracy_tolerance=std::numeric_limits<RealType>::epsilon()
        *RealType{1000.}*static_cast<RealType>(n)
        *std::max(RealType{1.},sigma_max);
    size_t begin{};
    while( begin<rank )
    {
        size_t end=begin+1;
        while( end<rank && std::abs(singular_values[end]-singular_values[begin])
              <=degeneracy_tolerance ) ++end;
        const size_t block_size=end-begin;
        ComplexDynamicMatrix correction(block_size,block_size,ComplexType{});
        if( block_size==1 )
        {
            ComplexType phase{};
            for( size_t k=0;k<n;++k ) phase+=V(begin,k)*std::conj(U(k,begin));
            if( std::abs(phase)<=rank_tolerance )
                throw std::runtime_error("SVD Takagi phase correction is singular");
            correction(0,0)=std::sqrt(phase/std::abs(phase));
        }
        else
        {
            ComplexDynamicMatrix phase_block(block_size,block_size,ComplexType{});
            for( size_t i=0;i<block_size;++i )
                for( size_t j=0;j<block_size;++j )
                    for( size_t k=0;k<n;++k )
                        phase_block(i,j)+=V(begin+i,k)*std::conj(U(k,begin+j));
            phase_block=RealType{0.5}*(phase_block+blaze::trans(phase_block));
            const TakagiFactor phase_factor=autonne_takagi(phase_block);
            if( phase_factor.numerical_rank!=block_size )
                throw std::runtime_error("SVD Takagi degenerate phase block lost rank");
            correction=phase_factor.L;
        }

        const RealType representative_sigma=std::accumulate(
            singular_values.begin()+begin,singular_values.begin()+end,RealType{})
            /static_cast<RealType>(block_size);
        const RealType root=std::sqrt(representative_sigma);
        for( size_t i=0;i<n;++i )
            for( size_t column=0;column<block_size;++column )
            {
                ComplexType value{};
                for( size_t j=0;j<block_size;++j )
                    value+=U(i,begin+j)*correction(j,column);
                result.L(i,begin+column)=root*value;
            }
        begin=end;
    }

    const ComplexDynamicMatrix reconstructed=result.L*blaze::trans(result.L);
    const RealType denominator=frobenius_norm(Gamma_input);
    result.reconstruction_error=denominator>RealType{0.}
        ?frobenius_norm(Gamma_input-reconstructed)/denominator
        :frobenius_norm(reconstructed);
    return result;
}

JointComplexGaussianSampler::ContourFieldSample
JointComplexGaussianSampler::contour_field_from_latent(
    const LatentVector& latent, const bool include_real_gauss_fields )
{
    // Physical-grid samplers do not have a frequency representation from
    // which to evaluate shifted grids. For dense CF4 the trajectory builder
    // obtains the internal nodes by cubic interpolation of this edge field.
    // FFT overrides this method and supplies spectral nodes.
    static_cast<void>(include_real_gauss_fields);
    ContourFieldSample result{};
    result.edge_field=field_from_latent(latent);
    return result;
}

JointComplexGaussianSampler::ContourFieldSample
JointComplexGaussianSampler::draw_contour_field(
    std::mt19937& engine, const bool include_real_gauss_fields )
{
    return contour_field_from_latent(
        draw_latent(engine),include_real_gauss_fields);
}

namespace
{
std::shared_ptr<GaussianBlockFactors> symmetry_factors(
    const ComplexDynamicMatrix& covariance,const bool use_svd )
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
            factor=std::make_shared<TakagiFactor>(use_svd?svd_takagi(block):autonne_takagi(block));
            for(const auto sigma:factor->singular_values)sigma_max=std::max(sigma_max,sigma);
            unique_covariances.push_back(std::move(block));unique_factors.push_back(factor);
        }
        result->largest=std::max(result->largest,indices.size());
        result->blocks.push_back({std::move(indices),factor});
    }
    // Local factorizations retain a superset of the whole-matrix numerical
    // rank. Apply its original global threshold to all blocks consistently.
    const RealType tolerance=std::numeric_limits<RealType>::epsilon()*RealType{100.}
        *static_cast<RealType>((use_svd?1:2)*n)*std::max(RealType{1.},sigma_max);
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
        }
        const auto& target=unique_covariances[k];
        const ComplexDynamicMatrix reconstructed=factor.L*blaze::trans(factor.L);
        const RealType scale=frobenius_norm(target);
        factor.reconstruction_error=scale>RealType{}?frobenius_norm(target-reconstructed)/scale:frobenius_norm(reconstructed);
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
    std::mt19937& engine,const size_t count,const bool gauss )
{
    std::vector<ContourFieldSample> fields;fields.reserve(count);
    for(size_t i=0;i<count;++i)fields.push_back(draw_contour_field(engine,gauss));
    return fields;
}

DenseComplexGaussianSampler::DenseComplexGaussianSampler( const ComplexDynamicMatrix& covariance )
    : m_factors(symmetry_factors(covariance,false)) {}
RealType DenseComplexGaussianSampler::reconstruction_error() const { return m_factors->error; }
size_t DenseComplexGaussianSampler::latent_dimension() const { return m_factors->rank; }
size_t DenseComplexGaussianSampler::size() const { return m_factors->size; }
size_t DenseComplexGaussianSampler::largest_factorization_dimension() const { return m_factors->largest; }
DenseComplexGaussianSampler::LatentVector DenseComplexGaussianSampler::draw_latent(std::mt19937& engine)
{
    LatentVector latent(latent_dimension());
    for(auto& value:latent)value=m_standard_normal(engine);
    return latent;
}
DenseComplexGaussianSampler::FieldVector DenseComplexGaussianSampler::field_from_latent(const LatentVector& latent)
{ return block_field(*m_factors,latent); }
DenseComplexGaussianSampler::FieldVector DenseComplexGaussianSampler::draw(std::mt19937& engine)
{ return field_from_latent(draw_latent(engine)); }
std::vector<JointComplexGaussianSampler::ContourFieldSample>
DenseComplexGaussianSampler::draw_contour_batch(std::mt19937& engine,size_t count,bool)
{ return block_batch(*this,*m_factors,engine,count); }

SVDComplexGaussianSampler::SVDComplexGaussianSampler( const ComplexDynamicMatrix& covariance )
    : m_factors(symmetry_factors(covariance,true)) {}
RealType SVDComplexGaussianSampler::reconstruction_error() const { return m_factors->error; }
size_t SVDComplexGaussianSampler::latent_dimension() const { return m_factors->rank; }
size_t SVDComplexGaussianSampler::size() const { return m_factors->size; }
size_t SVDComplexGaussianSampler::largest_factorization_dimension() const { return m_factors->largest; }
SVDComplexGaussianSampler::LatentVector SVDComplexGaussianSampler::draw_latent(std::mt19937& engine)
{
    LatentVector latent(latent_dimension());
    for(auto& value:latent)value=m_standard_normal(engine);
    return latent;
}
SVDComplexGaussianSampler::FieldVector SVDComplexGaussianSampler::field_from_latent(const LatentVector& latent)
{ return block_field(*m_factors,latent); }
SVDComplexGaussianSampler::FieldVector SVDComplexGaussianSampler::draw(std::mt19937& engine)
{ return field_from_latent(draw_latent(engine)); }
std::vector<JointComplexGaussianSampler::ContourFieldSample>
SVDComplexGaussianSampler::draw_contour_batch(std::mt19937& engine,size_t count,bool)
{ return block_batch(*this,*m_factors,engine,count); }

namespace
{
// Congruence by A*T (or its inverse), without assembling a dense T. Apply
// the same real transformation on rows and columns: never conjugate Gamma.
void transform_weighted_covariance(ComplexDynamicMatrix& matrix,
    const size_t nm,const size_t nr,const std::array<RealType,3>& roots,
    const bool inverse)
{
    const auto transform=[&](ComplexType& m, const size_t sector)
    { m=inverse?m/roots[sector]:m*roots[sector]; };
    const auto pair=[&](ComplexType& p,ComplexType& b)
    {
        const ComplexType first=p,second=b;
        if(inverse)
        {
            p=first/roots[1]+second/roots[2];
            b=first/roots[1]-second/roots[2];
        }
        else
        {
            p=(RealType{0.5}*first+RealType{0.5}*second)*roots[1];
            b=(RealType{0.5}*first-RealType{0.5}*second)*roots[2];
        }
    };
    for(size_t j=0;j<matrix.columns();++j)
    {
        for(size_t i=0;i<nm;++i)transform(matrix(i,j),0);
        for(size_t i=0;i<nr;++i)pair(matrix(nm+i,j),matrix(nm+nr+i,j));
    }
    for(size_t i=0;i<matrix.rows();++i)
    {
        for(size_t j=0;j<nm;++j)transform(matrix(i,j),0);
        for(size_t j=0;j<nr;++j)pair(matrix(i,nm+j),matrix(i,nm+nr+j));
    }
}
}

WeightedDenseComplexGaussianSampler::WeightedDenseComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance,const size_t num_matsubara_intervals,
    const size_t num_real_points,const std::array<RealType,3>& weights)
{
    const size_t n=covariance.rows(),points=n/3;
    if(n!=covariance.columns()||n%3!=0||num_matsubara_intervals>=points
       ||num_real_points==0||(points-num_matsubara_intervals-1)%2!=0
       ||(points-num_matsubara_intervals-1)/2!=num_real_points)
        throw std::invalid_argument("weighted-dense covariance does not match the physical contour grid");
    if(transpose_symmetry_error(covariance)>symmetry_tolerance(n))
        throw std::invalid_argument("weighted-dense covariance is not complex symmetric");
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)
        if(!std::isfinite(std::real(covariance(i,j)))||!std::isfinite(std::imag(covariance(i,j))))
            throw std::invalid_argument("weighted-dense covariance must be finite");
    for(const auto weight:weights)
        if(!std::isfinite(weight)||weight<=RealType{})
            throw std::invalid_argument("Gaussian noise weights must be finite and strictly positive");
    // A common rescaling of W must not change either the ensemble or the
    // absolute floor in the existing Takagi numerical-rank cutoff. Set the
    // largest weight to 2, so canonical (1,2,2) is an orthogonal change of
    // basis and retains the original dense singular-value scale as well.
    const RealType maximum=*std::max_element(weights.begin(),weights.end());
    for(size_t i=0;i<3;++i)
    {
        m_roots[i]=std::sqrt(weights[i]/maximum)*std::sqrt(RealType{2.});
        if(m_roots[i]==RealType{})
            throw std::invalid_argument("Gaussian noise weight ratio underflows working precision");
    }
    m_matsubara_size=3*(num_matsubara_intervals+1);
    m_real_size=3*num_real_points;
    ComplexDynamicMatrix scaled=covariance;
    transform_weighted_covariance(scaled,m_matsubara_size,m_real_size,m_roots,false);
    m_factors=symmetry_factors(scaled,false);

    // Measure rank-truncation and roundoff in the original physical basis.
    // Reconstruct blockwise, avoiding a second full physical factor in memory.
    reset(scaled);
    for(const auto& block:m_factors->blocks)
    {
        const auto& L=block.factor->L;
        const ComplexDynamicMatrix reconstructed=L*blaze::trans(L);
        for(size_t i=0;i<block.rows.size();++i)for(size_t j=0;j<block.rows.size();++j)
            scaled(block.rows[i],block.rows[j])=reconstructed(i,j);
    }
    transform_weighted_covariance(scaled,m_matsubara_size,m_real_size,m_roots,true);
    const RealType scale=frobenius_norm(covariance);
    m_reconstruction_error=scale>RealType{}?frobenius_norm(scaled-covariance)/scale
                                          :frobenius_norm(scaled);
    const RealType tolerance=std::max(RealType{1e-10},symmetry_tolerance(n));
    if(!std::isfinite(m_reconstruction_error)||m_reconstruction_error>tolerance)
        throw std::runtime_error("weighted-dense physical covariance reconstruction failed; reduce the noise weight contrast");
}

size_t WeightedDenseComplexGaussianSampler::latent_dimension() const { return m_factors->rank; }
size_t WeightedDenseComplexGaussianSampler::size() const { return m_factors->size; }
size_t WeightedDenseComplexGaussianSampler::largest_factorization_dimension() const { return m_factors->largest; }
WeightedDenseComplexGaussianSampler::LatentVector
WeightedDenseComplexGaussianSampler::draw_latent(std::mt19937& engine)
{
    LatentVector latent(latent_dimension());
    for(auto& value:latent)value=m_standard_normal(engine);
    return latent;
}
void WeightedDenseComplexGaussianSampler::to_physical_field(FieldVector& field) const
{
    for(size_t i=0;i<m_matsubara_size;++i)field[i]/=m_roots[0];
    for(size_t i=0;i<m_real_size;++i)
    {
        const auto eta=field[m_matsubara_size+i]/m_roots[1];
        const auto kappa=field[m_matsubara_size+m_real_size+i]/m_roots[2];
        field[m_matsubara_size+i]=eta+kappa;
        field[m_matsubara_size+m_real_size+i]=eta-kappa;
    }
}
WeightedDenseComplexGaussianSampler::FieldVector
WeightedDenseComplexGaussianSampler::field_from_latent(const LatentVector& latent)
{
    auto field=block_field(*m_factors,latent);
    to_physical_field(field);
    return field;
}
WeightedDenseComplexGaussianSampler::FieldVector
WeightedDenseComplexGaussianSampler::draw(std::mt19937& engine)
{ return field_from_latent(draw_latent(engine)); }
std::vector<JointComplexGaussianSampler::ContourFieldSample>
WeightedDenseComplexGaussianSampler::draw_contour_batch(std::mt19937& engine,size_t count,bool)
{
    auto fields=block_batch(*this,*m_factors,engine,count);
    for(auto& field:fields)to_physical_field(field.edge_field);
    return fields;
}

struct FFTDenseComplexGaussianSampler::FFTPlans
{
    FFTPlans( const size_t num_matsubara_intervals,
              const size_t num_matsubara_points,
              const size_t embedded_real, const size_t substeps )
        : matsubara(3*num_matsubara_points),real(6*embedded_real)
    {
        const RealType two_pi=RealType{2.}*std::acos(RealType{-1.});
        const RealType offset=std::sqrt(RealType{3.})/RealType{6.};
        const std::array<RealType,2> nodes{RealType{0.5}-offset,RealType{0.5}+offset};
        // Preserve the existing signed-frequency convention, including +Nyquist.
        // Each shifted inverse FFT evaluates the same frequency realization.
        phases.resize(2*substeps);
        for(size_t node=0;node<phases.size();++node)
        {
            const RealType fraction=(static_cast<RealType>(node/2)+nodes[node%2])
                                    /static_cast<RealType>(substeps);
            phases[node].resize(embedded_real);
            for(size_t mode=0;mode<embedded_real;++mode)
            {
                const auto signed_mode=mode<=embedded_real/2?static_cast<std::ptrdiff_t>(mode)
                    :static_cast<std::ptrdiff_t>(mode)-static_cast<std::ptrdiff_t>(embedded_real);
                const RealType angle=two_pi*static_cast<RealType>(signed_mode)*fraction
                    /static_cast<RealType>(embedded_real);
                phases[node][mode]=std::exp(ComplexType{RealType{},angle});
            }
        }
        matsubara_inverse=make_fft_plan(
            num_matsubara_intervals,3,matsubara.data(),FFTW_BACKWARD);
        real_inverse=make_fft_plan(
            embedded_real,6,real.data(),FFTW_BACKWARD);
        if( !matsubara_inverse||!real_inverse )
        {
            destroy_fft_plan(matsubara_inverse);
            destroy_fft_plan(real_inverse);
            matsubara_inverse={};
            real_inverse={};
            throw std::runtime_error("failed to create inverse FFTW plans");
        }
    }
    ~FFTPlans()
    {
        destroy_fft_plan(matsubara_inverse);
        destroy_fft_plan(real_inverse);
    }

    std::vector<ComplexType> matsubara{};
    std::vector<ComplexType> real{};
    std::vector<std::vector<ComplexType>> phases;
    FFTPlan matsubara_inverse{};
    FFTPlan real_inverse{};
};

struct FFTDenseComplexGaussianSampler::FrequencyFactors
{
    struct Block
    {
        std::vector<size_t> rows{};
        std::shared_ptr<RealFactorMatrix> factor{};
        RealFactorMatrix frequency_fields{};
    };

    std::vector<Block> blocks{};
};

FFTDenseComplexGaussianSampler::FFTDenseComplexGaussianSampler(
    const CovarianceSource& covariance,const size_t num_matsubara_intervals,
    const size_t num_real_points,const RealType delta_real_time,
    const RealType cross_frequency_cutoff, const size_t real_time_substeps )
    : m_num_matsubara_intervals(num_matsubara_intervals),
      m_num_matsubara_points(num_matsubara_intervals+1),
      m_num_real_points(num_real_points),
      m_embedded_real_points(2*num_real_points),
      m_real_time_substeps(real_time_substeps),
      m_physical_size(3*(m_num_matsubara_points+2*num_real_points))
{
    if( num_matsubara_intervals==0||num_real_points==0 )
        throw std::invalid_argument("FFT Gaussian sampler needs nonzero contour grids");
    if(num_real_points>=std::numeric_limits<size_t>::max()/6
        ||real_time_substeps>(std::numeric_limits<size_t>::max()/6-1)/num_real_points)
        throw std::invalid_argument("FFT real-time substep grid overflows size_t");
    const RealType marginal_error=physical_marginal_error(
        covariance,num_matsubara_intervals,num_real_points);
    FrequencyBlockFactorization block_factor=factor_frequency_blocks(
        covariance,num_matsubara_intervals,num_real_points,
        delta_real_time,cross_frequency_cutoff);
    m_latent_dimension=block_factor.total_rank;
    m_reconstruction_error=std::max(marginal_error,block_factor.reconstruction_error);
    m_covariance_approximation_error=block_factor.approximation_error;
    m_largest_factorization_dimension=block_factor.largest_dimension;
    m_frequency_factors=std::make_unique<FrequencyFactors>();
    m_frequency_factors->blocks.reserve(block_factor.blocks.size());
    std::unordered_map<const TakagiFactor*,std::shared_ptr<RealFactorMatrix>> packed_factors;
    for( auto& source:block_factor.blocks )
    {
        FrequencyFactors::Block block{};
        block.rows=std::move(source.rows);
        auto& packed=packed_factors[source.factor.get()];
        if(!packed)
        {
            const auto& L=source.factor->L;
            packed=std::make_shared<RealFactorMatrix>(2*L.rows(),L.columns());
            for(size_t j=0;j<L.columns();++j)for(size_t i=0;i<L.rows();++i)
            {
                (*packed)(2*i,j)=std::real(L(i,j));
                (*packed)(2*i+1,j)=std::imag(L(i,j));
            }
        }
        block.factor=packed;
        // Release the complex factors with block_factor after construction;
        // packing changes layout, not asymptotic factor storage.
        m_frequency_factors->blocks.push_back(std::move(block));
    }
    m_fft=std::make_unique<FFTPlans>(
        num_matsubara_intervals,m_num_matsubara_points,m_embedded_real_points,real_time_substeps);
}

FFTDenseComplexGaussianSampler::~FFTDenseComplexGaussianSampler()=default;

FFTDenseComplexGaussianSampler::LatentVector
FFTDenseComplexGaussianSampler::draw_latent( std::mt19937& engine )
{
    LatentVector latent(m_latent_dimension);
    for( auto& value:latent ) value=m_standard_normal(engine);
    return latent;
}

FFTDenseComplexGaussianSampler::FieldVector
FFTDenseComplexGaussianSampler::field_from_latent( const LatentVector& latent )
{
    return contour_field_from_latent(latent,false).edge_field;
}

FFTDenseComplexGaussianSampler::ContourFieldSample
FFTDenseComplexGaussianSampler::contour_field_from_latent(
    const LatentVector& latent, const bool include_real_gauss_fields )
{
    if( latent.size()!=m_latent_dimension )
        throw std::invalid_argument("FFT complex Gaussian latent-state rank mismatch");
    if(include_real_gauss_fields&&m_real_time_substeps==0)
        throw std::invalid_argument("q=0 FFT sampling provides native endpoints only");
    size_t latent_offset{};
    for(auto& block:m_frequency_factors->blocks)
    {
        // Avoid pointer arithmetic on an empty latent vector for zero-rank baths.
        const RealType* input=block.factor->columns()?latent.data()+latent_offset:nullptr;
        multiply_real_factor(*block.factor,input,1,latent.size(),block.frequency_fields);
        latent_offset+=block.factor->columns();
    }
    if(latent_offset!=latent.size())
        throw std::logic_error("FFT complex Gaussian latent layout mismatch");
    return contour_field_from_frequency(0,include_real_gauss_fields);
}

std::vector<FFTDenseComplexGaussianSampler::ContourFieldSample>
FFTDenseComplexGaussianSampler::draw_contour_batch(
    std::mt19937& engine,const size_t count,const bool include_real_gauss_fields)
{
    if(include_real_gauss_fields&&m_real_time_substeps==0)
        throw std::invalid_argument("q=0 FFT sampling provides native endpoints only");
    if(count==0)return {};
    RealFactorMatrix latent(m_latent_dimension,count);
    // Preserve sample-major normal draws, including distribution caching and
    // partial batches. Independent symmetry blocks still use disjoint rows.
    for(size_t sample=0;sample<count;++sample)
        for(size_t i=0;i<m_latent_dimension;++i)
            latent(i,sample)=m_standard_normal(engine);
    size_t offset{};
    for(auto& block:m_frequency_factors->blocks)
    {
        const RealType* input=block.factor->columns()?latent.data()+offset:nullptr;
        multiply_real_factor(*block.factor,input,count,latent.spacing(),block.frequency_fields);
        offset+=block.factor->columns();
    }
    if(offset!=m_latent_dimension)
        throw std::logic_error("FFT complex Gaussian batch latent layout mismatch");
    std::vector<ContourFieldSample> fields;
    fields.reserve(count);
    for(size_t sample=0;sample<count;++sample)
        fields.push_back(contour_field_from_frequency(sample,include_real_gauss_fields));
    return fields;
}

FFTDenseComplexGaussianSampler::ContourFieldSample
FFTDenseComplexGaussianSampler::contour_field_from_frequency(
    const size_t sample,const bool include_real_gauss_fields)
{
    const size_t matsubara_size=3*m_num_matsubara_points;
    const auto scatter_frequency_fields=[&]( const size_t phase_index )
    {
        std::fill(m_fft->matsubara.begin(),m_fft->matsubara.end(),ComplexType{});
        std::fill(m_fft->real.begin(),m_fft->real.end(),ComplexType{});
        for( const auto& block:m_frequency_factors->blocks )
            for( size_t i=0;i<block.rows.size();++i )
            {
                const size_t row=block.rows[i];
                const ComplexType value{block.frequency_fields(2*i,sample),
                                        block.frequency_fields(2*i+1,sample)};
                if( row<matsubara_size )
                {
                    m_fft->matsubara[row]=value;
                    continue;
                }
                const size_t real_row=row-matsubara_size;
                const size_t mode=real_row/6;
                m_fft->real[real_row]=phase_index==m_fft->phases.size()?value
                    :m_fft->phases[phase_index][mode]*value;
            }
    };

    scatter_frequency_fields(m_fft->phases.size());
    execute_fft(m_fft->matsubara_inverse);
    execute_fft(m_fft->real_inverse);
    const RealType matsubara_normalization=RealType{1.}/std::sqrt(
        static_cast<RealType>(m_num_matsubara_intervals));
    const RealType real_normalization=RealType{1.}/std::sqrt(
        static_cast<RealType>(m_embedded_real_points));

    ContourFieldSample result{};
    result.edge_field.resize(m_physical_size,false);
    reset(result.edge_field);
    const size_t transformed_matsubara_size=3*m_num_matsubara_intervals;
    for( size_t i=0;i<transformed_matsubara_size;++i )
        result.edge_field[i]=matsubara_normalization*m_fft->matsubara[i];
    for( size_t i=transformed_matsubara_size;i<matsubara_size;++i )
        result.edge_field[i]=m_fft->matsubara[i];
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t branch=0;branch<2;++branch )
            for( size_t component=0;component<3;++component )
                result.edge_field[
                    matsubara_size+3*(branch*m_num_real_points+t)+component]
                    =real_normalization*m_fft->real[6*t+3*branch+component];

    const size_t q=m_real_time_substeps;
    const size_t num_real_intervals=m_num_real_points-1;
    if( include_real_gauss_fields )
    {
        for(auto& nodes:result.real_gauss_fields)
            nodes.resize(6*num_real_intervals*q,false);
        for(size_t j=0;j<q;++j) for(size_t node=0;node<2;++node)
        {
            scatter_frequency_fields(2*j+node);
            execute_fft(m_fft->real_inverse);
            for(size_t t=0;t<num_real_intervals;++t)
                for(size_t c=0;c<6;++c)
                    result.real_gauss_fields[node][6*(t*q+j)+c]
                        =real_normalization*m_fft->real[6*t+c];
        }
    }
    return result;
}

FFTDenseComplexGaussianSampler::FieldVector
FFTDenseComplexGaussianSampler::draw( std::mt19937& engine )
{
    return field_from_latent(draw_latent(engine));
}

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const CovarianceSource& covariance,const size_t nm,const size_t nr,
    const RealType dt,const RealType cutoff,const size_t real_time_substeps )
{
    return std::make_unique<FFTDenseComplexGaussianSampler>(covariance,nm,nr,dt,cutoff,real_time_substeps);
}

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const std::string& algorithm,const ComplexDynamicMatrix& covariance,
    const size_t num_matsubara_intervals,const size_t num_real_points,
    const RealType delta_real_time,const RealType cross_frequency_cutoff,
    const std::array<RealType,3>& weights,const size_t real_time_substeps )
{
    if(real_time_substeps>1&&algorithm=="svd")
        throw std::invalid_argument("real-time substeps require dense, weighted-dense, or FFT sampling");
    if( algorithm=="dense" )
        return std::make_unique<DenseComplexGaussianSampler>(covariance);
    if( algorithm=="svd" )
        return std::make_unique<SVDComplexGaussianSampler>(covariance);
    if( algorithm=="weighted-dense" )
        return std::make_unique<WeightedDenseComplexGaussianSampler>(
            covariance,num_matsubara_intervals,num_real_points,weights);
    if( algorithm=="fft" )
        return std::make_unique<FFTDenseComplexGaussianSampler>(
            covariance,num_matsubara_intervals,num_real_points,
            delta_real_time,cross_frequency_cutoff,real_time_substeps);
    throw std::invalid_argument(
        "Unknown Gaussian factorization '"+algorithm+"'; use dense, svd, fft, or weighted-dense");
}

}
