#include"Complex_Gaussian.h"

#include<algorithm>
#include<cmath>
#include<limits>
#include<numeric>
#include<stdexcept>
#include<vector>

#include<fftw3.h>

namespace spinDMFT::Functions
{

namespace
{

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

void validate_fft_covariance( const ComplexDynamicMatrix& physical,
                              const size_t num_matsubara_intervals,
                              const size_t num_real )
{
    const size_t expected_size=3*(num_matsubara_intervals+1)+6*num_real;
    if( physical.rows()!=expected_size||physical.columns()!=expected_size )
        throw std::invalid_argument("FFT Gaussian covariance has incompatible contour dimensions");
}

ComplexType doubled_contour_value_unchecked(
    const ComplexDynamicMatrix& physical,const size_t num_matsubara_intervals,
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

RealType physical_marginal_error( const ComplexDynamicMatrix& physical,
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
        TakagiFactor factor{};
    };

    std::vector<Block> blocks{};
    size_t total_rank{};
    RealType approximation_error{};
    RealType reconstruction_error{};
    size_t largest_dimension{};
};

FrequencyBlockFactorization factor_frequency_blocks(
    const ComplexDynamicMatrix& physical,const size_t num_matsubara_intervals,
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
        TakagiFactor factor=svd_takagi(block);
        total_rank+=factor.numerical_rank;
        largest_dimension=std::max(largest_dimension,rows.size());
        factor_target_square+=block_norm*block_norm;
        factor_residual_square+=std::pow(factor.reconstruction_error*block_norm,2);
        factored.push_back({rows,std::move(factor)});
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

DenseComplexGaussianSampler::DenseComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance )
    : m_factor(autonne_takagi(covariance))
{}

DenseComplexGaussianSampler::LatentVector
DenseComplexGaussianSampler::draw_latent( std::mt19937& engine )
{
    LatentVector latent(m_factor.numerical_rank);
    for( auto& value:latent ) value=m_standard_normal(engine);
    return latent;
}

DenseComplexGaussianSampler::FieldVector
DenseComplexGaussianSampler::field_from_latent( const LatentVector& latent )
{
    if( latent.size()!=m_factor.numerical_rank )
        throw std::invalid_argument("Dense complex Gaussian latent-state rank mismatch");
    return m_factor.L*latent;
}

DenseComplexGaussianSampler::FieldVector
DenseComplexGaussianSampler::draw( std::mt19937& engine )
{
    return field_from_latent(draw_latent(engine));
}

SVDComplexGaussianSampler::SVDComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance )
    : m_factor(svd_takagi(covariance))
{}

SVDComplexGaussianSampler::LatentVector
SVDComplexGaussianSampler::draw_latent( std::mt19937& engine )
{
    LatentVector latent(m_factor.numerical_rank);
    for( auto& value:latent ) value=m_standard_normal(engine);
    return latent;
}

SVDComplexGaussianSampler::FieldVector
SVDComplexGaussianSampler::field_from_latent( const LatentVector& latent )
{
    if( latent.size()!=m_factor.numerical_rank )
        throw std::invalid_argument("SVD complex Gaussian latent-state rank mismatch");
    return m_factor.L*latent;
}

SVDComplexGaussianSampler::FieldVector
SVDComplexGaussianSampler::draw( std::mt19937& engine )
{
    return field_from_latent(draw_latent(engine));
}

struct FFTDenseComplexGaussianSampler::FFTPlans
{
    FFTPlans( const size_t num_matsubara_intervals,
              const size_t num_matsubara_points,
              const size_t embedded_real )
        : matsubara(3*num_matsubara_points),real(6*embedded_real)
    {
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
    FFTPlan matsubara_inverse{};
    FFTPlan real_inverse{};
};

struct FFTDenseComplexGaussianSampler::FrequencyFactors
{
    struct Block
    {
        std::vector<size_t> rows{};
        TakagiFactor factor{};
        LatentVector latent{};
        FieldVector frequency_field{};
    };

    std::vector<Block> blocks{};
};

FFTDenseComplexGaussianSampler::FFTDenseComplexGaussianSampler(
    const ComplexDynamicMatrix& covariance,const size_t num_matsubara_intervals,
    const size_t num_real_points,const RealType delta_real_time,
    const RealType cross_frequency_cutoff )
    : m_num_matsubara_intervals(num_matsubara_intervals),
      m_num_matsubara_points(num_matsubara_intervals+1),
      m_num_real_points(num_real_points),
      m_embedded_real_points(2*num_real_points),
      m_physical_size(3*(m_num_matsubara_points+2*num_real_points))
{
    if( num_matsubara_intervals==0||num_real_points==0 )
        throw std::invalid_argument("FFT Gaussian sampler needs nonzero contour grids");
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
    for( auto& source:block_factor.blocks )
    {
        FrequencyFactors::Block block{};
        block.rows=std::move(source.rows);
        block.factor=std::move(source.factor);
        block.latent.resize(block.factor.numerical_rank,false);
        block.frequency_field.resize(block.rows.size(),false);
        m_frequency_factors->blocks.push_back(std::move(block));
    }
    m_fft=std::make_unique<FFTPlans>(
        num_matsubara_intervals,m_num_matsubara_points,m_embedded_real_points);
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
    if( latent.size()!=m_latent_dimension )
        throw std::invalid_argument("FFT complex Gaussian latent-state rank mismatch");
    const size_t matsubara_size=3*m_num_matsubara_points;
    size_t latent_offset{};
    for( auto& block:m_frequency_factors->blocks )
    {
        for( auto& value:block.latent ) value=latent[latent_offset++];
        if( block.factor.numerical_rank==0 )
            reset(block.frequency_field);
        else
            block.frequency_field=block.factor.L*block.latent;
        for( size_t i=0;i<block.rows.size();++i )
        {
            const size_t row=block.rows[i];
            if( row<matsubara_size ) m_fft->matsubara[row]=block.frequency_field[i];
            else m_fft->real[row-matsubara_size]=block.frequency_field[i];
        }
    }
    if( latent_offset!=latent.size() )
        throw std::logic_error("FFT complex Gaussian latent layout mismatch");
    execute_fft(m_fft->matsubara_inverse);
    execute_fft(m_fft->real_inverse);
    const RealType matsubara_normalization=RealType{1.}/std::sqrt(
        static_cast<RealType>(m_num_matsubara_intervals));
    const RealType real_normalization=RealType{1.}/std::sqrt(
        static_cast<RealType>(m_embedded_real_points));

    FieldVector result(m_physical_size,ComplexType{});
    const size_t transformed_matsubara_size=3*m_num_matsubara_intervals;
    for( size_t i=0;i<transformed_matsubara_size;++i )
        result[i]=matsubara_normalization*m_fft->matsubara[i];
    for( size_t i=transformed_matsubara_size;i<matsubara_size;++i )
        result[i]=m_fft->matsubara[i];
    for( size_t t=0;t<m_num_real_points;++t )
        for( size_t branch=0;branch<2;++branch )
            for( size_t component=0;component<3;++component )
                result[matsubara_size+3*(branch*m_num_real_points+t)+component]
                    =real_normalization*m_fft->real[6*t+3*branch+component];
    return result;
}

FFTDenseComplexGaussianSampler::FieldVector
FFTDenseComplexGaussianSampler::draw( std::mt19937& engine )
{
    return field_from_latent(draw_latent(engine));
}

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const std::string& algorithm,const ComplexDynamicMatrix& covariance,
    const size_t num_matsubara_intervals,const size_t num_real_points,
    const RealType delta_real_time,const RealType cross_frequency_cutoff )
{
    if( algorithm=="dense" )
        return std::make_unique<DenseComplexGaussianSampler>(covariance);
    if( algorithm=="svd" )
        return std::make_unique<SVDComplexGaussianSampler>(covariance);
    if( algorithm=="fft" )
        return std::make_unique<FFTDenseComplexGaussianSampler>(
            covariance,num_matsubara_intervals,num_real_points,
            delta_real_time,cross_frequency_cutoff);
    throw std::invalid_argument(
        "Unknown Gaussian factorization '"+algorithm+"'; use dense, svd, or fft");
}

}
