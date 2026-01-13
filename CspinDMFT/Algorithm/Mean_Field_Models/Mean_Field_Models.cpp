#include"Mean_Field_Models.h"
#include<map>
#include<memory>
#include<numeric>
#include<iomanip>
#include<fstream>
#include<sstream>

#include<HDF5/HDF5_Routines.h>
namespace hdf5r = HDF5_Routines;

#include<Matrices/Matrices.h>
namespace mat = Matrices;

#include<Physics/Spin.h>
namespace sp = Physics::Spin;

#include"Error_Handling.h"
namespace error = Mean_Field_Models::Error_Handling;

namespace Mean_Field_Models
{
// ========================================================
// ============== CLASS MEAN FIELD MODEL ==================
// ========================================================
// constructor function
void MeanFieldModel::constructor_function( ps::ParameterSpace & pspace )
{
    // import general Parameters from pspace from boost_program_options
    spinmf_cmodel = pspace.spinmf_cmodel;
    rescale_meanfield = pspace.rescale_meanfield;
    rescale_spinspin = pspace.rescale_spinspin;
    config_file = pspace.config_file;

    // import physical parameters from file
    auto Jsq_linearized = read_configuration_from_file( pspace );
    interpret_model_specific_parameters( Jsq_linearized );
    num_HilbertSpaceDimension = sp::return_Hilbert_space_dimension( spin_float_list );

    // rescale if desired
    rescale();

    // hand back newly set parameters
    hand_back_general_parameters( pspace );
}

/* read the configuration data from an hdf5 file 
set num_Spins and J automatically
return the Jsq couplings which are model dependent */
std::vector<double> MeanFieldModel::read_configuration_from_file( const ps::ParameterSpace & pspace )
{
    // 1) open file and group:
    std::string total_filename = "Configuration_Data/";
    if( pspace.project_name != "" )
    {
        total_filename += pspace.project_name + "/";
    }
    total_filename += config_file + ".hdf5";
    hid_t file_id = H5Fopen( total_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( file_id < 0 )
    {
        error::CONFIG_FILE_NOT_FOUND( __PRETTY_FUNCTION__, total_filename );
    }
    
    // a) num_Spins:
    hdf5r::import_scalar( file_id, "num_Spins", num_Spins );

    // b) spin lengths:
    hid_t attr_exists = H5Aexists(file_id, "spin lengths");
    if(!(attr_exists <= 0))
    {
        hdf5r::import_list( file_id, "spin lengths", spin_float_list );
        std::transform( spin_float_list.cbegin(), spin_float_list.cend(), std::back_inserter( spin_list ), []( const double& spin )
        { 
            size_t two_spin = std::lround(2*spin);
            if( two_spin % 2 == 0 ) // integer spin
            {
                return std::to_string(std::lround(spin)); 
            }
            else // half integer spin
            {
                return std::to_string(two_spin) + "/2"; 
            }
        } );
    }
    else // spin lengths not contained in config file -> assume all spins to have S=1/2
    {
        spin_float_list.resize(num_Spins, RealType{0.5});
        spin_list.resize(num_Spins, "1/2");
    }

    // c) spin-spin couplings:
    std::vector<double> J_linearized{};
    hdf5r::import_ND_tensor_linearized( file_id, "spin-spin couplings", J_linearized );
    J.resize( num_Spins );
    for( size_t row = 0; row < num_Spins; ++row ) // write to Matrix
    {
        auto start = J_linearized.begin() + row*num_Spins;
        auto end = start + num_Spins;
        std::copy( std::make_move_iterator(start), std::make_move_iterator(end), J.begin(row) );
    }

    // d) mean-field couplings:
    std::vector<double> Jsq_linearized{};
    hdf5r::import_ND_tensor_linearized( file_id, "correlation weights", Jsq_linearized );

    // e) correlation_categories:
    attr_exists = H5Aexists(file_id, "num_Categories");
    if(!(attr_exists <= 0))
    {
        std::vector<uint> categories_linearized{};
        hdf5r::import_ND_tensor_linearized( file_id, "correlation_categories", categories_linearized );
        for( size_t c = 0; c < categories_linearized.size()/2; ++c ) // write to Matrix
        {
            correlation_categories.emplace_back(IndexPair{categories_linearized[2*c],categories_linearized[2*c+1]});
        }
    }

    // close file:
    H5Fclose( file_id );

    return Jsq_linearized;
}

// hand back the model-independent parameters to the ParameterSpace
void MeanFieldModel::hand_back_general_parameters( ps::ParameterSpace & pspace ) const
{
    pspace.num_Spins = num_Spins;
    pspace.num_HilbertSpaceDimension = num_HilbertSpaceDimension;
    pspace.spin_list = spin_list;
    pspace.spin_float_list = spin_float_list;
    pspace.J = J;
    pspace.correlation_categories = correlation_categories;
}

// ========================================================
// =================== CLASS MULTI MODEL ==================
// ========================================================
void MultiModel::compute_Hamiltonian( SparseObservable& new_Hamiltonian, const VecVec&& Vs_of_t, const char symmetry_type ) const
{
    new_Hamiltonian = H_CLUSTER;
    stda::for_n_each( [&new_Hamiltonian]( const auto& Sxi, const auto& Syi, const auto& Szi, const auto& Vi_of_t )
    {
        new_Hamiltonian += Vi_of_t[0]*Sxi + Vi_of_t[1]*Syi + Vi_of_t[2]*Szi;
    }, S_X_LIST.cbegin(), S_X_LIST.cend(), S_Y_LIST.cbegin(), S_Z_LIST.cbegin(), Vs_of_t.cbegin() );
    mat::truncate_sparse(new_Hamiltonian);
}

// ========================================================
// ========== CLASS CORRELATION REPLICA MODEL =============
// ========================================================
// interpret Jsq
void CorrelationReplicaModel::interpret_model_specific_parameters( const std::vector<double>& Jsq_linearized )
{
    J_mf_sq.resize( num_Spins ); 
    for( size_t i = 0; i < num_Spins; ++i ) // outer (symmetric) matrix : rows
    {
        for( size_t j = i; j < num_Spins; ++j ) // outer (symmetric) matrix : cols
        {
            Matrix J_mf_sq_ij(num_Spins, num_Spins);
            for( size_t k = 0; k < num_Spins; ++k ) // inner matrix : rows
            {
                for( size_t l = 0; l < num_Spins; ++l ) // inner matrix : cols
                {
                    size_t index = l + num_Spins*(k + num_Spins*(j + num_Spins*i));
                    J_mf_sq_ij(k,l) = Jsq_linearized[index]; 
                }
            }
            J_mf_sq(i,j) = std::move( J_mf_sq_ij );
        }
    }
}

// rescale the couplings
void CorrelationReplicaModel::rescale()
{
    // rescale if the factors are NOT very close to 1
    if( abs( rescale_meanfield - RealType{1.0} ) > RealType{pow(10,-6)} ) 
    {
        J_mf_sq = pow(rescale_meanfield,2) * J_mf_sq;
    }
    if( abs( rescale_spinspin - RealType{1.0} ) > RealType{pow(10,-6)} ) 
    {
        J = rescale_spinspin * J;
    }
}

/* this function computes the mean-field correlations from the given environment spin correlations according to 
<V^a_i(t)V^b_j(0)> = sum_{kl} sum_{cd} J_mf_sq_ij_kl * D^{ac} * D^{bd} <S^c_k(t) S^d_l(0)> */
mvgb::CovarianceMatrixBlocks CorrelationReplicaModel::self_consistency( const CluCorrTen& env_spin_corr ) const 
{
    char symmetry_type = env_spin_corr[0].get_symmetry();
    auto num_TimePoints = env_spin_corr[0][0].get_num_TimePoints();

    /* this scheme rotates the environment spin correlations according to the spin model
    R^ab_kl(t,0) := <D*S_k(t) (D*S_l(0))^T>_{ab} =  sum_{cd} D^ac * D^bd * <S^c_k(t) S^d_l(0)>
    where D is the spin-model rotation matrix
    spin_corr represents the correlation tensor of spin correlations <S^c(t) S^d(0)> 
    ipair_new is {a,b} and ipair is {c,d} */
    auto rotation_scheme = [&]( const CorrTen& spin_corr, const IndexPair& ipair_new ) -> Corr
    {
        IndexPairList dirs = ten::get_all_direction_pairs();
        return std::accumulate( dirs.cbegin(), dirs.cend(), Corr{num_TimePoints}, [&]( Corr sum, const IndexPair& ipair )
        {
            auto D = spinmf_cmodel.coupling_matrix;
            return sum + D(ipair_new[0],ipair[0]) * D(ipair_new[1],ipair[1]) * spin_corr(ipair[0],ipair[1]); // D^ac D^bd * <S^c(t) S^d(0)>
        } ); // +=
    };

    // 1.) create the rotated environment spin correlations using the defined transformation scheme
    CluCorrTen env_spin_corr_rotated{ correlation_categories };
    std::transform( env_spin_corr.cbegin(), env_spin_corr.cend(), env_spin_corr_rotated.begin(), [&]( const CorrTen& SkSl )
    {
        return CorrTen{ SkSl, rotation_scheme };
    } );

    // 2.) compute the mean-field correlations
    CluCorrTen mf_CCT{ num_Spins, symmetry_type, num_TimePoints };
    mf_CCT.iterate( [&]( auto& CT, const IndexPair& ij )
    {
        /* compute the mean-field correlations Vi(t)Vj(0) by weighted superposition of the rotated environment spin correlations tensor
        mfav( V^a_i(t) V^b_j(0) ) = sum_{kl} J_mf_sq_ij_kl * R^ab_kl(t,0) */
        CorrTen ViVj_CT{ symmetry_type, num_TimePoints }; 
        env_spin_corr_rotated.iterate( [&]( CorrTen& SkSl_rot, const IndexPair& kl )
        {
            ViVj_CT.iterate2( SkSl_rot, [&]( Corr& ViaVjb, Corr& SkSl_rot_ab, const IndexPair& ab )
            {
                ViaVjb += J_mf_sq(ij[0],ij[1])(kl[0],kl[1]) * SkSl_rot_ab;
            } );
        } );

        CT = std::move( ViVj_CT );
    } );

    // 3.) fill the covariance matrix from scheme
    auto mf_CCT_ptr = std::make_shared<CluCorrTen>(std::move(mf_CCT)); // pointer to mean-field correlations
    mss::CorrelationVectorTensorClusterFillingScheme<const CluCorrTen> scheme{ mf_CCT_ptr, symmetry_type, num_Spins };
    mvgb::CovarianceMatrixBlocks mf_Cov{};
    mf_Cov.fill_from_scheme( scheme );
    return mf_Cov;
}

};