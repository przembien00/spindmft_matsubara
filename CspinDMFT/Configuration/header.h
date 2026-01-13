/* Dim has to be defined as constexpr size_t to include this header */
#define USE_DOUBLE
#include<iostream>
#include<array>
#include<vector>
#include<string>
#include<HDF5/HDF5_Routines.h>
#include<File_Management/File_Management.h>
#include<Standard_Algorithms/Print_Routines.h>

namespace hdf5r = HDF5_Routines;
namespace fm = File_Management;
namespace print = Print_Routines;

using StringList = std::vector<std::string>;
using VMatrix = std::vector<std::vector<double>>;
using VMatrixOfMatrix = std::vector<std::vector<VMatrix>>;

#include<Physics/Couplings.h>
#include<Physics/Mean_Field_Couplings.h>
#include<Physics/Clusterization.h>
#include<Physics/Spin_Geometries.h>

namespace co = Physics::Couplings;
namespace mfco = Physics::Mean_Field_Couplings;
namespace cl = Physics::Clusterization;
namespace sg = Physics::Spin_Geometries;

using Spin = blaze::StaticVector<double, Dim>;
using Distance = blaze::StaticVector<double, Dim>;
using Direction = blaze::StaticVector<double, Dim>;
using CouplingFunction = std::function<double(const Spin&, const Spin&)>;
using Cluster = std::vector<Spin>;
using Ensemble = std::vector<Spin>;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;
using Matrix = blaze::DynamicMatrix<double,blaze::rowMajor>;
using SymmMatrix = blaze::SymmetricMatrix<Matrix>;
using DistanceFunction = std::function<Distance(const Spin&,const Spin&)>;
using RVector = blaze::DynamicVector<double, blaze::columnVector>;
using MatrixOfVector = blaze::DynamicMatrix<RVector, blaze::rowMajor>;
using SymmMatrixOfVector = blaze::SymmetricMatrix<MatrixOfVector>;
using MatrixOfMatrix = blaze::DynamicMatrix<Matrix,blaze::rowMajor>;
using SymmMatrixOfMatrix = blaze::SymmetricMatrix<MatrixOfMatrix>;
using Constraint = std::function<bool(const Spin&)>;

// create a list of all possible correlations based on num_Spins
IndexPairList create_whole_categories_list( const size_t& num_Spins )
{
    IndexPairList ipairs{};
    for( size_t i=0; i<num_Spins; ++i )
    {
        for( size_t j=i; j<num_Spins; ++j )
        {
            ipairs.emplace_back( IndexPair{i,j} );
        }
    }
    return ipairs;
}

// cast BlazeMatrix to VMatrix
template<typename BlazeMatrix>
VMatrix BlazeMatrix_to_VMatrix( const BlazeMatrix& symmMatrix )
{
    VMatrix vMatrix(symmMatrix.rows(), std::vector<double>(symmMatrix.columns()));
    for( size_t i = 0; i < symmMatrix.rows(); ++i )
    {
        for( size_t j = 0; j < symmMatrix.columns(); ++j )
        {
            vMatrix[i][j] = symmMatrix(i,j);
        }
    }
    return vMatrix;
}

// cast SymmMatrixOfMatrix to VMatrixOfMatrix
VMatrixOfMatrix SymmMatrixOfMatrix_to_VMatrixOfMatrix( const SymmMatrixOfMatrix& symmMatrixOfMatrix )
{
    VMatrixOfMatrix vMatrixOfMatrix(symmMatrixOfMatrix.rows(), std::vector<VMatrix>(symmMatrixOfMatrix.columns()));
    for( size_t i = 0; i < symmMatrixOfMatrix.rows(); ++i )
    {
        for( size_t j = 0; j < symmMatrixOfMatrix.columns(); ++j )
        {
            vMatrixOfMatrix[i][j] = BlazeMatrix_to_VMatrix(symmMatrixOfMatrix(i,j));
        }
    }
    return vMatrixOfMatrix;
}

// create the configuration file to store the computed couplings
void create_config_file( const std::string& project_name, const std::string& config_file, const SymmMatrix& J, 
    const SymmMatrixOfMatrix& J_mf_sq, IndexPairList correlation_categories = IndexPairList{} )
{
    // 1) create file
    std::string total_filename = "../Configuration_Data";
    if( project_name != "" )
    {
        total_filename += "/" + project_name;
        fm::create_folder(total_filename);
    }
    total_filename += "/" + config_file + ".hdf5";
    hid_t file_id = H5Fcreate( total_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // 2) add num_Spins
    hdf5r::store_scalar( file_id, "num_Spins", J.rows() );

    // 3) add J
    hdf5r::store_2D_tensor<double>( file_id, "spin-spin couplings", H5T_IEEE_F64LE, BlazeMatrix_to_VMatrix(J) );

    // 4) add J_mf_sq
    hdf5r::store_4D_tensor<double>( file_id, "correlation weights", H5T_IEEE_F64LE, SymmMatrixOfMatrix_to_VMatrixOfMatrix(J_mf_sq) );

    // 5) add correlation_categories
    if( correlation_categories.size() == 0 )
    {
        correlation_categories = create_whole_categories_list( J.rows() );
    }
    hdf5r::store_scalar( file_id, "num_Categories", correlation_categories.size() );
    hdf5r::store_2D_tensor<uint>( file_id, "correlation_categories", H5T_NATIVE_INT, correlation_categories );

    // 6) close file 
    H5Fclose(file_id);
}

