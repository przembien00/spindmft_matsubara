#pragma once

#include<hdf5.h>
#include<vector>
#include<string>
#include<numeric>
#include<iterator>
#include<type_traits>
#include <variant>

namespace HDF5_Routines
{

// ================================================
// =========== STORING ROUTINES: SIMPLE ===========
// ================================================
// bools need to be converted, e.g., to uint or string
template<typename T>
hid_t get_H5_Type( const T& value )
{
    if constexpr ( std::is_same<T,char>::value )
    {
        return H5T_C_S1;
    }
    else if constexpr ( std::is_same<T,unsigned int>::value ||  std::is_same<T,long unsigned int>::value || std::is_same<T,int>::value || std::is_same<T,long int>::value || std::is_same<T,size_t>::value )
    {
        return H5T_NATIVE_INT;
    }
    else if constexpr ( std::is_same<T,double>::value )
    {
        return H5T_IEEE_F64LE;
    }
    else if constexpr ( std::is_same<T,float>::value )
    {
        return H5T_IEEE_F32LE;
    }
    else
    {
        throw std::invalid_argument("Unsupported type in " + std::string(__PRETTY_FUNCTION__));
    }
}

template<typename T>
bool store_scalar( const hid_t ID, const std::string& name, const T& value )
{
    auto space_id = H5Screate( H5S_SCALAR );
    auto H5type = get_H5_Type(value);
    auto attr_id = H5Acreate( ID, name.c_str(), H5type, space_id, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Awrite( attr_id, H5type, &value );
    H5Aclose( attr_id );
    H5Sclose( space_id );
    return (status < 0) ? false : true;
}

// assumes the existence of begin() data() and size() for value_list
template<typename T>
bool store_list( const hid_t ID, const std::string& name, const T& value_list )
{
    hsize_t dims[1] = { value_list.size() }; 
    auto space_id = H5Screate_simple( 1, dims, NULL );
    auto H5type = get_H5_Type(*value_list.begin());
    auto attr_id = H5Acreate2( ID, name.c_str(), H5type, space_id, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Awrite(attr_id, H5type, value_list.data());
    H5Aclose( attr_id );
    H5Sclose( space_id );
    return (status < 0) ? false : true;
}

inline bool store_string( const hid_t ID, const std::string& name, const std::string& value )
{
    auto space_id = H5Screate( H5S_SCALAR );
    auto datatype_id = H5Tcopy( H5T_C_S1 );
    H5Tset_size( datatype_id, value.size() + 1 );
    auto attr_id = H5Acreate( ID, name.c_str(), datatype_id, space_id, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Awrite( attr_id, datatype_id, value.c_str() );
    H5Aclose( attr_id );
    H5Tclose( datatype_id );
    H5Sclose( space_id );
    return (status < 0) ? false : true;
}

inline bool store_string_list( const hid_t ID, const std::string& name, const std::vector<std::string>& values )
{
    hsize_t dims[1] = { values.size() }; 
    auto space_id = H5Screate_simple(1, dims, NULL);
    auto datatype_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype_id, H5T_VARIABLE);
    auto dataset_id = H5Dcreate2(ID, name.c_str(), datatype_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    std::vector<const char*> values_data{};
    for( auto& str : values ){ values_data.emplace_back( str.c_str() ); }
    auto status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, values_data.data());
    H5Dclose( dataset_id );
    H5Tclose( datatype_id );
    H5Sclose( space_id );
    return (status < 0) ? false : true;
}

inline std::string bool_to_string( const bool b )
{
    return b ? "yes" : "no";
}

inline std::string none_if_empty( const std::string s )
{
    return (s == "") ? "none" : s;
}


// ================================================
// =========== STORING ROUTINES: TENSORS ==========
// ================================================
/* the tensors are template variables and thus quite general
however, they need to contain the size() method and standard begin/end iterators at all tensor levels */

template <typename ElementType,typename Tensor>
std::vector<ElementType> linearize_2D_tensor_to_vector(const Tensor& tensor)
{
    std::vector<ElementType> linearized{};
    linearized.reserve(tensor.size() * tensor[0].size()); // assumes all subtensors have the same size
    for( auto it = tensor.cbegin(); it != tensor.cend(); ++it )
    {
        std::copy(it->cbegin(),it->cend(),std::back_inserter(linearized));
    }
    return linearized;
}

template <typename ElementType,typename Tensor>
inline bool store_2D_tensor( const hid_t ID, const std::string& name, const hid_t H5_BASE_TYPE, const Tensor& tensor, const std::string& info = "" )
{
    std::vector<hsize_t> tensor_sizes{ tensor.size(), tensor[0].size() };
    size_t tensor_rank = 2;
    auto linearized = linearize_2D_tensor_to_vector<ElementType>(tensor);
    auto dataspace = H5Screate_simple( tensor_rank, tensor_sizes.data(), NULL );
    auto dataset = H5Dcreate( ID, name.c_str(), H5_BASE_TYPE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Dwrite( dataset, H5_BASE_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized.data() );
    if( info != "" ){ store_string( dataset, "info", info ); }
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return (status < 0) ? false : true;
}

template <typename ElementType,typename Tensor>
std::vector<ElementType> linearize_3D_tensor_to_vector(const Tensor& tensor)
{
    std::vector<ElementType> linearized{};
    linearized.reserve(tensor.size() * tensor[0].size() * tensor[0][0].size());  // assumes all subtensors have the same size
    for( auto it1 = tensor.cbegin(); it1 != tensor.cend(); ++it1 )
    {
        for( auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2 )
        {
            std::copy(it2->cbegin(),it2->cend(),std::back_inserter(linearized));
        }
    }
    return linearized;
}

template <typename ElementType,typename Tensor>
inline bool store_3D_tensor( const hid_t ID, const std::string& name, const hid_t H5_BASE_TYPE, const Tensor& tensor, const std::string& info = "" )
{
    std::vector<hsize_t> tensor_sizes{ tensor.size(), tensor[0].size(), tensor[0][0].size() };
    size_t tensor_rank = 3;
    auto linearized = linearize_3D_tensor_to_vector<ElementType>(tensor);
    auto dataspace = H5Screate_simple( tensor_rank, tensor_sizes.data(), NULL );
    auto dataset = H5Dcreate( ID, name.c_str(), H5_BASE_TYPE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Dwrite( dataset, H5_BASE_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized.data() );
    if( info != "" ){ store_string( dataset, "info", info ); }
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return (status < 0) ? false : true;
}

template <typename ElementType,typename Tensor>
std::vector<ElementType> linearize_4D_tensor_to_vector(const Tensor& tensor)
{
    std::vector<ElementType> linearized{};
    linearized.reserve(tensor.size() * tensor[0].size() * tensor[0][0].size() * tensor[0][0][0].size());  // assumes all subtensors have the same size
    for( auto it1 = tensor.cbegin(); it1 != tensor.cend(); ++it1 )
    {
        for( auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2 )
        {
            for( auto it3 = it2->cbegin(); it3 != it2->cend(); ++it3 )
            {
                std::copy(it3->cbegin(),it3->cend(),std::back_inserter(linearized));
            }
        }
    }
    return linearized;
}

template <typename ElementType,typename Tensor>
inline bool store_4D_tensor( const hid_t ID, const std::string& name, const hid_t H5_BASE_TYPE, const Tensor& tensor, const std::string& info = "" )
{
    std::vector<hsize_t> tensor_sizes{ tensor.size(), tensor[0].size(), tensor[0][0].size(), tensor[0][0][0].size() };
    size_t tensor_rank = 4;
    auto linearized = linearize_4D_tensor_to_vector<ElementType>(tensor);
    auto dataspace = H5Screate_simple( tensor_rank, tensor_sizes.data(), NULL );
    auto dataset = H5Dcreate( ID, name.c_str(), H5_BASE_TYPE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    auto status = H5Dwrite( dataset, H5_BASE_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized.data() );
    if( info != "" ){ store_string( dataset, "info", info ); }
    H5Dclose( dataset );
    H5Sclose( dataspace );
    return (status < 0) ? false : true;
}


// ================================================
// ========== IMPORTING ROUTINES: SIMPLE ==========
// ================================================
// assumes value has matching type
template<typename T>
bool import_scalar( const hid_t ID, const std::string& name, T& value )
{
    hid_t attr = H5Aopen( ID, name.c_str(), H5P_DEFAULT );
    hid_t type = H5Aget_type( attr ); // datatype of the data
    auto status = H5Aread( attr, type, &value );
    H5Tclose( type );
    H5Aclose( attr );
    return (status < 0) ? false : true;
}

inline bool import_string( const hid_t ID, const std::string& name, std::string& str )
{
    hid_t attr = H5Aopen( ID, name.c_str(), H5P_DEFAULT );
    hid_t type = H5Aget_type( attr ); // datatype of the data
    hid_t type_mem = H5Tget_native_type(type, H5T_DIR_ASCEND);
    size_t size = static_cast<size_t>(H5Tget_size(type));
    std::vector<char> buffer(size + 1, 0); // ensure null termination
    auto status = H5Aread(attr, type_mem, buffer.data());
    if(status >= 0) { buffer[size] = '\0'; str = std::string(buffer.data()); }
    // char* buffer;
    // auto status = H5Aread(attr, type_mem, &buffer);
    // str = std::string(buffer);
    // H5free_memory(buffer);
    H5Tclose(type_mem);
    H5Tclose( type );
    H5Aclose( attr );
    return (status < 0) ? false : true;
}

// assumes value_list has matching type and methods resize() and data()
template<typename T>
bool import_list( const hid_t ID, const std::string& name, T& value_list )
{
    hid_t attr = H5Aopen( ID, name.c_str(), H5P_DEFAULT );
    hid_t space = H5Aget_space( attr );
    hid_t type = H5Aget_type( attr );
    hsize_t dims[1]; // 1D list
    H5Sget_simple_extent_dims(space, dims, nullptr);
    value_list.resize( dims[0] );
    auto status = H5Aread( attr, type, value_list.data() );
    H5Tclose( type );
    H5Sclose( space );
    H5Aclose( attr );
    return (status < 0) ? false : true;
}


// ================================================
// ========== IMPORTING ROUTINES: TENSORS =========
// ================================================
/* assumes linearized_data and tensor_dimensions has matching type and contains the methods resize() and data()
tensor_dimensions is an optional parameter! */
template <typename T, typename U = std::vector<hsize_t>>
bool import_ND_tensor_linearized( const hid_t ID, const std::string& name, T& linearized_data, U&& tensor_dimensions = U{} )
{
    hid_t dataset = H5Dopen2( ID, name.c_str(), H5P_DEFAULT );
    hid_t space = H5Dget_space( dataset );
    hid_t datatype = H5Dget_type( dataset ); // datatype of the dataset
    int rank = H5Sget_simple_extent_ndims( space );
    std::vector<hsize_t> dims( rank ); // Create a vector to hold the dimensions
    H5Sget_simple_extent_dims( space, dims.data(), nullptr);
    tensor_dimensions.resize( dims.size() );
    std::copy( dims.cbegin(), dims.cend(), tensor_dimensions.begin() );
    hsize_t size = std::accumulate( dims.cbegin(), dims.cend(), hsize_t{1}, []( hsize_t product, const hsize_t& factor ){ return std::move(product)*factor; } );
    linearized_data.resize(size);
    auto status = H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized_data.data() ); // read data
    H5Sclose(space);
    H5Dclose(dataset);
    return (status < 0) ? false : true;
}

};