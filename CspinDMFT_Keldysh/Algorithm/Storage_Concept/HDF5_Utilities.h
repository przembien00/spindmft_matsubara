#pragma once
#include <hdf5.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
namespace CspinDMFT::HDF5_Utilities {
struct Handle {
    hid_t id; herr_t (*close)(hid_t);
    Handle(hid_t v,herr_t (*c)(hid_t)):id(v),close(c) { if(v<0) throw std::runtime_error("HDF5 operation failed"); }
    ~Handle() { close(id); }
    operator hid_t() const { return id; }
    Handle(const Handle&)=delete;
};
inline void check(herr_t status) { if(status<0) throw std::runtime_error("HDF5 read/write failed"); }
inline std::vector<double> read(hid_t f,const char* name,const std::vector<hsize_t>& expected,bool attr=false) {
    Handle object(attr?H5Aopen(f,name,H5P_DEFAULT):H5Dopen2(f,name,H5P_DEFAULT),attr?H5Aclose:H5Dclose);
    Handle space(attr?H5Aget_space(object):H5Dget_space(object),H5Sclose);
    int rank=H5Sget_simple_extent_ndims(space);
    if(rank<0) throw std::runtime_error("invalid HDF5 dataspace");
    std::vector<hsize_t> shape(rank);
    H5Sget_simple_extent_dims(space,shape.data(),nullptr);
    if(shape!=expected) throw std::runtime_error(std::string("wrong shape for ")+name);
    size_t n=1;
    for(auto d:shape) n*=d;
    std::vector<double> values(n);
    check(attr?H5Aread(object,H5T_NATIVE_DOUBLE,values.data()):H5Dread(object,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()));
    for(double x:values) if(!std::isfinite(x)) throw std::runtime_error(std::string("nonfinite ")+name);
    return values;
}
}
