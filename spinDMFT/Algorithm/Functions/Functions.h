#pragma once

#include<Globals/Types.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Frequency_Multivariate_Gaussian/Frequency_Covariance_Matrix.h>
#include<Frequency_Multivariate_Gaussian/Frequency_Noise_Vectors.h>
#include"../matrices.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace spinDMFT::Functions
{

namespace stda = Standard_Algorithms;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace fmvg = Frequency_Multivariate_Gaussian;
namespace ps = spinDMFT::Parameter_Space;
namespace rtd = spinDMFT::Run_Time_Data;

using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;

void initialize_matrices( const ps::ParameterSpace& );

CorrTen generate_initial_spin_correlations( const ps::ParameterSpace& );

std::tuple<FieldVector, fmvg::FrequencyCovarianceMatrix> self_consistent_equations( const ps::ParameterSpace& pspace, const CorrTen& spin_correlations, FieldVector& spin_expval );

std::tuple<RealType, std::vector<Operator>, std::vector<Operator>> compute_propagators( const ps::ParameterSpace& pspace, std::vector<FieldVector>& noise_vectors, FieldVector& mf_mean );

void compute_S_of_t( const ps::ParameterSpace& pspace, const std::vector<Operator>& propagators, const std::vector<Operator>& propagators_inv, std::vector<Operator>& S_x_of_t, std::vector<Operator>& S_y_of_t, std::vector<Operator>& S_z_of_t );

void compute_spin_correlations( const ps::ParameterSpace& pspace, rtd::RunTimeData& rtdata,
    CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, FieldVector& spin_expval, RealType& Z,
    const std::vector<Operator>& S_x_of_t,
    const std::vector<Operator>& S_y_of_t, 
    const std::vector<Operator>& S_z_of_t );

void MPI_share_results( rtd::RunTimeData& rtdata, CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, FieldVector& spin_expval, RealType& partition_function );

void normalize( rtd::RunTimeData& rtdata, CorrTen& spin_correlations, RealType& partition_function );

};