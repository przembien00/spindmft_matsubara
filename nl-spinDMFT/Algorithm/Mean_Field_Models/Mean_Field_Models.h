#pragma once

#include<string>
#include<vector>
#include<iostream>
#include<HDF5/HDF5_Routines.h>

#include"../Parameter_Space/Parameter_Space.h"

#include<Multivariate_Gaussian/Multivariate_Gaussian_Blocks.h>
#include<Multivariate_Gaussian/Symmetry_Schemes.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Clusters.h>
#include<Physics/Physics.h>

namespace Mean_Field_Models
{

namespace ps = DMFT_parameter_space;
namespace mvgb = Multivariate_Gaussian::Blocks;
namespace mss = Multivariate_Gaussian::Blocks::Symmetry_Schemes;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace clu = Observables::Clusters;
namespace stda = Standard_Algorithms;
namespace ph = Physics;

using VecVec = std::vector<RealSpaceVector>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using CorrTenList = std::vector<CorrTen>;
using CluCorrTen = clu::CorrelationCluster<CorrTen>;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;

// ===============================================================
// ========== ABSTRACT BASE CLASS FOR MEAN FIELD MODELS ==========
// ===============================================================
class MeanFieldModel
{
 protected:
   size_t num_Spins{}; // number of spins
   size_t num_HilbertSpaceDimension{}; // Hilbert Space Dimension (deduced from num_Spins)
   std::vector<std::string> spin_list{}; // list of spin values as string
   std::vector<RealType> spin_float_list{}; // list of spin values as float
   SymmMatrix J{}; // local spin-spin couplings
   ph::SpinModel spinmf_cmodel{}; // coupling between spin and mean-field
   RealType rescale_meanfield{1.0}; // rescales spin-mean-field couplings by a factor
   RealType rescale_spinspin{1.0};  // rescales spin-spin couplings by a factor
   std::string config_file{}; // configuration name
   IndexPairList import_only_categories{}; // correlations to be imported
   IndexPairList compute_only_categories{}; // correlations to be simulated (optionally)
   ph::ChemicalShift chemical_shift{}; // chemical shift (optionally)

 public:
   // DESTRUCTOR
   virtual ~MeanFieldModel() = default;

   // CONSTRUCTOR FUNCTION
   void constructor_function( ps::ParameterSpace & pspace );

   // VIRTUAL AND NON-VIRTUAL METHODS AND OPERATORS
   std::vector<double> read_configuration_from_file( const ps::ParameterSpace & pspace ); // get parameters from hdf5 file
   virtual void interpret_model_specific_parameters( const std::vector<double>& Jsq_linearized ) = 0; // interpret Jsq (model specific)
   virtual void rescale() = 0; // rescale any Parameters
   void hand_back_general_parameters( ps::ParameterSpace & pspace ) const; // hand back the general Parameters to the ParameterSpace

   virtual mvgb::CovarianceMatrixBlocks self_consistency( const CorrTenList& environment_spin_corr ) const = 0;
   virtual void compute_Hamiltonian( SparseObservable& new_Hamiltonian, const VecVec&& Vs_of_t, const char symmetry_type ) const = 0;
};

// ===============================================================
// =========== ABSTRACT DERIVED CLASS FOR MULTI MODELS ===========
// ===============================================================
/* Each spin site includes an individual mean-field that has to be drawn */
class MultiModel : public MeanFieldModel
{
 public:
    // DESTRUCTOR
    virtual ~MultiModel() = default;

    // VIRTUAL AND NON-VIRTUAL METHODS AND OPERATORS
    virtual void compute_Hamiltonian( SparseObservable& new_Hamiltonian, const VecVec&& Vs_of_t, const char symmetry_type ) const override;
};

// ===============================================================
// ========= DERIVED CLASS FOR CORRELATION REPLICA MODEL =========
// ===============================================================
/* Each mean-field is set up by out-of-cluster correlations (that are imported and typically originate from a self-consistency) */
class CorrelationReplicaModel : public MultiModel 
{
 private:
    // PRIVATE MODEL SPECIFIC PARAMETERS
    SymmMatrixOfVector J_mf_sq{}; // mean-field correlation weights, symmetric matrix of vectors (3D tensor)

 public:
    // CONSTRUCTOR
    CorrelationReplicaModel() = default;

    // OVERRIDE METHODS AND OPERATORS
    virtual void interpret_model_specific_parameters( const std::vector<double>& Jsq_linearized ) override;
    virtual void rescale() override;
    virtual mvgb::CovarianceMatrixBlocks self_consistency( const CorrTenList& environment_spin_corr ) const override;
};

}