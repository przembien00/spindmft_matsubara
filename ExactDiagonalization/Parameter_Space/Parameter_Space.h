#pragma once

#include<string>
#include<vector>
#include<iostream>
#include<functional>
#include<hdf5.h>
#include"../Global/RealType_Global.h"
#include"../Global/Matrices_Global.h"


// FORWARD DECLARATION OF MODEL
namespace SpinED::Models
{
    class Model;
}


namespace SpinED::Parameter_Space
{

namespace mod = SpinED::Models;


// ================= CLASS DEFINITIONS =================
// MODEL PARAMETERS STRUCT
struct Model_Parameters
{
    RealType lambda{};
    RealType h_z{};
};

// PARAMETER SPACE CLASS
class ParameterSpace
{
 private:
  // PRIVATE MEMBERS
  int my_rank{};
  int world_size{};

 public:
  // CONSTRUCTORS
  ParameterSpace() = default;
  ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank );

  // PUBLIC METHODS
  uint get_num_Cores() const{ return static_cast<uint>(world_size); }
  std::string return_initial_correlations_src_folder( std::string src_folder_suggestion ) const;
  std::vector<std::string> create_folder_branch_list() const;
  std::string create_file_name( const bool with_extension ) const;
  std::string create_essentials_string() const;
  void read_SpinSystem();

  // PUBLIC MEMBERS
  // ========== model and physical parameters ==========
  std::string couplings_filename{};
  uint num_Spins{};
  uint num_HilbertSpaceDimension{};
  SpinSpinCouplings J{};
  Model_Parameters model_params{};
  std::shared_ptr<mod::Model> spin_model{nullptr};
  RealType rescale{};

  // ========== general numerical parameters ==========
  // ...concerning time discretization
  uint num_TimeSteps{};
  uint num_TimePoints{};
  RealType delta_t{};

  // ========== storing and naming ==========
  std::string information_text{};
  std::string project_name{};
  std::string filename_extension{};
  uint num_PrintDigits{};
  // bool save_first_n{};
};

/*
how to add a parameter:
1) add it in the header file Parameter_Space.h
2) add it in Parameter_Space.cpp, so that it can be adjusted with boost program options
3a) perhaps add it to the essential parameters print function
3b) add it to the filename string if its a physical parameter
4) add it to void HDF5_Storage::store_main in Storage_Concept.cpp
5) include the parameter in the code
*/

}
