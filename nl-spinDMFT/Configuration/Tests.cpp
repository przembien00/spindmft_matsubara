#include<cstddef>
constexpr size_t Dim = 3;
#include"header.h"

int main()
{
    std::string project_name = "Tests";
    
    std::string config_file = "Test_1";
    SymmMatrix J{{0.,0.},{0.,0.}};
    SymmMatrixOfVector J_mf_sq = blaze::declsym( MatrixOfVector{{{1.}, {0.}}, {{0.}, {1.}}} );
    IndexPairList import_only{IndexPair{0,0}}; // to be corrected
    create_config_file( project_name, config_file, J, J_mf_sq, import_only );

    config_file = "Test_2";
    J = SymmMatrix{{0.,1.},{1.,0.}};
    create_config_file( project_name, config_file, J, J_mf_sq, import_only );

    config_file = "Test_3";
    J = SymmMatrix{{0.,1.,1.,1.,1.},{1.,0.,1.,1.,1.},{1.,1.,0.,1.,1.},{1.,1.,1.,0.,1.},{1.,1.,1.,1.,0.}};
    J_mf_sq = blaze::declsym( MatrixOfVector{{{1.}, {0.}, {0.}, {0.}, {0.}}, {{0.}, {1.}, {0.}, {0.}, {0.}}, {{0.}, {0.}, {1.}, {0.}, {0.}}, {{0.}, {0.}, {0.}, {1.}, {0.}}, {{0.}, {0.}, {0.}, {0.}, {1.}}} );
    create_config_file( project_name, config_file, J, J_mf_sq, import_only );
}