#include<cstddef>
constexpr size_t Dim = 3;
#include"header.h"

int main()
{
    std::string project_name = "Tests";
    
    std::string config_file = "Test_1";
    SymmMatrix J{{0.,0.},{0.,0.}};
    Matrix ones{{1.,1.},{1.,1.}};
    SymmMatrixOfMatrix J_mf_sq = blaze::declsym( MatrixOfMatrix{{ones,ones},{ones,ones}} );
    create_config_file( project_name, config_file, J, J_mf_sq );

    config_file = "Test_2";
    J = SymmMatrix{{0.,1.,0.},{1.,0.,0.},{0.,0.,0.}};
    Matrix zero{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    Matrix third{{0.,0.,0.},{0.,0.,0.},{0.,0.,1.}};
    J_mf_sq = blaze::declsym( MatrixOfMatrix{{third, zero, zero}, {zero, third, zero}, {zero, zero, third}} );
    create_config_file( project_name, config_file, J, J_mf_sq );
}