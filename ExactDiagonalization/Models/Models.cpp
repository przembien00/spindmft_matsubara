#include"Models.h"

#include<map>
#include<memory>
#include<numeric>
#include<iomanip>
#include<fstream>
#include<sstream>
#include"Standard_Algorithms/Print_Routines.h"
#include"Standard_Algorithms/Numerics.h"
#include"../Global/InlineMatrices_Global.h"
#include"../Initialization/Initialization.h"
#include"../Functions/Functions.h"

namespace SpinED::Models
{

namespace print = Print_Routines;
namespace init = SpinED::Initialization;
namespace func = SpinED::Functions;

// ======================== CLASS IMPLEMENTATIONS ========================
#define SPARSE
#include"ISO_Disordered.cpp"
#include"ISO_Disordered_Blockwise.cpp"
#include"ISO_Polarized.cpp"
#include"XXZ_Disordered.cpp"
#include"DRF_Disordered.cpp"
//#include"DRF_Disordered_Blockwise.cpp" // needs to be implemented first


};