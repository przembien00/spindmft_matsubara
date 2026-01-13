#include"Global/RealType_Global.h"
#include"Global/MPI_Global.h"
#include"Global/Matrices_Global.h"
#include"Global/InlineMatrices_Global.h"


#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;
#include<Standard_Algorithms/Standard_Algorithms.h>

#include<spinDMFT_Tensors/Tensors.h>
namespace ten = spinDMFT::Tensors;

#include"Time_Measure/Time_Measure.h"
namespace tmm = Time_Measure;

#include"Parameter_Space/Parameter_Space.h"
namespace ps = SpinED::Parameter_Space;

#include"Run_Time_Data/Run_Time_Data.h"
namespace rtd = SpinED::Run_Time_Data;

#include"Initialization/Initialization.h"
namespace init = SpinED::Initialization;

#include"Functions/Functions.h"
namespace func = SpinED::Functions;

#include"Models/Models.h"
namespace mod = SpinED::Models;

#include"Storage_Concept/Storage_Concept.h"
namespace stoc = SpinED::Storage_Concept;