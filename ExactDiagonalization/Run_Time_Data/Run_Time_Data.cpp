#include"Run_Time_Data.h"

#include<fstream>
#include<iomanip>
#include<Standard_Algorithms/Print_Routines.h>
#include"RTD_Error_Handling.h"

namespace SpinED::Run_Time_Data
{

namespace print = Print_Routines;
namespace ps = SpinED::Parameter_Space;
namespace rtd = SpinED::Run_Time_Data;
namespace error = rtd::Error_Handling;

// ================= CLASS IMPLEMENTATIONS =================
// CONSTRUCTORS
RunTimeData::RunTimeData( const ps::ParameterSpace& pspace ):
    num_PrintDigits( pspace.num_PrintDigits )
{}

};