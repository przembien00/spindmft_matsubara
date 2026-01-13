#pragma once

#include<vector>
#include<string>
#include<functional>
#include<iostream>
#include<blaze/Math.h>
#include"../Global/RealType_Global.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace SpinED::Run_Time_Data
{

namespace ps = SpinED::Parameter_Space;
    
// ================= CLASS DEFINITIONS =================
// RUN TIME DATA CLASS
class RunTimeData
{ 
 public:
    // CONSTRUCTORS
    RunTimeData() = default;
    RunTimeData( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS : General
    uint get_num_SaveDigits() const;
   
 private:
    // IMPORTED FROM PARAMETER SPACE
    // CONST PRIVATE MEMBERS : General
    const uint num_PrintDigits{};
};

// HELPER FUNCTIONS AND OPERATORS
std::ostream& operator<<( std::ostream& os, const RunTimeData& rt_data );


}
