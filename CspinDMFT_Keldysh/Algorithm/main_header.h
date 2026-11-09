#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "matrices.h"
#include "Parameter_Space/Parameter_Space.h"
#include "Contour/Contour_Types.h"
#include "Contour/Contour_Kernel.h"
#include "Mean_Field_Models/Mean_Field_Models.h"
#include "Functions/Functions.h"
#include "Run_Time_Data/Run_Time_Data.h"
#include "Functions/Staged_Solver.h"
#include "Storage_Concept/Storage_Concept.h"
#include "Time_Measure/Time_Measure.h"

#include <Standard_Algorithms/Print_Routines.h>
#include "../../cpp_libs/Time_Measure/Time_Measure.h"
namespace ps = CspinDMFT::Parameter_Space;
namespace rtd = CspinDMFT::Run_Time_Data;
namespace stoc = CspinDMFT::Storage_Concept;
namespace tmm = CspinDMFT::Time_Measure;
namespace timing = ::Time_Measure;
namespace func = CspinDMFT::Functions;
namespace mfm = CspinDMFT::Mean_Field_Models;
namespace print = Print_Routines;
using Model = mfm::Model;
using MatsubaraCorrelations = CspinDMFT::Contour::MatsubaraCorrelations;
using ContourCorrelations = CspinDMFT::Contour::ContourCorrelations;
using ClusterMagnetization = CspinDMFT::Contour::ClusterMagnetization;
using TimeMagnetization = CspinDMFT::Contour::TimeMagnetization;

#include <memory>
namespace contour = CspinDMFT::Contour;
using C = CspinDMFT::C;
using Field = CspinDMFT::Field;
using Sampler = CspinDMFT::Sampler;
