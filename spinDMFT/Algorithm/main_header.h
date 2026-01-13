#include<Globals/Types.h>
#include<Globals/MPI_Types.h>
#include"matrices.h"

#include<Standard_Algorithms/Print_Routines.h>
#include<Standard_Algorithms/Standard_Algorithms.h>
namespace stda = Standard_Algorithms;
namespace print = Print_Routines;

#include<Frequency_Multivariate_Gaussian/Frequency_Covariance_Matrix.h>
#include<Frequency_Multivariate_Gaussian/Frequency_Noise_Vectors.h>
namespace fmvg = Frequency_Multivariate_Gaussian;

#include<Random/Random.h>
namespace rd = Random;

#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
using CorrTen = ten::CorrelationTensor<corr::CorrelationVector>;

#include"Time_Measure/Time_Measure.h"
namespace tmm = Time_Measure;

#include"Parameter_Space/Parameter_Space.h"
namespace ps = spinDMFT::Parameter_Space;

#include"Run_Time_Data/Run_Time_Data.h"
namespace rtd = spinDMFT::Run_Time_Data;

#include"Functions/Functions.h"
namespace func = spinDMFT::Functions;

#include"Storage_Concept/Storage_Concept.h"
namespace stoc = spinDMFT::Storage_Concept;