#pragma once
#include <Complex_Gaussian/Complex_Gaussian.h>
namespace CspinDMFT {
using C = ComplexType;
using Matrix = spinDMFT::Functions::ComplexDynamicMatrix;
using Sampler = spinDMFT::Functions::DenseComplexGaussianSampler;
using Field = Sampler::FieldVector;
}
