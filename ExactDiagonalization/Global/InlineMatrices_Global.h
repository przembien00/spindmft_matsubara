// Define global static matrices
#pragma once

#include"Matrices_Global.h"

// CONST VARIABLES BUILT AT COMPILE TIME
inline const DiagonalOperator SIGMA_0   { { ComplexType(1.,0.), ComplexType(0.,0.) }, 
                                          { ComplexType(0.,0.), ComplexType(1.,0.) } };

inline const Observable SIGMA_X         { { ComplexType(0.,0.), ComplexType(1.,0.) }, 
                                          { ComplexType(1.,0.), ComplexType(0.,0.) } };

inline const Observable SIGMA_Y         { { ComplexType(0.,0.), ComplexType(0.,-1.) }, 
                                          { ComplexType(0.,1.), ComplexType(0.,0.)  } };

inline const DiagonalOperator SIGMA_Z   { { ComplexType(1.,0.), ComplexType(0.,0.)  }, 
                                          { ComplexType(0.,0.), ComplexType(-1.,0.) } };

inline const RealDiagonalOperator SIGMA_0_REAL  { { RealType{1.}, RealType{0.}  }, 
                                                  { RealType{0.}, RealType{1.} } };

inline const RealObservable SIGMA_X_REAL        { { RealType{0.}, RealType{1.}  }, 
                                                  { RealType{1.}, RealType{0.}  } };

inline const RealOperator SIGMA_Y_REAL          { { RealType{0.}, RealType{-1.} }, 
                                                  { RealType{1.}, RealType{0.}  } };

inline const RealDiagonalOperator SIGMA_Z_REAL  { { RealType{1.}, RealType{0.}  }, 
                                                  { RealType{0.}, RealType{-1.} } };


// CONST VARIABLES BUILT EXTERNAL AT RUN TIME
// deprecated: allocation errors in case of large Hilbert spaces
/*
inline DiagonalOperator ZERO{};

inline RealDiagonalOperator ZERO_REAL{};

inline DiagonalOperator IDENTITY{};

inline RealDiagonalOperator IDENTITY_REAL{};
*/