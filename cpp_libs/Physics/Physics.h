#pragma once

#include<string>
#include<vector>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"P_Error_Handling.h"
#include"../Standard_Algorithms/Print_Routines.h"

namespace Physics
{

using FieldVector = blaze::StaticVector<RealType,3UL,blaze::columnVector>;
using FieldMatrix = blaze::StaticMatrix<RealType,3UL,3UL,blaze::rowMajor>;
using Rotation = blaze::StaticMatrix<RealType,3UL,3UL,blaze::rowMajor>;
using Observable = blaze::HermitianMatrix<blaze::DynamicMatrix<ComplexType,blaze::rowMajor>>;

namespace print = Print_Routines;
namespace error = Error_Handling;

// ============================================================
// ==================== PARAMETER CLASS =======================
// ============================================================
struct Parameter
{
    Parameter( const std::string& name, const RealType& value ) : m_name(name), m_value(value) {};
    bool operator==(const Parameter& p) const { return (m_name==p.m_name && m_value==p.m_value); };
    bool operator!=(const Parameter& p) const { return !(*this==p); };
    std::string m_name{};
    RealType m_value{};
};

// ============================================================
// ======== HEADER FOR DIAGONAL SPIN CORRELATION CLASS ========
// ============================================================
struct DiagonalSpinCorrelation
{
    DiagonalSpinCorrelation() = default;
    DiagonalSpinCorrelation( const std::string& name, const RealType& c_tscale = RealType{1.0}, const RealType& c_periodtime = RealType{1.0} );
    std::vector<Parameter> m_parameters{};
    std::function< RealType(const RealType&) > f; // returns the correlation in dependence of time and normalized to 1
    std::vector<RealType> create_discretization( const RealType& dt, const size_t num_points, const RealType& spin ) const; // returns a discretization of the correlation normalized to the spin length
    std::string m_name{};
};

// ============================================================
// ====== HEADER FOR NON DIAGONAL SPIN CORRELATION CLASS ======
// ============================================================
struct NonDiagonalSpinCorrelation
{
    NonDiagonalSpinCorrelation() = default;
    NonDiagonalSpinCorrelation( const std::string& name, const RealType& c_tscale = RealType{1.0}, const RealType& c_periodtime = RealType{1.0} );
    std::vector<Parameter> m_parameters{};
    std::function< RealType(const RealType&) > f; // returns the correlation in dependence of time and normalized to 1
    std::vector<RealType> create_discretization( const RealType& dt, const size_t num_points, const RealType& spin ) const; // returns a discretization of the correlation normalized to the spin length
    std::string m_name{};
};

// ============================================================
// =============== HEADER FOR SPIN MODEL CLASS ================
// ============================================================
/* the spin model refers to the specific coupling, vSi^T * D * vSj, between the spins such as Heisenberg XXZ and more */
struct SpinModel
{
    SpinModel() = default;
    SpinModel( const std::string& name, const RealType& lambda = RealType{1.0}, const RealType& rho = RealType{1.0} );
    bool operator==(const SpinModel& sm) const { return (std::equal(m_parameters.cbegin(),m_parameters.cend(),sm.m_parameters.cbegin()) && m_name==sm.m_name); };
    bool operator!=(const SpinModel& sm) const { return !(*this==sm); };
    std::vector<Parameter> m_parameters{};
    FieldMatrix coupling_matrix; // stores the coupling matrix D from vSi^T * D * vSj
    std::string compact_info( const size_t num_Digits ) const;
    std::string m_name{};
};

// ============================================================
// ============= HEADER FOR MAGNETIC FIELD CLASS ==============
// ============================================================
/* Static magnetic field, returns a single vector */
struct MagneticField
{
    MagneticField() = default;
    MagneticField( const std::string& name, const RealType& h_abs = RealType{0.0}, 
        const RealType& h_phi = RealType{0.0}, const RealType& h_theta = RealType{0.0} );
    std::vector<Parameter> m_parameters{};
    FieldVector m_h{}; // contains the global magnetic field
    std::string compact_info( const size_t num_Digits ) const;
    std::string m_name{};
};

// ============================================================
// ============= HEADER FOR CHEMICAL SHIFT CLASS ==============
// ============================================================
/* Local chemical shift: magnetic field in z-direction, may vary from site to site, returns a value at each site */
struct ChemicalShift
{
    ChemicalShift() = default;
    ChemicalShift( const std::string& name );
    ChemicalShift( const std::vector<RealType>& field );
    std::string compact_info( const size_t num_Digits ) const;
    std::function<RealType(const size_t)> at{}; // returns the magnetic field in dependence of the site
    std::string m_name{};
};

// ===============================================================
// =================== HEADER FOR NOISE CLASS ====================
// ===============================================================
/* Static Gaussian noise with variance tensor (3x3 dimensions), returns a variance in each direction */
struct Noise
{
  Noise() = default;
  Noise( const std::string& name, const RealType& C );
  std::vector<Parameter> m_parameters{};
  std::function<RealType(size_t,size_t)> m_variance_in{}; // returns the variance in dependence of alpha, beta (in x,y,z)
  std::string compact_info( const size_t num_Digits ) const;
  std::string m_name;
};

// ===============================================================
// ============= HEADER FOR EXTRA INTERACTION CLASS ==============
// ===============================================================
/* Extra interaction term build from the spin matrices Sx, Sy, Sz
an example for this is a quadrupolar interaction term HQ ~ Sz*Sz - 3 vS*vS
returns a single observable */
struct ExtraInteraction
{
  ExtraInteraction() = default;
  ExtraInteraction( const std::string& name, const RealType& strength );
  std::vector<Parameter> m_parameters{};
  std::function<Observable(const Observable&,const Observable&,const Observable&,const Observable&)> m_term{}; // returns the extra interaction term in dependence of Zero, Sx, Sy, Sz
  std::string compact_info( const size_t num_Digits ) const;
  std::string m_name;
};

// ===============================================================
// ============= HEADER FOR LOCAL EXTRA INTERACTION CLASS ==============
// ===============================================================
/* Extra interaction term build from the spin matrices Sx, Sy, Sz
an example for this is a quadrupolar interaction term HQ ~ Sz*Sz - 3 vS*vS
returns an observable in dependence of the provided site */
struct LocalExtraInteraction
{
  LocalExtraInteraction() = default;
  LocalExtraInteraction( const std::string& name, const RealType& strength, const size_t num_spins );
  LocalExtraInteraction( const std::string& name, const std::vector<RealType>& strengths );
  std::vector<Parameter> m_parameters{};
  std::function<Observable(const size_t,const Observable&,const Observable&,const Observable&,const Observable&)> m_term{}; // returns the extra interaction term in dependence of the site and the matrices Zero, Sx, Sy, Sz
  std::string compact_info( const size_t num_Digits ) const;
  std::string m_name;
};



// ============================================================
// ==== IMPLEMENTATION OF DIAGONAL SPIN CORRELATION CLASS =====
// ============================================================
// constructor
inline DiagonalSpinCorrelation::DiagonalSpinCorrelation( const std::string& name, const RealType& c_tscale, const RealType& c_periodtime ):
    m_name( name )
{
    if( m_name == "imagtime" )
    {
        m_parameters.emplace_back( Parameter{ "c_tscale", c_tscale } );
        f = [c_tscale]( const RealType& t ) -> RealType
        {
            return RealType{0.25} * ( RealType{0.25} * (t * t -  c_tscale * t ) + RealType{1.0});
        }; 
    }
    else if( m_name == "const" )
    {
        f = []( const RealType& t ) -> RealType
        {
            return RealType{1.};
        }; 
    }
    else if( m_name == "exponential" )
    {
        m_parameters.emplace_back( Parameter{ "c_tscale", c_tscale } );
        f = [c_tscale]( const RealType& t ) -> RealType
        {
            return std::exp( -RealType{4.0} * t/c_tscale );
        }; 
    }
    else if( m_name == "Gaussian" )
    {
        m_parameters.emplace_back( Parameter{ "c_tscale", c_tscale } );
        f = [c_tscale]( const RealType& t ) -> RealType
        {
            return std::exp( -RealType{8.0} * std::pow(t/c_tscale,2) );
        }; 
    }
    else if( m_name == "cos" )
    {
        m_parameters.emplace_back( Parameter{ "c_periodtime", c_periodtime } );
        f = [c_periodtime]( const RealType& t ) -> RealType
        {
            return std::cos( 2*M_PI/c_periodtime * t );
        };
    }
    else
    {
        error::CORRELATION( m_name, __PRETTY_FUNCTION__ );
    }
}

// evaluate the function at discrete time points
inline std::vector<RealType> DiagonalSpinCorrelation::create_discretization( const RealType& dt, const size_t num_points, const RealType& spin ) const
{
    std::vector<RealType> correlation{};
    RealType max_value = spin*(spin+1.) / 3.; // S(S+1)/3 is the value of av( S^alpha^2 )
    for( size_t t_index = 0; t_index < num_points; ++t_index )
    {
        correlation.emplace_back( max_value * f(static_cast<RealType>(t_index)*dt) );
    }
    return correlation;
}

// ============================================================
// == IMPLEMENTATION OF NON DIAGONAL SPIN CORRELATION CLASS ===
// ============================================================
// constructor
inline NonDiagonalSpinCorrelation::NonDiagonalSpinCorrelation( const std::string& name, const RealType& c_tscale, const RealType& c_periodtime ):
    m_name( name )
{
    if( m_name == "zero" )
    {
        f = []( const RealType& t ) -> RealType
        {
            return RealType{0.};
        }; 
    }
    else if( m_name == "sin" )
    {
        m_parameters.emplace_back( Parameter{ "c_periodtime", c_periodtime } );
        f = [c_periodtime]( const RealType& t ) -> RealType
        {
            return std::sin( 2*M_PI/c_periodtime * t );
        };
    }
    else
    {
        error::CORRELATION( m_name, __PRETTY_FUNCTION__ );
    }
}

// evaluate the function at discrete time points
inline std::vector<RealType> NonDiagonalSpinCorrelation::create_discretization( const RealType& dt, const size_t num_points, const RealType& spin ) const
{
    std::vector<RealType> correlation{};
    RealType max_value = spin*(spin+1.) / 3.; // S(S+1)/3 is the value of av( S^alpha^2 )
    for( size_t t_index = 0; t_index < num_points; ++t_index )
    {
        correlation.emplace_back( max_value * f(static_cast<RealType>(t_index)*dt) );
    }
    return correlation;
}

// ============================================================
// ============ IMPLEMENTATION OF SPIN MODEL CLASS ============
// ============================================================
// constructor
inline SpinModel::SpinModel( const std::string& name, const RealType& lambda, const RealType& rho ):
    m_name( name )
{
    if( name == "ISO" ) // Six Sjx + Siy Sjy + Siz Sjz
    {
        coupling_matrix = FieldMatrix{ {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} };
    }
    else if( name == "Ising" ) // Siz Sjz
    {
        coupling_matrix = FieldMatrix{ {0., 0., 0.}, {0., 0., 0.}, {0., 0., 1.} };
    }
    else if( name == "XXZ" ) // Six Sjx + Siy Sjy + lambda Siz Sjz
    {
        m_parameters.emplace_back( Parameter{ "lambda", lambda } );
        coupling_matrix = FieldMatrix{ {1., 0., 0.}, {0., 1., 0.}, {0., 0., lambda} };
    }
    else if( name == "DRF" ) // homonuclear Dipoles in the Rotating Frame:  -1/2 (Six Sjx + Siy Sjy) + Siz Sjz
    {
        coupling_matrix = FieldMatrix{ {-0.5, 0., 0.}, {0., -0.5, 0.}, {0., 0., 1.} };
    }
    else if( name == "XY" ) // Six Sjx + Siy Sjy + rho (Six Sjy + Siy Sjx) + lambda Siz Sjz  ...not 100% sure about this
    {
        m_parameters.emplace_back( Parameter{ "lambda", lambda } );
        m_parameters.emplace_back( Parameter{ "rho", rho } );
        coupling_matrix = FieldMatrix{ {1., rho, 0.}, {rho, 1., 0.}, {0., 0., lambda} };
    }
    else if( name == "het" ) // heteronuclear dipoles in the rotating frame: 2 Siz Sjz
    {
        coupling_matrix = FieldMatrix{ {0., 0., 0.}, {0., 0., 0.}, {0., 0., 2.} };
    }
    else if( name == "hom" ) // homonuclear dipoles in the rotating frame:  2 Siz Sjz - (Six Sjx + Siy Sjy) 
    {
        coupling_matrix = FieldMatrix{ {-1., 0., 0.}, {0., -1., 0.}, {0., 0., 2.} };
    }
    else
    {
        error::SPIN_MODEL( name, __PRETTY_FUNCTION__ );
    }
}

// create string with compact info
inline std::string SpinModel::compact_info( const size_t num_Digits ) const
{
    std::string info = "spinmodel=" + m_name;
    for( const auto& param : m_parameters )
    {
        info += "_" + param.m_name + "=" + print::remove_zeros(print::round_value_to_string( param.m_value, num_Digits ));
    }
    return info;
}   

// ============================================================
// ========= IMPLEMENTATION OF MAGNETIC FIELD CLASS ===========
// ============================================================
// constructor
inline MagneticField::MagneticField( const std::string& name, const RealType& h_abs, 
    const RealType& h_phi, const RealType& h_theta ):
    m_name( name )
{
    if( m_name == "none" )
    {
        m_h[0] = RealType{0.};
        m_h[1] = RealType{0.};
        m_h[2] = RealType{0.};
    }
    else if( m_name == "x" )
    {
        m_parameters.emplace_back( Parameter{ "h_abs", h_abs } );
        m_h[0] = h_abs;
        m_h[1] = RealType{0.};
        m_h[2] = RealType{0.};
    }
    else if( m_name == "y" )
    {
        m_parameters.emplace_back( Parameter{ "h_abs", h_abs } );
        m_h[0] = RealType{0.};
        m_h[1] = h_abs;
        m_h[2] = RealType{0.};
    } 
    else if( m_name == "z" )
    {
        m_parameters.emplace_back( Parameter{ "h_abs", h_abs } );
        m_h[0] = RealType{0.};
        m_h[1] = RealType{0.};
        m_h[2] = h_abs;
    }
    else if( m_name == "arb" )
    {
        m_parameters.emplace_back( Parameter{ "h_abs", h_abs } );
        m_parameters.emplace_back( Parameter{ "h_phi", h_phi } );
        m_parameters.emplace_back( Parameter{ "h_theta", h_theta } );
        m_h[0] = h_abs * std::sin(h_theta) * std::cos(h_phi);
        m_h[1] = h_abs * std::sin(h_theta) * std::sin(h_phi);
        m_h[2] = h_abs * std::cos(h_theta);
    }
    else
    {
        error::MAGNETIC_FIELD( m_name, __PRETTY_FUNCTION__ );
    }
};

// create string with compact info
inline std::string MagneticField::compact_info( const size_t num_Digits ) const
{
    std::string info = "h=" + m_name;
    for( auto& param : m_parameters )
    {
        info += "_" + param.m_name + "=" + print::remove_zeros(print::round_value_to_string( param.m_value, num_Digits ));
    }
    return info;
}   

// ============================================================
// ======== IMPLEMENTATION OF CHEMICAL SHIFT CLASS ============
// ============================================================
// constructor 1
inline ChemicalShift::ChemicalShift( const std::string& name ):
    m_name( name )
{
    if( m_name == "none" )
    {
        at = []( const size_t& site ){ return RealType{0.}; };
    }
    else
    {
        auto str_list = print::split_string_at_delimiter(m_name,',');
        std::vector<RealType> field{};
        std::transform( str_list.cbegin(), str_list.cend(), std::back_inserter(field), [](const auto& str){return static_cast<RealType>(std::stod(str));} );
        at = [field]( const size_t& site ){ return field[site]; };
    }
};

// constructor 2
inline ChemicalShift::ChemicalShift( const std::vector<RealType>& field ):
    m_name( "auto" )
{
    at = [field]( const size_t& site ){ return field[site]; };
};

// create string with compact info
inline std::string ChemicalShift::compact_info( const size_t num_Digits ) const
{
    std::string info = "chemshift=" + m_name;
    return info;
}

// ===============================================================
// ============== IMPLEMENTATION OF NOISE CLASS ==================
// ===============================================================
// constructor
inline Noise::Noise( const std::string& name, const RealType& C ):
    m_name( name )
{
    if( name == "none" )
    {
        m_variance_in = []( const size_t dir1, const size_t dir2 ){ return RealType{0.}; };
    }
    else if( name == "z" )
    {
        m_parameters.emplace_back( Parameter{ "C", C } );
        m_variance_in = [C]( const size_t dir1, const size_t dir2 )
        { 
            return ( dir1 == 2 && dir2 == 2 ) ? C : RealType{0.};
        };
    }
    else 
    {
        error::NOISE( m_name, __PRETTY_FUNCTION__ );
    }
}

// return compact info about the noise
inline std::string Noise::compact_info( const size_t num_Digits ) const
{
    std::string info = "noise=" + m_name;
    for( const auto& param : m_parameters )
    {
        info += "_" + param.m_name + "=" + print::remove_zeros(print::round_value_to_string( param.m_value, num_Digits ));
    }
    return info;
}

// ===============================================================
// ========= IMPLEMENTATION OF EXTRA INTERACTION CLASS ===========
// ===============================================================
// constructor
inline ExtraInteraction::ExtraInteraction( const std::string& name, const RealType& strength ):
    m_name( name )
{
    if( name == "none" )
    {
        m_term = []( const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz ){ return Zero; };
    }
    else if( name == "quadrupolar" )
    {
        m_parameters.emplace_back( Parameter{ "strength", strength } );
        m_term = [strength]( const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz )
        { 
            // return strength * ( RealType{2.}*Sz*Sz - Sx*Sx - Sy*Sy );
            return strength * RealType{3.} * Sz*Sz;
        };
    }
    else 
    {
        error::EXTRAINTERACTION( m_name, __PRETTY_FUNCTION__ );
    }
}

// return compact info about the extra interaction
inline std::string ExtraInteraction::compact_info( const size_t num_Digits ) const
{
    std::string info = "extraint=" + m_name;
    for( const auto& param : m_parameters )
    {
        info += "_" + param.m_name + "=" + print::remove_zeros(print::round_value_to_string( param.m_value, num_Digits ));
    }
    return info;
}

// ===============================================================
// ========= IMPLEMENTATION OF LOCAL EXTRA INTERACTION CLASS ===========
// ===============================================================
// constructor for same strength at all sites
inline LocalExtraInteraction::LocalExtraInteraction( const std::string& name, const RealType& strength, const size_t num_spins ):
    m_name( name )
{
    if( name == "none" )
    {
        m_term = []( const size_t i, const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz ){ return Zero; };
    }
    else if( name == "quadrupolar" )
    {
        m_parameters.emplace_back( Parameter{ "strengths (const)", strength } );
        m_term = [strength]( const size_t i, const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz )
        { 
            // return strength * ( RealType{2.}*Sz*Sz - Sx*Sx - Sy*Sy );
            return strength * RealType{3.} * Sz*Sz;
        };
    }
    else 
    {
        error::EXTRAINTERACTION( m_name, __PRETTY_FUNCTION__ );
    }
}

// constructor for varying strengths
inline LocalExtraInteraction::LocalExtraInteraction( const std::string& name, const std::vector<RealType>& strengths ):
    m_name( name )
{
    if( name == "none" )
    {
        m_term = []( const size_t i, const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz ){ return Zero; };
    }
    else if( name == "quadrupolar" )
    {
        // the strength values are not stored as parameters
        m_term = [strengths]( const size_t i, const Observable& Zero, const Observable& Sx, const Observable& Sy, const Observable& Sz )
        { 
            // return strength * ( RealType{2.}*Sz*Sz - Sx*Sx - Sy*Sy );
            return strengths[i] * RealType{3.} * Sz*Sz;
        };
    }
    else 
    {
        error::EXTRAINTERACTION( m_name, __PRETTY_FUNCTION__ );
    }
}

// return compact info about the extra interaction
inline std::string LocalExtraInteraction::compact_info( const size_t num_Digits ) const
{
    std::string info = "extraint=" + m_name;
    for( const auto& param : m_parameters )
    {
        info += "_" + param.m_name + "=" + print::remove_zeros(print::round_value_to_string( param.m_value, num_Digits ));
    }
    return info;
}

};