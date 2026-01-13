// =================== CLASS IMPLEMENTATION OF DERIVED CLASS ISO POLARIZED MODEL ===================
// CONSTRUCTOR
ISO_Polarized_Model::ISO_Polarized_Model( const ps::ParameterSpace& pspace ):
    Model( "ISO_Polarized", 'a', 'C', pspace, true ),
    h_z( pspace.model_params.h_z )
{} 

// generate string with compact info, e.g., for filenames
std::string ISO_Polarized_Model::compact_info( const uint num_Digits ) const
{
    return "ISO_Polarized_h_z=" + print::round_value_to_string( h_z, num_Digits );
}

// generate string with compact info about the initial state
std::string ISO_Polarized_Model::init_state_info( const uint num_Digits ) const
{
    return "Polarized_h_z=" + print::round_value_to_string( h_z, num_Digits );
}

// return vector of Parameters
std::vector<Parameter> ISO_Polarized_Model::return_params() const
{
    std::vector<Parameter> p{};
    p.emplace_back( Parameter{"h_z", h_z} );
    return p;
}

// compute the model specific Hamiltonian
void ISO_Polarized_Model::compute_Hamiltonian()
{
    H = blaze::ZeroMatrix<RealType,blaze::rowMajor>( num_HilbertSpaceDimension, num_HilbertSpaceDimension );
    for(uint spin_i = 0; spin_i < num_Spins-1; spin_i++)
    { 
        for(uint spin_j = spin_i+1; spin_j < num_Spins; spin_j++) // spin_i is always in front of spin_j
        {
            if( J( spin_i, spin_j ) != RealType{0.} ) // e.g. for nearest-neighbor interactions
            {
                std::vector<RealObservable> vx{ num_Spins-2, SIGMA_0_REAL };
                std::vector<RealOperator> vy{ num_Spins-2, SIGMA_0_REAL };          // the real version of sigma_y is not hermitian
                std::vector<RealDiagonalOperator> vz{ num_Spins-2, SIGMA_0_REAL };

                // Six * Sjx
                vx.insert( vx.begin()+spin_i, SIGMA_X_REAL );
                vx.insert( vx.begin()+spin_j, SIGMA_X_REAL );
                auto Vx{H};
                init::compute_multi_tensor_product_hermit( Vx, vx );

                // Siy * Sjy
                vy.insert( vy.begin()+spin_i, SIGMA_Y_REAL );
                vy.insert( vy.begin()+spin_j, SIGMA_Y_REAL );
                auto Vy{H};
                init::compute_multi_tensor_product_nonhermit( Vy, vy );

                // Siz * Sjz
                vz.insert( vz.begin()+spin_i, SIGMA_Z_REAL );
                vz.insert( vz.begin()+spin_j, SIGMA_Z_REAL ); 
                auto Vz{H};
                init::compute_multi_tensor_product_diag( Vz, vz );
                H += RealType{0.25} * J( spin_i, spin_j ) * ( Vx - Vy + Vz ); 
                /* the factor 0.25 is due to the conversion sigma->S
                the rotation is diagonal, so only three terms are occurring here 
                the minus sign is due to using the real version of sigma_y*/
            }
        }
    }
}

// compute the density operator (which is real)
void ISO_Polarized_Model::compute_density_operator()
{
    // compute Sz_tot:
    RealDiagonalOperator Sz_tot = RealDiagonalOperator{blaze::ZeroMatrix<RealType,blaze::rowMajor>( num_HilbertSpaceDimension, num_HilbertSpaceDimension )};
    for( uint spin_index = 0; spin_index < num_Spins; ++spin_index )
    {
        RealDiagonalOperator Sz_i{};
        init::build_linear_spin_operator_diag( Sz_i, num_Spins, spin_index, SIGMA_Z_REAL );
        Sz_tot += Sz_i;
    }
    RealDiagonalOperator Exponent{ RealType{0.5} * h_z * std::move(Sz_tot) };
    RealDiagonalOperator Exp = func::exponentiate( Exponent );
    rho = RealType{1.} / blaze::trace( Exp ) * Exp;
}

/* diagonalize H (which is real) */
void ISO_Polarized_Model::diagonalize()
{
    func::diagonalize_real( H, omega, K );
}

/* transform any relevant matrices to the diagonal basis 
all calculations can be done with real numbers */
void ISO_Polarized_Model::transform_to_diagonal_basis()
{
#ifndef SPARSE
    // compute the adjungate transformation
    auto K_T = blaze::trans( K );

    // construct Sx, Sy, Sz for the central spin
    RealObservable Sx{}; 
    RealOperator   Sy{}; // the real sigma_y is not hermitian
    RealDiagonalOperator Sz{};
    init::build_linear_spin_operator_hermit( Sx, num_Spins, 0, SIGMA_X_REAL );
    init::build_linear_spin_operator_nonhermit( Sy, num_Spins, 0, SIGMA_Y_REAL ); 
    init::build_linear_spin_operator_diag( Sz, num_Spins, 0, SIGMA_Z_REAL );

    // transform to the diagonal basis
    Sx_D = blaze::declsym(  K_T * Sx * K );
    Sy_D = K_T * Sy * K;                        // bottleneck 
    Sz_D = blaze::declsym(  K_T * Sz * K );
    func::write_to_diagonal( H_D, omega );

    rho_D = blaze::declsym( K_T * rho * K );    // this is quicker because rho is diagonal
    Sx_D_rho_D = Sx_D * rho_D;
    Sz_D_rho_D = Sz_D * rho_D;
#else
    // make K and K_T sparse
    SparseRealOperator K_SPARSE{K};
    RealType truncation = RealType{0.001}/std::sqrt( static_cast<RealType>( num_HilbertSpaceDimension ) );
    func::truncate_sparse( K_SPARSE, truncation );
    SparseRealOperator K_T_SPARSE = blaze::trans( K_SPARSE );

    // construct Sx, Sy, Sz for the central spin
    RealObservable Sx{}; 
    RealOperator   Sy{}; // the real sigma_y is not hermitian
    RealDiagonalOperator Sz_SPARSE{};
    init::build_linear_spin_operator_hermit( Sx, num_Spins, 0, SIGMA_X_REAL );
    init::build_linear_spin_operator_nonhermit( Sy, num_Spins, 0, SIGMA_Y_REAL ); 
    init::build_linear_spin_operator_diag( Sz_SPARSE, num_Spins, 0, SIGMA_Z_REAL );

    // make them sparse
    SparseRealObservable Sx_SPARSE{ Sx };
    SparseRealOperator   Sy_SPARSE{ Sy };
    func::truncate_sparse( Sx_SPARSE, RealType{0.1} );
    func::truncate_sparse( Sy_SPARSE, RealType{0.1} );

    // transform to the diagonal basis
    Sx_D = blaze::declsym(  K_T_SPARSE * Sx_SPARSE * K_SPARSE );
    Sy_D = K_T_SPARSE * Sy_SPARSE * K_SPARSE;                        // bottleneck 
    Sz_D = blaze::declsym(  K_T_SPARSE * Sz_SPARSE * K_SPARSE );
    func::write_to_diagonal( H_D, omega );

    rho_D = blaze::declsym( K_T_SPARSE * rho * K_SPARSE );    // this is quicker because rho is diagonal
    Sx_D_rho_D = Sx_D * rho_D;
    Sz_D_rho_D = Sz_D * rho_D;
#endif
}

/* compute the relevant autocorrelations at a specific time
the required matrices need to be complex here */
void ISO_Polarized_Model::compute_means_and_autocorrelations_at( const RealType& time )
{
    // AUTOCORRELATIONS:
    // compute 1/2 <{Sx(t),Sx(0)}>
    DiagonalOperator H_D_CPLX{H_D};
    DiagonalOperator Uplus  = func::exponentiate( ComplexType{0., time}, H_D_CPLX );
    DiagonalOperator Uminus = func::exponentiate( ComplexType{0.,-time}, H_D_CPLX );
    Operator L = Uplus  * Observable{Sx_D};
    Operator R = Uminus * Operator{Sx_D_rho_D};
    ComplexType c = blaze::trace( L * R ); // efficient trace
    RealType SxSx = std::real( c );
    m_correlations.get_xx().emplace_back( 0.25 * SxSx );

    // compute 1/2 <{Sx(t),Sy(0)}> = - 1/2 <{Sy(t),Sx(0)}> 
    L = Uplus  * Operator{ComplexType{0.,1.} * Sy_D}; // remember the factor i
    // R is above
    c = blaze::trace( L * R );
    RealType SxSy = std::real( c );
    m_correlations.get_xy().emplace_back( 0.25 * SxSy );

    // compute 1/2 <{Sy(t),Sx(0)}> = - 1/2 <{Sy(t),Sx(0)}>
    m_correlations.get_yx().emplace_back( 0. ); // NEEDS TO BE CHANGED, IS NOT ZERO !!!!!!!

    // compute 1/2 <{Sz(t),Sz(0)}>
    L = Uplus  * Observable{Sz_D};
    R = Uminus * Operator{Sz_D_rho_D};
    c = blaze::trace( L * R );
    RealType SzSz = std::real( c );
    m_correlations.get_zz().emplace_back( 0.25 * SzSz );


    // MEANS:
    // compute <Sz(t)>
    // L is above;
    // R = Uminus * Operator{rho_D};
    // c = blaze::trace( L * R );
    // RealType Sz = std::real( c );
    m_means.get_z().emplace_back( 0.5 * std::tanh( h_z/2 ) );
}
