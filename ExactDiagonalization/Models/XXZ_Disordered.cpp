// =================== CLASS IMPLEMENTATION OF DERIVED CLASS ISO POLARIZED MODEL ===================
// CONSTRUCTOR
XXZ_Disordered_Model::XXZ_Disordered_Model( const ps::ParameterSpace& pspace ):
    Model( "XXZ_Disordered", 'a', 'B', pspace, true ),
    lambda( pspace.model_params.lambda )
{} 

// generate string with compact info, e.g., for filenames
std::string XXZ_Disordered_Model::compact_info( const uint num_Digits ) const
{
    std::string info = "XXZ_lambda=" + print::remove_zeros(print::round_value_to_string(lambda,num_Digits)) + "_Disordered";
    return info;
}

// generate string with compact info about the initial state
std::string XXZ_Disordered_Model::init_state_info( const uint num_Digits ) const
{
    return "Disordered";
}

// return vector of Parameters
std::vector<Parameter> XXZ_Disordered_Model::return_params() const
{
    std::vector<Parameter> p{};
    p.emplace_back( Parameter{"lambda", lambda} );
    return p;
}

// compute the model specific Hamiltonian
void XXZ_Disordered_Model::compute_Hamiltonian()
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
                H += RealType{0.25} * J( spin_i, spin_j ) * ( Vx - Vy + lambda * Vz ); 
                /* the factor 0.25 is due to the conversion sigma->S
                the rotation is diagonal, so only three terms are occurring here 
                the minus sign is due to using the real version of sigma_y*/
            }
        }
    }
}

// compute the density operator rho = 1/d
void XXZ_Disordered_Model::compute_density_operator()
{
    rho = RealType{1.} / static_cast<RealType>( num_HilbertSpaceDimension );
}

/* diagonalize H (which is real) */
void XXZ_Disordered_Model::diagonalize()
{
    func::diagonalize_real( H, omega, K );
}

/* transform any relevant matrices to the diagonal basis 
all calculations can be done with real numbers */
void XXZ_Disordered_Model::transform_to_diagonal_basis()
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

#endif
}

/* compute the relevant autocorrelations at a specific time
the required matrices need to be complex here */
void XXZ_Disordered_Model::compute_means_and_autocorrelations_at( const RealType& time )
{
    // AUTOCORRELATIONS:
    // compute <Sx(t)Sx(0)>
    DiagonalOperator H_D_CPLX{H_D};
    DiagonalOperator Uplus  = func::exponentiate( ComplexType{0., time}, H_D_CPLX );
    DiagonalOperator Uminus = func::exponentiate( ComplexType{0.,-time}, H_D_CPLX );
    Operator L = Uplus  * Observable{Sx_D};
    Operator R = Uminus * Observable{Sx_D};
    RealType SxSx = rho * std::real( blaze::trace( L * R ) ); // efficient trace
    m_correlations.get_xx().emplace_back( 0.25 * SxSx );

    /*// compute <Sx(t)Sy(0)>
    L = Uplus  * Observable{Sx_D};
    R = Uminus * Operator{ComplexType{0.,1.} * Sy_D};
    RealType SxSy = rho * std::real( blaze::trace( L * R ) ); // efficient trace
    m_correlations.get_xy().emplace_back( 0.25 * SxSy );
    */

    // compute <Sz(t)Sz(0)>
    L = Uplus  * Observable{Sz_D};
    R = Uminus * Observable{Sz_D};
    RealType SzSz = rho * std::real( blaze::trace( L * R ) );
    m_correlations.get_zz().emplace_back( 0.25 * SzSz );

    // MEANS:
    // <Sz(t)> = 0
    m_means.get_z().emplace_back( RealType{0.} );
}
