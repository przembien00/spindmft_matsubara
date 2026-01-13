// =================== CLASS IMPLEMENTATION OF DERIVED CLASS ISO POLARIZED MODEL ===================
// CONSTRUCTOR
DRF_Disordered_Blockwise_Model::DRF_Disordered_Blockwise_Model( const ps::ParameterSpace& pspace ):
    Model( "DRF_Disordered_Blockwise", 'a', 'A', pspace, true ),
    num_blocks(pspace.num_Spins+1),
    half_blocks(std::ceil(static_cast<double>(num_blocks)/2.))
{} 

// generate string with compact info, e.g., for filenames
std::string DRF_Disordered_Blockwise_Model::compact_info( const uint num_Digits ) const
{
    return "DRF_Disordered_Blockwise";
}

// generate string with compact info about the initial state
std::string DRF_Disordered_Blockwise_Model::init_state_info( const uint num_Digits ) const
{
    return "Disordered";
}

// return vector of Parameters
std::vector<Parameter> DRF_Disordered_Blockwise_Model::return_params() const
{
    std::vector<Parameter> p{};
    return p;
}

// compute the model specific Hamiltonian
void DRF_Disordered_Blockwise_Model::compute_Hamiltonian()
{
    H.reserve(num_blocks);
    H.resize(num_blocks);
    states.resize(num_blocks);

    // create states:
    for(uint block=0; block<num_blocks; ++block)
    {
        std::vector<uint> sites(num_Spins);
        for(uint i=0; i<num_Spins; ++i){sites[i] = i;}
        states[block] =  Numerics::possibilities(sites,block); // block = num_up_spins
    }

    // create Hamiltonian (half of the blocks):
    for(uint block=0; block<half_blocks; ++block)
    {
        uint num_states_for_block = states[block].size();
        H[block] = blaze::ZeroMatrix<RealType,blaze::rowMajor>(num_states_for_block,num_states_for_block);
        
        auto it1 = states_for_block.cbegin();
        uint i=0;
        for(;it1 != states_for_block.cend(); ++it1)
        {
            // diagonal entry (SzSz):
            std::vector<RealType> state_srep(num_Spins,-0.5); // state in sign representation
            for(auto up: *it1){ state_srep[up] = 0.5; }
            for(uint s1=0; s1<num_Spins; ++s1) // sum over all Jij Siz Sjz
            {
                for(uint s2=s1+1; s2<num_Spins; ++s2) // (i<j)
                {
                    H[block](i,i) += state_srep[s1]*state_srep[s2]*J(s1,s2);
                }
            }

            // nondiagonal entries (S+S-):
            uint j=i+1;
            for(auto it2=it1+1; it2 != states_for_block.cend(); ++it2)
            {
                // compare states -> are they distinguished by flipping exactly one spin?
                auto psi_i = *it1;
                auto psi_j = *it2;
                uint up=0;
                while(true) // is safe because the states must be different 
                {
                    if(psi_i[up]!=psi_j[up]){ break; } // first difference between both states found
                    ++up;
                }
                bool flipflop = false;
                uint s1,s2;
                if(up==psi_i.size()-1) // last up spin is different => flipflop possible
                {
                    s1 = psi_i.back();
                    s2 = psi_j.back();
                    flipflop = true;
                }
                else if(psi_i[up]<psi_j[up])
                {
                    s1 = psi_i[up];
                    while(true)
                    {
                        if(up+1==psi_i.size()) // second difference between both states found at the end
                        {
                            s2 = psi_j.back();
                            flipflop = true;
                            break; 
                        } 
                        else if( psi_i[up+1]!=psi_j[up])  // second difference between both states found (not at the end)
                        {
                            if(psi_i[up+1]>psi_j[up]) // otherwise two flipflops would be required
                            {
                                s2 = psi_j[up];
                                flipflop = std::equal(psi_i.begin()+up+1,psi_i.end(),psi_j.begin()+up+1); // check, whether rest of the up-spin sequence matches
                            }
                            break; 
                        }
                        ++up;
                    }
                }
                else
                {
                    s1 = psi_j[up];
                    while(true)
                    {
                        if(up+1==psi_i.size()) // second difference between both states found at the end
                        {
                            s2 = psi_i.back();
                            flipflop = true;
                            break; 
                        } 
                        else if( psi_j[up+1]!=psi_i[up])  // second difference between both states found (not at the end)
                        {
                            if(psi_j[up+1]>psi_i[up]) // otherwise two flipflops would be required
                            {
                                s2 = psi_i[up];
                                flipflop = std::equal(psi_i.begin()+up+1,psi_i.end(),psi_j.begin()+up+1); // check, whether rest of the up-spin sequence matches
                            }
                            break; 
                        }
                        ++up;
                    }
                }
                if( flipflop ) // if the states differ by one flipflop process (Si+Sj-), add -Jij/4 to the Hamiltonian (Jij may be zero)
                {
                    H[block](i,j) -= 0.25 * J(s1,s2);
                }
                ++j;
            }
            ++i;
        }
    }
}

// compute the density operator rho = 1/d
void DRF_Disordered_Blockwise_Model::compute_density_operator()
{
    rho = RealType{1.} / static_cast<RealType>( num_HilbertSpaceDimension );
}

/* diagonalize H (which is real) */
void DRF_Disordered_Blockwise_Model::diagonalize()
{
    omega.resize(half_blocks);
    K.resize(half_blocks);
    for(uint block=0; block<half_blocks; ++block)
    {
        std::cout << "size of block no. " << block << ": " << H[block].rows() << "\n";
        std::cout << "starting diagonalization...\n";
        func::diagonalize_real(H[block], omega[block], K[block]);
    }
    // TODO : The other half of the blocks (for K and omega)
}

/* transform any relevant matrices to the diagonal basis 
all calculations can be done with real numbers */
void DRF_Disordered_Blockwise_Model::transform_to_diagonal_basis()
{
    // Sz
    std::vector<SparseRealOperator> K_SPARSE_blocks{};
    for(uint block=0; block<half_blocks; ++block)
    {
        // create Sz:
        SparseRealOperator Sz_block_SPARSE(H[block].rows(),H[block].rows());
        uint i=0;
        for(const auto& state : states[block]) // only diagonal entries (Sz)
        {
            if(state.size() == 0)
            {
                Sz_block_SPARSE(i,i) = -0.5;
            }
            else
            {
                Sz_block_SPARSE(i,i) = (state[0] == 0) ? 0.5 : -0.5 ;
            }
            ++i;
        }

        // make trafo sparse:
        SparseRealOperator K_block_SPARSE{K[block]};
        //RealType truncation = RealType{0.001}/std::sqrt( static_cast<RealType>( num_HilbertSpaceDimension ) );
        //func::truncate_sparse( K_block_SPARSE, truncation );

        // transform
        auto Sz_D_block = blaze::declsym(  blaze::trans( K_block_SPARSE ) * Sz_block_SPARSE * K_block_SPARSE );
        Sz_D.emplace_back(Sz_D_block);
        K_SPARSE_blocks.emplace_back(K_block_SPARSE);

        // create diagonal Hamiltonian
        RealDiagonalOperator H_D_block{};
        func::write_to_diagonal( H_D_block, omega[block] );
        H_D.emplace_back(H_D_block);
    }

    // make remaining blocks of K sparse:
    for(uint block=half_blocks; block<num_blocks; ++block)
    {
        K_SPARSE_blocks.emplace_back(SparseRealOperator{K[block]});
    }

    // copy all blocks to K_SPARSE:
    SparseRealOperator K_SPARSE(num_HilbertSpaceDimension, num_HilbertSpaceDimension);
    for(uint block=0; block<num_blocks; ++block)
    {
        uint bsize = K_SPARSE_blocks[block].rows();
        // TODO
    }

    // compute total index of first element of each block:
    std::vector<uint> first_index_of_blocks{};
    for(uint block=0; block<num_blocks; ++block)
    {
        first_index_of_blocks.emplace_back(first_index_of_blocks.back()+states[block].size());
    }

    // create S+:
    SparseRealOperator Sp(num_HilbertSpaceDimension, num_HilbertSpaceDimension);
    for(uint block1=0; block1<num_blocks; ++block1)
    {
        uint block2 = block1+1;
        auto b1 = states.cbegin()+block1, b2 = states.cbegin()+block2;
        uint b1size = b1->size(), b2size = b2->size();
        for(uint i=0; i<b1size; ++i)
        {
            for(uint j=0; i<b2size; ++j)
            {
                if( (*b1)[i][0] != 0 && (*b2)[j][0] == 0 ) // check whether first spin can be flipped to up
                {
                    if( std::equal((*b1)[i].cbegin(), (*b1)[i].cend(), (*b2)[j].cbegin()+1) // check, whether rest of the up-spin sequence matches
                    {
                        Sp(first_index_of_blocks[block1]+i, first_index_of_blocks[block2]+j) = 1.;
                    }
                }
            }
        }
    }

    // transform S+:
    Sp_D = blaze::trans( K_SPARSE ) * Sp * K_SPARSE

    // TODO : create H_D also in sparse
    // free K ?
    // unten: bestimme Uplus und Uminus einmalig und nicht in jedem Schritt?
    // Sz wirklich noch in Blockstruktur rechnen? -> macht alles komplizierter
}

/* compute the relevant autocorrelations at a specific time
the required matrices need to be complex here */
void DRF_Disordered_Blockwise_Model::compute_means_and_autocorrelations_at( const RealType& time )
{
    // compute <Sz(t)Sz(0)>
    RealType SzSz{};
    for(uint block=0; block<half_blocks; ++block)
    {
        SparseOperator H_D_CPLX_block{Operator{H_D[block]}};
        SparseOperator Uplus_block  = func::exponentiate_gen( ComplexType{0., time}, H_D_CPLX_block );
        SparseOperator Uminus_block = func::exponentiate_gen( ComplexType{0.,-time}, H_D_CPLX_block );
        SparseOperator L_block = Uplus_block  * SparseObservable{Sz_D[block]};
        SparseOperator R_block = Uminus_block * SparseObservable{Sz_D[block]};
        ComplexType c_block = blaze::trace( L_block * R_block );
        RealType SzSz_block = rho * std::real( c_block );

        if( num_Spins%2==0 && block==num_Spins/2 )
        {
            SzSz += SzSz_block;
        }
        else
        {
            SzSz += 2*SzSz_block;
        }
    }
    m_correlations.get_zz().emplace_back( SzSz );

    // compute <Sx(t)Sx(0)> = 1/2 <S+(t)S+>
    SparseOperator H_D_CPLX{Operator{H_D}};
    SparseOperator Uplus  = func::exponentiate_gen( ComplexType{0., time}, H_D_CPLX );
    SparseOperator Uminus = func::exponentiate_gen( ComplexType{0.,-time}, H_D_CPLX );
    SparseOperator L = Uplus * SparseOperator{Sp_D};
    SparseOperator R = Uminus * SparseOperator{Sp_D};
    ComplexType c_block = blaze::trace( L * R );
    RealType SxSx = rho * 0.5 * std::real( c_block );
    m_correlations.get_xx().emplace_back( SxSx );

    // insert <Sz(t)>
    m_means.get_z().emplace_back( 0. );
}
