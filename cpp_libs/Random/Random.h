#pragma once

#include<random>
#include<cmath>

namespace Random
{

// generates an integer seed for the Mersenne Twister Generator
inline size_t generate_seed( std::string seed_str, const size_t my_rank )
{

    if( seed_str.find("random") != std::string::npos ) // random seed
    {
        size_t seed{};
        if( my_rank == 0 ) // draw seed on rank 0 
        {
            std::random_device my_seed{};
            seed = my_seed();
        }
        MPI_Bcast( &seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); // broadcast seed to all the other ranks
        return seed;
    }  
    else // preset seed
    {
        return std::stoi( seed_str );
    }    
}


/* generates an iteration and rank dependent seed for the Mersenne Twister Generator
the local seed is build from the formula |seed - (iteration + 10^5) * rank| ensuring different seeds in each iteration step
equivalent seeds on certain cores can only happen for computations on more than a 10^5 cores (in this case the statistical error would be increased) */
inline size_t throw_seed( const size_t seed, const size_t num_Iteration, const size_t my_rank )
{
    return (size_t) std::abs( static_cast<int>(seed) - static_cast<int>((num_Iteration + pow(10,5)) * my_rank) );
}

};