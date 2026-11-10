# Complex Gaussian sampling

Shared dense Takagi factorization and Gaussian sampling, including the
Matsubara marginal and conditional joint latent sampling used by the cluster
and staged Keldysh solvers. Compile `Complex_Gaussian.cpp` and, when needed,
`Joint_Matsubara_PCN.cpp` with the caller's precision definition and LAPACK.

These files depend only on common libraries. The existing
`spinDMFT::Functions` namespace is retained for source compatibility; it does
not imply a dependency on a solver directory.
