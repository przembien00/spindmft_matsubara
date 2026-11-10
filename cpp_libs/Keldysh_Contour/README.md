# Single-site Keldysh contour primitives

Shared contour layouts, correlations, interpolation and staggered kernels.
The staged solver uses these directly; the cluster solver's tests use them
as the single-site reference. Compile `Contour_Kernel.cpp` when using the
kernel functions.

These files depend only on common libraries. The existing
`spinDMFT::Contour` namespace is retained for source compatibility.
