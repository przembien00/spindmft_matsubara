#pragma once

#include<fftw3.h>
#include<cstddef>
#include"../Globals/Types.h"

/* Thin wrapper around FFTW for the real inverse discrete Fourier transform used by the
   frequency-domain mean-field samplers. It maps FFTW's double/float API onto RealType so
   the build precision (USE_FLOAT) automatically selects the matching FFTW variant. */
namespace FFT
{

// ===================== PRECISION-DEPENDENT FFTW BINDINGS =====================
// RealType already equals the FFTW real type (double, or float under USE_FLOAT), so the
// fftw buffers can use RealType* directly. The macro selects the fftw_*/fftwf_* symbols.
#ifdef USE_FLOAT
    using FFTWPlan = fftwf_plan;
    #define SPINDMFT_FFTW( name ) fftwf_##name
#else
    using FFTWPlan = fftw_plan;
    #define SPINDMFT_FFTW( name ) fftw_##name
#endif


// =============================================================================
// ===================== INVERSE REAL DFT (HALF-COMPLEX) =======================
// =============================================================================
/* Inverse real DFT matching the orthonormal-style real Fourier convention used across the
   project (fmvg::FrequencyNoiseVectors::fourier_back_transform and
   spinDMFT::Functions::diag_freq_to_time):

     V_time[t] = V_full[0]/sqrt(N)
               +  sum_{n=1}^{N/2}     V_full[n] cos(2*pi*t*n/N)/sqrt(N/2)
               -  sum_{n=N/2+1}^{N-1} V_full[n] sin(2*pi*t*n/N)/sqrt(N/2)

   evaluated in O(N log N) via FFTW's half-complex -> real transform (FFTW_HC2R) instead of
   the explicit O(N^2) sums. The frequency vector V_full already uses the half-complex
   packing (index 0 = DC, 1..N/2 = cosine coefficients, N/2+1..N-1 = sine coefficients with
   frequency N-n), so the only work is a per-index rescaling into FFTW's HC2R input layout.

   FFTW_HC2R computes (unnormalized inverse), for half-complex input h:
     out[t] = h[0] + 2*sum_{k=1}^{(N-1)/2} [ h[k]cos(2*pi*t*k/N) - h[N-k]sin(2*pi*t*k/N) ]
              + (N even ? h[N/2]*(-1)^t : 0)
   so matching term by term gives the rescaling applied in transform():
     h[0]   = V_full[0]/sqrt(N)
     h[k]   =  V_full[k]  /(2*sqrt(N/2)),   h[N-k] = -V_full[N-k]/(2*sqrt(N/2))   (1<=k<=(N-1)/2)
     h[N/2] =  V_full[N/2]/sqrt(N/2)        (only for even N, the Nyquist mode)

   One plan is created per N and reused across calls (planning is expensive, execution is
   cheap). FFTW planning is not thread-safe and the instance owns its in/out buffers, so use
   one instance per thread/process -- which matches the MPI-per-process sampling here. */
class InverseRealDFT
{
public:
    explicit InverseRealDFT( const size_t N )
        : m_N( N ),
          m_norm0( std::sqrt( static_cast<RealType>( N ) ) ),
          m_normNZ( std::sqrt( static_cast<RealType>( N ) / RealType{2.} ) )
    {
        m_in  = static_cast<RealType*>( SPINDMFT_FFTW( malloc )( sizeof( RealType ) * N ) );
        m_out = static_cast<RealType*>( SPINDMFT_FFTW( malloc )( sizeof( RealType ) * N ) );
        // FFTW_ESTIMATE: cheap deterministic planning that does not clobber the buffers, so
        // the plan can be created once and executed repeatedly on its own m_in/m_out.
        m_plan = SPINDMFT_FFTW( plan_r2r_1d )( static_cast<int>( N ), m_in, m_out,
                                               FFTW_HC2R, FFTW_ESTIMATE );
    }

    ~InverseRealDFT()
    {
        SPINDMFT_FFTW( destroy_plan )( m_plan );
        SPINDMFT_FFTW( free )( m_in );
        SPINDMFT_FFTW( free )( m_out );
    }

    // owns FFTW resources -> non-copyable, non-movable
    InverseRealDFT( const InverseRealDFT& ) = delete;
    InverseRealDFT& operator=( const InverseRealDFT& ) = delete;
    InverseRealDFT( InverseRealDFT&& ) = delete;
    InverseRealDFT& operator=( InverseRealDFT&& ) = delete;

    // transform a single length-N real frequency sequence V_full into the time sequence
    // V_time (also length N); the buffers may not alias.
    void transform( const RealType* V_full, RealType* V_time )
    {
        // rescale V_full into FFTW's half-complex input layout
        m_in[0] = V_full[0] / m_norm0;
        const RealType inv_2normNZ = RealType{1.} / ( RealType{2.} * m_normNZ );
        for( size_t k = 1; k <= ( m_N - 1 ) / 2; ++k )
        {
            m_in[k]         =  V_full[k]         * inv_2normNZ;
            m_in[m_N - k]   = -V_full[m_N - k]   * inv_2normNZ;
        }
        if( m_N % 2 == 0 ) // Nyquist mode exists only for even N
        {
            m_in[m_N / 2] = V_full[m_N / 2] / m_normNZ;
        }

        SPINDMFT_FFTW( execute )( m_plan );

        for( size_t t = 0; t < m_N; ++t ) V_time[t] = m_out[t];
    }

    size_t size() const { return m_N; }

private:
    size_t m_N;
    RealType m_norm0;   // sqrt(N)   -- DC normalization
    RealType m_normNZ;  // sqrt(N/2) -- non-zero-mode normalization
    RealType* m_in{ nullptr };
    RealType* m_out{ nullptr };
    FFTWPlan m_plan{};
};

#undef SPINDMFT_FFTW

}
