#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <rtseis/utilities/transforms/utilities.hpp>
#include <rtseis/utilities/transforms/dftRealToComplex.hpp>
#include "acousticGreens2D.hpp"
#include "rickerWavelet.hpp"
#include "sourceTimeFunction.hpp"

namespace DFTUtilities = RTSeis::Utilities::Transforms::DFTUtilities;

class AcousticGreens2D::AcousticGreens2DImpl
{
public:
    /// Source time function
    std::unique_ptr<ISourceTimeFunction> mSTF;
    /// Receiver (x, y) position in meters
    std::vector<std::pair<double,double>> mReceiver;
    /// Source (x, y) position in meters
    std::vector<std::pair<double,double>> mSource;
    /// Green's functions
    std::vector<std::vector<double>> mGreens;
    /// Velocity (m/s)
    double mVelocity = 0;
    /// Density (kg/m**3)
    double mDensity = 0;
    /// Quality factor
    double mQ = 9999;
};

/// Constructor
AcousticGreens2D::AcousticGreens2D() :
    pImpl(std::make_unique<AcousticGreens2DImpl> ())
{
}

/// Copy c'tor
AcousticGreens2D::AcousticGreens2D(const AcousticGreens2D &greens)
{
    *this = greens;
}

/// Move c'tor
AcousticGreens2D::AcousticGreens2D(AcousticGreens2D &&greens) noexcept
{
    *this = std::move(greens);
}

/// Copy assignment
AcousticGreens2D& AcousticGreens2D::operator=(const AcousticGreens2D &greens)
{
    if (&greens == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl->mSTF = greens.pImpl->mSTF->clone();
    pImpl->mReceiver = greens.pImpl->mReceiver;
    pImpl->mSource = greens.pImpl->mSource;
    pImpl->mVelocity = greens.pImpl->mVelocity;
    pImpl->mDensity = greens.pImpl->mDensity;
    pImpl->mQ = greens.pImpl->mQ;
    return *this;
}

/// Move assignment
AcousticGreens2D&
AcousticGreens2D::operator=(AcousticGreens2D &&greens) noexcept
{
    if (&greens == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(greens.pImpl);
    return *this;
}

/// Destructor
AcousticGreens2D::~AcousticGreens2D() = default;

/// Set the source time function
void AcousticGreens2D::setSourceTimeFunction(const ISourceTimeFunction &stf)
{
    pImpl->mGreens.resize(0);
    if (pImpl->mSTF){pImpl->mSTF.reset();}
    if (!stf.isValid())
    {
        throw std::invalid_argument("Wavelet is invalid\n");
    }
    pImpl->mSTF = stf.clone();
}

/// Have the source time function?
bool AcousticGreens2D::haveSourceTimeFunction() const noexcept
{
    if (pImpl->mSTF){return true;}
    return false;
}

void AcousticGreens2D::setVelocity(const double vel)
{
    pImpl->mGreens.resize(0);
    pImpl->mVelocity = 0;
    if (vel <= 0)
    {
        throw std::invalid_argument("Velocity must be positive\n");
    }
    pImpl->mVelocity = vel;
}

void AcousticGreens2D::setDensity(const double density)
{
    pImpl->mGreens.resize(0);
    pImpl->mDensity = 0;
    if (density <= 0)
    {
        throw std::invalid_argument("Density must be positive\n");
    }
    pImpl->mDensity = density;
}

void AcousticGreens2D::setQualityFactor(const double q)
{
    pImpl->mGreens.resize(0);
    pImpl->mQ = 9999;
    if (q <= 0)
    {
        throw std::invalid_argument("Quality factor must be positive\n");
    }
    pImpl->mQ = q;
}

void AcousticGreens2D::setReceiverPositions(
    const std::vector<std::pair<double, double>> &positions)
{
    pImpl->mGreens.resize(0);
    pImpl->mReceiver.resize(0);
    if (positions.empty())
    {
        throw std::invalid_argument("No receivers!\n");
    }
    pImpl->mReceiver = positions;
}

void AcousticGreens2D::setSourcePosition(
    const std::pair<double, double> &position) noexcept
{
    pImpl->mGreens.resize(0);
    pImpl->mSource.resize(1);
    pImpl->mSource[0] = position;
}

/// Compute Green's functions
void AcousticGreens2D::compute()
{
    pImpl->mGreens.resize(0);
    if (pImpl->mVelocity <= 0){throw std::runtime_error("Velocity not set\n");}
    if (pImpl->mDensity <= 0){throw std::runtime_error("Density not set\n");}
    if (pImpl->mQ <= 0){throw std::runtime_error("Q not set\n");}
    // Compute the DFT frequencies
    if (!haveSourceTimeFunction())
    {
        throw std::runtime_error("Source time function not set\n");
    }
    // Work out the source time function and Fourier transform frequencies
    auto npts = pImpl->mSTF->getNumberOfSamples();
    auto wstf = pImpl->mSTF->getWaveletFourierTransform();
    auto samplingRate = pImpl->mSTF->getSamplingRate();
    auto freqs = DFTUtilities::realToComplexDFTFrequencies(npts,
                                                           1./samplingRate);
    // Initialize the inverse Fourier transform
    RTSeis::Utilities::Transforms::DFTRealToComplex rdft;
    rdft.initialize(npts);
    int lenft = rdft.getTransformLength();
    assert(lenft == static_cast<int> (freqs.size()));
    // Set some factors used in the greens function computation
    const std::complex<double> expipi4(std::sqrt(2.0)/2.0, sqrt(2.0)/2.0);
    const double vel = pImpl->mVelocity;
    const double Q = pImpl->mQ;
    const double vq2 = 2*vel*Q;
    const double fact =-1/(4*pImpl->mDensity*vel*vel);
    auto coeff = std::complex<double> (0, fact)*expipi4;
    int nsrc = static_cast<int> (pImpl->mSource.size());
    int nrec = static_cast<int> (pImpl->mReceiver.size());
    int nfreqs = static_cast<int> (freqs.size());
    pImpl->mGreens.resize(nrec);
    // Initialize the greens functions
    std::vector<double> zeros(npts, 0);
    for (int irec=0; irec<nrec; ++irec)
    {
        pImpl->mGreens[irec] = zeros;
    }
    // Loop on the receivers
    for (int irec=0; irec<nrec; ++irec)
    {
        for (int isrc=0; isrc<nsrc; ++isrc)
        {
            auto xs = pImpl->mSource[isrc].first;
            auto ys = pImpl->mSource[isrc].second;
            auto xr = pImpl->mReceiver[irec].first;
            auto yr = pImpl->mReceiver[irec].second;
            auto r = std::sqrt(std::pow(xr - xs, 2) + std::pow(yr - ys, 2));
            constexpr double twopi = 2*M_PI;
            std::vector<std::complex<double>> G(nfreqs);
            // Tabulate Green's functions for all frequencies
            #pragma omp simd
            for (int i=0; i<nfreqs; ++i)
            {
                double xscal = 0;
                if (freqs[i] > 0)
                {
                    xscal = std::sqrt(vel/(M_PI*M_PI*freqs[i]*r));
                }
                auto omega = twopi*freqs[i];
                auto omegar = omega*r;
                std::complex<double> carg(-omegar/vq2, -omegar/vel);
                G[i] = (coeff*xscal)*std::exp(carg);
            } // Loop on frequencies
            // Convolve the STF
            #pragma omp simd
            for (int i=0; i<nfreqs; ++i)
            {
                G[i] = wstf[i]*G[i];
            }
            // Inverse Fourier transform
            std::vector<double> timeSeries(npts);
            double *tsPtr = timeSeries.data();
            rdft.inverseTransform(lenft, G.data(), npts, &tsPtr); 
            // Stack the Green's functions
            std::transform(pImpl->mGreens[irec].begin(),
                           pImpl->mGreens[irec].end(),
                           timeSeries.begin(),
                           pImpl->mGreens[irec].begin(),
                           std::plus<double>());
        } // Loop on sources
    } // Loop on receivers
}

std::vector<double> AcousticGreens2D::getGreensFunction(const int irec) const
{
    if (pImpl->mGreens.empty())
    {
        throw std::runtime_error("Greens functions not yet computed\n");
    }
    int ngrns = static_cast<int> (pImpl->mGreens.size());
    if (irec < 0 || irec >= ngrns)
    {
        throw std::invalid_argument("irec = " + std::to_string(irec)
                                  + " must be in range [0,"
                                  + std::to_string(ngrns-1) + "]\n");
    }
    return pImpl->mGreens[irec];
}
