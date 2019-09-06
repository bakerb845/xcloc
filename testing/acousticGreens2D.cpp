#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "acousticGreens2D.hpp"
#include "rickerWavelet.hpp"

class AcousticGreens2D::AcousticGreens2DImpl
{
public:
    /// Source time function
    RickerWavelet mWavelet;
    /// Receiver (x, y, z) position in meters
    std::vector<double> mReceiverX;
    std::vector<double> mReceiverY;
    std::vector<double> mReceiverZ;
    /// Source (x, y, z) position in meters
    std::vector<double> mSourceX;
    std::vector<double> mSourceY;
    std::vector<double> mSourceZ;
    /// Velocity (m/s)
    double mVelocity = 0;
    /// Density (kg/m**3)
    double mDensity = 0;
    /// Quality factor
    double mQ = 0;
};

/// Constructor
AcousticGreens2D::AcousticGreens2D() :
    pImpl(std::make_unique<AcousticGreens2DImpl> ())
{
}

/// Destructor
AcousticGreens2D::~AcousticGreens2D() = default;

/// Set ricker wavelet
void AcousticGreens2D::setRickerSourceTimeFunction(
    const RickerWavelet &wavelet)
{
    pImpl->mWavelet.clear();
    if (!wavelet.isValid())
    {
        throw std::invalid_argument("Wavelet is invalid\n");
    }
    pImpl->mWavelet = wavelet;
}
