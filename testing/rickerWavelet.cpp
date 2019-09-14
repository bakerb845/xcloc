#include <cstdio>
#include <cstdlib>
#include <string>
#include <complex>
#include <cmath>
#include <cstring>
#include <ipps.h>
#include "rickerWavelet.hpp"
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"

class RickerWavelet::RickerWaveletImpl
{
public:
    std::vector<double> mRicker;
    std::vector<std::complex<double>> mRickerDFT;
    double mSamplingRate = 0;
    double mCenterFrequency = 0;
    double mTol = 0.005;
    int mSamples = 0;
    bool mShift = false;
    bool mNormalize = true;
    bool mHaveWavelet = false;
    bool mHaveDFT = false;
};

RickerWavelet::RickerWavelet() :
    pImpl(std::make_unique<RickerWaveletImpl>())
{
}

RickerWavelet::RickerWavelet(const RickerWavelet &wavelet)
{
    *this = wavelet;
}

RickerWavelet::RickerWavelet(RickerWavelet &&wavelet) noexcept
{
    *this = std::move(wavelet);
}

/// Copy assignment
RickerWavelet& RickerWavelet::operator=(const RickerWavelet &wavelet)
{
    if (&wavelet == this){return *this;}
    pImpl = std::make_unique<RickerWaveletImpl> (*wavelet.pImpl);
    return *this;
}

/// Deep copy
std::unique_ptr<ISourceTimeFunction> RickerWavelet::clone() const
{
    auto res = std::make_unique<RickerWavelet> (*this);
    return res;
}

/// Move assignment
RickerWavelet& RickerWavelet::operator=(RickerWavelet &&wavelet) noexcept
{
    if (&wavelet == this){return *this;}
    pImpl = std::move(wavelet.pImpl);
    return *this;
}

/// Destructor
RickerWavelet::~RickerWavelet() = default;

/// Clear
void RickerWavelet::clear() noexcept
{
    pImpl->mRicker.clear();
    pImpl->mRickerDFT.clear();
    pImpl->mSamplingRate = 0;
    pImpl->mCenterFrequency = 0;
    pImpl->mTol = 0.005;
    pImpl->mSamples = 0;
    pImpl->mShift = false;
    pImpl->mNormalize = true;
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
}

/// Sampling rate
void RickerWavelet::setSamplingRate(const double samplingRate)
{
    pImpl->mSamplingRate = 0;
    pImpl->mCenterFrequency = 0;
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("Sampling rate = " 
                                  + std::to_string(samplingRate)
                                  + " must be postiive\n");
    }
    pImpl->mSamplingRate = samplingRate;
}

double RickerWavelet::getSamplingRate() const
{
    if (!haveSamplingRate())
    {
        throw std::runtime_error("Sampling rate not set\n");
    }
    return pImpl->mSamplingRate;
}

bool RickerWavelet::haveSamplingRate() const noexcept
{
    if (pImpl->mSamplingRate <= 0){return false;}
    return true;
}

/// Center frequency
void RickerWavelet::setCenterFrequency(const double centerFrequency)
{
    pImpl->mCenterFrequency = 0;
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
    double fnyq = getSamplingRate()/2;
    if (centerFrequency <= 0)
    {
        throw std::invalid_argument("Center frequency = "
                                  + std::to_string(centerFrequency)
                                  + " must be postiive\n");
    }
    if (centerFrequency > fnyq)
    {
        throw std::invalid_argument("Center frequency = "
                                  + std::to_string(centerFrequency)
                                  + " cannot exceed Nyquist = "
                                  + std::to_string(fnyq) + "\n");
    }
    pImpl->mCenterFrequency = centerFrequency;
}

double RickerWavelet::getCenterFrequency() const
{
    if (!haveCenterFrequency())
    {
        throw std::runtime_error("Center frequency not set\n");
    }
    return pImpl->mCenterFrequency;
}

bool RickerWavelet::haveCenterFrequency() const noexcept
{
    if (pImpl->mCenterFrequency <= 0){return false;}
    return true;
}

/// Number of samples
void RickerWavelet::setNumberOfSamples(const int nSamples)
{
    pImpl->mSamples = 0;
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
    if (nSamples < 1)
    {
        throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                  + " must be positive\n");
    }
    pImpl->mSamples = nSamples;
}


int RickerWavelet::getNumberOfSamples() const
{
    if (!haveNumberOfSamples())
    {
        throw std::runtime_error("nSamples not set\n");
    }
    return pImpl->mSamples;
}

bool RickerWavelet::haveNumberOfSamples() const noexcept
{
    if (pImpl->mSamples < 1){return false;}
    return true;
}

/// Normalize? 
void RickerWavelet::setNormalizeByEnergy(const bool lnorm) noexcept
{
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
    pImpl->mNormalize = lnorm;
}

bool RickerWavelet::getNormalizeByEnergy() const noexcept
{
     return pImpl->mNormalize;
} 

/// Shift?
void RickerWavelet::setShiftWaveletToTraceStart(const bool lshift,
                                                const double tol) noexcept
{
    pImpl->mHaveWavelet = false;
    pImpl->mHaveDFT = false;
    if (tol <= 0 || tol >= 1)
    {
        pImpl->mShift = false;
        pImpl->mTol = 0;
    }
    else
    {
        pImpl->mShift = lshift;
        pImpl->mTol = tol;
    }
}

bool RickerWavelet::getShiftWaveletToTraceStart() const noexcept
{
    return pImpl->mShift;
}

/// Compute the ricker wavelet
std::vector<double> RickerWavelet::getWavelet() const
{
    if (!pImpl->mHaveWavelet)
    {
        pImpl->mHaveDFT = false;
        // Compute a Ricker wavelet whose peak is at the middle of the time series
        const double peakFreq = pImpl->mCenterFrequency; 
        const double pi2f2 = std::pow(M_PI*peakFreq, 2);
        const double dt = 1./getSamplingRate();
        int npts = getNumberOfSamples();
        int npts2 = npts/2;
        double xmax = 0;
        pImpl->mRicker.resize(npts);
        for (int i=0; i<npts; ++i)
        {
            double t = static_cast<double> (-npts2 + i)*dt;
            double t2 = t*t;
            double pi2f2t2 = pi2f2*t2;
            // Make this negative so first motion is up
            pImpl->mRicker[i] =-(1.0 - 2.0*pi2f2t2)*std::exp(-pi2f2t2);
            xmax = std::max(std::abs(pImpl->mRicker[i]), xmax);
        }
        // Shift the wavelet to the start of the trace - easier to work in
        // pct of 1 b/c the wavelet has max value of 1
        double *ricker = pImpl->mRicker.data();
        if (pImpl->mShift)
        {
            double tolXmax = pImpl->mTol*xmax;
            std::vector<double> work(npts, 0); // Initialize to zeros
            for (int i=0; i<npts; ++i)
            {
                if (std::abs(ricker[i]) > tolXmax)
                {
                    int ncopy = npts - i;
                    std::copy(ricker+i, ricker+i+ncopy, work.begin());
                    break;
                }
            }
            std::fill(ricker, ricker+static_cast<size_t> (npts), 0);
            std::copy(work.begin(), work.end(), ricker);
        }
        // Normalize the area energy in the signal
        if (pImpl->mNormalize)
        {
            double area;
            ippsNorm_L2_64f(ricker, npts, &area);
            ippsDivC_64f_I(area, ricker, npts);
        }
        pImpl->mHaveWavelet = true;
    }
    return pImpl->mRicker; 
}

/// Compute Fourier transform of wavelet
std::vector<std::complex<double>>
RickerWavelet::getWaveletFourierTransform() const
{
    RTSeis::Utilities::Transforms::DFTRealToComplex<double> dft;
    auto npts = getNumberOfSamples();
    if (!pImpl->mHaveWavelet)
    {
        auto wavelet = getWavelet();
        pImpl->mRickerDFT.resize(npts/2+1); 
        dft.initialize(npts);
        std::complex<double> *ptr = pImpl->mRickerDFT.data();
        dft.forwardTransform(npts, pImpl->mRicker.data(),
                             npts/2+1, &ptr);
    }
    else
    {
        if (!pImpl->mHaveDFT)
        {
            pImpl->mRickerDFT.resize(npts/2+1);
            dft.initialize(npts);
            std::complex<double> *ptr = pImpl->mRickerDFT.data();
            dft.forwardTransform(npts, pImpl->mRicker.data(),
                                 npts/2+1, &ptr);
            pImpl->mHaveDFT = true;
        }
    }
    return pImpl->mRickerDFT; 
} 

/// Valid?
bool RickerWavelet::isValid() const noexcept
{
    if (!haveSamplingRate()){return false;}
    if (!haveCenterFrequency()){return false;}
    if (!haveNumberOfSamples()){return false;}
    return true;
}
