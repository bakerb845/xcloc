#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cassert>
#include "xcloc/correlogramParameters.hpp"
#include "xcloc/enums.hpp"

using namespace XCLoc;

class CorrelogramParameters::CorrelogramParametersImpl
{
public:
    /// Cross-correlation pairs
    std::vector<std::pair<int,int>> mXCPairs;
    /// Number of samples in each input signal
    int mSamples = 0;
    /// The number of samples to which to pad
    int mPadSamples = 0;
    /// FIR filter length for envelope
    int mEnvelopeLength = 0;
    /// Filtering post-processing strategy
    CorrelogramFilteringType mFilterType
         = CorrelogramFilteringType::NO_FILTERING;
    /// MKL accuracy
    MKLFloatingPointAccuracy mMKLAccuracy
         = MKLFloatingPointAccuracy::HIGH_ACCURACY;
};

/// Constructors
CorrelogramParameters::CorrelogramParameters() :
    pImpl(std::make_unique<CorrelogramParametersImpl> ())
{
}

CorrelogramParameters::CorrelogramParameters(
    const CorrelogramParameters &parameters)
{
    *this = parameters;
}

CorrelogramParameters::CorrelogramParameters(
    CorrelogramParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Operators
CorrelogramParameters&
CorrelogramParameters::operator=(
    const CorrelogramParameters &parms)
{
    if (&parms == this){return *this;}
    pImpl = std::make_unique<CorrelogramParametersImpl> (*parms.pImpl);
    return *this; 
}

CorrelogramParameters&
CorrelogramParameters::operator=(
   CorrelogramParameters &&parms) noexcept
{
    if (&parms == this){return *this;}
    pImpl = std::move(parms.pImpl);
    return *this;
}

/// Destructors
CorrelogramParameters::~CorrelogramParameters() = default;

void CorrelogramParameters::clear() noexcept
{
    pImpl->mXCPairs.clear();
    pImpl->mSamples = 0;
    pImpl->mPadSamples = 0;
    pImpl->mEnvelopeLength = 0;
    pImpl->mFilterType = CorrelogramFilteringType::NO_FILTERING;
    pImpl->mMKLAccuracy = MKLFloatingPointAccuracy::HIGH_ACCURACY;
}

/// Number of samples
void CorrelogramParameters::setNumberOfSamples(const int nSamples)
{
    pImpl->mSamples = 0;
    pImpl->mPadSamples = 0;
    if (nSamples < 1)
    {
        throw std::invalid_argument("Number of samples = "
                                  + std::to_string(nSamples)
                                  + " must be positive\n");
    }
    pImpl->mSamples = nSamples;
    pImpl->mPadSamples = nSamples;
}

int CorrelogramParameters::getNumberOfSamples() const
{
    if (pImpl->mSamples < 1)
    {
        throw std::runtime_error("Number of samples was never specified\n");
    }
    return pImpl->mSamples;
}

/// Creat a default correlation pairs table
void CorrelogramParameters::setCorrelationPairs(
    const int nSignals, const bool ldoAutoCorrelations)
{
    pImpl->mXCPairs.clear();
    if (nSignals < 2)
    {
        throw std::invalid_argument("nSignals = " + std::to_string(nSignals)
                                  + " must be at least 2\n");
    }
    std::vector<std::pair<int, int>> xcPairs;
    // Only superdiagonal
    if (!ldoAutoCorrelations)
    {
        int nxcPairs = (nSignals*(nSignals - 1))/2;
        xcPairs.reserve(nxcPairs);
        for (int i=0; i<nSignals; ++i)
        {
            for (int j=i+1; j<nSignals; ++j)
            {
                xcPairs.push_back(std::make_pair(i, j));
            }
        }
#ifdef DEBUG
        cassert(static_cast<int> (xcPairs.size()) == nxcPairs);
#endif
    }
    // Upper triangle (includes diagonal)
    else
    {
        int nxcPairs = (nSignals*(nSignals + 1))/2;
        xcPairs.reserve(nxcPairs);
        for (int i=0; i<nSignals; ++i)
        {
            for (int j=i; j<nSignals; ++j)
            {
                xcPairs.push_back(std::make_pair(i,j));
            }
        }
#ifdef DEBUG
        cassert(static_cast<int> (xcPairs.size()) == nxcPairs);
#endif
    }
    setCorrelationPairs(xcPairs);
}

void CorrelogramParameters::setCorrelationPairs(
    const std::vector<std::pair<int, int>> &xcPairs)
{
    pImpl->mXCPairs.clear();
    if (xcPairs.size() < 1)
    {
        throw std::invalid_argument("xcPairs.size() must be at least 1\n");
    }
    pImpl->mXCPairs = xcPairs; 
}

std::vector<std::pair<int, int>>
CorrelogramParameters::getCorrelationPairs() const
{
    if (pImpl->mXCPairs.empty())
    {
        throw std::runtime_error("Correlation pairs never set\n");
    }
    return pImpl->mXCPairs;
}

std::pair<int, int>
CorrelogramParameters::getCorrelationPair(const int ixc) const
{
    auto nXCPairs = getNumberOfCorrelationPairs(); // throws
    if (ixc < 0 || ixc >= nXCPairs)
    {
        throw std::invalid_argument("ixc = " + std::to_string(ixc)
                                  + " must be in range [0,"
                                  + std::to_string(nXCPairs-1) + "]\n");
    }
    std::pair<int, int> result = pImpl->mXCPairs[ixc];
    return result;
}

/// Number of correlation pairs
int CorrelogramParameters::getNumberOfCorrelationPairs() const
{
    if (pImpl->mXCPairs.empty())
    {
        throw std::runtime_error("Correlation pairs never set\n");
    }
    return static_cast<int> (pImpl->mXCPairs.size());
}


/// Number of padded samples
void CorrelogramParameters::setNumberOfPaddedSamples(
    const int nPadSamples)
{
    pImpl->mPadSamples = pImpl->mSamples;
    int nSamples = getNumberOfSamples(); /// Throws
    if (nPadSamples < nSamples)
    {
        throw std::invalid_argument("nPadSamples = " 
                                  + std::to_string(nPadSamples)
                                  + " must be at least "
                                  + std::to_string(nSamples) + "\n");
    }
    pImpl->mPadSamples = nPadSamples; 
}

int CorrelogramParameters::getNumberOfPaddedSamples() const
{
    if (pImpl->mPadSamples < 1)
    {
        throw std::runtime_error("Number of samples never set\n");
    }
    return pImpl->mPadSamples;
}

int CorrelogramParameters::getCorrelogramLength() const
{
    int lenxc = 2*getNumberOfPaddedSamples() - 1; // throws
    return lenxc;
}

/// MKL accuracy
void CorrelogramParameters::setMKLFloatingPointAccuracy(
    const MKLFloatingPointAccuracy accuracy) noexcept
{
    pImpl->mMKLAccuracy = accuracy;
}

MKLFloatingPointAccuracy 
CorrelogramParameters::getMKLFloatingPointAccuracy() const noexcept
{
    return pImpl->mMKLAccuracy;
}

/// No filtering - compute raw correlograms
void CorrelogramParameters::setNoFiltering() noexcept
{
    pImpl->mFilterType = CorrelogramFilteringType::NO_FILTERING;
}

/// FIR envelope
void CorrelogramParameters::setFIREnvelopeFiltering(
    const int envelopeLength)
{
    pImpl->mFilterType = CorrelogramFilteringType::NO_FILTERING;
    if (envelopeLength < 1)
    {
        throw std::invalid_argument("Envelope length = "
                                  + std::to_string(envelopeLength)
                                  + " must be positive\n");
    }
    int lenxc = getCorrelogramLength(); // Throws
    if (envelopeLength%2 == 1)
    {
        if (envelopeLength > lenxc)
        {
            throw std::invalid_argument("Envelope length = "
                                      + std::to_string(envelopeLength)
                                      + " cannot exceed "
                                      + std::to_string(lenxc) + "\n");
        }
        pImpl->mEnvelopeLength = envelopeLength;
    }
    else
    {
        if (envelopeLength + 1 > lenxc)
        {
            throw std::invalid_argument("Envelope length = "
                                      + std::to_string(envelopeLength + 1)
                                      + " cannot exceed "
                                      + std::to_string(lenxc) + "\n");
        }
        pImpl->mEnvelopeLength = envelopeLength + 1;
    }
    pImpl->mFilterType = CorrelogramFilteringType::FIR_ENVELOPE_FILTERING;
}

int CorrelogramParameters::getFIREnvelopeFilterLength() const
{
    if (getFilteringType() != CorrelogramFilteringType::FIR_ENVELOPE_FILTERING)
    {
        throw std::runtime_error("Must set an FIR envelope filter\n");
    }
    return pImpl->mEnvelopeLength;
}

/// Filtering strategy 
CorrelogramFilteringType
CorrelogramParameters::getFilteringType() const noexcept
{
    return pImpl->mFilterType;
}

// Valid?
bool CorrelogramParameters::isValid() const noexcept
{
    if (pImpl->mSamples < 1){return false;}
    if (pImpl->mXCPairs.empty()){return false;}
    if (pImpl->mFilterType == CorrelogramFilteringType::FIR_ENVELOPE_FILTERING)
    {
        if (pImpl->mEnvelopeLength%2 != 1){return false;}
    }
    return true;
}
