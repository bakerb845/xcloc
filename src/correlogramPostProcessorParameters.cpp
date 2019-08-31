#include <cstdio>
#include <cstdlib>
#include <string>
#include "xcloc/correlogramPostProcessorParameters.hpp"

using namespace XCLoc;

class CorrelogramPostProcessorParameters::CorrelogramPostProcessorImpl
{
public:
    int mEnvelopeLength = 301; 
    CorrelogramFilteringType mFilterType = CorrelogramFilteringType::UNKNOWN;
};

/// Constructor
CorrelogramPostProcessorParameters::CorrelogramPostProcessorParameters() :
    pImpl(std::make_unique<CorrelogramPostProcessorImpl>())
{
}

/// Copy constructor
CorrelogramPostProcessorParameters::CorrelogramPostProcessorParameters(
    const CorrelogramPostProcessorParameters &parameters)
{
    *this = parameters;
}

/// Move constructor
CorrelogramPostProcessorParameters::CorrelogramPostProcessorParameters(
    CorrelogramPostProcessorParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy operator
CorrelogramPostProcessorParameters&
CorrelogramPostProcessorParameters::operator=(
    const CorrelogramPostProcessorParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<CorrelogramPostProcessorImpl> (*parameters.pImpl);
    return *this;
}

/// Move operator
CorrelogramPostProcessorParameters&
CorrelogramPostProcessorParameters::operator=(
    CorrelogramPostProcessorParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
CorrelogramPostProcessorParameters::~CorrelogramPostProcessorParameters() = default;

/// Clears class
void CorrelogramPostProcessorParameters::clear() noexcept
{
    pImpl->mEnvelopeLength = 301;
    pImpl->mFilterType = CorrelogramFilteringType::UNKNOWN;
}

/// No filtering - compute raw correlograms
void CorrelogramPostProcessorParameters::setNoFiltering() noexcept
{
    pImpl->mFilterType = CorrelogramFilteringType::NO_FILTERING;
}

/// FIR envelope
void CorrelogramPostProcessorParameters::setFIREnvelopeFiltering(
    const int envelopeLength)
{
    pImpl->mFilterType = CorrelogramFilteringType::UNKNOWN;
    if (envelopeLength < 1)
    {
        throw std::invalid_argument("Envelope length = " 
                                  + std::to_string(envelopeLength)
                                  + " must be positive\n");
    }
    if (envelopeLength%2 == 1)
    {
        pImpl->mEnvelopeLength = envelopeLength;
    }
    else
    {
        pImpl->mEnvelopeLength = envelopeLength + 1;
    }
    pImpl->mFilterType = CorrelogramFilteringType::FIR_ENVELOPE_FILTERING;
}

int CorrelogramPostProcessorParameters::getFIREnvelopeFilterLength() const
{
    if (getFilteringType() != CorrelogramFilteringType::FIR_ENVELOPE_FILTERING)
    {
        throw std::runtime_error("Must set an FIR envelope filter\n");
    }
    return pImpl->mEnvelopeLength;
}

CorrelogramFilteringType
CorrelogramPostProcessorParameters::getFilteringType() const
{
    if (pImpl->mFilterType == CorrelogramFilteringType::UNKNOWN)
    {
        throw std::runtime_error("Filter strategy not yet set\n");
    }
    return pImpl->mFilterType;
}

/// Check parameters
bool CorrelogramPostProcessorParameters::isValid() const noexcept
{
    if (pImpl->mFilterType == CorrelogramFilteringType::UNKNOWN){return false;}
    return true;
}
