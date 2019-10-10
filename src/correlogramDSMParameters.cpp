#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include "xcloc/correlogramDSMParameters.hpp"

using namespace XCLoc;

class CorrelogramDSMParameters::CorrelogramDSMParametersImpl
{
public:
    /// Correlogram chunk size
    size_t mCacheBlockSize = 4096;
};

/// Constructor
CorrelogramDSMParameters::CorrelogramDSMParameters() :
    pImpl(std::make_unique<CorrelogramDSMParametersImpl> ())
{
}

CorrelogramDSMParameters::CorrelogramDSMParameters(
    const CorrelogramDSMParameters &parameters)
{
    *this = parameters;
}

CorrelogramDSMParameters::CorrelogramDSMParameters(
    CorrelogramDSMParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment operator
CorrelogramDSMParameters&
CorrelogramDSMParameters::operator=(const CorrelogramDSMParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<CorrelogramDSMParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment operator
CorrelogramDSMParameters&
CorrelogramDSMParameters::operator=(
    CorrelogramDSMParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
CorrelogramDSMParameters::~CorrelogramDSMParameters() = default;

/// Reset class
void CorrelogramDSMParameters::clear() noexcept
{
    pImpl->mCacheBlockSize = 4096;
}

/// Block size
void CorrelogramDSMParameters::setCacheBlockSize(const size_t cacheBlockSize)
{
    // Check minimum size
    if (cacheBlockSize < 64)
    {
        throw std::invalid_argument("cache block size = " 
                                  + std::to_string(cacheBlockSize)
                                  + " must be at least 64\n");
    }
    // Check this is a power of 2: block = 2**log2(block)
    auto p = static_cast<int> (std::log2(static_cast<double> (cacheBlockSize)));
    auto result = static_cast<size_t> (std::exp2(p));
    if (result != cacheBlockSize)
    {
        throw std::invalid_argument("cache block size = "
                                  + std::to_string(cacheBlockSize)
                                  + " is not a power of 2\n"); 
    }
    pImpl->mCacheBlockSize = cacheBlockSize;
}

size_t CorrelogramDSMParameters::getCacheBlockSize() const noexcept
{
    return pImpl->mCacheBlockSize;
}

bool CorrelogramDSMParameters::isValid() const noexcept
{
    return true;
}
