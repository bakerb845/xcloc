#include <cstdio>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#else
namespace
{
//int omp_get_thread_num(){return 0;}
int omp_get_num_threads(){return 1;}
}
#endif
#include <vector>
#include <rtseis/utilities/transforms/firEnvelope.hpp>
#include "xcloc/correlogramPostProcessor.hpp"
#include "xcloc/correlogramPostProcessorParameters.hpp"

using namespace XCLoc;

template<class T>
class CorrelogramPostProcessor<T>::CorrelogramPostProcessorImpl
{
public:
    CorrelogramPostProcessorImpl()
    {
        mNumberOfThreads = omp_get_num_threads();
        mEnvelopes.resize(mNumberOfThreads);
    }
    void clear() noexcept
    {
        mEnvelopes.clear();
        mNumberOfThreads = 1;
    }
    std::vector<RTSeis::Utilities::Transforms::FIREnvelope<T>> mEnvelopes;
    std::shared_ptr<const XCLoc::CorrelationEngine<T>> mCorrelationEngine;
    int mNumberOfThreads = 1;
};


/// Constructor
template<class T>
CorrelogramPostProcessor<T>::CorrelogramPostProcessor() :
    pImpl(std::make_unique<CorrelogramPostProcessorImpl>())
{
}

/*
/// Copy constructor 
template<class T>
CorrelogramPostProcessor<T>::CorrelogramPostProcessor(
    const CorrelgramPostProcessor &processor)
{
    *this = processor;
}

/// Move constructor 
template<class T>
CorrelogramPostProcessor<T>::CorrelogramPostProcessor(
    CorrelogramPostProcessor &&processor) noexcept
{
    *this = std::move(processor);
}

template<class T>
CorrelogramPostProcessor& CorrelogramPostProcessor<T>::operator=(
    const CorrelogramPostProcessor &processor)
{
    if (&processor = this){return *this;}
    if (pImpl){pImpl.clear();}
    pImpl = std::make_unique<CorrelogramPostProcessor> (*processor.pImpl);
    return *this;
}

template<class T>
CorrelogramPostProcessor& CorrelogramPostProcessor<T>::operator=(
    CorrelogramPostProcessor &&processor) noexcept
{
    if (&processor = this){return *this;}
    if (pImpl){pImpl.clear();}
    pImpl = std::move(processor.pImpl);
    return *this; 
}
*/

/// Destructor
template<class T>
CorrelogramPostProcessor<T>::~CorrelogramPostProcessor() = default;

template<class T>
void CorrelogramPostProcessor<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initialize the class
template<class T>
void CorrelogramPostProcessor<T>::initialize(
    const CorrelogramPostProcessorParameters &parameters,
    std::shared_ptr<const XCLoc::CorrelationEngine<T>> engine)
{
    if (!parameters.isValid())
    {
        throw std::invalid_argument("parameters are not valid\n");
    }
    pImpl->mCorrelationEngine.reset();
    pImpl->mCorrelationEngine = engine;
}

/// Template class instantiation
template class XCLoc::CorrelogramPostProcessor<double>;
template class XCLoc::CorrelogramPostProcessor<float>;
