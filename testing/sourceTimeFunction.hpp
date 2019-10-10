#ifndef XCLOC_SOURCETIMEFUNCTION_HPP
#define XCLOC_SOURCETIMEFUNCTION_HPP
#include <complex>
#include <vector>
/*!
 * @brief Abstract base class for source time function.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class ISourceTimeFunction
{
public:
    //ISourceTimeFunction(const ISourceTimeFunction &stf);
    //ISourceTimeFunction& operator=(const ISourceTimeFunction &stf);
    virtual std::unique_ptr<ISourceTimeFunction> clone() const = 0;
    virtual ~ISourceTimeFunction() = default;
    virtual void clear() = 0;
    virtual void setSamplingRate(double samplingRate) = 0;
    virtual double getSamplingRate() const = 0;
    virtual bool haveSamplingRate() const noexcept = 0;
    virtual void setNumberOfSamples(int nSamples) = 0;
    virtual int getNumberOfSamples() const = 0;
    virtual bool haveNumberOfSamples() const noexcept = 0;
    virtual bool isValid() const noexcept = 0;
    virtual void setShiftWaveletToTraceStart(bool lshift, double tol) noexcept = 0;
    virtual bool getShiftWaveletToTraceStart() const noexcept = 0;
    virtual void setNormalizeByEnergy(bool lnorm) noexcept = 0;
    virtual bool getNormalizeByEnergy() const noexcept = 0;
    virtual std::vector<double> getWavelet() const = 0;
    virtual std::vector<std::complex<double>> getWaveletFourierTransform() const = 0;
};
#endif
