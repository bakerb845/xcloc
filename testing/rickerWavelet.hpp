#ifndef XCLOC_RICKERWAVELET_HPP
#define XCLOC_RICKERWAVELET_HPP
#include <complex>
#include <vector>
#include <memory>
#include "sourceTimeFunction.hpp"
/*!
 * @name RickerWavelet "gaussianWavelet.hpp"
 * @brief A class for computing Gaussian wavelet source time functions.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class RickerWavelet : public ISourceTimeFunction
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    RickerWavelet();
    /*!
     * @brief Copy constructor
     * @param[in] wavelet  The wavelet class from which to initialize this.
     */
    RickerWavelet(const RickerWavelet &wavelet);
    /*!
     * @brief Move constructor.
     * @param[in,out] wavelet  On input this contains the wavelet class from
     *                         which to initialize this.  On exit, wavelet's
     *                         behavior is undefined.
     */
    RickerWavelet(RickerWavelet &&wavelet) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] wavelet  The wavelet to copy.
     * @result A deep copy of the input wavelet class. 
     */
    RickerWavelet& operator=(const RickerWavelet &wavelet);
    /*!
     * @brief Makes a copy the source time function.
     * @result A deep copy of the source time function.
     */
    std::unique_ptr<ISourceTimeFunction> clone() const override;
    /*!
     * @brief Move assignment operator.
     * @param[in,out] wavelet  The wavelet whose memory will be moved to this.
     *                         On exit wavelet's behavior is undefined.
     */
    RickerWavelet& operator=(RickerWavelet &&wavelet) noexcept;
    /*! } */

    /*! @name Destructors
     * @{
     */
    virtual ~RickerWavelet();
    /*!
     * @brief Resets class and releases memory.
     */
    void clear() noexcept override;
    /*! @} */

    /*! @name Sampling Rate
     * @{
     */
    /*!
     * @brief The trace sampling rate.
     * @param[in] samplingRate  The sampling rate in Hz.
     * @note This will invalidate the center frequency.
     */
    void setSamplingRate(double samplingRate) override;
    /*!
     * @brief Gets the sampling rate.
     * @result The sampling rate.
     * @throws std::runtime_error if this was not set.
     */
    double getSamplingRate() const override;
    /*!
     * @brief Determines if the sampling rate was set.
     * @result True indicates that the sampling rate was set.
     */ 
    bool haveSamplingRate() const noexcept override;
    /*! @} */

    /*! @name Number of Samples
     *{
     */
    /*!
     * @brief Sets the number of samples in the wavelet time series.
     * @throws std::invalid_argument if this is not positive.
     */
    void setNumberOfSamples(int nSamples) override;
    /*!
     * @brief Gets the number of samples in the wavelet time series.
     * @throws std::runtime_error if this was not set.
     */
    int getNumberOfSamples() const override;
    /*!
     * @brief Determines whether or not the number of samples was set.
     * @result True indicates that the number of samples was set.
     */
    bool haveNumberOfSamples() const noexcept override;
    /*! @} */

    /*!
     * @brief Checks if this is usable source time function class.
     * @result True indicates that this is a usable source time function.
     */
    bool isValid() const noexcept override;

    /*! @name Center Frequency
     * @{
     */
    /*!
     * @brief Sets the center frequency of the Ricker wavelet.
     * @param[in] centerFrequency  The center frequency in Hz.
     * @throws std::invalid_argument if this is not postive.
     * @throws std::runtime_error if this exceeds the 1/2 samplingRate.
     * @sa \c getSamplingRate()
     */
    void setCenterFrequency(double centerFrequency);
    /*!
     * @brief Gets the center frequency.
     * @result The wavelet's center frequency in Hz.
     * @sa \c haveCenterFrequency()
     */
    double getCenterFrequency() const;
    /*!
     * @brief Determines whether or not the center frequency is set.
     * @result True indicates that center frequency is set.
     */
    bool haveCenterFrequency() const noexcept;
    /*! @} */

    /*!
     * @brief Enables/disables shifting the wavelet start to the beginning of
     *        the trace.
     * @param[in] lshift  If true then attempt to shift the wavelet to
     *                    the start of the trace. 
     * @param[in] tol     This is a tolerance that when the absolute value of
     *                    the i'th  wavelet exceeds tol*max(abs(wavelet)) then
     *                    the i'th sample will be declared the trace start.
     *                    This is only relevant if lshift is true.   
     *                    If this is <= 0 or >= 1 then lshift will be set to
     *                    false.
     */
    void setShiftWaveletToTraceStart(bool lshift, double tol = 0.0001) noexcept override;
    /*!
     * @brief Determines whether or not to shift the wavelet to the start.
     */
    bool getShiftWaveletToTraceStart() const noexcept override;
    /*!
     * @brief Enables/disabling normalizing the wavelet by the square root of
     *        the energy.
     * @param[in] lnorm  If true then the wavelet will be normalized by the
     *                   square root of the energy.
     */
    void setNormalizeByEnergy(bool lnorm) noexcept override;
    /*!
     * @brief Determines whether or not to normalize the wavelet.
     * @result If true then the wavelet is normalized by the energy. 
     */
    bool getNormalizeByEnergy() const noexcept override;

    /*! @name Wavelets
     * @{
     */
    /*!
     * @brief Computes the Ricker wavelet.
     * @result The Ricker wavelet.
     * @throws std::runtime_error if the number of points, center frequency, or
     *         sampling rate was not set.
     */
    std::vector<double> getWavelet() const override;
    /*!
     * @brief Gets DFT of Ricker wavelet.
     * @result The DFT of the Ricker wavelet.  This goes from 0 to the Nyquist
     *         frequency if \c getNumberOfSamples() is even.  When 
     *         \c getNumberOfSamples() is odd then the last sample is just
     *         just shy of the Nyquist frequency.
     */
    std::vector<std::complex<double>> getWaveletFourierTransform() const override;
    /*! @} */
private:
    class RickerWaveletImpl;
    std::unique_ptr<RickerWaveletImpl> pImpl;
};

#endif
