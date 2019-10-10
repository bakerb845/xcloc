#ifndef XCLOC_CROSSCORRELATIONLOCATOR_HPP
#define XCLOC_CROSSCORRELATIONLOCATOR_HPP
#include <memory>
template<class T>
class CrossCorrelationLocator
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    CrossCorrelationLocator();
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~CrossCorrelationLocator();
    /*!
     * @brief Clears all memory from module and resets class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Initialization
     * @{
     */
    void setCorrelationParameters( ); 
    /*!
     * @brief True indicates that the travel time tables and cross-correlation
     *        parameters set set.
     * @result True indicates that the class is ready to be applied to data.
     */
    bool isInitialized() const noexcept;
    /*! @} */

    /*! @name Input Signals
     * @{
     */
    /*!
     * @brief Sets the input signal.
     * @param[in] waveid   The waveform identifier.
     * @throws std::invalid_argument if the waveform ID cannot be found,
     *         the number of points is unexpected, or signal is NULL. 
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void setInputSignal(const WaveformIdentifier &waveid,
                        int npts,
                        const double signal[]);
    /*! @copydoc setInputSignal */ 
    void setInputSignal(const WaveformIdentifier &waveid,
                        int npts,
                        const float signal[]);
    /*!
     * @brief Zeros all elements of the given waveform.  This i suseful
     *        if the station is late or has a gap.
     * @param[in] waveid   The waveform identifier.
     * @throws std::invalid_argument if the waveform ID cannot be found.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void zeroInputSignal(const WaveformIdentifier &waveid);
    /*! @} */

    /*! @name Compute
     * @{
     */
    /*!
     * @brief Performs the correlations then migrates the correlograms.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void compute();
    /*! @} */

private:
    class CrossCorrelationLocatorImpl;
    std::unique_ptr<CrossCorrelationLocatorImpl> pImpl;
};

#endif
