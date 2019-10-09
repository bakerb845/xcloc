#ifndef XCLOC_CORRELOGRAMS_HPP
#define XCLOC_CORRELOGRAMS_HPP 1
#include <memory>
#include "xcloc/enums.hpp"

namespace XCLoc
{
class WaveformIdentifier;
class CorrelogramParameters;
/*!
 * @class Correlograms "correlationEngine.hpp" "xcloc/correlationEngine.hpp"
 * @brief This computes the correlograms.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T=double>
class Correlograms
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    Correlograms();
    /*!
     * @brief Constructs the correlation engine from the parameters.
     * @param[in] parameters  The correlation engine parameters.
     * @throws std::invalid_argument if parameters is not valid.
     */
    explicit Correlograms(const CorrelogramParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] engine  The correlation engine from which to initialize
     *                        this class.  On exit, engine's behavior is
     *                        undefined.
     */
    Correlograms(Correlograms &&engine) noexcept;
    /*! @} */

    /*! @brief Operators
     * @{
     */
    /*!
     * @brief Move assignment operator.
     * @param[in,out] engine  The correlation engine whose memory will be moved
     *                        to this.  On exit, engine's behavior is undefined.
     * @result The memory moved from engine to this.
     */
    Correlograms& operator=(Correlograms &&engine) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor 
     */
    ~Correlograms();
    /*!
     * @brief Clears the memory and resets the module.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Initialization
     * @{
     */
    /*!
     * @brief Initializes the cross-correlation engine.
     * @throws std::invalid_argument if parameters are not valid.
     */
    void initialize(const CorrelogramParameters &parameters);
    /*!
     * @brief Sets the iw'th waveform.
     * @param[in] waveid    The waveform identifier.  This refers to the signal
     *                      ID in the correlation pair table.
     * @param[in] nSamples  The number of samples in the waveform.  This must
     *                      match \c getInputSignalLength().
     * @param[in] x         The signal to set.  This is an array of dimension
     *                      [nSamples].
     * @throws std::runtime_error if the class is not intialized.
     * @throws std::invalid_argument if nSamples is incorrect, x is NULL, or 
     *         the waveform identifier does not exist.
     * @sa \c getInputSignalLength() \c isInitialized()
     */
    void setInputSignal(int64_t waveid, int nSamples, const double x[]);
    /*! @copydoc setInputSignal */
    void setInputSignal(int64_t waveid, int nSamples, const float x[]);
    /*!
     * @brief Sets all elements of the iw'th waveform to zero.  This is useful
     *        if the station is late or has a gap.
     * @param[in] waveid   Sets all elements of the waveid'th signal to zero.
     */
    void zeroInputSignal(int waveid); 
    /*! @} */

    /*! @brief Compute
     * @{
     */
    /*!
     * @brief Computes the cross-correlograms.
     * @throws std::runtime_error 
     */
    void computeCrossCorrelograms();
    /*!
     * @brief Computes the phase-correlograms.
     * @throws std::runtime_error if the class is not initialized and all the
     *         signals are not set.
     */
    void computePhaseCorrelograms();
    /*! @} */

    /*! @name Correlograms
     * @{
     */
    /*!
     * @brief Gets a pointer to ixc'th unprocessed correlogram.
     * @param[in] ixc   The correlogram index.  This must be in the range of
     *                  [0, \c getNumberOfCorrelograms() - 1].
     * @result A pointer to the ixc'th correlogram.  This is an array whose
     *         length is [\c getCorrelogramLength()].
     * @throws std::runtime_error if the correlograms are not yet computed. 
     * @sa \haveCorrelograms(), \c getCorrelogramLength()
     */
    const T* getRawCorrelogramPointer(int ixc) const;
    /*!
     * @brief Gets a pointer to ixc'th processed correlogram.
     * @param[in] ixc   The correlogram index.  This must be in the range of
     *                  [0, \c getNumberOfCorrelograms() - 1].
     * @result A pointer to the ixc'th correlogram.  This is an array whose
     *         length is [\c getCorrelogramLength()].
     * @throws std::runtime_error if the correlograms are not yet computed.
     * @sa \haveCorrelograms(), \c getCorrelogramLength()
     */
    const T* getProcessedCorrelogramPointer(int ixc) const;
    /*!
     * @brief Gets the ixc'th unprocessed correlogram.
     * @param[in] ixc    The correlogram index.  This must be in the range of
     *                   [0, \c getNumberOfCorrelograms() - 1].
     * @param[in] nwork  The length of xc.  This must be at least
     *                   \c getCorrelogramLength().
     * @param[out] xc    The ixc'th correlogram.  This is an array of dimension
     *                   [nwork] however only the first
     *                   \c getCorrelogramLength() points are accessed.
     * @throws std::invalid_argument if nwork is too small or xc is NULL.
     * @throws std::runtime_error if the correlograms are not yet computed.
     * @sa \haveCorrelograms(), \c getCorrelogramLength()
     */
    void getRawCorrelogram(int ixc, int nwork, T *xc[]) const;
    /*!
     * @brief Gets the ixc'th processed correlogram.
     * @param[in] ixc    The correlogram index.  This must be in the range of
     *                   [0, \c getNumberOfCorrelograms() - 1].
     * @param[in] nwork  The length of xc.  This must be at least
     *                   \c getCorrelogramLength().
     * @param[out] xc    The ixc'th correlogram.  This is an array of dimension
     *                   [nwork] however only the first
     *                   \c getCorrelogramLength() points are accessed.
     * @throws std::invalid_argument if nwork is too small or xc is NULL.
     * @throws std::runtime_error if the correlograms are not yet computed.
     * @sa \haveCorrelograms(), \c getCorrelogramLength()
     */
    void getProcessedCorrelogram(int ixc, int nwork, T *xc[]) const;
    /*!
     * @brief Get the waveform identifiers comprising the ixc'th correlation
     *        pair.
     * @param[in] ixc  The correlogram index.  This must be in the range of
     *                 [0, \c getNumberOfCorrelograms() - 1].
     * @result The waveform identifiers comprising the ixc'th correlation.
     * @throws std::runtime_error if the class is not initialized.
     */
    std::pair<WaveformIdentifier,WaveformIdentifier> getCorrelationPair(int ixc) const;
    /*! @} */

    /*! @name Properties
     * @{
     */
    /*!
     * @brief Gets the expected length of the input signals.
     * @result The expected length of the input signals.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    int getInputSignalLength() const;
    /*!
     * @brief Gets the number of correlograms.
     * @result The number of correlograms.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    int getNumberOfCorrelograms() const;
    /*!
     * @brief Gets the number of samples in the correlgoram.
     * @result The number of samples in each correlogram.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    int getCorrelogramLength() const;
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Flag indicating that all signals have been set.
     * @result True indicates that all signals have been set.
     */
    bool haveAllInputSignals() const noexcept;
    /*!
     * @brief Flag indicating that the correlograms have been computed.
     * @result True indicates that the correlograms have been computed.
     */
    bool haveCorrelograms() const noexcept;
    /*! @} */
private:
    class CorrelogramsImpl;
    std::unique_ptr<CorrelogramsImpl> pImpl;   
};
}
#endif
