#ifndef XCLOC_CORRELATIONENGINEPARAMETERS_HPP
#define XCLOC_CORRELATIONENGINEPARAMETERS_HPP 1
#include <vector>
#include <memory>
#include "xcloc/enums.hpp"

namespace XCLoc
{
class CorrelationEngineParameters;
/*!
 * @class CorrelationEngineParameters "correlationEngineParameters.hpp" "xcloc/correlationEngineParameters.hpp"
 * @brief This defines the parameters for the cross-correlation engine.
 */
class CorrelationEngineParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    CorrelationEngineParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  The parameters from which to initialize 
     *                        this class.
     */
    CorrelationEngineParameters(const CorrelationEngineParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The parameters class from which to initialize
     *                            this class.  On exit, parameters' behavior is
     *                            undefined.
     */
    CorrelationEngineParameters(CorrelationEngineParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in,out] parameters  The parameters to copy.
     * @result A deep copy of the parameters.
     */
    CorrelationEngineParameters& operator=(const CorrelationEngineParameters &parameters); 
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters  The parameters memory to be moved to this.
     *                            On exit, parameters' behavior is undefined.
     * @result The moved memory from parameters.
     */
    CorrelationEngineParameters& operator=(CorrelationEngineParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~CorrelationEngineParameters();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Number of samples
     * @{
     */
    /*!
     * @brief Sets the number of samples in each signal.
     * @param[in] nSamples  The number of samples in each signal.
     *                      In this instance the correlation lengths will be
     *                      2*nSamples - 1.  This can be overridden by calling
     *                      \c setNumberOfPaddedSamples().
     * @throws std::invalid_argument if this is not positive. 
     * @note This must be specified.
     */
    void setNumberOfSamples(int nSamples);
    /*!
     * @brief Returns the number of samples in each signal.
     * @throws std::runtime_error if this was never specified.
     */
    int getNumberOfSamples() const;
    /*! @} */

    /*! @name Correlation pairs
     * @{
     */
    /*!
     * @brief Sets the default cross-correlation pairs.
     *        By default this will generate the (nSignals*(nSignals-1))/2 table
     *        of correlolograms.
     * @param[in] nSignals            The number of input signals to correlate.
     * @param[in] doAutoCorrelations  If true then the cross-correlation pair
     *                                will include auto-correlations.
     * @throws std::invalid_argument If nSignals is less than 2.
     * @note This must be specified.
     */
    void setCorrelationPairs(const int nSignals,
                             const bool doAutoCorrelations = false);
    /*!
     * @brief Sets the correlation pairs.
     * @param[in] xcPairs   Each index defines the signal identifier
     *                      pairing defining a correlation pair.
     * @throws std::invalid_argument if xcPairs.size() is less than 1.
     * @note This must be specified.
     */
    void setCorrelationPairs(const std::vector<std::pair<int,int>> &xcPairs);
    /*!
     * @brief Gets the number of correlation pairs.
     * @result The number of correlation pairs.
     * @throws std::runtime_error if \c setCorrelationPairs was not called.
     */
    int getNumberOfCorrelationPairs() const;
    /*!
     * @brief Gets all of the correlation pairs.
     * @result The correlation pair table.
     * @throws std::invalid_argument if \c setCorrelationPairs was not called.
     */
    std::vector<std::pair<int, int>> getCorrelationPairs() const;
    /*!
     * @brief Gets the correlation pair for the ixc'th correlogram.
     * @param[in] ixc   The ixc'th correlogram.  This must be in the range
     *                  [0, \c getNumberOfCorrelationPairs() - 1].
     * @result The waveform ID's comprising the ixc'th correlation pair.
     */
    std::pair<int, int> getCorrelationPair(const int ixc) const;
    /*!
     * @brief Gets the unique signal IDs as defined in the correlation table.
     * @result A sorted list of unique signal IDs that are accounted for 
     *         in the cross-correlation table.  The length is defined by
     *         \c getNumberOfCorrelationPairs().
     * @throws std::runtime_error if \c setCorrelationPairs() was not called.
     */
    std::vector<int> getUniqueSignalIDs() const;
    /*! @} */

    /*! @name DFT zero padding length
     * @{
     */
    /*!
     * @brief Sets the number of zero padded samples in the input signals prior
     *        to transforming.
     * @param[in] nPadSamples   The number of samples to which to zero pad.
     *                          In this case, case the cross-correlation lengths
     *                          will be 2*(nPadSamples - 1).
     * @throws std::invalid_argument if this is less than 
     *         \c getNumberOfSamples().
     * @note This will be reset to nSamples if \c setNumberOfSamples() is 
     *       called after calling this.
     */
    void setNumberOfPaddedSamples(const int nPadSamples);
    /*!
     * @brief Gets the number zero padded samples in the input signals.
     * @result The number of samples to which to zero pad the input signals
     *         prior to transforming.
     * @throws std::runtime_error if \c setNumberOfSamples() was not called.
     */
    int getNumberOfPaddedSamples() const;
    /*! @} */

    /*! @name MKL Accuracy
     * @{
     */
    /*!
     * @brief Sets the accuracy of some MKL floating point computations.
     * @param[in] accuracy  The desired accuracy of floating point computations.
     */ 
    void setMKLFloatingPointAccuracy(const MKLFloatingPointAccuracy accuracy) noexcept;
    /*!
     * @brief Gets the accuracy of some MKL floating point computations.
     * @result The accuracy of some MKL floating point computations.
     */
    MKLFloatingPointAccuracy getMKLFloatingPointAccuracy() const noexcept; 
    /*! @} */

    /*!
     * @brief Checks if the parameters are valid.
     * @result True indicates that the parameters are valid.
     */
    bool isValid() const noexcept;
private:
    class CorrelationEngineParametersImpl;
    std::unique_ptr<CorrelationEngineParametersImpl> pImpl;
};
}
#endif
