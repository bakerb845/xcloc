#ifndef XCLOC_CORRELOGRAM_POSTPROCESSORPARAMETERS_HPP
#define XCLOC_CORRELOGRAM_POSTPROCESSORPARAMETERS_HPP
#include <memory>
#include "xcloc/enums.hpp"
namespace XCLoc
{
class CorrelogramPostProcessorParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    CorrelogramPostProcessorParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  The correlogram post-processor parameters from
     *                        which to initialize this class.
     */
    CorrelogramPostProcessorParameters(
        const CorrelogramPostProcessorParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The correlogram post-processor parameters from
     *                            which to initialize this class.  On exit
     *                            parameters' behavior is undefined.
     */
    CorrelogramPostProcessorParameters(
        CorrelogramPostProcessorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters  The parameters to copy.
     * @result A deep copy of the input parameters.
     */
    CorrelogramPostProcessorParameters&
        operator=(const CorrelogramPostProcessorParameters &parameters);
    /*!
     * @name Move assignment operator.
     * @param[in,out] parameters  The parameters whose memory will be moved to
     *                            this.  On exit, parameters' behavior is
     *                            undefined. 
     * @result The memory from parameters moved to this.
     */
    CorrelogramPostProcessorParameters&
        operator=(CorrelogramPostProcessorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~CorrelogramPostProcessorParameters();
    /*!
     * @brief Resets all variables and reverts to no filtering of correlograms.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Filtering strategies
     * @{
     */
    /*!
     * @brief Disables filtering.  This should only be when travel times are
     *        computed in very accurate velocity models.  Otherwise, 
     *        destructive interference can severely degrade the imaging.
     */
    void setNoFiltering() noexcept;
    /*!
     * @brief The correlogram post-processor will compute the
     * @param[in] envelopeLength  The FIR envelope length.  Since the FIR
     *                            delay is removed analytically the length
     *                            must be an odd number.  If envelopeLength
     *                            is even then envelopeLength will internally
     *                            be set to envelopeLength + 1.
     * @throws std::invalid_argument if the envelope length is not positive. 
     */
    void setFIREnvelopeFiltering(int envelopeLength);
    /*!
     * @brief Gets the FIR-based envelope filter length.
     * @result The FIR filter length.
     * @throws std::runtime_error if the filtering strategy is no
     *         CorrelogramFilteringType::FIR_ENVELOPE_FILTERING.
     * @sa \c getFilteringType().
     */
    int getFIREnvelopeFilterLength() const;
    /*! @} */

    /*!
     * @brief Gets the correlogram filtering strategy.
     * @result The filtering strategy.
     * @throws std::runtime_error if the filtering strategy was not set.
     * @sa \c isValid()
     */
    CorrelogramFilteringType getFilteringType() const;
    /*!
     * @brief Flag indicating that the post-processing strategy is valid.
     */
    bool isValid() const noexcept;
private:
    class CorrelogramPostProcessorImpl;
    std::unique_ptr<CorrelogramPostProcessorImpl> pImpl;
};
}
#endif
