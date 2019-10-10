#ifndef XCLOC_CORRELOGRAMDSMPARAMETERS_HPP
#define XCLOC_CORRELOGRAMDSMPARAMETERS_HPP
#include <memory>
namespace XCLoc
{
/*!
 * @class CorrelogramDSMParameters "correlogramDSMParameters.hpp" "xcloc/correlogramDSMParameters.hpp"
 * @brief Defines the parameters for the diffraction stack migration of the
 *        correlograms.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class CorrelogramDSMParameters
{
public:
    /*!
     * @brief Constructor.
     */
    CorrelogramDSMParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  Parameters class from which to initialize
     *                        this class.
     */
    CorrelogramDSMParameters(const CorrelogramDSMParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  Parameters class from which to initialize
     *                            this class.  On exit, parameters' behavior
     *                            is undefined.
     */
    CorrelogramDSMParameters(CorrelogramDSMParameters &&parameters) noexcept;

    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters  The parameters class to copy.
     * @result A deep copy of the parameters class.
     */
    CorrelogramDSMParameters& operator=(const CorrelogramDSMParameters &parameters);
    /*!
     * @brief Move assignment operator.
     * @param[in] parameters  The parameters class whose memory is will be
     *                        moved to this.  On exit, parameters' behavior is
     *                        undefined.
     * @result The moved memory from parameters to this.
     */
    CorrelogramDSMParameters& operator=(CorrelogramDSMParameters &&parameters) noexcept;

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~CorrelogramDSMParameters();
    /*!
     * @brief Erases and resets the parameters.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Cache-block Size
     * @{
     */
    /*!
     * @brief To ease memory pressure the migration updates a subset of the 
     *        migration image.
     * @param[in] cacheBlockSize  The number of bytes in the cache-block.
     *                            You may consider starting with 4096.
     * @throws std::invalid_argument if the cacheBlockSize is less than 64
     *         or the blockSize is not a power of 2.
     */
    void setCacheBlockSize(size_t cacheBlockSize);
    /*!
     * @brief Gets the number of bytes in the cache block.
     * @result The cache block size.
     */
    size_t getCacheBlockSize() const noexcept; 
    /*! @} */

    /*!
     * @brief Determines whether or not the parameters class is valid.
     */
    bool isValid() const noexcept;
private:
    class CorrelogramDSMParametersImpl;
    std::unique_ptr<CorrelogramDSMParametersImpl> pImpl;
};
}
#endif
