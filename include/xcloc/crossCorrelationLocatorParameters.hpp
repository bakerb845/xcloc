#ifndef XCLOC_CROSSCORRELATIONLOCATORPARAMETERS_HPP
#define XCLOC_CROSSCORRELATIONLOCATORPARAMETERS_HPP 1
#include <memory>
namespace XCLoc
{
class WaveformIdentifier;
/*!
 * @brief Defines the parameters for the cross-correlation-based location.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class CrossCorrelationLocationParameters
{
public:
    CrossCorrelationLocationParameters();

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~CrossCorrelationLocationParameters();
    /*!
     * @brief Clears the memory then resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Sets the cross-correlation pairs.
     * @param[in] xcPairs
     * @throws std::invalid_argument if xcPairs has no elements.
     * @note The waveform identifiers given here correspond to the travel time
     *       table names.
     */
    void setCorrelationTable(const std::vector<std::pair<WaveformIdentifier, WaveformIdentifier>> &xcPairs); 

    void addTravelTimeTable(const WaveformIdentifier &waveid,
                            const TravelTimeTable<T> &table);
    void clearTravelTimeTables() noexcept;

    /*!
     * @brief Determines whether or not.
     * @result True indicates that this is a valid parameters class.
     */
    bool isValid() const noexcept;
private:
    class CrossCorrelationLocationParametersImpl;
    std::unique_ptr<CrossCorrelationLocationParametersImpl> pImpl;
};
}
#endif
