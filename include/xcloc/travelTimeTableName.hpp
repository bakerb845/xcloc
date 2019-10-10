#ifndef XCLOC_TRAVELTIMETABLENAME_HPP
#define XCLOC_TRAVELTIMETABLENAME_HPP
#include <string>
#include <memory>
namespace XCLoc
{
/*!
 * @class TravelTimeTableName "travelTimeTableName.hpp" "xcloc/travelTimeTableName.hpp"
 * @brief Sets a travel time table name which is defined by the network code,
 *        the name of the station in the network, and the name of the seismic
 *        phase that the table predicts.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class TravelTimeTableName
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    TravelTimeTableName();
    /*!
     * @brief Copy constructor.
     * @param[in] tableName  The travel time table name from which to initialize
     *                       this class.
     */
    TravelTimeTableName(const TravelTimeTableName &tableName);
    /*!
     * @brief Move constructor.
     * @param[in,out] tableName  The travel time table name class whose memory
     *                           will be moved to this.  On exit, tableName's
     *                           behavior is undefined.
     */
    TravelTimeTableName(TravelTimeTableName &&tableName) noexcept;
    /*!
     * @brief Constructs a travel time table name by specifying the network,
     *        station, and phase.
     * @param[in] network  The network code.
     * @param[in] station  The station code.
     * @param[in] phase    The phase name.
     * @param[in] polarization  The polarization.  This can be "P" or "S".
     */
    TravelTimeTableName(const std::string &network,
                        const std::string &station,
                        const std::string &phase,
                        const std::string &polarization); 
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] tableName  The travel time table name class to copy.
     * @result A deep copy of tableName.
     */
    TravelTimeTableName& operator=(const TravelTimeTableName &tableName);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] tableName  The travel time table class whose memory will
     *                           be moved to this.  On exit, tableName's 
     *                           behavior is undefined.
     * @result The memory from tableName moved to this.
     */
    TravelTimeTableName& operator=(TravelTimeTableName &tableName) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~TravelTimeTableName();
    /*!
     * @brief Resets the class and clears all memory.
     */ 
    void clear() noexcept;
    /*! @} */

    /*! @name Network Name
     * @{
     */
    /*!
     * @brief Sets the name of the network.
     * @param[in] network  The name of the network to which the station belongs.
     */
    void setNetwork(const std::string &network) noexcept;
    /*!
     * @brief Gets the name of the network.
     * @result The name of the network.
     * @throws std::runtime_error if the network name was not set.
     * @sa \c haveNetwork()
     */
    std::string getNetwork() const;
    /*!
     * @brief Determines if the network was set.
     * @result True indicates that the network was set.
     */
    bool haveNetwork() const noexcept;
    /*! @} */

    /*! @name Station Name
     * @{
     */
    /*!
     * @brief Sets the name of the station.
     * @param[in] network  The name of the station in the seismic network.
     */
    void setStation(const std::string &station) noexcept;
    /*!
     * @brief Gets the name of the station.
     * @result The name of the station.
     * @throws std::runtime_error if the station name was not set.
     * @sa \c haveStation()
     */
    std::string getStation() const;
    /*!
     * @brief Determines if the station was set.
     * @result True indicates that the station was set.
     */
    bool haveStation() const noexcept;
    /*! @} */

    /*! @name Phase Name
     * @{
     */
    /*!
     * @brief Sets the name of the seismic phase.
     * @param[in] phase  The name of the seismic phase.
     */
    void setPhase(const std::string &phase) noexcept;
    /*!
     * @brief Gets the name of the seismic phase.
     * @result The name of the seismic phase.
     * @throws std::runtime_error if the phase name was not set.
     * @sa \c havePhase()
     */
    std::string getPhase() const;
    /*!
     * @brief Determines if the seismic phase was set.
     * @result True indicates that the phase name was set.
     */
    bool havePhase() const noexcept;
    /*! @} */

    /*! @name Polarization
     * @{
     */
    /*!
     * @brief The polarization so that this travel time table can be tied to
     *        the appropriate correlogram.  This can be P or S.
     * @param[in] polarization  The polarization.
     */
    void setPolarization(const std::string &polarization) noexcept;
    /*!
     * @brief Gets the polarization.
     * @result The polarization.
     * @throws std::runtime_error if the polarization was not set.
     * @sa \c havePhase()
     */
    std::string getPolarization() const;
    /*!
     * @brief Determines if polarizatin was set.
     * @result True indicates that the polarization was set.
     */
    bool havePolarization() const noexcept;
    /*! @} */

    /*!
     * @brief Determines if the network, station, and phase name are set.
     * @result True indicates that the network, station, the phase name,
     *         and polarization have been set and this is therefore a valid 
     *         travel time table name.
     */
    bool isValid() const noexcept;
private:
    class TravelTimeTableNameImpl;
    std::unique_ptr<TravelTimeTableNameImpl> pImpl;    
};
/*!
 * @brief Determines if two table names are equivalent, i.e., lhs == rhs.
 * @param[in] lhs  The left hand side in the equality test.
 * @param[in] rhs  The right hand side in the equality test.
 * @result True indicates that the two names are equivalent.
 */
bool operator==(const TravelTimeTableName &lhs, 
                const TravelTimeTableName &rhs) noexcept;
/*!
 * @brief Determines if two table names are not equivalent, i.e., lhs != rhs.
 * @param[in] lhs  The left hand side in the equality test.
 * @param[in] rhs  The right hand side in the equality test.
 * @result True indicates that the two names are not equivalent.
 */
bool operator!=(const TravelTimeTableName &lhs,
                const TravelTimeTableName &rhs)  noexcept;
}
#endif
