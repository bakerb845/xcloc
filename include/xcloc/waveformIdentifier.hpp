#ifndef XCLOC_WAVEFORMIDENTIFIER_HPP
#define XCLOC_WAVEFORMIDENTIFIER_HPP
#include <string>
#include <memory>
namespace XCLoc
{
class WaveformIdentifier
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    WaveformIdentifier();
    /*!
     * @brief Constructor.
     * @param[in] network       The network code.
     * @param[in] station       The station name.
     * @param[in] polarization  The signal polarization, e.g. "P" or "S".
     */
    WaveformIdentifier(const std::string &network,
                       const std::string &station,
                       const std::string &phase);
    /*!
     * @brief Copy constructor.
     * @param[in] waveid  The waveform identifier class from which to
     *                    initialize this class.
     */
    WaveformIdentifier(const WaveformIdentifier &waveid);
    /*!
     * @brief Move constructor.
     * @param[in,out] waveid  The waveform identifier class from which to
     *                        initialize this class.  On exit, waveid's
     *                        behavior is undefined.
     */
    WaveformIdentifier(WaveformIdentifier &&waveid) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] waveid  The wave identifier class to copy to this.
     * @result A deep copy of waveid.
     */ 
    WaveformIdentifier& operator=(const WaveformIdentifier &waveid);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] waveid  The wave identifier class whose memory will be
     *                        moved to this.  On exit, waveid's behavior is
     *                        undefined.
     * @result The memory from waveid moved to this.
     */
    WaveformIdentifier& operator=(WaveformIdentifier &&waveid) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~WaveformIdentifier();
    /*!
     * @brief Clears the memory from the class and resets the values.
     */
    void clear() noexcept;
    /*! @} */

    void setNetwork(const std::string &network);
    std::string getNetwork() const noexcept;

    void setStation(const std::string &station);
    std::string getStation() const noexcept;

    void setPolarization(const std::string &polarization);
    std::string getPolarization() const noexcept;

    int64_t getWaveformIdentifier() const;
private:
    class WaveformIdentifierImpl;
    std::unique_ptr<WaveformIdentifierImpl> pImpl;
};

 /*!
  * @brief Tests if this lhs is greater than rhs.
  * @result True indicates that lhs is greater than rhs.
  * @note Comparisons are done with hash values.
  */
 bool operator>(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs);
 /*!
  * @brief Tests if this lhs is less than rhs.
  * @result True indicates that lhs is less than rhs.
  * @note Comparisons are done with hash values.
  */
bool operator<(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs);
 /*!
  * @brief Tests if lhs is equal to rhs.
  * @result True indicates that lhs is equal to rhs.
  */
bool operator==(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs);
 /*!
  * @brief Tests if lhs is not equal to rhs.
  * @result True indicates that lhs is not equal to rhs.
  */
bool operator!=(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs);
}
#endif
