#ifndef XCLOC_VERSION_HPP 
#define XCLOC_VERSION_HPP
namespace XCLoc
{
/*!
 * @class Version version.hpp "xcloc/version.hpp"
 * @brief Defines XCLoc version information.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class Version
{
public:
    /*!
     * @brief Returns the API major version number.
     * @result The major version number.
     */
    static int getMajor() noexcept;
    /*!
     * @brief Returns the API minor version number.
     * @result The minor version number.
     */
    static int getMinor() noexcept;
    /*!
     * @brief Returns the API patch version number.
     * @result The patch version number.
     */
    static int getPatch() noexcept;
    /*!
     * @brief Returns the full version number as a string.
     * @result The full version number, e.g., "1.2.3".
     */
    static std::string getVersion() noexcept;
    /*!
     * @brief Determines if the version is greater than or equal to
     *        the current (major, minor, patch).
     * @param[in] major  The major version number.
     * @param[in] minor  The minor version number.
     * @param[in] patch  The patch number.
     * @result True indicates taht the version is at least equal to the
     *         given major, minor, patch.
     */
    static bool isAtLeast(int major, int minor, int patch) noexcept;
};
}
#endif
