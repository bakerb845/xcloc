#ifndef XCLOC_TRAVELTIMETABLE_HPP
#define XCLOC_TRAVELTIMETABLE_HPP
#include <memory>
#include "xcloc/mesh/enums.hpp"
namespace XCLoc
{
namespace Mesh
{
template<class T>class RegularMesh2D;
template<class T>class RegularMesh3D;
template<class T>class IMesh;
}
/*!
 * @class TravelTimeTable "travelTimeTable.hpp" "xcloc/travelTimeTable.hpp"
 * @brief Defines a travel time table which is a tabulation of travel times
 *        from many points to a given point (i.e., candidate source points
 *        to a receiver point.).
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class TravelTimeTable
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    TravelTimeTable();
    /*!
     * @brief Copy constructor.
     * @param[in] table  The travel time table class from which to initialize
     *                   this class. 
     */ 
    TravelTimeTable(const TravelTimeTable &table);
    /*!
     * @brief Move constructor.
     * @param[in,out] table  The travel time table class from which to
     *                       initialize this class.  On exit, table's 
     *                       behavior is undefined.
     */ 
    TravelTimeTable(TravelTimeTable &&table) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] table  The travel time table to copy.
     * @result A deep coyp of the travel time table.
     */
    TravelTimeTable& operator=(const TravelTimeTable &table);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] table  The travel time table's whose memory will be
     *                       copied to this.  On exit, table's behavior is
     *                       undefined.
     * @result The memory moved from table to this.
     */
    TravelTimeTable& operator=(TravelTimeTable &&table) noexcept;
    /*! @} */
    /*!
     * @name Destructor
     * @{
     */
    ~TravelTimeTable();
    /*!
     * @brief Removes the travel time table currently in memory and resets the
     *        class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Travel Time Table
     * @{
     */
    /*!
     * @brief Sets the travel time table.
     * @param[in] fieldName  The name of the travel time field.
     * @param[in] mesh       The mesh on which the travel time field is stored. 
     * @throws std::invalid_argument if the travel time field cannot be
     *         found or has no points.
     */
    void setTravelTimeTable(const std::string &fieldName,
                            const Mesh::IMesh<T> &mesh);
    /*!
     * @brief Determines if the nodal field was set or not.
     * @result True indicates that the travel time field is set.
     */
    bool haveTravelTimeTable() const noexcept;
    /*!
     * @brief Gets a pointer to the travel time field.
     * @result A pointer to the travel time field.  This is an array whose
     *         dimension is defined by \c getNumberOfPoints().
     */
    const T *getTravelTimeTablePointer() const;
    /*!
     * @brief Gets the number of points in the travel time field. 
     * @result The number of points in the travel time field.
     * @throws std::runtime_error if the travel time field was never set.
     * @sa \c haveTravelTimeField()
     */
    int getNumberOfPoints() const;
    /*!
     * @brief Determines whether or not the field node-based or cell-based.
     * @result True indicates the field is node-based. 
     * @throws std::runtime_error if the travel time field was never set.
     * @sa \c haveTravelTimeField()
     */
    bool isNodalField() const;
    /*! @} */

    /*! @name Geometry
     * @{
     */
    /*!
     * @brief Determines the underlying mesh type.
     * @result The underlying mesh type.
     * @throws std::runtime_error if the travel time field was never set.
     * @sa \c haveTravelTimeField()
     */
    XCLoc::Mesh::MeshType getMeshType() const;
    /*! @} */
private:
    class TravelTimeTableImpl;
    std::unique_ptr<TravelTimeTableImpl> pImpl;
};
}
#endif
