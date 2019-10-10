#ifndef XCLOC_TRAVELTIMETABLES_HPP
#define XCLOC_TRAVELTIMETABLES_HPP
#include <vector>
#include <memory>
#include "xcloc/mesh/enums.hpp"
namespace XCLoc
{
template<class T> class TravelTimeTable; 
class TravelTimeTableName;
/*!
 * @class TravelTimeTables travelTimeTables.hpp "xcloc/travelTimeTables.hpp"
 * @brief Defines a collection of travel time tables.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class TravelTimeTables
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor 
     */
    TravelTimeTables();
    /*!
     * @brief Copy constructor.
     * @param[in] tables  The travel time tables from which to initialize
     *                    this class.
     */
    TravelTimeTables(const TravelTimeTables &tables);
    /*!
     * @brief Move constructor.
     * @param[in,out] tables  The travel time tables from which to initialize
     *                        this class.  On exit, tables's behavior is
     *                        undefined.
     */
    TravelTimeTables(TravelTimeTables &&tables) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] tables  The tables class to copy to this.
     * @result A deep copy of tables.
     */
    TravelTimeTables& operator=(const TravelTimeTables &tables);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] tables  The tables class whose memory will be moved to
     *                        this.  On exit, tables's behavior is undefined.
     * @result The memory from tables moved to this. 
     */
    TravelTimeTables& operator=(TravelTimeTables &&tables) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~TravelTimeTables();
    /*!
     * @brief Releases all tables and clears memory.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Tables
     * @{
     */
    /*!
     * @brief Adds a travel time table that predicts the arrival time of a
     *        given phase to the station. 
     * @param[in] tableName  The name of the table.
     * @throws std::invalid_argument if the table name is not set, the
     *         travel time table is not set on table, or if the table's
     *         geometry does not match the geometry of the first table.
     * @note The first table set will define the geometry.
     * @note If a table corresponding to the network, station, and phase 
     *       exists then it will be overwritten.
     */
    void addTable(const TravelTimeTableName &tableName,
                  const TravelTimeTable<T> &table);
    /*!
     * @brief Gets the travel time table for the given phase and station.
     * @param[in] tableName  The name of the travel time table to extract.
     * @throws std::invalid_argument if the table name does not exist.
     */
    TravelTimeTable<T> getTable(const TravelTimeTableName &tableName) const;
    /*!
     * @brief Gets a pointer to the travel time field.
     * @param[in] tableName  The name of the travel time table to find.
     * @result A pointer to the travel time field.  This is an array whose
     *         dimension is defined by \c getNumberOfPoints().
     * @throws std::invalid_argument if the table name cannot be found.
     */
    const T *getTravelTimeTablePointer(const TravelTimeTableName &tableName) const;
    /*!
     * @brief Gets the name of all the travel time table names.
     * @result A vector with all the travel time table names.
     */
    std::vector<TravelTimeTableName> getTableNames() const;
    /*!
     * @brief Determines if the table for the given network, station, and
     *        phase exists.
     * @param[in] tableName  The name of the table.
     * @throws std::invalid_argument if tableName is not valid.
     * @result True indicates that the table exists.
     */
    bool haveTable(const TravelTimeTableName &tableName) const;
    /*!
     * @brief Returns the number of tables.
     * @result The number of travel time tables.
     */
    int getNumberOfTables() const noexcept;
    /*! @} */

    /*! @name Geometry
     * @{
     */
    /*!
     * @brief Determines the underlying mesh type.
     * @result The underlying mesh type.
     * @throws std::runtime_error if there are no travel time tables set.
     * @sa \c getNumberOfTables()
     */
    XCLoc::Mesh::MeshType getMeshType() const;
    /*!
     * @brief Determines if the underlying travel time field is node-based
     *        or cell-based.
     * @result True indicates that the underlying field is of nodal type.
     * @throws std::runtime_error if there are no travel time tables set.
     * @sa \c getNumberOfTables()
     */ 
    bool isNodalField() const;
    /*!
     * @brief Gets the number of points in the travel time field.
     * @result The number of points in the travel time field.
     * @throws std::runtime_error if the travel time field was never set.
     * @sa \c haveTravelTimeField()
     */
    int getNumberOfPoints() const;
    /*! @} */
private:
    class TravelTimeTablesImpl;
    std::unique_ptr<TravelTimeTablesImpl> pImpl;
};
}
#endif
