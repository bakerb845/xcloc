#ifndef XCLOC_MESH_MESH_HPP
#define XCLOC_MESH_MESH_HPP
#include <memory>
#include "xcloc/mesh/enums.hpp"
namespace XCLoc::Mesh
{
/*!
 * @class IMesh "mesh.hpp" "xcloc/geometry/mesh.hpp"
 * @brief Abstract base class for defining a model geometry.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class IMesh
{
public:
    /*!
     * @brief Ensures abstract base class can copy itself.
     */
    virtual std::unique_ptr<IMesh<T>> clone() const = 0;
    /*!
     * @brief Default destructor.
     */
    virtual ~IMesh() = default;
    /*!
     * @brief Determines if the nodal scalar field is set.
     * @param[in] fieldName  The name of the scalar field.
     * @result True indicates that the field exists.
     */
    virtual bool haveNodalScalarField(const std::string &fieldName) const noexcept = 0;
    /*!
     * @brief Determines if the cell-based scalar field is set.
     * @param[in] fieldName  The name of the scalar field.
     * @result True indicates that the field exists.
     */
    virtual bool haveCellularScalarField(const std::string &fieldName) const noexcept = 0;
    /*!
     * @brief Gets a pointer to a nodal scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result A pointer to the nodal-based data.  This is an array whose
     *         dimension is [\c getNumberOfGridPoints()].
     */
    virtual const T *getNodalScalarFieldPointer(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets a pointer to a cellular scalar field.
     * @param[in] fieldName  The name of the cellular field.
     * @result A pointer to the cell-based data.  This is an array whose
     *         dimension is [\c getNumberOfCells()].
     */ 
    virtual const T *getCellularScalarFieldPointer(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets the number of cells.
     * @result The number of cells in the mesh.
     */
    virtual int getNumberOfCells() const = 0;
    /*!
     * @brief Gets the number of grid points.
     * @result The number of grid points in the mesh.
     */
    virtual int getNumberOfGridPoints() const = 0;
    /*!
     * @brief Gets the max index and value of a cell-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The maximum value (result.first) and the cell index
     *         (result.second) in the field at which the max was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    virtual std::pair<T, int> getCellularScalarFieldMaxValueAndIndex(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets the min index and value of a cell-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The minimum value (result.first) and the cell index
     *         (result.second) in the field at which the max was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    virtual std::pair<T, int> getCellularScalarFieldMinValueAndIndex(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets the min index and value of a nodal-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The minimum value (result.first) and the nodal index
     *         (result.second) in the field at which the min was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    virtual std::pair<T, int> getNodalScalarFieldMinValueAndIndex(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets the max index and value of a nodal-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The maximum value (result.first) and the nodal index
     *         (result.second) in the field at which the max was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    virtual std::pair<T, int> getNodalScalarFieldMaxValueAndIndex(
        const std::string &fieldName) const = 0;
    /*!
     * @brief Gets the mesh type.
     * @result The mesh type.
     */
    virtual XCLoc::Mesh::MeshType getMeshType() const = 0;
};
}
#endif
