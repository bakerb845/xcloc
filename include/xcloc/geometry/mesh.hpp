#ifndef XCLOC_GEOMETRY_MESH_HPP
#define XCLOC_GEOMETRY_MESH_HPP
#include "xcloc/geometry/enums.hpp"
namespace XCLoc::Geometry
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
     * @brief Default destructor.
     */
    virtual ~IMesh() = default;
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
     * @brief Gets the mesh type.
     * @result The mesh type.
     */
    virtual XCLoc::Geometry::MeshType getMeshType() const = 0;
};
}
#endif
