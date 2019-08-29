#ifndef XCLOC_GEOMETRY_REGULARMESH3D_HPP
#define XCLOC_GEOMETRY_REGULARMESH3D_HPP
#include <memory>
#include "xcloc/geometry/mesh.hpp"
namespace XCLoc::Geometry
{
/*!
 * @brief Defines a structured grid.
 */
template<class T>
class RegularMesh3D : public IMesh<T>
{
public:
    /*!
     * @brief Defines the storage ordering of the field in the regular mesh.
     */
    enum class FieldOrdering
    {
        NX_NY_NZ, /*!< The field is packed in row major order [nx, ny, nz]
                       where nz changes most rapidly and nx changes slowest. */
        NZ_NY_NX  /*!< The field is packed in row major order [nz, ny, nx]
                       where nx changes most rapidly and nz changes slowest. */
    };
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    RegularMesh3D(); 
    /*!
     * @brief Copy constructor.
     * @param[in] mesh  The mesh from which to initialize this class.
     */
    RegularMesh3D(const RegularMesh3D &mesh);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] mesh  The mesh to copy.
     * @result A deep copy of the mesh.
     */
    RegularMesh3D& operator=(const RegularMesh3D &mesh);
    /*! @} */
    /*! @name Destructors 
     * @{
     */
    /*!
     * @brief Destructor.
     */
    virtual ~RegularMesh3D();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Grid spacing
     * @{
     */
    /*!
     * @brief Sets the grid spacing in x.
     * @param[in] dx   The grid spacing in x in meters.
     * @throws std::invalid_argument if the grid spacing is 0.
     */
    void setGridSpacingInX(double dx);
    /*!
     * @brief Gets the grid spacing in x.
     * @result The grid spacing in x in meters.
     * @throws std::runtime_error if the grid spacing in x was not set. 
     */
    double getGridSpacingInX() const;
    /*!
     * @brief Sets the grid spacing in y.
     * @param[in] dy   The grid spacing in y in meters.
     * @throws std::invalid_argument if the grid spacing is 0.
     */
    void setGridSpacingInY(double dy);
    /*!
     * @brief Gets the grid spacing in y.
     * @result The grid spacing in y in meters.
     * @throws std::runtime_error if the grid spacing iy x was not set.
     */
    double getGridSpacingInY() const;
    /*!
     * @brief Sets the grid spacing in z.
     * @param[in] dz   The grid spacing in z in meters.
     * @throws std::invalid_argument if the grid spacing is 0.
     */
    void setGridSpacingInZ(double dz);
    /*!
     * @brief Gets the grid spacing in z.
     * @result The grid spacing in z in meters.
     * @throws std::runtime_error if the grid spacing in z was not set.
     */
    double getGridSpacingInZ() const;
    /*!@ } */

    /*! @name Grid size
     * @{
     */
    /*!
     * @brief Sets the number of grid points in x.
     * @param[in] nx  The number of grid points in x.
     * @throws std::invalid_argument if this is not at least 2.
     * @note This will release all internal arrays as this is a resize
     *       operation.
     */
    void setNumberOfGridPointsInX(int nx);
    /*!
     * @brief Gets the number of of grid points in x.
     * @throws std::runtime_error if this was not set.
     * @sa \c haveGridDimensions()
     */
    int getNumberOfGridPointsInX() const;
    /*!
     * @brief Sets the number of grid points in y.
     * @param[in] ny   The number of grid points in y.
     * @throws std::invalid_argument if this is not at least 2.
     * @note This will release all internal arrays as this is a resize
     *       operation.
     */
    void setNumberOfGridPointsInY(int ny);
    /*!
     * @brief Gets the number of of grid points in y.
     * @throws std::runtime_error if this was not set.
     * @sa \c haveGridDimensions()
     */
    int getNumberOfGridPointsInY() const;
    /*!
     * @brief Sets the number of grid points in z.
     * @param[in] nz   The number of grid points in z.
     * @throws std::invalid_argument if this is not at least 2.
     * @note This will release all internal arrays as this is a resize
     *       operation.
     */
    void setNumberOfGridPointsInZ(int nz);
    /*!
     * @brief Gets the number of of grid points in z.
     * @throws std::runtime_error if this was not set.
     * @sa \c haveGridDimensions()
     */
    int getNumberOfGridPointsInZ() const;
    /*!
     * @brief Gets the number of cells in the mesh.
     * @result The number of cells in the mesh.
     * @throws std::runtime_error if \c setNumberOfGridPointsInX(),
     *         \c setNumberOfGridPointsInY(), or \c setNumberOfGridPointsInZ()
     *         was not called.
     * @sa \c haveGridDimensions()
     */
    int getNumberOfCells() const override;
    /*!
     * @brief Gets the number of grid points in the mesh.
     * @result The number of grid points in the mesh.
     * @throws std::runtime_error if \c setNumberOfGridPointsInX(),
     *         \c setNumberOfGridPointsInY(), or \c setNumberOfGridPointsInZ()
     *         was not called.
     * @sa \c haveGridDimensions()
     */ 
    int getNumberOfGridPoints() const override;
    /*!
     * @brief Determines if the grid dimensions are set.
     * @result True indicates that it is safe to get the number of grid points
     *         in any or all directions.
     */
    bool haveGridDimensions() const noexcept;
    /*! @} */

    /*! @name Cell-based Scalar Field
     * @{
     */
    /*!
     * @brief Sets a cellular scalar field.
     * @param[in] fieldName  The name of the field.  If this field exists
     *                       then it will be overwritten.
     * @param[in] ncell      The number of cells.  This must match 
     *                       \c getNumberOfCells().
     * @param[in] field      The scalar field's value at each node.  This
     *                       is an array whose dimension [ncell].
     * @param[in] order      Defines the ordering of field.
     * @throws std::runtime_error if the number of grid points in x, y, and z
     *         is not set.
     * @throws std::invalid_argument if field is NULL.
     * @note If the field already exists then it will be overwritten.
     * @sa \c haveGridDimensions(), \c haveCellularScalarField()
     */
    void setCellularScalarField(const std::string &fieldName,
                                int ncell,
                                const T field[],
                                const FieldOrdering order);
    /*!
     * @brief Gets a pointer to the cell-based scalar field data.
     * @param[in] fieldName  The name of the scalar field.
     * @result A pointer to the cell-based scalar field.  This is an array whose
     *         dimensions is [getNumberOfCells()].
     * @throws std::invalid_argument if fieldName's data has not been set.
     * @sa \c haveCellularScalarField(), \c getNumberOfCells()
     */
    const T *getCellularScalarFieldPointer(
        const std::string &fieldName) const override;
    /*!
     * @brief Determines whether or not the cell-based field was set.
     * @param[in] fieldName  The name of the cell-based scalar field.
     * @result True indicates that the cell-based field exists.
     */
    bool haveCellularScalarField(const std::string &fieldName) const noexcept;
    /*! @} */

    /*! @name Nodal-based Scalar Field
     * @{
     */
    /*!
     * @brief Sets a nodal scalar field.
     * @param[in] fieldName  The name of the field.  If this field exists
     *                       then it will be overwritten.
     * @param[in] ngrd       The number of grid points.  This must match
     *                       \c getNumberOfGridPoints().
     * @param[in] field      The scalar field's value at each node.  This
     *                       is an array whose dimension [ngrd].
     * @param[in] order      Defines the ordering of field.
     * @throws std::runtime_error if the number of grid points in x, y, and z
     *         is not set.
     * @throws std::invalid_argument if field is NULL.
     * @note If the field already exists then it will be overwritten.
     * @sa \c haveGridDimensions(), \c haveNodalScalarField()
     */
    void setNodalScalarField(const std::string &fieldName,
                             int ngrd,
                             const T field[],
                             const FieldOrdering order);
    /*!
     * @brief Gets a pointer to the node-based scalar field data.
     * @param[in] fieldName  The name of the scalar field.
     * @result A pointer to the node-based scalar field.  This is an array whose
     *         dimensions is [getNumberOfGridPoints()].
     * @throws std::invalid_argument if fieldName's data has not been set.
     * @sa \c haveNodalScalarField(), \c getNumberOfGridPoints()
     */
    const T *getNodalScalarFieldPointer(
        const std::string &fieldName) const override;
    /*!
     * @brief Determines whether or not the nodal-based field was set.
     * @param[in] fieldName  The name of the nodal-based scalar field.
     * @result True indicates that the nodal-based field exists.
     */
    bool haveNodalScalarField(const std::string &fieldName) const noexcept;
    /*! @} */
    
    /*!
     * @brief Gets the mesh type.
     * @result The mesh type.
     */
    XCLoc::Geometry::MeshType getMeshType() const noexcept override;
private:
    class RegularMeshImpl;
    std::unique_ptr<RegularMeshImpl> pImpl;
};
}
#endif
