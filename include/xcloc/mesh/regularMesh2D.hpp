#ifndef XCLOC_MESH_REGULARMESH2D_HPP
#define XCLOC_MESH_REGULARMESH2D_HPP
#include <memory>
#include "xcloc/mesh/mesh.hpp"
namespace XCLoc::Mesh
{
/*!
 * @brief Defines a structured grid.
 */
template<class T = double>
class RegularMesh2D : public IMesh<T>
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    RegularMesh2D(); 
    /*!
     * @brief Copy constructor.
     * @param[in] mesh  The mesh from which to initialize this class.
     */
    RegularMesh2D(const RegularMesh2D &mesh);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] mesh  The mesh to copy.
     * @result A deep copy of the mesh.
     */
    RegularMesh2D& operator=(const RegularMesh2D &mesh);
    /*!
     * @brief Clone self.
     * @result A deep copy of this.
     */
    std::unique_ptr<IMesh<T>> clone() const override;
    /*! @} */
    /*! @name Destructors 
     * @{
     */
    /*!
     * @brief Destructor.
     */
    virtual ~RegularMesh2D();
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
     * @throws std::runtime_error if \c setNumberOfGridPointsInX()
     *         or \c setNumberOfGridPointsInZ() was not called.
     * @sa \c haveGridDimensions()
     */
    int getNumberOfCells() const override;
    /*!
     * @brief Gets the number of grid points in the mesh.
     * @result The number of grid points in the mesh.
     * @throws std::runtime_error if \c setNumberOfGridPointsInX() or
     *         \c setNumberOfGridPointsInZ() was not called.
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

    /*! @name Origin
     * @{
     */
    /*!
     * @brief Sets the x origin.
     * @param[in] x0  The x origin in meters.
     */
    void setOriginInX(double x0) noexcept;
    /*!
     * @brief Gets the x origin.
     * @result The x origin in meters.
     */
    double getOriginInX() const noexcept;
    /*!
     * @brief Sets the z origin.
     * @param[in] z0  The z origin in meters.
     */
    void setOriginInZ(double z0) noexcept;
    /*!
     * @brief Gets the z origin.
     * @result The z origin in meters.
     */
    double getOriginInZ() const noexcept;
    /*!@ } */

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
     * @throws std::runtime_error if the number of grid points in x and z
     *         is not set.
     * @throws std::invalid_argument if field is NULL.
     * @note If the field already exists then it will be overwritten.
     * @sa \c haveGridDimensions(), \c haveCellularScalarField()
     */
    void setCellularScalarField(const std::string &fieldName,
                                int ncell,
                                const T field[],
                                const RegularMesh2DOrderingType order);
    /*!
     * @brief Gets a pointer to the cell-based scalar field data.
     * @param[in] fieldName  The name of the scalar field.
     * @result A pointer to the cell-based scalar field.  This is an array whose
     *         dimensions is [getNumberOfCells()].  Its ordering is defined
     *         by \c getCellularScalarFieldOrdering().
     * @throws std::invalid_argument if fieldName's data has not been set.
     * @sa \c haveCellularScalarField(), \c getNumberOfCells(),
     *     \c getCellularScalarFieldOrdering()
     */
    const T *getCellularScalarFieldPointer(
        const std::string &fieldName) const override;
    /*!
     * @brief Determines whether or not the cell-based field was set.
     * @param[in] fieldName  The name of the cell-based scalar field.
     * @result True indicates that the cell-based field exists.
     */
    bool haveCellularScalarField(const std::string &fieldName) const noexcept override;
    /*!
     * @brief Gets the ordering of the cell-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @note This may not match your input format.
     * @sa \c haveCellularScalarField()
     */
    RegularMesh2DOrderingType getCellularScalarFieldOrdering(
        const std::string &fieldName) const;
    /*!
     * @brief Gets the max index and value of a cell-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The maximum value (result.first) and the cell index
     *         (result.second) in the field at which the max was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    std::pair<T, int> getCellularScalarFieldMaxValueAndIndex(
        const std::string &fieldName) const override;
    /*!
     * @brief Gets the min index and value of a cell-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The minimum value (result.first) and the cell index
     *         (result.second) in the field at which the min was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    std::pair<T, int> getCellularScalarFieldMinValueAndIndex(
        const std::string &fieldName) const override;
    /*!
     * @brief Converts a cell index to the cell-based pair (icellx, icellz).
     * @param[in] index    The cell index.  This must be in the range
     *                     [0, \c getNumberOfCells() - 1].
     * @param[out] icellx  The cell index in x.
     * @param[out] icellz  The cell index in z.
     * @throws std::invalid_argument if the index is out of range. 
     * @throws std::runtime_error if the number of grid points in x,
     *         y, or z was not set.
     * @sa \c getNumberOfCells()
     */
    void convertCellIndexToGrid(int index, int *icellx, int *icellz) const;
    /*!
     * @brief Converts a cell index to position pair (x,z)
     * @param[in] index    The cell index.  This must be in the range
     *                     [0, \c getNumberOfCells() - 1].
     * @param[out] x       The cell's position in x (meters).
     * @param[out] z       The cell'z position in z (meters).
     * @throws std::invalid_argument if the index is out of range.
     * @throws std::runtime_error if the number of grid points in x
     *         or z was not set or the grid spacing in x or z was not set.
     */
    void convertCellIndexToPosition(int index, double *x, double *z) const;
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
     * @throws std::runtime_error if the number of grid points in x and z
     *         is not set.
     * @throws std::invalid_argument if field is NULL.
     * @note If the field already exists then it will be overwritten.
     * @sa \c haveGridDimensions(), \c haveNodalScalarField()
     */
    void setNodalScalarField(const std::string &fieldName,
                             int ngrd,
                             const T field[],
                             const RegularMesh2DOrderingType order);
    /*!
     * @brief Gets a pointer to the node-based scalar field data.
     * @param[in] fieldName  The name of the scalar field.
     * @result A pointer to the node-based scalar field.  This is an array whose
     *         dimensions is [getNumberOfGridPoints()].  Its ordering is 
     *         defined by \c getNodalScalarFieldOrdering().
     * @throws std::invalid_argument if fieldName's data has not been set.
     * @sa \c haveNodalScalarField(), \c getNumberOfGridPoints(),
     *     \c getNodalScalarFieldOrdering()
     */
    const T *getNodalScalarFieldPointer(
        const std::string &fieldName) const override;
    /*!
     * @brief Determines whether or not the nodal-based field was set.
     * @param[in] fieldName  The name of the nodal-based scalar field.
     * @result True indicates that the nodal-based field exists.
     */
    bool haveNodalScalarField(const std::string &fieldName) const noexcept override;
    /*!
     * @brief Gets the ordering of the node-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @note This may not match your input format.
     * @sa \c haveNodalScalarField() 
     */
    RegularMesh2DOrderingType getNodalScalarFieldOrdering(
        const std::string &fieldName) const;
    /*!
     * @brief Gets the min index and value of a node-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The minimum value (result.first) and node index (result.second)
     *         in the field at which the min was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    std::pair<T, int> getNodalScalarFieldMinValueAndIndex(
        const std::string &fieldName) const override;
    /*!
     * @brief Gets the max index and value of a node-based scalar field.
     * @param[in] fieldName  The name of the scalar field.
     * @result The maximum value (result.first) and node index (result.second)
     *         in the field at which the max was found.
     * @throws std::invalid_argument if the scalar field does not exist.
     * @sa \c haveCellularScalarField().
     */
    std::pair<T, int> getNodalScalarFieldMaxValueAndIndex(
        const std::string &fieldName) const override;
    /*!
     * @brief Converts a grid index to the grid-based pair (ix,iz).
     * @param[in] index  The grid index.  This must be in the range
     *                   [0, \c getNumberOfGridPoints() - 1].
     * @param[out] ix    The grid index in x.
     * @param[out] iz    The grid index in z.
     * @throws std::invalid_argument if the index is out of range.
     * @throws std::runtime_error if the number of grid points in x,
     *         y, or z was not set.
     * @sa \c getNumberOfGridPoints()
     */
    void convertNodeIndexToGrid(int index, int *ix, int *iz) const;
    /*!
     * @brief Converts a grid index to position pair (x,z)
     * @param[in] index    The cell index.  This must be in the range
     *                     [0, \c getNumberOfGridPoints() - 1].
     * @param[out] x       The grid index's position in x (meters).
     * @param[out] z       The grid index's position in z (meters).
     * @throws std::invalid_argument if the index is out of range.
     * @throws std::runtime_error if the number of grid points in x
     *         or z was not set or the grid spacing in x or z
     *         was not set.
     */
    void convertGridIndexToPosition(int index, double *x, double *z) const;
    /*! @} */
    
    /*!
     * @brief Gets the mesh type.
     * @result The mesh type.
     */
    XCLoc::Mesh::MeshType getMeshType() const noexcept override;
private:
    class RegularMeshImpl;
    std::unique_ptr<RegularMeshImpl> pImpl;
};
}
#endif
