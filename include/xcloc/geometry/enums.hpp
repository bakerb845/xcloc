#ifndef XCLOC_GEOMETRY_ENUMS_HPP
#define XCLOC_GEOMETRY_ENUMS_HPP
namespace XCLoc::Geometry
{
/*!
 * @brief Defines the mesh type.
 */
enum class MeshType
{
    REGULAR_MESH,         /*!< A regular mesh with fixed a grid spacing in x,
                               a fixed grid spacing in y, and a fixed grid
                               spacing in z. */
    UNSTRUCTURED_POINTS,  /*!< Unstructured points in (X,Y,Z) points. */
    UNSTRUCTURED_MESH     /*!< General unstructured mesh (X,Y,Z). */
};

/*!
 * @brief Defines the mesh orientation.
 */
enum class GeometryOrientationType
{
    EAST_NORTH_DOWN,   /*!< This is a common seismological coordinate system
                            where +x is positive east, +y is positive north,
                            and +z is positive down. */
    EAST_NORTH_UP,     /*!< Like EAST_NORTH_DOWN except +z is positive up.
                            This is a right handed coordinate system which
                            is useful for plotting software. */
};
}
#endif
