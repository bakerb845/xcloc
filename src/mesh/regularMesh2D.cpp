#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/align/aligned_allocator.hpp>
#include "xcloc/mesh/regularMesh2D.hpp"
#include "xcloc/mesh/enums.hpp"
#include "xcloc/enums.hpp"

using namespace XCLoc::Mesh;

template<class T>
class RegularMesh2D<T>::RegularMeshImpl
{
public:
    /// Destructor - class clear to actually perform deallocation
    ~RegularMeshImpl()
    {
        clear();
    }
    /// Releases memory on the map
    void clearMap()
    {
        mNodalScalarFields.clear();
        mCellularScalarFields.clear();
    }
    /// Releases memory on module
    void clear() noexcept
    {
        clearMap();
        mDeltaX = 0;
        mDeltaZ = 0;
        mGridPointsInX = 0;
        mGridPointsInZ = 0;
    }
//private:
    /// Holds the nodal-based scalar fields
    std::map<std::string,
             std::vector<T, boost::alignment::aligned_allocator<T, 64>>
            > mNodalScalarFields;
    /// Holds the cell-based scalar fields
    std::map<std::string,
             std::vector<T, boost::alignment::aligned_allocator<T, 64>>
            > mCellularScalarFields;
    /// Grid spacing in x (meters)
    double mDeltaX = 0;
    /// Grid spacing in z (meters)
    double mDeltaZ = 0;
    /// x origin (meters)
    double mOriginX = 0;
    /// z origin (meters)
    double mOriginZ = 0;
    /// Number of grid points in x.
    int mGridPointsInX = 0;
    /// Number of grid points in z.
    int mGridPointsInZ = 0;
    /// Defines the geometry orientation
    GeometryOrientationType mGeometryOrientation
        = GeometryOrientationType::EAST_NORTH_DOWN;
    /// Indicates the type of mesh that this is.
    const XCLoc::Mesh::MeshType mMeshType = XCLoc::Mesh::MeshType::REGULAR_MESH;
    /// Defines the mesh ordering
    const RegularMesh2DOrderingType mOrder = RegularMesh2DOrderingType::NX_NZ;
};

/// Constructor
template<class T>
RegularMesh2D<T>::RegularMesh2D() :
    pImpl(std::make_unique<RegularMeshImpl>())
{
}

/// Copy constructor
template<class T>
RegularMesh2D<T>::RegularMesh2D(const RegularMesh2D &mesh)
{
    *this = mesh;
}

/// Destructor
template<class T>
RegularMesh2D<T>::~RegularMesh2D() = default;

template<class T>
void RegularMesh2D<T>::clear() noexcept
{
    pImpl->clear();
}

/// Copy assignment
template<class T>
RegularMesh2D<T>& RegularMesh2D<T>::operator=(const RegularMesh2D &mesh)
{
    if (&mesh == this){return *this;}
    if (pImpl){clear();}
    pImpl = std::make_unique<RegularMeshImpl> (*mesh.pImpl);
    return *this;
}

/// Deep copy
template<class T>
std::unique_ptr<IMesh<T>> RegularMesh2D<T>::clone() const
{
    auto res = std::make_unique<RegularMesh2D> (*this);
    return res;
}

/// Number of grid points in x
template<class T>
void RegularMesh2D<T>::setNumberOfGridPointsInX(const int nx)
{
    pImpl->mGridPointsInX = 0;
    pImpl->clearMap();
    if (nx < 2)
    {
        throw std::invalid_argument("nx = " + std::to_string(nx)
                                  + " must be at least 2\n");
    }
    pImpl->mGridPointsInX = nx;
}

template<class T>
int RegularMesh2D<T>::getNumberOfGridPointsInX() const
{
    if (pImpl->mGridPointsInX < 1)
    {
        throw std::runtime_error("Number of x grid points not yet set\n");
    }
    return pImpl->mGridPointsInX;
}

/// Number of grid points in z
template<class T>
void RegularMesh2D<T>::setNumberOfGridPointsInZ(const int nz)
{
    pImpl->mGridPointsInZ = 0;
    pImpl->clearMap();
    if (nz < 2)
    {
        throw std::invalid_argument("nz = " + std::to_string(nz)
                                  + " must be at least 2\n");
    }
    pImpl->mGridPointsInZ = nz;
}

template<class T>
int RegularMesh2D<T>::getNumberOfGridPointsInZ() const
{
    if (pImpl->mGridPointsInZ < 1)
    {
        throw std::runtime_error("Number of z grid points not yet set\n");
    }
    return pImpl->mGridPointsInZ;
}

/// Indicates whether or not the grid dimensions are set
template<class T>
bool RegularMesh2D<T>::haveGridDimensions() const noexcept
{
    if (pImpl->mGridPointsInX < 1){return false;}
    if (pImpl->mGridPointsInZ < 1){return false;}
    return true;
}

template<class T>
int RegularMesh2D<T>::getNumberOfCells() const
{
    if (!haveGridDimensions())
    {
        throw std::runtime_error("Grid points not yet set\n");
    }
    int ncells = (pImpl->mGridPointsInX - 1) 
                *(pImpl->mGridPointsInZ - 1);
    return ncells;
}

template<class T>
int RegularMesh2D<T>::getNumberOfGridPoints() const
{
    if (!haveGridDimensions())
    {
        throw std::runtime_error("Grid points not yet set\n");
    }
    int ngrid = pImpl->mGridPointsInX
               *pImpl->mGridPointsInZ;
    return ngrid;
}

/// Grid spacing in x
template<class T>
void RegularMesh2D<T>::setGridSpacingInX(const double dx)
{
    pImpl->mDeltaX = 0;
    if (dx == 0)
    {
        throw std::invalid_argument("dx cannot be 0\n");
    }
    pImpl->mDeltaX = dx;
}

template<class T>
double RegularMesh2D<T>::getGridSpacingInX() const
{
    if (pImpl->mDeltaX == 0)
    {
        throw std::runtime_error("dx never set\n");
    }
    return pImpl->mDeltaX;
}

/// Grid spacing in z
template<class T>
void RegularMesh2D<T>::setGridSpacingInZ(const double dz)
{
    pImpl->mDeltaZ = 0;
    if (dz == 0)
    {
        throw std::invalid_argument("dz cannot be 0\n");
    }
    pImpl->mDeltaZ = dz;
}

template<class T>
double RegularMesh2D<T>::getGridSpacingInZ() const
{
    if (pImpl->mDeltaZ == 0)
    {
        throw std::runtime_error("dz never set\n");
    }
    return pImpl->mDeltaZ;
}

/// Grid origin 
template<class T>
double RegularMesh2D<T>::getOriginInX() const noexcept
{
    return pImpl->mOriginX;
}

template<class T>
void RegularMesh2D<T>::setOriginInX(const double x0) noexcept
{
    pImpl->mOriginX = x0;
}

template<class T>
double RegularMesh2D<T>::getOriginInZ() const noexcept
{
    return pImpl->mOriginZ;
}

template<class T>
void RegularMesh2D<T>::setOriginInZ(const double z0) noexcept
{
    pImpl->mOriginZ = z0;
}

/// Cell-based scalar field
template<class T>
void RegularMesh2D<T>::setCellularScalarField(
    const std::string &fieldName,
    const int ncell,
    const T field[],
    const RegularMesh2DOrderingType order)
{
    int ncellRef = getNumberOfCells(); // Throws
    if (ncell != ncellRef)
    {
        throw std::invalid_argument("ncell = " + std::to_string(ncell)
                                  + " must equal "
                                  + std::to_string(ncellRef) + "\n");
    }
    if (field == nullptr){throw std::invalid_argument("field is NULL\n");}
    std::vector<T, boost::alignment::aligned_allocator<T, 64>> fvec(ncell);
    // Need to permute to row major (nz, ny, nx)
    if (order == RegularMesh2DOrderingType::NZ_NX)
    {
        auto ncellx = getNumberOfGridPointsInX() - 1;
        auto ncellz = getNumberOfGridPointsInZ() - 1;
        for (int k=0; k<ncellz; ++k)
        {
            #pragma omp simd
            for (int i=0; i<ncellx; ++i)
            {
                auto igrd = k*ncellx + i; 
                auto jgrd = i*ncellz + k; 
                fvec[igrd] = field[jgrd];
            }
        }
    }
    else
    {
        std::copy(field, field+ncell, fvec.begin()); 
    }
    pImpl->mCellularScalarFields.insert_or_assign(fieldName, fvec);
}

template<class T>
const T* RegularMesh2D<T>::getCellularScalarFieldPointer(
    const std::string &fieldName) const
{
    if (!haveCellularScalarField(fieldName))
    {
        throw std::invalid_argument("cellular scalar field name = " + fieldName
                                  + " does not exist\n");
    }
    auto result = pImpl->mCellularScalarFields.find(fieldName);
    return result->second.data();
}

template<class T>
bool RegularMesh2D<T>::haveCellularScalarField(
    const std::string &fieldName) const noexcept
{
    //return pImpl->mCellularScalarFields.contains(fieldName);
    auto result = pImpl->mCellularScalarFields.find(fieldName);
    if (result != pImpl->mCellularScalarFields.end()){return true;}
    return false;
}

/// Get the cellular ordering in the scalar field
template<class T> XCLoc::Mesh::RegularMesh2DOrderingType
RegularMesh2D<T>::getCellularScalarFieldOrdering(
    const std::string &fieldName) const
{
    if (!haveCellularScalarField(fieldName))
    {
        throw std::invalid_argument("field " + fieldName
                                  + " does not exist\n");
    }
    return RegularMesh2DOrderingType::NX_NZ;
}

template<class T>
std::pair<T, int> RegularMesh2D<T>::getCellularScalarFieldMinValueAndIndex(
    const std::string &fieldName) const
{
    auto field = getCellularScalarFieldPointer(fieldName); // This throws
    int ncell = getNumberOfCells(); // This throws
    auto element = std::min_element(field, field+ncell);
    int index = std::distance(field, element);
    std::pair<T, int> result(*element, index);
    return result;
}

template<class T>
std::pair<T, int> RegularMesh2D<T>::getCellularScalarFieldMaxValueAndIndex(
    const std::string &fieldName) const
{
    auto field = getCellularScalarFieldPointer(fieldName); // This throws
    int ncell = getNumberOfCells(); // This throws
    auto element = std::max_element(field, field+ncell);
    int index = std::distance(field, element);
    std::pair<T, int> result(*element, index);
    return result;
}

/// Grid conversions
template<class T>
void RegularMesh2D<T>::convertCellIndexToGrid(
    const int index, int *icellx, int *icellz) const
{
    int ncell = getNumberOfCells(); // Throws
    if (index < 0 || index >= ncell)
    {
        throw std::invalid_argument("index = " + std::to_string(index)
                                  + " must be in range [0,"
                                  + std::to_string(ncell-1) + "]\n");
    }
    auto nx = pImpl->mGridPointsInX - 1;
    // Solve for slowest changing index (z) then fastest changing index (x)
    *icellz = index/nx;
    *icellx = index - *icellz*nx;
}

template<class T>
void RegularMesh2D<T>::convertCellIndexToPosition(
    const int index, double *x, double *z) const
{
    int ix, iz;
    convertCellIndexToGrid(index, &ix, &iz); // Throws
    auto dx = getGridSpacingInX(); // Throws
    auto dz = getGridSpacingInZ(); // Throws
    // Shift origins so that we finish half way between grid points
    auto x0 = getOriginInX() + dx/2;
    auto z0 = getOriginInZ() + dz/2;
    *x = x0 + static_cast<double> (ix)*dx;
    *z = z0 + static_cast<double> (iz)*dz;
}

/// Nodal-based scalar field
template<class T>
void RegularMesh2D<T>::setNodalScalarField(
    const std::string &fieldName,
    const int ngrd,
    const T field[],
    const RegularMesh2DOrderingType order)
{
    int ngrdRef = getNumberOfGridPoints(); // Throws
    if (ngrd != ngrdRef)
    {
        throw std::invalid_argument("ngrd = " + std::to_string(ngrd)
                                  + " must equal "
                                  + std::to_string(ngrdRef) + "\n");
    }
    if (field == nullptr){throw std::invalid_argument("field is NULL\n");}
    std::vector<T, boost::alignment::aligned_allocator<T, 64>> fvec(ngrd);
    // Need to permute to row major (nz, ny, nx)
    if (order == RegularMesh2DOrderingType::NZ_NX)
    {
        auto ngrdx = getNumberOfGridPointsInX();
        auto ngrdz = getNumberOfGridPointsInZ();
        for (int k=0; k<ngrdz; ++k)
        {
            #pragma omp simd
            for (int i=0; i<ngrdx; ++i)
            {
                auto igrd = k*ngrdx + i;
                auto jgrd = i*ngrdz + k;
                fvec[igrd] = field[jgrd];
            }
        }
    }
    else
    {
        std::copy(field, field+ngrd, fvec.begin());
    }
    pImpl->mNodalScalarFields.insert_or_assign(fieldName, fvec);
}

template<class T>
const T* RegularMesh2D<T>::getNodalScalarFieldPointer(
    const std::string &fieldName) const
{
    if (!haveNodalScalarField(fieldName))
    {
        throw std::invalid_argument("nodal scalar field name = " + fieldName
                                  + " does not exist\n");
    }
    auto result = pImpl->mNodalScalarFields.find(fieldName);
    return result->second.data();
}

template<class T>
bool RegularMesh2D<T>::haveNodalScalarField(
    const std::string &fieldName) const noexcept 
{
    //return pImpl->mNodalScalarFields.contains(fieldName);
    auto result = pImpl->mNodalScalarFields.find(fieldName);
    if (result != pImpl->mNodalScalarFields.end()){return true;}
    return false;
}

template<class T>
std::pair<T, int> RegularMesh2D<T>::getNodalScalarFieldMinValueAndIndex(
    const std::string &fieldName) const
{
    auto field = getNodalScalarFieldPointer(fieldName); // This throws
    int ngrd = getNumberOfGridPoints(); // This throws
    auto element = std::min_element(field, field+ngrd);
    int index = std::distance(field, element);
    std::pair<T, int> result(*element, index);
    return result;
}

template<class T>
std::pair<T, int> RegularMesh2D<T>::getNodalScalarFieldMaxValueAndIndex(
    const std::string &fieldName) const
{
    auto field = getNodalScalarFieldPointer(fieldName); // This throws
    int ngrd = getNumberOfGridPoints(); // This throws
    auto element = std::max_element(field, field+ngrd);
    int index = std::distance(field, element);
    std::pair<T, int> result(*element, index);
    return result;
}

/// Get the node ordering in the scalar field
template<class T> XCLoc::Mesh::RegularMesh2DOrderingType
RegularMesh2D<T>::getNodalScalarFieldOrdering(
    const std::string &fieldName) const
{
    if (!haveNodalScalarField(fieldName))
    {
        throw std::invalid_argument("field " + fieldName + " does not exist\n");
    }
    return RegularMesh2DOrderingType::NX_NZ;
}

/// Grid conversions
template<class T>
void RegularMesh2D<T>::convertNodeIndexToGrid(const int index,
                                              int *ix, int *iz) const
{
    int ngrd = getNumberOfGridPoints(); // Throws
    if (index < 0 || index >= ngrd)
    {
        throw std::invalid_argument("index = " + std::to_string(index)
                                  + " must be in range [0,"
                                  + std::to_string(ngrd-1) + "]\n");
    }
    auto nx = pImpl->mGridPointsInX;
    // Solve for slowest changing index (z) then fastest changing index (x)
    *iz = index/nx;
    *ix = index - *iz*nx;
}

template<class T>
void RegularMesh2D<T>::convertGridIndexToPosition(
    const int index, double *x, double *z) const
{
    int ix, iz;
    convertNodeIndexToGrid(index, &ix, &iz); // Throws
    auto dx = getGridSpacingInX(); // Throws
    auto dz = getGridSpacingInZ(); // Throws
    auto x0 = getOriginInX();
    auto z0 = getOriginInZ();
    *x = x0 + static_cast<double> (ix)*dx;
    *z = z0 + static_cast<double> (iz)*dz;
}

/// Gets the geometry type
template<class T>
XCLoc::Mesh::MeshType RegularMesh2D<T>::getMeshType() const noexcept
{
    return pImpl->mMeshType;
}

/// Class instantiation
template class XCLoc::Mesh::RegularMesh2D<double>;
template class XCLoc::Mesh::RegularMesh2D<float>;
template class XCLoc::Mesh::RegularMesh2D<int>;