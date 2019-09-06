#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/align/aligned_allocator.hpp>
#include "xcloc/mesh/regularMesh3D.hpp"
#include "xcloc/mesh/enums.hpp"
#include "xcloc/enums.hpp"

using namespace XCLoc::Mesh;

template<class T>
class RegularMesh3D<T>::RegularMeshImpl
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
        mDeltaY = 0;
        mDeltaZ = 0;
        mGridPointsInX = 0;
        mGridPointsInY = 0;
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
    /// Grid spacing in y (meters)
    double mDeltaY = 0;
    /// Grid spacing in z (meters)
    double mDeltaZ = 0;
    /// x origin (meters)
    double mOriginX = 0;
    /// y origin (meters)
    double mOriginY = 0;
    /// z origin (meters)
    double mOriginZ = 0;
    /// Number of grid points in x.
    int mGridPointsInX = 0;
    /// Number of grid points in y.
    int mGridPointsInY = 0;
    /// Number of grid points in z.
    int mGridPointsInZ = 0;
    /// Defines the geometry orientation
    GeometryOrientationType mGeometryOrientation
        = GeometryOrientationType::EAST_NORTH_DOWN;
    /// Indicates the type of mesh that this is.
    const XCLoc::Mesh::MeshType mMeshType = XCLoc::Mesh::MeshType::REGULAR_MESH;
    /// Defines the mesh ordering
    const RegularMesh3DOrderingType mOrder = RegularMesh3DOrderingType::NX_NY_NZ;
};

/// Constructor
template<class T>
RegularMesh3D<T>::RegularMesh3D() :
    pImpl(std::make_unique<RegularMeshImpl>())
{
}

/// Copy constructor
template<class T>
RegularMesh3D<T>::RegularMesh3D(const RegularMesh3D &mesh)
{
    *this = mesh;
}

/// Destructor
template<class T>
RegularMesh3D<T>::~RegularMesh3D() = default;

template<class T>
void RegularMesh3D<T>::clear() noexcept
{
    pImpl->clear();
}

/// Copy assignment
template<class T>
RegularMesh3D<T>& RegularMesh3D<T>::operator=(const RegularMesh3D &mesh)
{
    if (&mesh == this){return *this;}
    if (pImpl){clear();}
    pImpl = std::make_unique<RegularMeshImpl> (*mesh.pImpl);
    return *this;
}

/// Number of grid points in x
template<class T>
void RegularMesh3D<T>::setNumberOfGridPointsInX(const int nx)
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
int RegularMesh3D<T>::getNumberOfGridPointsInX() const
{
    if (pImpl->mGridPointsInX < 1)
    {
        throw std::runtime_error("Number of x grid points not yet set\n");
    }
    return pImpl->mGridPointsInX;
}

/// Number of grid points in y
template<class T>
void RegularMesh3D<T>::setNumberOfGridPointsInY(const int ny)
{
    pImpl->mGridPointsInY = 0;
    pImpl->clearMap();
    if (ny < 2)
    {
        throw std::invalid_argument("ny = " + std::to_string(ny)
                                  + " must be at least 2\n");
    }
    pImpl->mGridPointsInY = ny;
}

template<class T>
int RegularMesh3D<T>::getNumberOfGridPointsInY() const
{
    if (pImpl->mGridPointsInY < 1)
    {
        throw std::runtime_error("Number of y grid points not yet set\n");
    }
    return pImpl->mGridPointsInY;
}

/// Number of grid points in z
template<class T>
void RegularMesh3D<T>::setNumberOfGridPointsInZ(const int nz)
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
int RegularMesh3D<T>::getNumberOfGridPointsInZ() const
{
    if (pImpl->mGridPointsInZ < 1)
    {
        throw std::runtime_error("Number of z grid points not yet set\n");
    }
    return pImpl->mGridPointsInZ;
}

/// Indicates whether or not the grid dimensions are set
template<class T>
bool RegularMesh3D<T>::haveGridDimensions() const noexcept
{
    if (pImpl->mGridPointsInX < 1){return false;}
    if (pImpl->mGridPointsInY < 1){return false;}
    if (pImpl->mGridPointsInZ < 1){return false;}
    return true;
}

template<class T>
int RegularMesh3D<T>::getNumberOfCells() const
{
    if (!haveGridDimensions())
    {
        throw std::runtime_error("Grid points not yet set\n");
    }
    int ncells = (pImpl->mGridPointsInX - 1) 
                *(pImpl->mGridPointsInY - 1)
                *(pImpl->mGridPointsInZ - 1);
    return ncells;
}

template<class T>
int RegularMesh3D<T>::getNumberOfGridPoints() const
{
    if (!haveGridDimensions())
    {
        throw std::runtime_error("Grid points not yet set\n");
    }
    int ngrid = pImpl->mGridPointsInX
               *pImpl->mGridPointsInY
               *pImpl->mGridPointsInZ;
    return ngrid;
}

/// Grid spacing in x
template<class T>
void RegularMesh3D<T>::setGridSpacingInX(const double dx)
{
    pImpl->mDeltaX = 0;
    if (dx == 0)
    {
        throw std::invalid_argument("dx cannot be 0\n");
    }
    pImpl->mDeltaX = dx;
}

template<class T>
double RegularMesh3D<T>::getGridSpacingInX() const
{
    if (pImpl->mDeltaX == 0)
    {
        throw std::runtime_error("dx never set\n");
    }
    return pImpl->mDeltaX;
}

/// Grid spacing in y 
template<class T>
void RegularMesh3D<T>::setGridSpacingInY(const double dy)
{
    pImpl->mDeltaY = 0;
    if (dy == 0)
    {
        throw std::invalid_argument("dy cannot be 0\n");
    }
    pImpl->mDeltaY = dy;
}

template<class T>
double RegularMesh3D<T>::getGridSpacingInY() const
{
    if (pImpl->mDeltaY == 0)
    {
        throw std::runtime_error("dy never set\n");
    }
    return pImpl->mDeltaY;
}

/// Grid spacing in z
template<class T>
void RegularMesh3D<T>::setGridSpacingInZ(const double dz)
{
    pImpl->mDeltaZ = 0;
    if (dz == 0)
    {
        throw std::invalid_argument("dz cannot be 0\n");
    }
    pImpl->mDeltaZ = dz;
}

template<class T>
double RegularMesh3D<T>::getGridSpacingInZ() const
{
    if (pImpl->mDeltaZ == 0)
    {
        throw std::runtime_error("dz never set\n");
    }
    return pImpl->mDeltaZ;
}

/// Grid origin 
template<class T>
double RegularMesh3D<T>::getOriginInX() const noexcept
{
    return pImpl->mOriginX;
}

template<class T>
void RegularMesh3D<T>::setOriginInX(const double x0) noexcept
{
    pImpl->mOriginX = x0;
}

template<class T>
double RegularMesh3D<T>::getOriginInY() const noexcept
{
    return pImpl->mOriginY;
}

template<class T>
void RegularMesh3D<T>::setOriginInY(const double y0) noexcept
{
    pImpl->mOriginY = y0;
}

template<class T>
double RegularMesh3D<T>::getOriginInZ() const noexcept
{
    return pImpl->mOriginZ;
}

template<class T>
void RegularMesh3D<T>::setOriginInZ(const double z0) noexcept
{
    pImpl->mOriginZ = z0;
}

/// Cell-based scalar field
template<class T>
void RegularMesh3D<T>::setCellularScalarField(
    const std::string &fieldName,
    const int ncell,
    const T field[],
    const RegularMesh3DOrderingType order)
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
    if (order == RegularMesh3DOrderingType::NZ_NY_NX)
    {
        auto ncellx = getNumberOfGridPointsInX() - 1;
        auto ncelly = getNumberOfGridPointsInY() - 1;
        auto ncellz = getNumberOfGridPointsInZ() - 1;
        for (int k=0; k<ncellz; ++k)
        {
            for (int j=0; j<ncelly; ++j)
            {
                #pragma omp simd
                for (int i=0; i<ncellx; ++i)
                {
                    auto igrd = k*ncelly*ncellx + j*ncellx + i; 
                    auto jgrd = i*ncelly*ncellz + j*ncellz + k; 
                    fvec[igrd] = field[jgrd];
                }
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
const T* RegularMesh3D<T>::getCellularScalarFieldPointer(
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
bool RegularMesh3D<T>::haveCellularScalarField(
    const std::string &fieldName) const noexcept
{
    //return pImpl->mCellularScalarFields.contains(fieldName);
    auto result = pImpl->mCellularScalarFields.find(fieldName);
    if (result != pImpl->mCellularScalarFields.end()){return true;}
    return false;
}

/// Get the cellular ordering in the scalar field
template<class T> XCLoc::Mesh::RegularMesh3DOrderingType
RegularMesh3D<T>::getCellularScalarFieldOrdering(
    const std::string &fieldName) const
{
    if (!haveCellularScalarField(fieldName))
    {
        throw std::invalid_argument("field " + fieldName
                                  + " does not exist\n");
    }
    return RegularMesh3DOrderingType::NX_NY_NZ;
}

template<class T>
std::pair<T, int> RegularMesh3D<T>::getCellularScalarFieldMinValueAndIndex(
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
std::pair<T, int> RegularMesh3D<T>::getCellularScalarFieldMaxValueAndIndex(
    const std::string &fieldName) const
{
    auto field = getCellularScalarFieldPointer(fieldName); // This throws
    int ncell = getNumberOfCells(); // This throws
    auto element = std::max_element(field, field+ncell);
    int index = std::distance(field, element);
    std::pair<T, int> result(*element, index);
    return result;
}

/// Nodal-based scalar field
template<class T>
void RegularMesh3D<T>::setNodalScalarField(
    const std::string &fieldName,
    const int ngrd,
    const T field[],
    const RegularMesh3DOrderingType order)
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
    if (order == RegularMesh3DOrderingType::NZ_NY_NX)
    {
        auto ngrdx = getNumberOfGridPointsInX();
        auto ngrdy = getNumberOfGridPointsInY();
        auto ngrdz = getNumberOfGridPointsInZ();
        for (int k=0; k<ngrdz; ++k)
        {
            for (int j=0; j<ngrdy; ++j)
            {
                #pragma omp simd
                for (int i=0; i<ngrdx; ++i)
                {
                    auto igrd = k*ngrdy*ngrdx + j*ngrdx + i;
                    auto jgrd = i*ngrdy*ngrdz + j*ngrdz + k;
                    fvec[igrd] = field[jgrd];
                }
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
const T* RegularMesh3D<T>::getNodalScalarFieldPointer(
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
bool RegularMesh3D<T>::haveNodalScalarField(
    const std::string &fieldName) const noexcept 
{
    //return pImpl->mNodalScalarFields.contains(fieldName);
    auto result = pImpl->mNodalScalarFields.find(fieldName);
    if (result != pImpl->mNodalScalarFields.end()){return true;}
    return false;
}

template<class T>
std::pair<T, int> RegularMesh3D<T>::getNodalScalarFieldMinValueAndIndex(
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
std::pair<T, int> RegularMesh3D<T>::getNodalScalarFieldMaxValueAndIndex(
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
template<class T> XCLoc::Mesh::RegularMesh3DOrderingType
RegularMesh3D<T>::getNodalScalarFieldOrdering(
    const std::string &fieldName) const
{
    if (!haveNodalScalarField(fieldName))
    {
        throw std::invalid_argument("field " + fieldName + " does not exist\n");
    }
    return RegularMesh3DOrderingType::NX_NY_NZ;
}

/// Gets the geometry type
template<class T>
XCLoc::Mesh::MeshType RegularMesh3D<T>::getMeshType() const noexcept
{
    return pImpl->mMeshType;
}

/// Class instantiation
template class XCLoc::Mesh::RegularMesh3D<double>;
template class XCLoc::Mesh::RegularMesh3D<float>;
template class XCLoc::Mesh::RegularMesh3D<int>;
