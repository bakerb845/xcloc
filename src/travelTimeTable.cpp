#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "xcloc/travelTimeTable.hpp"
#include "xcloc/mesh/mesh.hpp"
#include "xcloc/mesh/regularMesh2D.hpp"
#include "xcloc/mesh/regularMesh3D.hpp"

using namespace XCLoc;

template<class T>
class TravelTimeTable<T>::TravelTimeTableImpl
{
public:
    TravelTimeTableImpl() = default;
    ~TravelTimeTableImpl() = default;
    TravelTimeTableImpl(const TravelTimeTable &table)
    {
        *this = table;
    }
    TravelTimeTableImpl& operator=(const TravelTimeTableImpl &table)
    {
        mTable = mTable->clone(); 
        return *this;
    }
//    Mesh::RegularMesh2D<T> mRegular2DMesh;
//    Mesh::RegularMesh3D<T> mRegular3DMesh;
    std::unique_ptr<Mesh::IMesh<T>> mTable = NULL;
    std::string mTableName;
    bool mIsNodalField = false;
};

/// Constructor
template<class T>
TravelTimeTable<T>::TravelTimeTable() :
    pImpl(std::make_unique<TravelTimeTableImpl> ())
{
}

/// Copy constructor
template<class T>
TravelTimeTable<T>::TravelTimeTable(const TravelTimeTable &table)
{
    *this = table;
}

/// Move constructor
template<class T>
TravelTimeTable<T>::TravelTimeTable(TravelTimeTable &&table) noexcept
{
    *this = std::move(table);
}

/// Copy assignment operator
template<class T>
TravelTimeTable<T>& TravelTimeTable<T>::operator=(const TravelTimeTable &table)
{
    if (&table == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<TravelTimeTableImpl> ();
    pImpl->mTable = table.pImpl->mTable->clone();
    return *this;
}

/// Move assignment operator
template<class T>
TravelTimeTable<T>& TravelTimeTable<T>::operator=(TravelTimeTable &&table) noexcept
{
    if (&table == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(table.pImpl);
    return *this;
}

/// Destructor
template<class T>
TravelTimeTable<T>::~TravelTimeTable() = default;

/// Releases memory on the module
template<class T>
void TravelTimeTable<T>::clear() noexcept
{
    if (pImpl->mTable)
    {
        pImpl->mTable.reset();
    }
    pImpl->mTableName.clear();
}

/// Sets the travel time field
template<class T>
void TravelTimeTable<T>::setTravelTimeTable(
    const std::string &fieldName, const Mesh::IMesh<T> &mesh)
{
    clear();
    // Determine if the field makes sense
    if (mesh.haveNodalScalarField(fieldName))
    {
        auto ptr = mesh.getNodalScalarFieldPointer(fieldName);
        auto npts = mesh.getNumberOfGridPoints();
        if (npts < 1)
        {
            throw std::invalid_argument("No points in field: " + fieldName);
        }
        auto minIndex = std::distance(ptr, std::min_element(ptr, ptr+npts));
        if (ptr[minIndex] < 0)
        {
            throw std::invalid_argument("Negative travel time detected!\n");
        }
        pImpl->mTable = mesh.clone();
    }
    else if (mesh.haveCellularScalarField(fieldName))
    {
        auto ptr = mesh.getCellularScalarFieldPointer(fieldName);
        auto npts = mesh.getNumberOfCells();
        if (npts < 1)
        {
            throw std::invalid_argument("No points in field: " + fieldName);
        }
        auto minIndex = std::distance(ptr, std::min_element(ptr, ptr+npts));
        if (ptr[minIndex] < 0)
        {
            throw std::invalid_argument("Negative travel time detected!\n");
        }
        pImpl->mTable = mesh.clone();
    }
    else
    {
        throw std::invalid_argument("Field name " + fieldName
                                  + " doesn't exist in mesh\n");
    }
    pImpl->mTableName = fieldName;
}

template<class T>
const T* TravelTimeTable<T>::getTravelTimeTablePointer() const
{
    if (!haveTravelTimeTable())
    {
        throw std::invalid_argument("Travel time table never set\n");
    }
    if (isNodalField())
    {
        return pImpl->mTable->getNodalScalarFieldPointer(pImpl->mTableName);
    }
    else
    {
        return pImpl->mTable->getCellularScalarFieldPointer(pImpl->mTableName);
    }
}

template<class T>
bool TravelTimeTable<T>::haveTravelTimeTable() const noexcept
{
    if (pImpl->mTable){return true;}
    return false;
}

/// Gets the mesh type 
template<class T>
XCLoc::Mesh::MeshType TravelTimeTable<T>::getMeshType() const
{
    if (!haveTravelTimeTable())
    {
        throw std::runtime_error("Travel time table not yet set\n");
    }
    return pImpl->mTable->getMeshType(); 
}

/// Determines if the underlying travel time field is nodal or cell-based
template<class T> bool TravelTimeTable<T>::isNodalField() const
{
    if (!haveTravelTimeTable())
    {
        throw std::runtime_error("Travel time table not yet set\n");
    }
    return pImpl->mTable->haveNodalScalarField(pImpl->mTableName);
}

template<class T> int TravelTimeTable<T>::getNumberOfPoints() const
{
    auto isNodal = isNodalField(); // Throws 
    if (isNodal)
    {
        return pImpl->mTable->getNumberOfGridPoints();
    }
    else
    {
        return pImpl->mTable->getNumberOfCells();
    }
}

/// Template instantiation
template class XCLoc::TravelTimeTable<float>;
template class XCLoc::TravelTimeTable<double>;
