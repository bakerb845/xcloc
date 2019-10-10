#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include "xcloc/travelTimeTables.hpp"
#include "xcloc/travelTimeTable.hpp"
#include "xcloc/travelTimeTableName.hpp"
#include "xcloc/waveformIdentifier.hpp"

using namespace XCLoc;

namespace
{

int getIndex(const TravelTimeTableName &tableName,
             const std::vector<TravelTimeTableName> mTableNames)
{
    int index = -1;
    auto it = std::find(mTableNames.begin(), mTableNames.end(),
                        tableName);
    if (it != mTableNames.end())
    {
        index = std::distance(mTableNames.begin(), it);
    }
    return index;
}
 
}

template<class T>
class TravelTimeTables<T>::TravelTimeTablesImpl
{
public:
    std::map<WaveformIdentifier, TravelTimeTable<T>> mwork;//, compareWaveformIdentifiers> mwork;
    std::vector<TravelTimeTableName> mTableNames;
    std::vector<TravelTimeTable<T>> mTables;
};

/// Constructor
template<class T>
TravelTimeTables<T>::TravelTimeTables() :
    pImpl(std::make_unique<TravelTimeTablesImpl> ())
{
}

/// Copy constructor
template<class T>
TravelTimeTables<T>::TravelTimeTables(const TravelTimeTables &tables)
{
    *this = tables;
}

/// Move constructor
template<class T>
TravelTimeTables<T>::TravelTimeTables(TravelTimeTables &&tables) noexcept
{
    *this = std::move(tables);
}

/// Move assignmment operator
template<class T>
TravelTimeTables<T>& TravelTimeTables<T>::operator=(
    TravelTimeTables &&tables) noexcept
{
    if (&tables == this){return *this;}
    pImpl = std::move(tables.pImpl);
    return *this;
}

/// Copy assignment operator
template<class T>
TravelTimeTables<T>& TravelTimeTables<T>::operator=(
    const TravelTimeTables &tables)
{
    if (&tables == this){return *this;}
    pImpl = std::make_unique<TravelTimeTablesImpl> (*tables.pImpl);
    return *this;
}

/// Destructor
template<class T>
TravelTimeTables<T>::~TravelTimeTables() = default;

/// Clear the class
template<class T>
void TravelTimeTables<T>::clear() noexcept
{
    pImpl->mTableNames.clear();
    pImpl->mTables.clear();
}

/// Add or append a travel time table
template<class T>
void TravelTimeTables<T>::addTable(const TravelTimeTableName &tableName,
                                   const TravelTimeTable<T> &table)
{
    if (!tableName.isValid())
    {
        throw std::invalid_argument("tableName is not valid\n");
    }
    // Check meshes are consistent
    if (getNumberOfTables() > 0)
    {
        if (getMeshType() != table.getMeshType())
        {
            throw std::invalid_argument("Inconsistent mesh types\n");
        }
        if (isNodalField() != table.isNodalField())
        {
            throw std::invalid_argument("Cant mix node and cell-based field\n");
        }
        if (getNumberOfPoints() != table.getNumberOfPoints())
        {
            throw std::invalid_argument("Inconsistent number of points\n");
        }
    }
    // Look for the table in the list
    auto index = getIndex(tableName, pImpl->mTableNames);
    if (index < 0)
    {
        pImpl->mTableNames.push_back(tableName);
        pImpl->mTables.push_back(table);
    }
    else
    {
        pImpl->mTables[index] = table;
    }
}

/// Get a travel time table
template<class T> TravelTimeTable<T>
TravelTimeTables<T>::getTable(const TravelTimeTableName &tableName) const
{
    auto index = getIndex(tableName, pImpl->mTableNames);
    if (index < 0)
    {
        throw std::invalid_argument("Could not find table\n");
    }
    return pImpl->mTables[index];
}

/// Gets the travel time table names
template<class T> std::vector<TravelTimeTableName>
TravelTimeTables<T>::getTableNames() const
{
    return pImpl->mTableNames;
}

/// Get a pointer to the travel time table
template<class T>
const T *TravelTimeTables<T>::getTravelTimeTablePointer(
    const TravelTimeTableName &tableName) const
{
    auto index = getIndex(tableName, pImpl->mTableNames);
    if (index < 0)
    {
        throw std::invalid_argument("Could not find table\n");
    }
    return pImpl->mTables[index].getTravelTimeTablePointer();
}

/// Check if a table exists
template<class T>
bool TravelTimeTables<T>::haveTable(const TravelTimeTableName &tableName) const
{
    // Look for the table in the list
    auto it = std::find(pImpl->mTableNames.begin(), pImpl->mTableNames.end(),
                        tableName);
    if (it == pImpl->mTableNames.end()){return false;}
    return true;
}

/// Get the travel time table type
template<class T>
XCLoc::Mesh::MeshType TravelTimeTables<T>::getMeshType() const
{
    if (getNumberOfTables() < 1)
    {
        throw std::runtime_error("No tables set\n");
    }
    return pImpl->mTables[0].getMeshType();
}

/// Is the underlying travel time field a nodal or cell-based?
template<class T>
bool TravelTimeTables<T>::isNodalField() const
{
    if (getNumberOfTables() < 1)
    {
        throw std::runtime_error("No tables set\n");
    }
    return pImpl->mTables[0].isNodalField();
}

/// Gets the number of tables
template<class T>
int TravelTimeTables<T>::getNumberOfTables() const noexcept
{
    return static_cast<int> (pImpl->mTables.size());
}

/// Gets the number of points in the field
template<class T>
int TravelTimeTables<T>::getNumberOfPoints() const
{
    if (getNumberOfTables() < 1)
    {
        throw std::runtime_error("No tables set\n");
    }
    return pImpl->mTables[0].getNumberOfPoints(); 
}

/// Template instantiation
template class XCLoc::TravelTimeTables<float>;
template class XCLoc::TravelTimeTables<double>;
