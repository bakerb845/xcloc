#include <cstdio>
#include <cstdlib>
#include <string>
#include "xcloc/travelTimeTableName.hpp"

using namespace XCLoc;

class TravelTimeTableName::TravelTimeTableNameImpl
{
public:
    std::string mNetwork;
    std::string mStation;
    std::string mPhase;
    std::string mPolarization;
    bool mHaveNetwork = false;
    bool mHaveStation = false;
    bool mHavePhase = false;
    bool mHavePolarization = false;
};

/// Constructor
TravelTimeTableName::TravelTimeTableName() :
    pImpl(std::make_unique<TravelTimeTableNameImpl> ())
{
}

/// Constructor
TravelTimeTableName::TravelTimeTableName(
    const std::string &network,
    const std::string &station,
    const std::string &phase,
    const std::string &polarization) :
    pImpl(std::make_unique<TravelTimeTableNameImpl> ())
{
    setNetwork(network);
    setStation(station);
    setPhase(phase);
    setPolarization(polarization);
}

/// Copy constructor
TravelTimeTableName::TravelTimeTableName(
    const TravelTimeTableName &table)
{
    *this = table;
}

/// Move constructor
TravelTimeTableName::TravelTimeTableName(
    TravelTimeTableName &&table) noexcept
{
    *this = std::move(table);
}

/// Copy assignment
TravelTimeTableName& 
TravelTimeTableName::operator=(const TravelTimeTableName &tableName)
{
    if (&tableName == this){return *this;}
    pImpl = std::make_unique<TravelTimeTableNameImpl> (*tableName.pImpl);
    return *this;
}

/// Move assignment
TravelTimeTableName&
TravelTimeTableName::operator=(TravelTimeTableName &tableName) noexcept
{
    if (&tableName == this){return *this;}
    pImpl = std::move(tableName.pImpl);
    return *this;
}

/*
/// Equality operator
bool TravelTimeTableName::operator==(
    const TravelTimeTableName &tableName) const noexcept
{
    // Network
    if (haveNetwork() && tableName.haveNetwork())
    {
        if (getNetwork() != tableName.getNetwork()){return false;}
    }
    else
    {
         return false;
    }
    // Station
    if (haveStation() && tableName.haveStation())
    {
        if (getStation() != tableName.getStation()){return false;}
    }
    else
    {
         return false;
    }
    // Phase
    if (havePhase() && tableName.havePhase())
    {
        if (getPhase() != tableName.getPhase()){return false;}
    }
    else
    {
        return false;
    }
    return true;
}

/// Inequality operator
bool TravelTimeTableName::operator!=(
    const TravelTimeTableName &tableName) const noexcept
{
    return !(*this == tableName);
}
*/

/// Destructor
TravelTimeTableName::~TravelTimeTableName() = default;

/// Resets the class
void TravelTimeTableName::clear() noexcept
{
    pImpl->mNetwork.clear();
    pImpl->mStation.clear();
    pImpl->mPhase.clear();
    pImpl->mPolarization.clear();
    pImpl->mHaveNetwork = false;
    pImpl->mHaveStation = false;
    pImpl->mHavePhase = false;
    pImpl->mHavePolarization = false;
}

/// Network
void TravelTimeTableName::setNetwork(const std::string &network) noexcept
{
    pImpl->mNetwork = network;
    pImpl->mHaveNetwork = true;
}

std::string TravelTimeTableName::getNetwork() const
{
    if (!haveNetwork())
    {
        throw std::runtime_error("Network name not set\n");
    }
    return pImpl->mNetwork;
}

bool TravelTimeTableName::haveNetwork() const noexcept
{
    return pImpl->mHaveNetwork;
}

/// Station
void TravelTimeTableName::setStation(const std::string &station) noexcept
{
    pImpl->mStation = station;
    pImpl->mHaveStation = true;
}

std::string TravelTimeTableName::getStation() const
{
    if (!haveStation())
    {
        throw std::runtime_error("Station name not set\n");
    }
    return pImpl->mStation;
}

bool TravelTimeTableName::haveStation() const noexcept
{
    return pImpl->mHaveStation;
}

/// Phase name
void TravelTimeTableName::setPhase(const std::string &phase) noexcept
{
    pImpl->mPhase = phase;
    pImpl->mHavePhase = true;
}

std::string TravelTimeTableName::getPhase() const
{
    if (!havePhase())
    {
        throw std::runtime_error("Phase name not set\n");
    }
    return pImpl->mPhase;
}

bool TravelTimeTableName::havePhase() const noexcept
{
    return pImpl->mHavePhase;
}

/// Polarization
void TravelTimeTableName::setPolarization(
    const std::string &polarization) noexcept
{
    pImpl->mPolarization = polarization;
    for (auto i=0; i<pImpl->mPolarization.size(); ++i)
    {
        pImpl->mPolarization[i] = std::toupper(pImpl->mPolarization[i]);
    }
    pImpl->mHavePolarization = true;
}

std::string TravelTimeTableName::getPolarization() const
{
    if (!havePolarization())
    {
        throw std::runtime_error("Polarization name not set\n");
    }
    return pImpl->mPolarization;
}

bool TravelTimeTableName::havePolarization() const noexcept
{
    return pImpl->mHavePolarization;
}

/// Check if network, station, phase were set
bool TravelTimeTableName::isValid() const noexcept
{
    if (!haveNetwork()){return false;}
    if (!haveStation()){return false;}
    if (!havePhase()){return false;}
    if (!havePolarization()){return false;}
    return true;
}

/// Equality operator
bool XCLoc::operator==(const TravelTimeTableName &lhs,
                       const TravelTimeTableName &rhs) noexcept
{
    // Network
    if (lhs.haveNetwork() && rhs.haveNetwork())
    {
        if (lhs.getNetwork() != rhs.getNetwork()){return false;}
    }
    else
    {
         return false;
    }
    // Station
    if (lhs.haveStation() && rhs.haveStation())
    {
        if (lhs.getStation() != rhs.getStation()){return false;}
    }
    else
    {
         return false;
    }
    // Phase
    if (lhs.havePhase() && rhs.havePhase())
    {
        if (lhs.getPhase() != rhs.getPhase()){return false;}
    }
    else
    {
        return false;
    }
    // Polarization
    if (lhs.havePolarization() && rhs.havePolarization())
    {
        if (lhs.getPolarization() != rhs.getPolarization()){return false;}
    }
    else
    {
        return false;
    }
    return true;
}

/// Inequality operator
bool XCLoc::operator!=(const TravelTimeTableName &lhs,
                       const TravelTimeTableName &rhs) noexcept
{
   return !(lhs == rhs);
}
