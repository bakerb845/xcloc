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
    bool mHaveNetwork = false;
    bool mHaveStation = false;
    bool mHavePhase = false;
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
    const std::string &phase) :
    pImpl(std::make_unique<TravelTimeTableNameImpl> ())
{
    setNetwork(network);
    setStation(station);
    setPhase(phase);
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

/// Destructor
TravelTimeTableName::~TravelTimeTableName() = default;

/// Resets the class
void TravelTimeTableName::clear() noexcept
{
    pImpl->mNetwork.clear();
    pImpl->mStation.clear();
    pImpl->mPhase.clear();
    pImpl->mHaveNetwork = false;
    pImpl->mHaveStation = false;
    pImpl->mHavePhase = false;
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

/// Check if network, station, phase were set
bool TravelTimeTableName::isValid() const noexcept
{
    if (!haveNetwork()){return false;}
    if (!haveStation()){return false;}
    if (!havePhase()){return false;}
    return true;
}
