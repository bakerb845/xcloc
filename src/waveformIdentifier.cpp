#include <cstdio>
#include <cstdlib>
#include <string>
#include <functional>
#include "xcloc/waveformIdentifier.hpp"

using namespace XCLoc;
namespace
{
std::string makeHashString(const std::string &network,
                           const std::string &station,
                           const std::string &polarization)
{
    auto result = network + "." + station + "." + polarization;
    return result;
}
}

class WaveformIdentifier::WaveformIdentifierImpl
{
public:
    WaveformIdentifierImpl() :
        mNetwork(""),
        mStation(""),
        mPolarization("") 
    {
        mHash = std::hash<std::string>{}(
           makeHashString(mNetwork, mStation, mPolarization));
    }
    std::string mNetwork = "";
    std::string mStation = "";
    std::string mPolarization = "";
    int64_t mHash = 0;
};

/// Constructor
WaveformIdentifier::WaveformIdentifier() :
    pImpl(std::make_unique<WaveformIdentifierImpl> ())
{
}

/// Constructor
WaveformIdentifier::WaveformIdentifier(
    const std::string &network,
    const std::string &station,
    const std::string &polarization) :
    pImpl(std::make_unique<WaveformIdentifierImpl> ())
{
    setNetwork(network);
    setStation(station);
    setPolarization(polarization);
}
/// Copy constructor
WaveformIdentifier::WaveformIdentifier(const WaveformIdentifier &waveid)
{
    *this = waveid;
}

/// Move constructor
WaveformIdentifier::WaveformIdentifier(WaveformIdentifier &&waveid) noexcept
{
    *this = std::move(waveid);
}

/// Copy assignment 
WaveformIdentifier&
WaveformIdentifier::operator=(const WaveformIdentifier &waveid)
{
    if (&waveid == this){return *this;}
    pImpl = std::make_unique<WaveformIdentifierImpl> (*waveid.pImpl);
    return *this;
}

/// Move assignment
WaveformIdentifier&
WaveformIdentifier::operator=(WaveformIdentifier &&waveid) noexcept
{
    if (&waveid == this){return *this;}
    pImpl = std::move(waveid.pImpl);
    return *this;
}

/*
/// Comparitors
bool WaveformIdentifier::operator>(const WaveformIdentifier &waveid)
{
    return (pImpl->mHash > waveid.pImpl->mHash);
}

bool WaveformIdentifier::operator<(const WaveformIdentifier &waveid)
{
    return (pImpl->mHash < waveid.pImpl->mHash);
}

bool WaveformIdentifier::operator==(const WaveformIdentifier &waveid)
{
    return (pImpl->mHash == waveid.pImpl->mHash);
}

bool WaveformIdentifier::operator!=(const WaveformIdentifier &waveid)
{
    return !(*this == waveid);
}
*/

/// Destructor
WaveformIdentifier::~WaveformIdentifier() = default;

/// Clears module
void WaveformIdentifier::clear() noexcept
{
    pImpl->mNetwork = "";
    pImpl->mStation = "";
    pImpl->mPolarization = "";
    pImpl->mHash = std::hash<std::string>{}(
        makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPolarization));
}

/// Set/get network
void WaveformIdentifier::setNetwork(const std::string &network)
{
    pImpl->mNetwork = network;
    pImpl->mHash = std::hash<std::string>{}(
        makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPolarization));
}

std::string WaveformIdentifier::getNetwork() const noexcept
{
    return pImpl->mNetwork;
}

/// Set/get station
void WaveformIdentifier::setStation(const std::string &station)
{
    pImpl->mStation = station;
    pImpl->mHash = std::hash<std::string>{}(
        makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPolarization));
}

std::string WaveformIdentifier::getStation() const noexcept
{
    return pImpl->mStation;
}

/// Set/get polarization
void WaveformIdentifier::setPolarization(const std::string &polarization)
{
    pImpl->mPolarization = polarization;
    for (auto i=0; i<pImpl->mPolarization.size(); ++i)
    {
        pImpl->mPolarization[i] = std::toupper(pImpl->mPolarization[i]);
    }
    pImpl->mHash = std::hash<std::string>{}(
        makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPolarization));
}

std::string WaveformIdentifier::getPolarization() const noexcept
{
    return pImpl->mPolarization;
}


/// Gets the waveform identifier
int64_t WaveformIdentifier::getWaveformIdentifier() const
{
    return pImpl->mHash;
}

/// Comparitors
bool XCLoc::operator > (
    const WaveformIdentifier &lhs, const WaveformIdentifier &rhs)
{
    return (lhs.getWaveformIdentifier() > rhs.getWaveformIdentifier());
}

bool XCLoc::operator<(
    const WaveformIdentifier &lhs, const WaveformIdentifier &rhs)
{
    return (lhs.getWaveformIdentifier() < rhs.getWaveformIdentifier());
}

bool XCLoc::operator==(
    const WaveformIdentifier &lhs, const WaveformIdentifier &rhs)
{
    return (lhs.getWaveformIdentifier() == rhs.getWaveformIdentifier());
}

bool XCLoc::operator!=(
    const WaveformIdentifier &lhs, const WaveformIdentifier &rhs)
{
    return !(lhs == rhs);
}
