#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include "xcloc/correlogramDSM.hpp"
#include "xcloc/correlograms.hpp"
#include "xcloc/mesh/regularMesh2D.hpp"
#include "xcloc/mesh/regularMesh3D.hpp"

using namespace XCLoc;

#define IMAGE_NAME "dsmImage"

template<class T>
class CorrelogramDiffractionStackMigration<T>::DSMImpl
{
public:
    /// Holds the input travel time tables
    std::vector<XCLoc::Mesh::IMesh<T>> mTravelTimeTables;
    /// Holds the tables for performing the DSM
    std::vector<XCLoc::Mesh::IMesh<int>> mDSMTables;
    /// Holds the DSM image
    XCLoc::Mesh::IMesh<T> mImage;
    /// Contains the correlograms
    std::shared_ptr<const Correlograms<T>> mCorrelograms;
    /// Sampling period of seismograms
    double mSamplingPeriod = 0;
    /// Chunk size to emphasize locality in DSM for read/write operations
    int mChunkSize = 2048;
    /// Flag indicating this is nodal-based imaging or cell-based imaging.
    bool mNodalBasedImaging = true;
    /// Flag indicating that the correlograms are set
    bool mHaveCorrelogramPointer = false;
};

/// Constructor
template<class T>
CorrelogramDiffractionStackMigration<T>::CorrelogramDiffractionStackMigration() :
    pImpl(std::make_unique<DSMImpl> ())
{
}

/// Destructor
template<class T>
CorrelogramDiffractionStackMigration<T>::~CorrelogramDiffractionStackMigration()
    = default;

template<class T>
void CorrelogramDiffractionStackMigration<T>::clear() noexcept
{
    pImpl->mTravelTimeTables.clear();
    pImpl->mCorrelograms.reset();
    pImpl->mImage.clear();
    pImpl->mSamplingPeriod = 0;
    pImpl->mChunkSize = 2048;
    pImpl->mNodalBasedImaging = true;
    pImpl->mHaveCorrelogramPointer = false; 
}

/// Sets the correlograms
template<class T>
void CorrelogramDiffractionStackMigration<T>::setCorrelograms(
    std::shared_ptr<const Correlograms<T>> correlograms)
{
    if (pImpl->mCorrelograms){pImpl->mCorrelograms.release();}
    pImpl->mCorrelograms
        = std::make_shared<const Correlograms<T>> (correlograms);
    pImpl->mHaveCorrelogramPointer = true;
}

/// Determines if the correlation engine was set
template<class T>
bool CorrelogramDiffractionStackMigration<T>::haveCorrelogramPointer() const noexcept
{
    return pImpl->mHaveCorrelogramPointer;
}

/// Sets the sampling rate of the correlograms
template<class T>
void CorrelogramDiffractionStackMigration<T>::setCorrelogramSamplingRate(
    const double samplingRate)
{
    pImpl->mSamplingPeriod = 0;
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("Sampling rate = "
                                  + std::to_string(samplingRate)
                                  + " must be positive\n");
    }
    pImpl->mSamplingPeriod = 1/samplingRate;
}

/// Gets sampling rate of correlograms
template<class T>
double CorrelogramDiffractionStackMigration<T>::getCorrelogramSamplingRate() const
{
    if (!haveCorrelogramSamplingRate())
    {
        throw std::runtime_error("Sampling rate not set\n");
    }
    return 1.0/pImpl->mSamplingPeriod;
}

/// Checks if correlogram sampling rate was set
template<class T>
bool CorrelogramDiffractionStackMigration<T>::haveCorrelogramSamplingRate() const noexcept
{
    if (pImpl->mSamplingPeriod <= 0){return false;}
    return true;
}

/// Computes the diffraction stack migration
template<class T>
void CorrelogramDiffractionStackMigration<T>::compute()
{
    // Get the number of grid points
    int nGridPoints = pImpl->mDSMTables[0].getNumberOfGridPoints();
    if (!pImpl->mNodalBasedImaging)
    {
        nGridPoints = pImpl->mDSMTables[0].getNumberOfCells();
    }
    // Get a pointer to the image
    T __attribute__((aligned(64))) *imagePtr = nullptr;
    if (pImpl->mNodalImaging)
    {
        imagePtr = pImpl->mImage.getNodalPointer(IMAGE_NAME);
    }
    else
    {
        imagePtr = pImpl->mImage.getCellularPointer(IMAGE_NAME);
    }
    // Determine the correlogram length
    int nxc = pImpl->mCorrelograms.getNumberOfCorrelograms();
    int lxc2 = nxc/2;
    // Loop on chunks
    for (int igrd=0; igrd<nGridPoints; igrd=igrd+pImpl->mChunkSize)
    {
        // Figure out the number of local grid points to process 
        int nLocalGridPoints = std::min(pImpl->mChunkSize, nGridPoints - igrd);
        // Get the migration image pointer
        T __attribute__((aligned(64))) *localImagePtr =  imagePtr + igrd;
        // Loop on correlograms
        for (int ixc=0; ixc<nxc; ++ixc)
        {
            // Get the ixc'th correlogram
            auto __attribute((aligned(64)))
                xc = pImpl->mCorrelograms.getProcessedCorrelogram(ixc);
            auto xcPair = pImpl->mCorrelograms.getCorreloationPair(ixc);
            // Extract the appropriate travel time tables
            int it1 = xcPair.first;
            int it2 = xcPair.second;
            int *tt1 = nullptr;
            int *tt2 = nullptr;
            if (pImpl->mNodalImaging)
            {
                tt1 = pImpl->mDSMTables[it1].getNodalPointer("something\n");
                tt2 = pImpl->mDSMTables[it2].getNodalPointer("somethign\n");
            }
            else 
            {
                tt1 = pImpl->mDSMTables[it1].getCellularPointer("somethign\n");
                tt2 = pImpl->mDSMTables[it2].getCelluarlPointer("something\n");
            } 
            // Loop on the local grid 
            #pragma omp simd
            for (int i=0; i<nLocalGridPoints; ++i)
            {
                int correlogramIndex = lxc2 + tt1[i] - tt2[i];
                localImagePtr[i] = localImagePtr[i] + xc[correlogramIndex];
            }
        } // Loop on correlograms
    } // Loop on chunks
}
