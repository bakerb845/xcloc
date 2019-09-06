#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include "xcloc/diffractionStackMigrationEngine.hpp"
#include "xcloc/correlograms.hpp"
#include "xcloc/mesh/regularMesh3D.hpp"

using namespace XCLoc;

template<class T>
class DiffractionStackMigrationEngine<T>::DSMImpl
{
public:
    std::vector<XCLoc::Mesh::IMesh<T>> mTravelTimeTables;
    /// Contains the correlograms
    std::shared_ptr<const Correlograms<T>> mCorrelograms;
    /// Sampling period of seismograms
    double mSamplingPeriod = 1;
    /// Chunk size to emphasize locality in DSM for read/write operations
    int mChunkSize = 2048;
    /// Flag indicating this is nodal-based imaging or cell-based imaging.
    bool mNodalBasedImaging = true;
    /// Flag indicating that the correlograms are set
    bool mHaveCorrelograms = false;
};

/// Constructor
template<class T>
DiffractionStackMigrationEngine<T>::DiffractionStackMigrationEngine() :
    pImpl(std::make_unique<DSMImpl> ())
{
XCLoc::Mesh::RegularMesh3D<T> mGrid;
pImpl->mTravelTimeTables.push_back(mGrid);
}

/// Destructor
template<class T>
DiffractionStackMigrationEngine<T>::~DiffractionStackMigrationEngine()
    = default;

/// Sets the correlation engine
/*
template<class T>
void DiffractionStackMigrationEngine<T>::setCorrelograms(
    std::shared_ptr<const Correlograms<T>> correlograms)
{
    if (pImpl->mCorrelograms){pImpl->mCorrelograms.release();}
    pImpl->mCorrelograms
        = std::make_shared<const Correlograms<T>> (correlograms);
    pImpl->mHaveCorrelograms = true;
}
*/

/// Determines if the correlation engine was set
template<class T>
bool DiffractionStackMigrationEngine<T>::haveCorrelograms() const noexcept
{
    return pImpl->mHaveCorrelograms;
}

/// Computes the diffraction stack migration
template<class T>
void DiffractionStackMigrationEngine<T>::compute()
{
    // Get the number of grid points
    int nGridPoints = pImpl->mTravelTimeTables[0].getNumberOfGridPoints();
    if (!pImpl->mNodalBasedImaging)
    {
        nGridPoints = pImpl->mTravelTimeTables[0].getNumberOfCells();
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
        T __attribute__((aligned(64))) *imagePtr = nullptr;
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
            // Loop on the local grid 
            #pragma omp simd
            for (int i=0; i<nLocalGridPoints; ++i)
            {
                int correlogramIndex = lxc2 + tt1[i] - tt2[i];
                imagePtr[i] = imagePtr[i] + xc[correlogramIndex];
            }
        } // Loop on correlograms
    } // Loop on chunks
}
