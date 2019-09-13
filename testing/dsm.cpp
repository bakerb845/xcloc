#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <complex>
#include <random>
#include <rtseis/utilities/transforms/utilities.hpp>
#include "xcloc/correlogramDSM.hpp"
#include "xcloc/correlogramDSMParameters.hpp"
#include "xcloc/correlograms.hpp"
#include "xcloc/mesh/regularMesh2D.hpp"
#include "xcloc/travelTimeTable.hpp"
#include "acousticGreens2D.hpp"
#include "rickerWavelet.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace XCLoc;

AcousticGreens2D makeSimulation();
void createTravelTimeTable(const int nx, const double x0, const double x1,
                           const int nz, const double z0, const double z1,
                           const double vp, const double vs,
                           const double xRec, const double zRec);

TEST(CorrelogramDiffractionStackMigration, parameters)
{
    CorrelogramDSMParameters parameters;
    EXPECT_THROW(parameters.setCacheBlockSize(63),
                 std::invalid_argument);
    EXPECT_THROW(parameters.setCacheBlockSize(65),
                 std::invalid_argument);
    for (int i=6; i<32; ++i)
    {
        auto result = static_cast<size_t>(std::exp2(i));
        EXPECT_NO_THROW(parameters.setCacheBlockSize(result));
        EXPECT_EQ(parameters.getCacheBlockSize(), result);
    }
    parameters.setCacheBlockSize(128);
    // Test copy c'tor
    CorrelogramDSMParameters parmsCopy(parameters);
    EXPECT_EQ(parmsCopy.getCacheBlockSize(), 128);
}

TEST(CorrelogramDiffractionStackMigration, dsm)
{
    // Set DSM parameters
    CorrelogramDSMParameters parameters;
    parameters.setCacheBlockSize(2048); // Optional
    // Set correlograms parameters
    Correlograms<double> correlograms;
    // Create the diffraction stack migration
    auto sim = makeSimulation(); 
}

/// Make a simulation consisting of 10 stations randomly distributed on a 2D
/// grid.
AcousticGreens2D makeSimulation()
{
    int nrec = 10;
    double vp = 4500;    // Vp velocity
    double Qp = 600;     // Quality factor
    double rho = 2700;   // Density
    // Model width
    int nx = 101;
    int nz = 501;
    double x0 = 0;
    double x1 = 5000;
    double z0 = 0;
    double z1 = 10000;
    // Receiver locations

    // Source time function
    double ts = 10;      // 10 seconds source time function
    double fcent = 8.5;  // Center frequency of ricker wavelet 
    double samplingRate = 100; // Sampling rate 
    RickerWavelet wavelet;
    int npts = static_cast<int> (ts*samplingRate + 0.5);
    wavelet.setSamplingRate(samplingRate);
    wavelet.setCenterFrequency(fcent);
    wavelet.setNumberOfSamples(npts);
    wavelet.setShiftWaveletToTraceStart(true, 1.e-5);
    EXPECT_EQ(npts,  wavelet.getNumberOfSamples());
    // Green's functions
    AcousticGreens2D greens;
    greens.setSourceTimeFunction(wavelet);
    greens.setVelocity(vp);
    greens.setDensity(rho);
    greens.setQualityFactor(Qp);

auto stf = wavelet.getWavelet();
FILE *fout = fopen("stf.txt", "w");
int k =0;
for (auto s : stf)
{
 fprintf(fout, "%lf, %lf\n", k/samplingRate, s);
 k= k+1;
}
fclose(fout);
EXPECT_EQ(stf.size(), wavelet.getNumberOfSamples());
auto wstf = wavelet.getWaveletFourierTransform();
auto wfreqs = RTSeis::Utilities::Transforms::DFTUtilities::realToComplexDFTFrequencies(stf.size(), 1./samplingRate);
fout = fopen("wstf.txt", "w");
k=0;
EXPECT_EQ(wstf.size(), wfreqs.size());
for (auto s : wstf)
{
 fprintf(fout, "%lf, %lf\n", wfreqs[k], std::abs(s));
  k = k + 1;
}
fclose(fout);

    return greens;
}

void createTravelTimeTable(const int nx, const double x0, const double x1,
                           const int nz, const double z0, const double z1,
                           const double vp, const double vs,
                           const double xRec, const double zRec)
{
    auto xwidth = x1 - x0;
    auto zwidth = z1 - z0;
    if (nx < 2){throw std::invalid_argument("nx must be at least 2\n");}
    if (nz < 2){throw std::invalid_argument("nz must be at least 2\n");}
    /// Create a travel time table
    auto dx = xwidth/static_cast<double> (nx - 1);
    auto dz = zwidth/static_cast<double> (nz - 1); 
    XCLoc::Mesh::RegularMesh2D<double> mesh;
    mesh.setNumberOfGridPointsInX(nx);
    mesh.setNumberOfGridPointsInZ(nz);
    mesh.setGridSpacingInX(dx);
    mesh.setGridSpacingInZ(dz);
    int ngrd = mesh.getNumberOfGridPoints();
    std::vector<double> tp(ngrd, 0);
    std::vector<double> ts(ngrd, 0);
    for (int ix=0; ix<nx; ++ix)
    {
        for (int iz=0; iz<nz; ++iz)
        {
            auto xGrd = dx*ix;
            auto zGrd = dz*iz;
            auto diffx = xRec - xGrd;
            auto diffz = zRec - zGrd;
            auto dist = std::sqrt(diffx*diffx + diffz*diffz);
            tp[ix*nz + iz] = dist/vp;
            ts[ix*nz + iz] = dist/vs;
        }
    }
    mesh.setNodalScalarField("pTravelTimes", ngrd, tp.data(),
                             XCLoc::Mesh::RegularMesh2DOrderingType::NX_NZ);
    mesh.setNodalScalarField("sTravelTimes", ngrd, ts.data(),
                             XCLoc::Mesh::RegularMesh2DOrderingType::NX_NZ);
    // Create the travel time field
    TravelTimeTable<double> pTable, sTable;
    pTable.setTravelTimeTable("pTravelTimes", mesh);
    sTable.setTravelTimeTable("sTravelTimes", mesh);
}

void createRandomReceiverLocations2D(const int nrec,
                                     const double x0, const double x1,
                                     const double y0, const double y1,
                                     std::vector<double> *xrec,
                                     std::vector<double> *yrec)
{
    xrec->resize(0);
    yrec->resize(0);
    if (nrec < 0)
    {
        throw std::invalid_argument("nrec cannot be negative\n");
    }
    std::random_device device;
    std::mt19937 rng(device());
    std::uniform_real_distribution<double> uniform_dist1(x0, x1);
    std::uniform_real_distribution<double> uniform_dist2(y0, y1); 
    std::vector<double> x(nrec), y(nrec);
    for (int i=0; i<nrec; ++i)
    {
        x[i] = uniform_dist1(rng);
        y[i] = uniform_dist2(rng);
    }
    *xrec = x;
    *yrec = y;
}

}
