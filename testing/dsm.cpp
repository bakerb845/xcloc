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
std::pair<TravelTimeTable<double>, TravelTimeTable<double>>
createTravelTimeTable(const int nx, const double x0, const double x1,
                      const int ny, const double y0, const double y1,
                      const double vp, const double vs,
                      const double xRec, const double yRec);
std::vector<std::pair<double,double>>
createRandomLocations2D(const int nrec,
                        const double x0, const double x1,
                        const double y0, const double y1);

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
    int ny = 501;
    double x0 =-2500;
    double x1 = 2500;
    double y0 =-5000;
    double y1 = 5000;
    auto dx = (x1 - x0)/static_cast<double> (nx - 1);
    auto dy = (y1 - y0)/static_cast<double> (ny - 1);
    // Source location
    auto sourcePosition
        = std::make_pair<double, double> (x0 + 41*dx, y0 + 265*dy);
    // Receiver locations
    auto receiverPositions = createRandomLocations2D(nrec, x0, x1, y0, y1);
    // Source time function
    double ts = 6;       // 6 seconds source time function
    double fcent = 8.5;  // Center frequency of ricker wavelet 
    double samplingRate = 100; // Sampling rate 
    RickerWavelet wavelet;
    int npts = static_cast<int> (ts*samplingRate + 0.5);
    wavelet.setSamplingRate(samplingRate);
    wavelet.setCenterFrequency(fcent);
    wavelet.setNumberOfSamples(npts);
    wavelet.setShiftWaveletToTraceStart(true, 1.e-5);
    wavelet.setNormalizeByEnergy(false);
    EXPECT_EQ(npts,  wavelet.getNumberOfSamples());
    // Green's functions
    AcousticGreens2D greens;
    greens.setSourceTimeFunction(wavelet);
    greens.setVelocity(vp);
    greens.setDensity(rho);
    greens.setQualityFactor(Qp);
    greens.setSourcePosition(sourcePosition);
    greens.setReceiverPositions(receiverPositions); 
    greens.compute();
    // Travel time tables
    for (int i=0; i<static_cast<int> (receiverPositions.size()); ++i)
    {
        auto xRec = receiverPositions[i].first;
        auto yRec = receiverPositions[i].second;
        auto pstable = createTravelTimeTable(nx, x0, x1,
                                             ny, y0, y1,
                                             vp, vp/1.73,
                                             xRec, yRec);
    }

FILE *fout = fopen("greens.txt", "w");
for (int irec=0; irec<nrec; ++irec)
{
 auto g = greens.getGreensFunction(irec);
 for (int i=0; i<g.size(); ++i)
 {
  fprintf(fout, "%lf, %e\n", i/samplingRate, irec*5.e-13 + g[i]);
 }
 fprintf(fout, "\n");
}
fclose(fout);
/*
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
*/

    return greens;
}

std::pair<TravelTimeTable<double>, TravelTimeTable<double>>
createTravelTimeTable(const int nx, const double x0, const double x1,
                      const int ny, const double y0, const double y1,
                      const double vp, const double vs,
                      const double xRec, const double yRec)
{
    auto xwidth = x1 - x0;
    auto ywidth = y1 - y0;
    if (nx < 2){throw std::invalid_argument("nx must be at least 2\n");}
    if (ny < 2){throw std::invalid_argument("ny must be at least 2\n");}
    /// Create a travel time table
    auto dx = xwidth/static_cast<double> (nx - 1);
    auto dy = ywidth/static_cast<double> (ny - 1); 
    XCLoc::Mesh::RegularMesh2D<double> mesh;
    mesh.setNumberOfGridPointsInX(nx);
    mesh.setNumberOfGridPointsInZ(ny);
    mesh.setGridSpacingInX(dx);
    mesh.setGridSpacingInZ(dy);
    int ngrd = mesh.getNumberOfGridPoints();
    std::vector<double> tp(ngrd, 0);
    std::vector<double> ts(ngrd, 0);
    for (int ix=0; ix<nx; ++ix)
    {
        for (int iy=0; iy<ny; ++iy)
        {
            auto xGrd = dx*ix;
            auto yGrd = dy*iy;
            auto diffx = xRec - xGrd;
            auto diffy = yRec - yGrd;
            auto dist = std::sqrt(diffx*diffx + diffy*diffy);
            tp[ix*ny + iy] = dist/vp;
            ts[ix*ny + iy] = dist/vs;
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
    return std::make_pair(pTable, sTable);
}

/// Generates uniformly distributed random positions on
std::vector<std::pair<double,double>>
createRandomLocations2D(const int nrec,
                        const double x0, const double x1,
                        const double y0, const double y1)
{
    if (nrec < 0)
    {
        throw std::invalid_argument("nrec cannot be negative\n");
    }
    //std::random_device device;
    std::uint_least32_t seed = 3348;
    std::mt19937 rng(seed); //device());
    std::uniform_real_distribution<double> uniform_dist1(x0, x1);
    std::uniform_real_distribution<double> uniform_dist2(y0, y1); 
    std::vector<std::pair<double, double>> positions(nrec);
    for (int i=0; i<nrec; ++i)
    {
        positions[i] = std::make_pair(uniform_dist1(rng), uniform_dist2(rng));
    }
    return positions;
}

}
