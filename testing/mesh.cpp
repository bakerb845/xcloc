#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <ipps.h>
#include "xcloc/mesh/regularMesh3D.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace XCLoc::Mesh;
TEST(testMesh, regularMesh3D)
{
    RegularMesh3D mesh;
    int nx = 9;
    int ny = 10;
    int nz = 11;
    double dx = 3.0;
    double dy = 4.0;
    double dz = 5.0;
    double x0 =-1;
    double y0 =-2;
    double z0 =-3;
    // Grid size 
    EXPECT_FALSE(mesh.haveGridDimensions());

    EXPECT_NO_THROW(mesh.setNumberOfGridPointsInX(nx));
    EXPECT_FALSE(mesh.haveGridDimensions());

    EXPECT_NO_THROW(mesh.setNumberOfGridPointsInY(ny));
    EXPECT_FALSE(mesh.haveGridDimensions());

    EXPECT_NO_THROW(mesh.setNumberOfGridPointsInZ(nz));
    EXPECT_TRUE(mesh.haveGridDimensions());

    EXPECT_EQ(mesh.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(mesh.getNumberOfGridPointsInY(), ny);
    EXPECT_EQ(mesh.getNumberOfGridPointsInZ(), nz);

    EXPECT_EQ(mesh.getNumberOfCells(), (nx-1)*(ny-1)*(nz-1));
    EXPECT_EQ(mesh.getNumberOfGridPoints(), nx*ny*nz);

    // Grid spacing
    EXPECT_NO_THROW(mesh.setGridSpacingInX(dx));
    EXPECT_NO_THROW(mesh.setGridSpacingInY(dy));
    EXPECT_NO_THROW(mesh.setGridSpacingInZ(dz));

    EXPECT_EQ(mesh.getGridSpacingInX(), dx);
    EXPECT_EQ(mesh.getGridSpacingInY(), dy);
    EXPECT_EQ(mesh.getGridSpacingInZ(), dz);

    // Origin
    mesh.setOriginInX(x0);
    mesh.setOriginInY(y0);
    mesh.setOriginInZ(z0);
    EXPECT_EQ(mesh.getOriginInX(), x0);
    EXPECT_EQ(mesh.getOriginInY(), y0);
    EXPECT_EQ(mesh.getOriginInZ(), z0);

    std::vector<double> cellScalar1(mesh.getNumberOfCells());
    std::vector<double> cellScalar2(mesh.getNumberOfCells());
    std::vector<double> nodeScalar1(mesh.getNumberOfGridPoints());
    std::vector<double> nodeScalar2(mesh.getNumberOfGridPoints());
    int iter = 1;
    for (int k=0; k<nz-1; ++k)
    {
        for (int j=0; j<ny-1; ++j)
        {
            for (int i=0; i<nx-1; ++i)
            {
                int icell = i*(nz-1)*(ny-1) + j*(nz-1) + k;
                int jcell = k*(nx-1)*(ny-1) + j*(nx-1) + i;
                cellScalar1[icell] = -iter;
                cellScalar2[jcell] = +iter;
                int icellx, icelly, icellz;
                mesh.convertCellIndexToGrid(jcell, &icellx, &icelly, &icellz);
                EXPECT_EQ(icellx, i);
                EXPECT_EQ(icelly, j);
                EXPECT_EQ(icellz, k);
                double x, y, z;
                double xref = x0 + i*dx + dx/2;
                double yref = y0 + j*dy + dy/2;
                double zref = z0 + k*dz + dz/2;
                mesh.convertCellIndexToPosition(jcell, &x, &y, &z);
                EXPECT_NEAR(x, xref, 1.e-14);
                EXPECT_NEAR(y, yref, 1.e-14);
                EXPECT_NEAR(z, zref, 1.e-14);
                iter = iter + 1;
            }
        }
    }
    iter = 1;
    for (int k=0; k<nz; ++k)
    {
        for (int j=0; j<ny; ++j)
        {
            for (int i=0; i<nx; ++i)
            {
                int igrd = i*nz*ny + j*nz + k;
                int jgrd = k*nx*ny + j*nx + i;
                int ix, iy, iz;
                mesh.convertNodeIndexToGrid(jgrd, &ix, &iy, &iz);
                EXPECT_EQ(i, ix);
                EXPECT_EQ(j, iy);
                EXPECT_EQ(k, iz);
                double x, y, z;
                double xref = x0 + i*dx;
                double yref = y0 + j*dy;
                double zref = z0 + k*dz;
                mesh.convertGridIndexToPosition(jgrd, &x, &y, &z);
                EXPECT_NEAR(x, xref, 1.e-14);
                EXPECT_NEAR(y, yref, 1.e-14);
                EXPECT_NEAR(z, zref, 1.e-14);
                nodeScalar1[igrd] =-iter; //(igrd + 1);
                nodeScalar2[jgrd] = iter; //(igrd + 1);
                iter = iter + 1;
            }
        }
    }
    // Set cell data
    EXPECT_FALSE(mesh.haveCellularScalarField("cellTest1"));
    EXPECT_FALSE(mesh.haveCellularScalarField("cellTest2"));
    mesh.setCellularScalarField("cellTest1",
                                cellScalar1.size(),
                                cellScalar1.data(),
                                RegularMesh3DOrderingType::NX_NY_NZ); 
    mesh.setCellularScalarField("cellTest2",
                                cellScalar2.size(),
                                cellScalar2.data(),
                                RegularMesh3DOrderingType::NZ_NY_NX);
    EXPECT_TRUE(mesh.haveCellularScalarField("cellTest1"));
    EXPECT_TRUE(mesh.haveCellularScalarField("cellTest2"));
    // Set nodal data
    EXPECT_FALSE(mesh.haveNodalScalarField("nodalTest1"));
    EXPECT_FALSE(mesh.haveNodalScalarField("nodalTest2"));
    mesh.setNodalScalarField("nodalTest1",
                             nodeScalar1.size(),
                             nodeScalar1.data(),
                             RegularMesh3DOrderingType::NX_NY_NZ);
    mesh.setNodalScalarField("nodalTest2",
                             nodeScalar2.size(),
                             nodeScalar2.data(),
                             RegularMesh3DOrderingType::NZ_NY_NX);
    EXPECT_TRUE(mesh.haveNodalScalarField("nodalTest1"));
    EXPECT_TRUE(mesh.haveNodalScalarField("nodalTest2"));

    // Verify mins/maxes are right
    auto cellMin1 = mesh.getCellularScalarFieldMinValueAndIndex("cellTest1");
    auto cellMax1 = mesh.getCellularScalarFieldMaxValueAndIndex("cellTest1");
    auto nodeMin1 = mesh.getNodalScalarFieldMinValueAndIndex("nodalTest1");
    auto nodeMax1 = mesh.getNodalScalarFieldMaxValueAndIndex("nodalTest1");
    EXPECT_EQ(static_cast<int> (nodeMin1.first),  -mesh.getNumberOfGridPoints());
    EXPECT_EQ(static_cast<int> (nodeMax1.first),  -1);
    EXPECT_EQ(static_cast<int> (cellMin1.first),  -mesh.getNumberOfCells());
    EXPECT_EQ(static_cast<int> (cellMax1.first),  -1);
    EXPECT_EQ(static_cast<int> (nodeMin1.second), mesh.getNumberOfGridPoints()-1);
    EXPECT_EQ(static_cast<int> (nodeMax1.second), 0);
    EXPECT_EQ(static_cast<int> (cellMin1.second), mesh.getNumberOfCells()-1);
    EXPECT_EQ(static_cast<int> (cellMax1.second), 0);

    auto cellMin2 = mesh.getCellularScalarFieldMinValueAndIndex("cellTest2");
    auto cellMax2 = mesh.getCellularScalarFieldMaxValueAndIndex("cellTest2");
    auto nodeMin2 = mesh.getNodalScalarFieldMinValueAndIndex("nodalTest2");
    auto nodeMax2 = mesh.getNodalScalarFieldMaxValueAndIndex("nodalTest2");
    EXPECT_EQ(static_cast<int> (nodeMin2.first),  1);
    EXPECT_EQ(static_cast<int> (nodeMax2.first),  mesh.getNumberOfGridPoints());
    EXPECT_EQ(static_cast<int> (cellMin2.first),  1);
    EXPECT_EQ(static_cast<int> (cellMax2.first),  mesh.getNumberOfCells());

    EXPECT_EQ(static_cast<int> (nodeMin2.second), 0);
    EXPECT_EQ(static_cast<int> (nodeMax2.second), mesh.getNumberOfGridPoints()-1);
    EXPECT_EQ(static_cast<int> (cellMin2.second), 0);
    EXPECT_EQ(static_cast<int> (cellMax2.second), mesh.getNumberOfCells()-1);

    // Test copy
    RegularMesh3D meshCopy(mesh);
    EXPECT_TRUE(mesh.haveGridDimensions());
    EXPECT_EQ(meshCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(meshCopy.getNumberOfGridPointsInY(), ny);
    EXPECT_EQ(meshCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_EQ(meshCopy.getNumberOfCells(), (nx-1)*(ny-1)*(nz-1));
    EXPECT_EQ(meshCopy.getNumberOfGridPoints(), nx*ny*nz);
    EXPECT_EQ(meshCopy.getOriginInX(), x0);
    EXPECT_EQ(meshCopy.getOriginInY(), y0);
    EXPECT_EQ(meshCopy.getOriginInZ(), z0);


}

}
