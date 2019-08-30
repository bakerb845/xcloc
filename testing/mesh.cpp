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
    int iter = 0;
    for (int k=0; k<nz-1; ++k)
    {
        for (int j=0; j<ny-1; ++j)
        {
            for (int i=0; i<nx-1; ++i)
            {
                cellScalar1[i*(nz-1)*(ny-1) + j*(nz-1) + k] = -iter;
                cellScalar2[k*(nx-1)*(ny-1) + j*(nx-1) + i] = +iter;
                iter = iter + 1;
            }
        }
    }
    for (int i=0; i<mesh.getNumberOfGridPoints(); ++i)
    {
        nodeScalar1[i] =-i;
        nodeScalar2[i] = i;
    }
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
    // Test copy
}

}
