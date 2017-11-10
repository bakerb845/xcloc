# Introduction

xcloc is a [multilateration](https://en.wikipedia.org/wiki/Multilateration)-based earthquake location tool.  The algorithmic details are summarized in [Dales et. al., 2017](https://github.com/bakerb845/xcloc/blob/master/docs/DalesGJI-2017.pdf).  To gain a much deeper understanding of the methodology it is recommended to read [Ficthner et. al., 2016](  https://github.com/bakerb845/xcloc/blob/master/docs/Fichtner2016.pdf).

# Building

Before continuing, understand that this code is intended to run on Intel hardware.  I believe many portions of the algorithm could be ported to GPU's via OpenCL however this task would require a substantial amount of kernel code - hence the reason it doesn't exist.  Given this the Intel build instructions are discussed.

## Prerequisites

To build it is required that one have 

- [CMake](https://cmake.org/) version at least 2.6
- Message Passing Interface version at least 3.  Recommended free versions are [MPICH](https://www.mpich.org/) and [OpenMPI](https://www.open-mpi.org/).
- [Intel MKL](https://software.intel.com/en-us/mkl)
- [Intel Performance Primitives](https://software.intel.com/en-us/intel-ipp)
- [HDF5](https://support.hdfgroup.org/HDF5/)
- [libxml2](http://xmlsoft.org/) 

Optionally, for compiling the unit tests, one much also build

- [libiscl](https://github.com/bakerb845/libiscl)
- [sacio](https://github.com/bakerb845/sacio)

However, these will eventually become project submodules and wrangled into the CMake.

## Configuring CMake

    #!/bin/sh
    /usr/bin/cmake ./ \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_C_FLAGS="-g3 -O2 -xCORE-AVX2 -std=c11 -qopenmp -Wall -Wcomment -Wunused -Wcheck -qopt-report=5" \
    -DXCLOC_USE_MPI=TRUE \
    -DXCLOC_USE_INTEL=TRUE \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include \
    -DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.so;/opt/intel/mkl/lib/intel64_lin/libmkl_core.so" \
    -DIPP_INCLUDE_DIR=/opt/intel/ipp/include \
    -DIPP_LIBRARY="/opt/intel/ipp/lib/intel64_lin/libipps.so;/opt/intel/ipp/lib/intel64_lin/libippvm.so;/opt/intel/ipp/lib/intel64_lin/libippcore.so" \
    -DMPI_C_INCLUDE_PATH=/opt/intel/impi/2018.0.128/include64 \
    -DMPI_C_LIBRARIES=/opt/intel/impi/2018.0.128/lib64/libmpi.so \
    -DH5_C_INCLUDE_DIR=/home/bakerb25/C/hdf5-1.10.1_intel/include \
    -DH5_LIBRARY=/home/bakerb25/C/hdf5-1.10.1_intel/lib/libhdf5.so \
    -DSACIO_INCLUDE_DIR=/home/bakerb25/C/sacio/include \
    -DSACIO_LIBRARY=/home/bakerb25/C/sacio/lib/libsacio_shared.so \
    -DISCL_LIBRARY=/home/bakerb25/C/libiscl/lib/libiscl_shared.so \
    -DISCL_INCLUDE_DIR=/home/bakerb25/C/libiscl/include \
    -DXML2_INCLUDE_DIR=/usr/include/libxml2 \
    -DXML2_LIBRARY=/usr/lib/x86_64-linux-gnu/libxml2.so

