# Introduction

xcloc is a [multilateration](https://en.wikipedia.org/wiki/Multilateration)-based earthquake location tool.  The algorithmic details are summarized in [Dales et. al., 2017](https://github.com/bakerb845/xcloc/blob/master/docs/DalesGJI-2017.pdf).  To gain a much deeper understanding of the methodology it is recommended to read [Ficthner et. al., 2016](  https://github.com/bakerb845/xcloc/blob/master/docs/Fichtner2016.pdf).

# Building

Before continuing, understand that this code is intended to run on Intel hardware.  It likely could see a performance boost with coprocessors.  However, the idea of porting to a GPU is daunting and I instead would prefer to run the code on a Xeon Phi first to confirm my suspicions before writing that mountain of kernel code.  Given this the Intel build instructions are discussed.

## Prerequisites

To build it is required that one have 

- [CMake](https://cmake.org/) version at least 2.6
- A fully C11 compliant C compiler.
- A fully Fortran2008 compliant Fortran compiler.
- Message Passing Interface version at least 3.  Recommended free versions are [MPICH](https://www.mpich.org/) and [OpenMPI](https://www.open-mpi.org/).  The Fortran 2008 bindings must be available.
- [Intel MKL](https://software.intel.com/en-us/mkl)
- [Intel Performance Primitives](https://software.intel.com/en-us/intel-ipp)
- [HDF5](https://support.hdfgroup.org/HDF5/)
- [libxml2](http://xmlsoft.org/) 
- [iniparser](https://github.com/ndevilla/iniparser)

Optionally, for compiling the unit tests, one must also build

- [ISCL](https://gitlab.isti.com/bbaker/iscl)
- [SacIO](https://gitlab.isti.com/bbaker/sacio)

However, these will eventually become project submodules and wrangled into the CMake.

## Configuring CMake

I typically find it useful to generate a configuration script like the one below

    #!/bin/sh
    export CC=/opt/intel/bin/icc
    export CXX=/opt/intel/bin/icpc
    export FC=/opt/intel/bin/ifort
    export F90=/opt/intel/bin/ifort
    export MKL_LIB_ROOT=/opt/intel/mkl/lib/intel64
    export IPP_LIB_ROOT=/opt/intel/ipp/lib/intel64
    if [ -f Makefile]; then
       make clean
    fi
    if [ -f CMakeCache.txt ]; then
       echo "Removing CMakeCache.txt"
       rm CMakeCache.txt
    fi
    cmake ./ \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DCMAKE_C_FLAGS="-g -O2 -xCORE-AVX2 -std=c11 -qopenmp -Wall -Wcomment -Wunused -Wcheck -qopt-report=5" \
    -DCMAKE_CXX_FLAGS="-g -O2 -qopenmp -std=c++11 -Wall -Wcomment -Wunused -Wcheck -qopt-report=5" \
    -DCMAKE_Fortran_FLAGS="-g -O2 -qopenmp -W 1 -warn unused -align array64byte -qopt-report=4 -xHOST -nofor-main" \
    -DXCLOC_USE_MPI=TRUE \
    -DXCLOC_USE_INTEL=TRUE \
    -DXCLOC_PROFILE=TRUE \
    -DMKL_INCLUDE_DIR=/opt/intel/mkl/include \
    -DMKL_LIBRARY="${MKL_LIB_ROOT}/libmkl_intel_lp64.so;${MKL_LIB_ROOT}/libmkl_sequential.so;${MKL_LIB_ROOT}/libmkl_core.so;${MKL_LIB_ROOT}/libmkl_vv
ml_avx2.so;${MKL_LIB_ROOT}/libmkl_avx2.so" \
    -DIPP_INCLUDE_DIR=/opt/intel/ipp/include \
    -DIPP_LIBRARY="${IPP_LIB_ROOT}/libipps.so;${IPP_LIB_ROOT}/libippvm.so;${IPP_LIB_ROOT}/libippcore.so" \
    -DMPI_C_INCLUDE_PATH=/opt/intel/impi/2018.0.128/include64 \
    -DMPI_C_LIBRARIES="/opt/intel/impi/2018.0.128/lib64/libmpi.so" \
    -DMPI_Fortran_INCLUDE_PATH=/opt/intel/impi/2018.1.163/include64 \
    -DMPI_Fortran_LIBRARIES="/opt/intel/impi/2018.0.128/lib64/libmpifort.so" \
    -DH5_C_INCLUDE_DIR=/home/bakerb25/C/hdf5-1.10.1_intel/include \
    -DH5_LIBRARY=/home/bakerb25/C/hdf5-1.10.1_intel/lib/libhdf5.so \
    -DSACIO_INCLUDE_DIR=/home/bakerb25/C/sacio/include \
    -DSACIO_LIBRARY=/home/bakerb25/C/sacio/lib/libsacio_shared.so \
    -DISCL_LIBRARY=/home/bakerb25/C/iscl/lib/libiscl_shared.so \
    -DISCL_INCLUDE_DIR=/home/bakerb25/C/iscl/include \
    -DXML2_INCLUDE_DIR=/usr/include/libxml2 \
    -DXML2_LIBRARY=/usr/lib/x86_64-linux-gnu/libxml2.so \
    -DINIPARSER_INCLUDE_DIR=/home/bakerb25/C/iniparser/src \
    -DINIPARSER_LIBRARY=/home/bakerb25/C/iniparser/libiniparser.so.1 \
    -DADVISOR_INCLUDE_DIR=/opt/intel/advisor/include

and run it in the root source directory.  After CMake has been successfully configured one then simply types

    make

in the root source directory and

    make install

to install the library and executables.  If [Doxygen](http://www.stack.nl/~dimitri/doxygen/) is available then the developer level documentation can be built by specifying

    doxygen Doxyfile

in the root source directory.