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

Optionally, it will help to have

- [libiscl](https://github.com/bakerb845/libiscl)
- [sacio](https://github.com/bakerb845/sacio)

which will eventually be incorporated as submodules.


