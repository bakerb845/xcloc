# Introduction

xcloc is a [multilateration](https://en.wikipedia.org/wiki/Multilateration)-based earthquake location tool.  The algorithmic details are summarized in [Dales et. al., 2017](https://github.com/bakerb845/xcloc/blob/master/docs/DalesGJI-2017.pdf).  To gain a much deeper understanding of the methodology it is recommended to read [Ficthner et. al., 2016](  https://github.com/bakerb845/xcloc/blob/master/docs/Fichtner2016.pdf).

# Building

Before continuing, understand that this code is intended to run on Intel hardware.  I believe many portions of the algorithm could be ported to GPU's via OpenCL however this task would require a substantial amount of kernel code - hence the reason it doesn't exist.  Given this the Intel build instructions are discussed.

## Prerequisites

- CMake
- MPI
- Intel MKL
- Intel Performance Primitives
- HDF5
- libxml2 

Optionally, it will help to have
- libiscl
- sacio
to build the unit tests.

