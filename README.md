# libpolycomp

A library which implements a polynomial compression, as well as other
well-known compression techniques.

## Building the library

You need the following dependencies:

1. CMake
2. GNU GSL
3. FFTW

Enter the directory containing the `libpolycomp` sources and run the
following commands:

    mkdir build
    cd build
    cmake ..
    make

The source code is written using C98, so it should compile almost
everywhere.
