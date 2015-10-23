# libpolycomp

A library which implements a polynomial compression, as well as other
well-known compression techniques.

## Building the library

You need the following dependencies:

1. CMake (https://cmake.org/)
2. GNU Scientific Library (http://www.gnu.org/software/gsl/)
3. FFTW (http://www.fftw.org/)

Enter the directory containing the `libpolycomp` sources and run the
following commands:

    mkdir build
    cd build
    cmake ..
    make

By default, both a static and dynamic library will be compiled, with
names "libpolycomp_static.a" and "libpolycomp.so". To install them
system-wide, run the following command as root:

    make install

The source code is written using C98, so it should compile almost
everywhere.

## Documentation

The documentation is kept using Doxygen, an online copy is available
at the site http://ziotom78.github.io/libpolycomp .
