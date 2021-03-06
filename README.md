# libpolycomp

[![ascl:1604.002](https://img.shields.io/badge/ascl-1604.002-blue.svg?colorB=262255)](http://ascl.net/1604.002)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A C library which implements a polynomial compression, as well as other
well-known compression techniques.

## Building the library

You need the following dependencies:

1. A C compiler (tested compilers: GCC 4.8, Clang 3.4, Intel C compiler 16.0)
2. CMake (https://cmake.org/)
3. GNU Scientific Library (http://www.gnu.org/software/gsl/)
4. FFTW (http://www.fftw.org/)

Enter the directory containing the `libpolycomp` sources and run the
following commands:

    mkdir build
    cd build
    cmake ..
    make

If you need to specify custom locations for the include/library files,
use CMake's variables `CMAKE_INCLUDE_PATH` and `CMAKE_LIBRARY_PATH`:

    cmake -DCMAKE_INCLUDE_PATH=/myinclude -DCMAKE_LIBRARY_PATH=/mylib ..
    
By default, both a static and dynamic library will be compiled, with
names "libpolycomp_static.a" and "libpolycomp.so". To install them
system-wide, run the following command as root:

    make install

If you do *not* want to install the library system-wide but prefer to
specify a custom location, use the variable
`CMAKE_INSTALL_PREFIX:PATH`:

    cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/mypolycomp ..

The source code is written using C89, so it should compile almost
everywhere.

### OpenMP

By default, OpenMP will be used if available. To disable OpenMP,
change the ``cmake`` invocation above as follows:

    # snip
    cmake -DENABLE_OPENMP=OFF ..
    # snip

## Documentation

The documentation is kept using Doxygen, an online copy is available
at the site http://ziotom78.github.io/libpolycomp .

## Python bindings

A set of Python 3 bindings, as well as a stand-alone executable which
implements compression/decompression of FITS files, is available here:
https://github.com/ziotom78/polycomp

## How to cite this library

If you use this library in your work, please cite the paper _Polycomp:
efficient and configurable compression of astronomical timelines_
(http://dx.doi.org/10.1016/j.ascom.2016.04.004):
`````
@article{Tomasi201688,
	title = "Polycomp: Efficient and configurable compression of astronomical timelines",
	journal = "Astronomy and Computing",
	volume = "16",
	number = "",
	pages = "88-98",
	year = "2016",
	issn = "2213-1337",
	doi = "http://dx.doi.org/10.1016/j.ascom.2016.04.004",
	author = "M. Tomasi",
}
`````

## License

The code is released under the MIT license.
