# Copyright (c) 2015 Maurizio Tomasi
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required (VERSION 3.0.2)
project (libpolycomp
  VERSION 1.0
  LANGUAGES C)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
#    "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake_modules")

# This works with GCC and clang
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -g")
  
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CFITSIO_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

add_library (polycomp
  diff_rle.c
  quant.c
  rle.c
  version.c)

target_include_directories (polycomp
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

########################################################################
# Test harness

# add_pcomp_test(<target> <sources>...)
#
# Build an executable <target> from the specified list of source files
# and and register it as a test with name <target>.

function(add_pcomp_test target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} polycomp)
  target_link_libraries(${target} m)

  add_test(${target} ${target})
endfunction()

enable_testing()

add_pcomp_test(test_version test_version.c)
add_pcomp_test(test_rle test_rle.c)
add_pcomp_test(test_diff_rle test_diff_rle.c)
add_pcomp_test(test_quant test_quant.c)
