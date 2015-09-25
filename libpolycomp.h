/* libpolycomp.h - interface to "polycomp", a compression library
 *                 aimed to smooth, noise-free data
 *
 * Copyright (c) 2015 Maurizio Tomasi
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef LIBPOLYCOMP_H_GUARD
#define LIBPOLYCOMP_H_GUARD

#include <stddef.h>
#include <stdint.h>

/***********************************************************************
 * Error codes
 */
#define PCOMP_STAT_SUCCESS 0 /* All ok */
#define PCOMP_STAT_INVALID_ENCODING 1 /* Decompression error */
#define PCOMP_STAT_INVALID_BUFFER 2 /* Output buffer too small */
#define PCOMP_STAT_INVALID_FIT 3 /* Least-square fit error */

/***********************************************************************
 * Version information
 */

/* Return the major/minor version of the library */
void pcomp_version(int* major, int* minor);

/***********************************************************************
 * RLE compression functions
 *
 * All of the following function perform a direct/inverse run-length
 * encoding of the input data. The buffer "output_buf" must be
 * preallocated with a number of elements which is at least equal to
 * the value returned by the function "pcomp_rle_bufsize".
 *
 * The inverse encoding (decompression) routines require the user to
 * know how large "output_buf" it will be. This can be achieved by
 * saving this information somewhere during compression.
 */

/* Return the minimum number of elements necessary for the buffer that
 * will contain the compressed data, if the input data has
 * "input_size" elements. */
size_t pcomp_rle_bufsize(size_t input_size);

/* Compression of integer sequences */
int pcomp_compress_rle_int8(int8_t* output_buf, size_t* output_size,
                            const int8_t* input_buf, size_t input_size);
int pcomp_compress_rle_int16(int16_t* output_buf, size_t* output_size,
                             const int16_t* input_buf,
                             size_t input_size);
int pcomp_compress_rle_int32(int32_t* output_buf, size_t* output_size,
                             const int32_t* input_buf,
                             size_t input_size);
int pcomp_compress_rle_int64(int64_t* output_buf, size_t* output_size,
                             const int64_t* input_buf,
                             size_t input_size);

/* Compression of unsigned integer sequences */
int pcomp_compress_rle_uint8(uint8_t* output_buf, size_t* output_size,
                             const uint8_t* input_buf,
                             size_t input_size);
int pcomp_compress_rle_uint16(uint16_t* output_buf, size_t* output_size,
                              const uint16_t* input_buf,
                              size_t input_size);
int pcomp_compress_rle_uint32(uint32_t* output_buf, size_t* output_size,
                              const uint32_t* input_buf,
                              size_t input_size);
int pcomp_compress_rle_uint64(uint64_t* output_buf, size_t* output_size,
                              const uint64_t* input_buf,
                              size_t input_size);

/* Decompression of integer sequences */
int pcomp_decompress_rle_int8(int8_t* output_buf, size_t output_size,
                              const int8_t* input_buf,
                              size_t input_size);
int pcomp_decompress_rle_int16(int16_t* output_buf, size_t output_size,
                               const int16_t* input_buf,
                               size_t input_size);
int pcomp_decompress_rle_int32(int32_t* output_buf, size_t output_size,
                               const int32_t* input_buf,
                               size_t input_size);
int pcomp_decompress_rle_int64(int64_t* output_buf, size_t output_size,
                               const int64_t* input_buf,
                               size_t input_size);

/* Decompression of unsigned integer sequences */
int pcomp_decompress_rle_uint8(uint8_t* output_buf, size_t output_size,
                               const uint8_t* input_buf,
                               size_t input_size);
int pcomp_decompress_rle_uint16(uint16_t* output_buf,
                                size_t output_size,
                                const uint16_t* input_buf,
                                size_t input_size);
int pcomp_decompress_rle_uint32(uint32_t* output_buf,
                                size_t output_size,
                                const uint32_t* input_buf,
                                size_t input_size);
int pcomp_decompress_rle_uint64(uint64_t* output_buf,
                                size_t output_size,
                                const uint64_t* input_buf,
                                size_t input_size);

/***********************************************************************
 * Differenced RLE compression functions
 *
 * All of the following function perform a direct/inverse run-length
 * encoding of consecutive differences in the samples of the input
 * data. The buffer "output_buf" must be preallocated with a number of
 * elements which is at least equal to the value returned by the
 * function "pcomp_diffrle_bufsize".
 *
 * The inverse encoding (decompression) routines require the user to
 * know how large "output_buf" it will be. This can be achieved by
 * saving this information somewhere during compression.
 */

/* Return the minimum number of elements necessary for the buffer that
 * will contain the compressed data. */
size_t pcomp_diffrle_bufsize(size_t input_size);

/* Compression of integer sequences */
int pcomp_compress_diffrle_int8(int8_t* output_buf, size_t* output_size,
                                const int8_t* input_buf,
                                size_t input_size);
int pcomp_compress_diffrle_int16(int16_t* output_buf,
                                 size_t* output_size,
                                 const int16_t* input_buf,
                                 size_t input_size);
int pcomp_compress_diffrle_int32(int32_t* output_buf,
                                 size_t* output_size,
                                 const int32_t* input_buf,
                                 size_t input_size);
int pcomp_compress_diffrle_int64(int64_t* output_buf,
                                 size_t* output_size,
                                 const int64_t* input_buf,
                                 size_t input_size);

/* Compression of unsigned integer sequences */
int pcomp_compress_diffrle_uint8(uint8_t* output_buf,
                                 size_t* output_size,
                                 const uint8_t* input_buf,
                                 size_t input_size);
int pcomp_compress_diffrle_uint16(uint16_t* output_buf,
                                  size_t* output_size,
                                  const uint16_t* input_buf,
                                  size_t input_size);
int pcomp_compress_diffrle_uint32(uint32_t* output_buf,
                                  size_t* output_size,
                                  const uint32_t* input_buf,
                                  size_t input_size);
int pcomp_compress_diffrle_uint64(uint64_t* output_buf,
                                  size_t* output_size,
                                  const uint64_t* input_buf,
                                  size_t input_size);

/* Decompression of integer sequences */
int pcomp_decompress_diffrle_int8(int8_t* output_buf,
                                  size_t output_size,
                                  const int8_t* input_buf,
                                  size_t input_size);
int pcomp_decompress_diffrle_int16(int16_t* output_buf,
                                   size_t output_size,
                                   const int16_t* input_buf,
                                   size_t input_size);
int pcomp_decompress_diffrle_int32(int32_t* output_buf,
                                   size_t output_size,
                                   const int32_t* input_buf,
                                   size_t input_size);
int pcomp_decompress_diffrle_int64(int64_t* output_buf,
                                   size_t output_size,
                                   const int64_t* input_buf,
                                   size_t input_size);

/* Decompression of unsigned integer sequences */
int pcomp_decompress_diffrle_uint8(uint8_t* output_buf,
                                   size_t output_size,
                                   const uint8_t* input_buf,
                                   size_t input_size);
int pcomp_decompress_diffrle_uint16(uint16_t* output_buf,
                                    size_t output_size,
                                    const uint16_t* input_buf,
                                    size_t input_size);
int pcomp_decompress_diffrle_uint32(uint32_t* output_buf,
                                    size_t output_size,
                                    const uint32_t* input_buf,
                                    size_t input_size);
int pcomp_decompress_diffrle_uint64(uint64_t* output_buf,
                                    size_t output_size,
                                    const uint64_t* input_buf,
                                    size_t input_size);

/***********************************************************************
 * Quantization
 *
 * The following functions perform a quantization/dequantization of
 * the input data. The buffer "output_buf" must have been preallocated
 * with a number of bytes which is at least equal to the value of the
 * function.
 *
 * The inverse encoding (decompression) routines require the user to
 * know how large "output_buf" it will be. This can be achieved by
 * saving this information somewhere during compression.
 */

/* Parameters used in the quantization. */
typedef struct {
    size_t element_size;
    size_t bits_per_sample;
    double min_value;
    double normalization;
} pcomp_quant_params_t;

/* Return the minimum number of bytes necessary for the buffer that
 * will contain the compressed data, if the input data has
 * "input_size" elements, each requiring "element_size" bytes, and if
 * quantization will use a number of bits per sample equal to
 * "bits_per_sample". Return zero if the input is invalid. */
size_t pcomp_quant_bufsize(size_t input_size,
                           const pcomp_quant_params_t* params);

int pcomp_compress_quant_float(void* output_buf, size_t* output_size,
                               const float* input_buf,
                               size_t input_size,
                               pcomp_quant_params_t* params);
int pcomp_compress_quant_double(void* output_buf, size_t* output_size,
                                const double* input_buf,
                                size_t input_size,
                                pcomp_quant_params_t* params);

int pcomp_decompress_quant_float(float* output_buf, size_t output_size,
                                 const void* input_buf,
                                 size_t input_size,
                                 const pcomp_quant_params_t* params);
int pcomp_decompress_quant_double(double* output_buf,
                                  size_t output_size,
                                  const void* input_buf,
                                  size_t input_size,
                                  const pcomp_quant_params_t* params);

/***********************************************************************
 * Polynomial fitting routines
 */

struct __pcomp_poly_fit_data_t;
typedef struct __pcomp_poly_fit_data_t pcomp_poly_fit_data_t;

pcomp_poly_fit_data_t* pcomp_init_poly_fit(size_t num_of_points,
                                           size_t num_of_coeffs);
void pcomp_free_poly_fit(pcomp_poly_fit_data_t* poly_fit);
int pcomp_run_poly_fit(pcomp_poly_fit_data_t* poly_fit, double* coeffs,
                       const double* points);

/***********************************************************************
 * Chebyshev transform routines
 */

struct __pcomp_chebyshev_t;
typedef struct __pcomp_chebyshev_t pcomp_chebyshev_t;
typedef enum {
    PCOMP_TD_DIRECT,
    PCOMP_TD_INVERSE
} pcomp_transform_direction_t;

pcomp_chebyshev_t*
pcomp_init_chebyshev(size_t num_of_elements,
                     pcomp_transform_direction_t dir);
void pcomp_free_chebyshev(pcomp_chebyshev_t* plan);
int pcomp_run_chebyshev(pcomp_chebyshev_t* plan,
                        pcomp_transform_direction_t dir, double* output,
                        const double* input);

/***********************************************************************
 * Polynomial compression (low-level functions)
 */

struct __pcomp_polycomp_data_t;
typedef struct __pcomp_polycomp_data_t pcomp_polycomp_data_t;

pcomp_polycomp_data_t* pcomp_init_polycomp_data(size_t num_of_samples,
                                                size_t num_of_coeffs);
void pcomp_free_polycomp_data(pcomp_polycomp_data_t* params);

/***********************************************************************
 * Polynomial compression (high-level functions)
 *
 * The following functions implement the "polynomial compression",
 * which is based on a mixture of polynomial least-square fitting and
 * Chebyshev transform techniques. Unlike the other functions defined
 * above, this group handles memory allocation autonomously. This
 * means that "output_buf" must not be preallocated before calling one
 * of the *_compress_* functions, and the pre-existing value of
 * "output_size" is ignored in the call (it is only set on exit).
 */

/* Structure used to hold information about a chunk of data compressed
 * using the polynomial compression */
typedef struct {
    /* Number of samples in this chunk */
    size_t num_of_elements;

    /* Is this chunk compressed using polynomial/Chebyshev
     * coefficients? */
    int is_compressed;
    /* If the chunk is not compressed (is_compressed == 0), this
     * points to a buffer which holds "num_of_elements" uncompressed
     * samples */
    double* uncompressed;

    /* Polynomial coefficients, from the lowest-order to the
     * highest-order */
    size_t num_of_poly_coeffs;
    double* poly_coeffs;

    /* Chebyshev coefficients */
    size_t num_of_cheby_coeffs; /* This is always less than
                                 * num_of_elements, as the Chebyshev
                                 * series is truncated. */
    double* cheby_coeffs;
} pcomp_poly_chunk_t;

/* Parameters used for the polynomial compression */
typedef struct {
    size_t chunk_size;
    size_t num_of_poly_coeffs;
    double max_error;
} pcomp_poly_parameters;

int pcomp_compress_poly_float(pcomp_poly_chunk_t** output_buf,
                              size_t* num_of_chunks,
                              const float* input_buf, size_t input_size,
                              const pcomp_poly_parameters* params);
int pcomp_compress_poly_double(pcomp_poly_chunk_t** output_buf,
                               size_t* num_of_chunks,
                               const double* input_buf,
                               size_t input_size,
                               const pcomp_poly_parameters* params);

int pcomp_decompress_poly_float(float* output_buf, size_t* output_size,
                                const pcomp_poly_chunk_t* chunk_array,
                                size_t num_of_chunks,
                                const pcomp_poly_parameters* params);
int pcomp_decompress_poly_double(double* output_buf,
                                 size_t* output_size,
                                 const pcomp_poly_chunk_t* chunk_array,
                                 size_t num_of_chunks,
                                 const pcomp_poly_parameters* params);

/* Free the heap memory allocated for this chunk */
int pcomp_free_chunk(pcomp_poly_chunk_t* chunk);

/* Free the heap memory allocated for an array chunk */
int pcomp_free_chunks(pcomp_poly_chunk_t* chunk_array,
                      size_t num_of_chunks);

#endif /* LIBPOLYCOMP_H_GUARD */
