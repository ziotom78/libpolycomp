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

/** \file libpolycomp.h
 * \brief Header file for Libpolycomp
 *
 * This file is the only header the user has to include in order to
 * use all the facilities offered by the library. Just use
 * \code
 * #include <libpolycomp.h>
 * \endcode
 * at the beginning of your source files to have the functions and the
 * types implemented by Libpolycomp available.
 *
 * There are several groups of compression functions available;
 * currently, they are:
 *   - Run-Length Encoding (RLE);
 *   - Differenced Run-Length Encoding (diffRLE);
 *   - Quantization;
 *   - Polynomial compression.
 */

#ifndef LIBPOLYCOMP_H_GUARD
#define LIBPOLYCOMP_H_GUARD

#include <stddef.h>
#include <stdint.h>

/***********************************************************************
 * Error codes
 */
#define PCOMP_STAT_SUCCESS 0 /** \brief All ok */
#define PCOMP_STAT_INVALID_ENCODING                                    \
    1 /** \brief Decompression error                                   \
         */
#define PCOMP_STAT_INVALID_BUFFER                                      \
    2 /** \brief Output buffer too small */
#define PCOMP_STAT_INVALID_FIT 3 /** \brief Least-square fit error */

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
struct __pcomp_quant_params_t;
typedef struct __pcomp_quant_params_t pcomp_quant_params_t;

pcomp_quant_params_t* pcomp_init_quant_params(size_t element_size,
                                              size_t bits_per_sample);
void pcomp_free_quant_params(pcomp_quant_params_t* params);

size_t pcomp_quant_element_size(const pcomp_quant_params_t* params);
size_t pcomp_quant_bits_per_sample(const pcomp_quant_params_t* params);
double pcomp_quant_normalization(const pcomp_quant_params_t* params);
double pcomp_quant_offset(const pcomp_quant_params_t* params);

void pcomp_quant_set_normalization(pcomp_quant_params_t* params,
                                   double normalization, double offset);

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

pcomp_poly_fit_data_t* pcomp_init_poly_fit(size_t num_of_samples,
                                           size_t num_of_coeffs);
void pcomp_free_poly_fit(pcomp_poly_fit_data_t* poly_fit);
size_t
pcomp_poly_fit_num_of_samples(const pcomp_poly_fit_data_t* poly_fit);
size_t
pcomp_poly_fit_num_of_coeffs(const pcomp_poly_fit_data_t* poly_fit);

int pcomp_run_poly_fit(pcomp_poly_fit_data_t* poly_fit, double* coeffs,
                       const double* points);

/***********************************************************************
 * Chebyshev transform routines
 */

struct __pcomp_chebyshev_t;
typedef struct __pcomp_chebyshev_t pcomp_chebyshev_t;

/** \ingroup poly
 *
 * \brief Direction of a Chebyshev transform
 */
typedef enum {
    /** \brief Compute a forward Chebyshev transform, with a
     * normalization factor 1/(N + 1) */
    PCOMP_TD_DIRECT = 0,
    /** \brief Compute a backward Chebyshev transform, with a
        normalization factor equal to one. */
    PCOMP_TD_INVERSE = 1
} pcomp_transform_direction_t;

pcomp_chebyshev_t*
pcomp_init_chebyshev(size_t num_of_samples,
                     pcomp_transform_direction_t dir);
void pcomp_free_chebyshev(pcomp_chebyshev_t* plan);
size_t pcomp_chebyshev_num_of_samples(const pcomp_chebyshev_t* plan);
pcomp_transform_direction_t
pcomp_chebyshev_direction(const pcomp_chebyshev_t* plan);
int pcomp_run_chebyshev(pcomp_chebyshev_t* plan,
                        pcomp_transform_direction_t dir, double* output,
                        const double* input);
const double* pcomp_chebyshev_input(const pcomp_chebyshev_t* plan);
const double* pcomp_chebyshev_output(const pcomp_chebyshev_t* plan);

/***********************************************************************
 * Polynomial compression (low-level functions)
 */

typedef uint8_t pcomp_poly_size_t;
typedef uint16_t pcomp_chunk_size_t;

/** \ingroup poly
 *
 * \brief Kind of algorithm used for the polynomial compression
 *
 * See the discussion in the section \ref poly for more information.
 */
typedef enum {
    /** \brief When needed, apply the Chebyshev transform to the
     * residuals of the polynomial fit */
    PCOMP_ALG_USE_CHEBYSHEV = 0,
    /** \brief If the absolute value of the residuals of a polynomial
     * fit are too large, store the data in uncompressed form */
    PCOMP_ALG_NO_CHEBYSHEV = 1
} pcomp_polycomp_algorithm_t;

struct __pcomp_polycomp_t;
typedef struct __pcomp_polycomp_t pcomp_polycomp_t;

struct __pcomp_polycomp_chunk_t;
typedef struct __pcomp_polycomp_chunk_t pcomp_polycomp_chunk_t;

pcomp_polycomp_chunk_t*
pcomp_init_chunk(pcomp_chunk_size_t num_of_samples);
pcomp_polycomp_chunk_t*
pcomp_init_uncompressed_chunk(pcomp_chunk_size_t num_of_samples,
                              const double* samples);
pcomp_polycomp_chunk_t* pcomp_init_compressed_chunk(
    pcomp_chunk_size_t num_of_samples,
    pcomp_poly_size_t num_of_poly_coeffs, const double* poly_coeffs,
    pcomp_chunk_size_t num_of_cheby_coeffs,
    const uint8_t* chebyshev_mask, const double* cheby_coeffs);

void pcomp_free_chunk(pcomp_polycomp_chunk_t* chunk);

pcomp_chunk_size_t
pcomp_chunk_num_of_samples(const pcomp_polycomp_chunk_t* chunk);
size_t pcomp_chunk_num_of_bytes(const pcomp_polycomp_chunk_t* chunk);
int pcomp_chunk_is_compressed(const pcomp_polycomp_chunk_t* chunk);
const double*
pcomp_chunk_uncompressed_data(const pcomp_polycomp_chunk_t* chunk);

pcomp_poly_size_t
pcomp_chunk_num_of_poly_coeffs(const pcomp_polycomp_chunk_t* chunk);
const double*
pcomp_chunk_poly_coeffs(const pcomp_polycomp_chunk_t* chunk);

pcomp_chunk_size_t
pcomp_chunk_num_of_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk);
const double*
pcomp_chunk_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk);
size_t pcomp_chunk_cheby_mask_size(pcomp_chunk_size_t chunk_size);
const uint8_t*
pcomp_chunk_cheby_mask(const pcomp_polycomp_chunk_t* chunk);

void pcomp_straighten(double* output, const double* input,
                      size_t num_of_samples, double period);

pcomp_polycomp_t*
pcomp_init_polycomp(pcomp_chunk_size_t samples_per_chunk,
                    pcomp_poly_size_t num_of_coeffs,
                    double max_allowable_error,
                    pcomp_polycomp_algorithm_t algorithm);
void pcomp_free_polycomp(pcomp_polycomp_t* params);

int pcomp_run_polycomp_on_chunk(pcomp_polycomp_t* params,
                                const double* input,
                                pcomp_chunk_size_t num_of_samples,
                                pcomp_polycomp_chunk_t* chunk,
                                double* max_error);

int pcomp_decompress_polycomp_chunk(double* output,
                                    const pcomp_polycomp_chunk_t* chunk,
                                    pcomp_chebyshev_t* inv_chebyshev);

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

pcomp_chunk_size_t
pcomp_polycomp_samples_per_chunk(const pcomp_polycomp_t* params);
pcomp_poly_size_t
pcomp_polycomp_num_of_poly_coeffs(const pcomp_polycomp_t* params);
double pcomp_polycomp_max_error(const pcomp_polycomp_t* params);
pcomp_polycomp_algorithm_t
pcomp_polycomp_algorithm(const pcomp_polycomp_t* params);
pcomp_chebyshev_t*
pcomp_polycomp_forward_cheby(const pcomp_polycomp_t* params);
pcomp_chebyshev_t*
pcomp_polycomp_backward_cheby(const pcomp_polycomp_t* params);
double pcomp_polycomp_period(const pcomp_polycomp_t* params);

void pcomp_polycomp_set_period(pcomp_polycomp_t* params, double period);

int pcomp_polyfit_and_chebyshev(pcomp_polycomp_t* params,
                                double* coeffs, double* cheby_residuals,
                                const double* input,
                                double* max_residual);

int pcomp_mask_get_bit(const uint8_t* mask, size_t pos);
void pcomp_mask_set_bit(uint8_t* mask, size_t pos, int value);
size_t pcomp_find_chebyshev_mask(pcomp_chebyshev_t* chebyshev,
                                 pcomp_chebyshev_t* inv_chebyshev,
                                 double max_allowable_error,
                                 uint8_t* mask, double* max_error);

int pcomp_compress_polycomp(pcomp_polycomp_chunk_t** chunk_array[],
                            size_t* num_of_chunks,
                            const double* input_buf, size_t input_size,
                            const pcomp_polycomp_t* params);

size_t
pcomp_total_num_of_samples(pcomp_polycomp_chunk_t* const chunk_array[],
                           size_t num_of_chunks);

/* The buffer "output_buf" must have room for a number of "double"
 * values which is at least the number returned by
 * "pcomp_total_num_of_samples". */
int pcomp_decompress_polycomp(
    double* output_buf, pcomp_polycomp_chunk_t* const chunk_array[],
    size_t num_of_chunks);

/* Free the heap memory allocated for an array of chunks */
void pcomp_free_chunks(pcomp_polycomp_chunk_t* chunk_array[],
                       size_t num_of_chunks);

size_t pcomp_chunks_num_of_bytes(pcomp_polycomp_chunk_t* const chunks[],
                                 size_t num_of_chunks);

int pcomp_encode_chunks(void* buf, size_t* buf_size,
                        pcomp_polycomp_chunk_t* const chunk_array[],
                        size_t num_of_chunks);

int pcomp_decode_chunks(pcomp_polycomp_chunk_t** chunk_array[],
                        size_t* num_of_chunks, const void* buf);

#endif /* LIBPOLYCOMP_H_GUARD */
