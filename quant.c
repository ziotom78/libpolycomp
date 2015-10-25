/* quant.c - Quantization of floating-point numbers
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

#include "libpolycomp.h"
#include <limits.h> /* CHAR_BIT */
#include <math.h>
#include <stdlib.h>

/** \defgroup quant Floating-point quantization
 *
 * ### The algorithm and its applicability
 *
 * The quantization compression is applicable only to series of
 * floating-point numbers. It works by converting the series into a
 * sequence of integer numbers, whose bit size is less than the number
 * of bits required for the input numbers.
 *
 * This compression scheme works well for data acquired by means of
 * some digital process. Typically, the bitsize of such samples is
 * smaller than 32 or 64 bits: in such cases, this encoding allows to
 * achieve good compression ratios with negligible loss of
 * information.
 *
 * This kind of compression is lossy, because the decompressed stream
 * is not identical to the stream before the compression. It is
 * possible to prove that, under quite general assumptions, the
 * difference between the two is a stream of zero-average random
 * numbers. This stream of residuals follows a symmetric, non-Gaussian
 * distribution whose RMS is \f$q^2/12\f$, where \f$q\f$ is the
 * quantization step (i.e., \f$q = 2^{-N}\f$, where \f$N\f$ is the
 * number of bits used in the quantization).
 *
 * ### Implementation details
 *
 * Before compressing the data, the caller must initialize an opaque
 * structure, \ref pcomp_quant_params_t, that contains information
 * about the quantization process. Such structure is created via a
 * call to \ref pcomp_init_quant_params and freed by \ref
 * pcomp_free_quant_params. Access to the fields of the structure is
 * only possible through the following functions:
 *
 * - \ref pcomp_quant_element_size
 * - \ref pcomp_quant_bits_per_sample
 * - \ref pcomp_quant_min_value
 * - \ref pcomp_quant_normalization
 *
 * The caller must pre-allocate the buffer that will hold the
 * quantized stream.
 */

/***********************************************************************
 * Creation/destruction/getter functions for "pcomp_quant_params_t"
 */

struct __pcomp_quant_params_t {
    size_t element_size;
    size_t bits_per_sample;
    double min_value;
    double normalization;
};

/** \ingroup quant
 *
 * \brief Initialize a \ref pcomp_quant_params_t structure
 *
 * The function returns a pointer to a heap-allocated structure. It
 * must be freed by the caller once is no longer used via a call to
 * \ref pcomp_free_quant_params.
 *
 * \param[in] element_size Width (in bytes) of the values to be
 * quantized. It can either be 4 (32-bit floating points) or 8 (64-bit
 * floating points), but the function does not enforce this.
 *
 * \param[in] bits_per_sample Number of bits to be used for each
 * sample after quantization. This value should be less than the
 * number of bits used in the input floating-point numbers, in order
 * to achieve a compression ratio greater than 1.
 *
 * \returns A pointer to the newly created \ref pcomp_quant_params_t
 * structure.
 */
pcomp_quant_params_t* pcomp_init_quant_params(size_t element_size,
                                              size_t bits_per_sample)
{
    pcomp_quant_params_t* params = malloc(sizeof(pcomp_quant_params_t));
    params->element_size = element_size;
    params->bits_per_sample = bits_per_sample;
    params->min_value = 0.0;
    params->normalization = 1.0;

    return params;
}

/** \ingroup quant
 *
 * \brief Free a \ref pcomp_quant_params_t structure
 *
 * Free the memory associated with the structured pointed by \a
 * params. Such structure must have been allocated via a call to \ref
 * pcomp_init_quant_params.
 *
 * If \a params is \c NULL, the function does nothing.
 *
 * \params[in] params Pointer to the structure to free.
 */
void pcomp_free_quant_params(pcomp_quant_params_t* params)
{
    if (params == NULL)
        return;
    free(params);
}

/** \ingroup quant
 *
 * \brief Return the size (in bytes) of the elements to be quantized
 *
 * \returns Either 4 (32-bit floating points) or 8 (64-bit floating
 * points).
 */
size_t pcomp_quant_element_size(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->element_size;
}

/** \ingroup quant
 *
 * \brief Return the number of bits that must be used for each
 * quantized sample
 */
size_t pcomp_quant_bits_per_sample(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->bits_per_sample;
}

/** \ingroup quant
 *
 * \brief Return the normalization constant used for converting
 * floating-point numbers into quantized integers
 *
 * See \ref pcomp_quant_set_normalization for further information.
 */
double pcomp_quant_normalization(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->normalization;
}

/** \ingroup quant
 *
 * \brief Return the additive constant used for converting
 * floating-point numbers into quantized integers
 *
 * See \ref pcomp_quant_set_normalization for further information.
 */
double pcomp_quant_offset(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->min_value;
}

/** \ingroup quant
 *
 * \brief Set the normalization constants (multiplicative and additive)
 * used to quantize floating-point numbers
 *
 * This function allows to change the quantization transform described
 * by \a params via the two additive values \f$S\f$ (\a normalization)
 * and \f$\Delta\f$ (\a offset) used in the formula \f[x_\mathrm{out}
 * = [S \times x_\mathrm{in} + \Delta]\f], where square brackets
 * denote a rounding operation.
 *
 * \param[in] params Pointer to the \ref pcomp_quant_params_t structure
 *to modify
 *
 * \param[in] normalization New value for the multiplicative constant
 *\f$S\f$
 *
 * \param[in] offset New value for the additive constant \f$\Delta\f$
 */
void pcomp_quant_set_normalization(pcomp_quant_params_t* params,
                                   double normalization, double offset)
{
    if (params == NULL)
        abort();

    params->normalization = normalization;
    params->min_value = offset;
}

/***********************************************************************
 * Estimate the size of the buffer needed to store quantized data
 */

/** \ingroup quant
 *
 * \brief Return the size (in bytes) of the buffer that will contain a
 * quantized stream of \a input_size floating point values.
 *
 * Unlike functions like \ref pcomp_rle_bufsize, the value returned by
 * this function is exact, not an upper bound. It represents the
 * number of *bytes* needed to store all the compressed floating point
 * values.
 *
 * \param[in] input_size Number of elements to compress
 *
 * \param[in] params Parameters describing the quantization process
 * (created by \ref pcomp_init_quant_params)
 *
 * \returns The size in bytes of the buffer needed to hold all the \a
 * input_size compressed values.
 */
size_t pcomp_quant_bufsize(size_t input_size,
                           const pcomp_quant_params_t* params)
{
    size_t num_of_bits;
    size_t result;

    if (params == NULL)
        abort();

    if (params->element_size == 0 || params->bits_per_sample == 0)
        return 0;

    num_of_bits = input_size * params->element_size * CHAR_BIT;
    result = num_of_bits / params->bits_per_sample;
    if (num_of_bits % params->bits_per_sample > 0)
        result++;

    return result;
}

/***********************************************************************
 * Quantization compression functions
 */

/** \ingroup quant
 *
 * \fn int pcomp_compress_quant_double(void* output_buf,
 *                                     size_t* output_size,
 *                                     const double* input_buf,
 *                                     size_t input_size,
 *                                     pcomp_quant_params_t* params)
 *
 * \brief Quantize a stream of 64-bit floating point numbers.
 *
 * This function applies a quantization formula to all the numbers in
 * the array \a input_buf. It saves the result as a stream of raw
 * bytes in \a output_buf.
 *
 * An usage example of the function is the following, which uses 5
 * bits per every 64-bit input value:
 * \code{.c}
 * double input_buf[] = { 3.06, 5.31, 2.25, 7.92, 4.86 };
 * size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);
 * pcomp_quant_params_t* params;
 * size_t output_size;
 * const size_t bits_per_element = 5;
 *
 * params = pcomp_init_quant_params(sizeof(input_buf[0]),
 *                                  bits_per_element);
 * output_size = pcomp_quant_bufsize(input_size, params);
 * pcomp_compress_quant_double(output_buf, &output_size, input_buf,
 *                             input_size, params);
 * \endcode
 *
 * \param[out] output_buf Pointer to the buffer that will contain the
 * quantized values
 *
 * \param[inout] output_size Pointer to a variable containing the
 * number of bytes allocated for \a output_buf. This number must be at
 * least equal to the return value of the function \ref
 * pcomp_quant_bufsize. On exit, this will contain the actual number
 * of bytes written to \a output_buf.
 *
 * \param[in] input_buf Pointer to the array of 64-bit floating points
 * to quantize
 *
 * \param[in] input_size Number of 64-bit floating point *elements*
 * (not bytes) in the array \a input_buf that must be quantized.
 *
 * \param[in] params Structure defining the details of the
 * quantization. It should be created via a call to \ref
 * pcomp_init_quant_params.
 *
 * \returns \ref PCOMP_STAT_SUCCESS if the encoding completed
 * successfully. Otherwise, the error code specifies the kind of error
 * occurred during the call.
 */

/** \ingroup quant
 *
 * \fn int pcomp_compress_quant_float(void* output_buf,
 *                                    size_t* output_size,
 *                                    const float* input_buf,
 *                                    size_t input_size,
 *                                    pcomp_quant_params_t* params)
 *
 * \brief Quantize a stream of 32-bit floating point numbers.
 *
 * \param[out] output_buf Pointer to the buffer that will contain the
 * quantized values
 *
 * \param[inout] output_size Pointer to a variable containing the
 * number of bytes allocated for \a output_buf. This number must be at
 * least equal to the return value of the function \ref
 * pcomp_quant_bufsize. On exit, this will contain the actual number
 * of bytes written to \a output_buf.
 *
 * \param[in] input_buf Pointer to the array of 32-bit floating points
 * to quantize
 *
 * \param[in] input_size Number of 32-bit floating point *elements*
 * (not bytes) in the array \a input_buf that must be quantized.
 *
 * \param[in] params Structure defining the details of the
 * quantization. It should be created via a call to \ref
 * pcomp_init_quant_params.
 *
 * \returns \ref PCOMP_STAT_SUCCESS if the encoding completed
 * successfully. Otherwise, the error code specifies the kind of error
 * occurred during the call.
 */

#define DEFINE_FIND_BOUNDS_FN(name, datatype_t)                        \
    static void name(const datatype_t* values, size_t num,             \
                     double* min, double* max)                         \
    {                                                                  \
        size_t idx;                                                    \
                                                                       \
        if (values == NULL || num == 0 || min == NULL || max == NULL)  \
            abort();                                                   \
                                                                       \
        *min = *max = values[0];                                       \
        for (idx = 1; idx < num; ++idx) {                              \
            if (values[idx] < *min)                                    \
                *min = values[idx];                                    \
            else if (values[idx] > *max)                               \
                *max = values[idx];                                    \
        }                                                              \
    }

DEFINE_FIND_BOUNDS_FN(find_bounds_float, float)
DEFINE_FIND_BOUNDS_FN(find_bounds_double, double)

static double two_to(int exponent)
{
    double result = 2.0;
    int idx;
    for (idx = 2; idx <= exponent; ++idx)
        result *= 2.0;

    return result;
}

#define DEFINE_COMPRESS_QUANT_FN(name, find_bounds_fn, datatype_t)     \
    int name(void* output_buf, size_t* output_size,                    \
             const datatype_t* input_buf, size_t input_size,           \
             pcomp_quant_params_t* params)                             \
    {                                                                  \
        double max;                                                    \
        uint8_t* byte_buf                                              \
            = output_buf; /* Casting to "uint8_t" is far easier... */  \
        uint8_t cur_byte_buf = 0;                                      \
        size_t bits_in_buffer = 0;                                     \
        size_t idx;                                                    \
                                                                       \
        if (output_buf == NULL || output_size == NULL                  \
            || input_buf == NULL || params == NULL)                    \
            abort();                                                   \
                                                                       \
        find_bounds_fn(input_buf, input_size, &params->min_value,      \
                       &max);                                          \
        params->normalization = (two_to(params->bits_per_sample)       \
                                 - 1.0) / (max - params->min_value);   \
                                                                       \
        for (idx = 0; idx < input_size; ++idx) {                       \
            size_t bit_idx;                                            \
            double scaled_value                                        \
                = floor((input_buf[idx] - params->min_value)           \
                            * params->normalization                    \
                        + 0.5);                                        \
                                                                       \
            for (bit_idx = 0; bit_idx < params->bits_per_sample;       \
                 ++bit_idx) {                                          \
                uint8_t bit = (uint8_t)floor(fmod(scaled_value, 2.0)); \
                cur_byte_buf = (cur_byte_buf << 1) + bit;              \
                bits_in_buffer++;                                      \
                scaled_value /= 2.0;                                   \
                                                                       \
                if (bits_in_buffer >= CHAR_BIT) {                      \
                    /* Is there enough room for more bytes? */         \
                    if (byte_buf - (uint8_t*)output_buf                \
                        > *output_size)                                \
                        return PCOMP_STAT_INVALID_BUFFER;              \
                                                                       \
                    /* Flush "cur_byte_buf" into "output_buf" (this is \
                     * where byte_buf points to) */                    \
                    *byte_buf++ = cur_byte_buf;                        \
                                                                       \
                    bits_in_buffer = 0;                                \
                    cur_byte_buf = 0;                                  \
                }                                                      \
            }                                                          \
        }                                                              \
                                                                       \
        /* If there are still bytes in the cache, zero-fill the byte   \
         * and                                                         \
         * flush it */                                                 \
        if (bits_in_buffer > 0) {                                      \
            *byte_buf++ = cur_byte_buf;                                \
        }                                                              \
                                                                       \
        *output_size = byte_buf - (uint8_t*)output_buf;                \
        return PCOMP_STAT_SUCCESS;                                     \
    }

DEFINE_COMPRESS_QUANT_FN(pcomp_compress_quant_float, find_bounds_float,
                         float)
DEFINE_COMPRESS_QUANT_FN(pcomp_compress_quant_double,
                         find_bounds_double, double)

/***********************************************************************
 * Quantization decompression functions
 */

#define DEFINE_DECOMPRESS_QUANT_FN(name, datatype_t)                   \
    int name(datatype_t* output_buf, size_t output_size,               \
             const void* input_buf, size_t input_size,                 \
             const pcomp_quant_params_t* params)                       \
    {                                                                  \
        size_t output_idx;                                             \
        uint8_t byte;                                                  \
        size_t bits_in_cache;                                          \
        const uint8_t* byte_buf = input_buf;                           \
                                                                       \
        if (output_buf == NULL || input_buf == NULL || params == NULL) \
            abort();                                                   \
                                                                       \
        if (input_size == 0)                                           \
            return PCOMP_STAT_SUCCESS;                                 \
                                                                       \
        if (output_size == 0)                                          \
            return PCOMP_STAT_INVALID_BUFFER;                          \
                                                                       \
        byte = *byte_buf++;                                            \
        bits_in_cache = CHAR_BIT;                                      \
                                                                       \
        output_idx = 0;                                                \
        while (output_idx < output_size) {                             \
            size_t bit_idx;                                            \
            uint64_t buffer;                                           \
                                                                       \
            buffer = 0;                                                \
            for (bit_idx = 0; bit_idx < params->bits_per_sample;       \
                 ++bit_idx) {                                          \
                buffer = (buffer << 1)                                 \
                         + ((uint64_t)byte >> (CHAR_BIT - 1));         \
                byte <<= 1;                                            \
                bits_in_cache--;                                       \
                if (bits_in_cache == 0 && output_idx < output_size     \
                    && (byte_buf - (uint8_t*)input_buf)                \
                           < input_size) {                             \
                    byte = *byte_buf++;                                \
                    bits_in_cache = CHAR_BIT;                          \
                }                                                      \
            }                                                          \
                                                                       \
            output_buf[output_idx] = buffer / params->normalization    \
                                     + params->min_value;              \
            output_idx++;                                              \
        }                                                              \
                                                                       \
        return PCOMP_STAT_SUCCESS;                                     \
    }

DEFINE_DECOMPRESS_QUANT_FN(pcomp_decompress_quant_float, float)
DEFINE_DECOMPRESS_QUANT_FN(pcomp_decompress_quant_double, double)
