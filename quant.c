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

/***********************************************************************
 * Creation/destruction/getter functions for "pcomp_quant_params_t"
 */

struct __pcomp_quant_params_t {
    size_t element_size;
    size_t bits_per_sample;
    double min_value;
    double normalization;
};

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

void pcomp_free_quant_params(pcomp_quant_params_t* params)
{
    if (params == NULL)
        return;
    free(params);
}

size_t pcomp_quant_element_size(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->element_size;
}

size_t pcomp_quant_bits_per_sample(const pcomp_quant_params_t* params)
{
    if (params == NULL)
        abort();

    return params->bits_per_sample;
}

/***********************************************************************
 * Estimate the size of the buffer needed to store quantized data
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
