/* diff_rle.c - Differenced run-length encoding
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
#include <stdlib.h>

size_t pcomp_diffrle_bufsize(size_t input_size)
{
    return 1 + 2 * input_size;
}

/***********************************************************************
 * Differenced run-length compression routines
 */

#define IMPLEMENT_DIFFRLE_COMPR_FN(name, rle_compr_fn, datatype_t)     \
    int name(datatype_t* output_buf, size_t* output_size,              \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t rle_output_size;                                        \
        datatype_t* diff_buf = NULL;                                   \
        int rle_compr_result;                                          \
        size_t idx;                                                    \
                                                                       \
        if (output_buf == NULL || output_size == NULL                  \
            || input_buf == NULL)                                      \
            abort();                                                   \
                                                                       \
        if (input_size == 0) {                                         \
            *output_size = 0;                                          \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        output_buf[0] = input_buf[0];                                  \
        if (input_size == 1) {                                         \
            *output_size = 1;                                          \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        diff_buf = malloc((input_size - 1) * sizeof(datatype_t));      \
        for (idx = 0; idx < input_size - 1; ++idx) {                   \
            diff_buf[idx] = input_buf[idx + 1] - input_buf[idx];       \
        }                                                              \
                                                                       \
        rle_output_size = *output_size - 1;                            \
        rle_compr_result                                               \
            = rle_compr_fn(output_buf + 1, &rle_output_size, diff_buf, \
                           input_size - 1);                            \
        free(diff_buf);                                                \
        if (rle_compr_result != PCOMP_STAT_SUCCESS)                    \
            return rle_compr_result;                                   \
                                                                       \
        *output_size = rle_output_size + 1;                            \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int8,
                           pcomp_compress_rle_int8, int8_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int16,
                           pcomp_compress_rle_int16, int16_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int32,
                           pcomp_compress_rle_int32, int32_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int64,
                           pcomp_compress_rle_int64, int64_t)

IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint8,
                           pcomp_compress_rle_uint8, uint8_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint16,
                           pcomp_compress_rle_uint16, uint16_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint32,
                           pcomp_compress_rle_uint32, uint32_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint64,
                           pcomp_compress_rle_uint64, uint64_t)

/***********************************************************************
 * Differenced run-length decompression routines
 */

#define IMPLEMENT_DIFFRLE_DECOMPR_FN(name, datatype_t)                 \
    int name(datatype_t* output_buf, size_t output_size,               \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t input_idx = 0;                                          \
        size_t output_idx = 0;                                         \
                                                                       \
        if (output_buf == NULL || input_buf == NULL)                   \
            abort();                                                   \
                                                                       \
        if (input_size == 0) {                                         \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        if (input_size % 2 != 1) {                                     \
            return PCOMP_STAT_INVALID_ENCODING;                        \
        }                                                              \
                                                                       \
        output_buf[output_idx++] = input_buf[input_idx++];             \
                                                                       \
        while (output_idx < output_size                                \
               && input_idx < input_size - 1) {                        \
            datatype_t count = input_buf[input_idx];                   \
            datatype_t incr = input_buf[input_idx + 1];                \
            datatype_t idx;                                            \
                                                                       \
            for (idx = 0; idx < count; ++idx) {                        \
                output_buf[output_idx] = output_buf[output_idx - 1]    \
                                         + incr;                       \
                output_idx++;                                          \
            }                                                          \
                                                                       \
            input_idx += 2;                                            \
        }                                                              \
                                                                       \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int8, int8_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int16, int16_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int32, int32_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int64, int64_t)

IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint8, uint8_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint16, uint16_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint32, uint32_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint64, uint64_t)
