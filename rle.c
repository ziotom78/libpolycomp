/* rle.c - Run-Length encoding
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

size_t pcomp_rle_bufsize(size_t input_size) { return 2 * input_size; }

/***********************************************************************
 * Run-length compression routines
 */

#define IMPLEMENT_RLE_COMPR_FN(name, datatype_t, max_value)            \
    int name(datatype_t* output_buf, size_t* output_size,              \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t true_output_size = 0;                                   \
        size_t input_idx = 0;                                          \
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
        if (*output_size < 2) {                                        \
            return PCOMP_STAT_INVALID_BUFFER;                          \
        }                                                              \
                                                                       \
        while (input_idx < input_size) {                               \
            datatype_t first_datum_in_the_run = input_buf[input_idx];  \
            datatype_t count = 0;                                      \
            while (count < max_value && input_idx < input_size         \
                   && input_buf[input_idx]                             \
                          == first_datum_in_the_run) {                 \
                input_idx++;                                           \
                count++;                                               \
            }                                                          \
                                                                       \
            if (true_output_size >= *output_size - 2) {                \
                return PCOMP_STAT_INVALID_BUFFER;                      \
            }                                                          \
                                                                       \
            output_buf[true_output_size++] = count;                    \
            output_buf[true_output_size++] = first_datum_in_the_run;   \
        }                                                              \
                                                                       \
        *output_size = true_output_size;                               \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int8, int8_t, INT8_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int16, int16_t, INT16_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int32, int32_t, INT32_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int64, int64_t, INT64_MAX)

IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint8, uint8_t, UINT8_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint16, uint16_t, UINT16_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint32, uint32_t, UINT32_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint64, uint64_t, UINT64_MAX)

/***********************************************************************
 * Run-Length decompression routines
 */

#define IMPLEMENT_RLE_DECOMPR_FN(name, datatype_t)                     \
    int name(datatype_t* output_buf, size_t* output_size,              \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t input_idx = 0;                                          \
        size_t output_idx = 0;                                         \
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
        if (input_size % 2 != 0) {                                     \
            return PCOMP_STAT_INVALID_ENCODING;                        \
        }                                                              \
                                                                       \
        while (output_idx < *output_size                               \
               && input_idx < input_size - 1) {                        \
            datatype_t count = input_buf[input_idx];                   \
            datatype_t value = input_buf[input_idx + 1];               \
            datatype_t idx;                                            \
                                                                       \
            for (idx = 0; idx < count; ++idx) {                        \
                output_buf[output_idx++] = value;                      \
            }                                                          \
                                                                       \
            input_idx += 2;                                            \
        }                                                              \
                                                                       \
        *output_size = output_idx;                                     \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int8, int8_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int16, int16_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int32, int32_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int64, int64_t)

IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint8, uint8_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint16, uint16_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint32, uint32_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint64, uint64_t)
