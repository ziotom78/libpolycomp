/* test_rle.c - Tests for RLE functions
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

#include <libpolycomp.h>
#include <assert.h>
#include <stdlib.h>

/***********************************************************************
 * Check that the RLE routines encode the output sequence in the
 * required format */

#define DEFINE_DIFFRLE_BINARY_FORMAT_TEST(fn_name, pcomp_fn,           \
                                          datatype_t)                  \
    void fn_name(void)                                                 \
    {                                                                  \
        datatype_t input_buf[] = { 10, 12, 14, 16, 17, 18 };           \
        size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);  \
        size_t output_size = pcomp_diffrle_bufsize(input_size)         \
                             * sizeof(datatype_t);                     \
        datatype_t* output_buf = malloc(output_size);                  \
                                                                       \
        assert(pcomp_fn(output_buf, &output_size, input_buf,           \
                        input_size) == PCOMP_STAT_SUCCESS);            \
                                                                       \
        assert(output_size == 5);                                      \
        assert(output_buf[0] == 10); /* Starting value */              \
        assert(output_buf[1] == 3); /* Count */                        \
        assert(output_buf[2] == 2); /* Increment */                    \
        assert(output_buf[3] == 2); /* Count */                        \
        assert(output_buf[4] == 1); /* Increment */                    \
                                                                       \
        free(output_buf);                                              \
    }

DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_int8,
                                  pcomp_compress_diffrle_int8, int8_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_int16,
                                  pcomp_compress_diffrle_int16, int16_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_int32,
                                  pcomp_compress_diffrle_int32, int32_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_int64,
                                  pcomp_compress_diffrle_int64, int64_t)

DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_uint8,
                                  pcomp_compress_diffrle_uint8, uint8_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_uint16,
                                  pcomp_compress_diffrle_uint16,
                                  uint16_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_uint32,
                                  pcomp_compress_diffrle_uint32,
                                  uint32_t)
DEFINE_DIFFRLE_BINARY_FORMAT_TEST(test_diffrle_binary_format_uint64,
                                  pcomp_compress_diffrle_uint64,
                                  uint64_t)

void test_diffrle_binary_format(void)
{
    test_diffrle_binary_format_int8();
    test_diffrle_binary_format_int16();
    test_diffrle_binary_format_int32();
    test_diffrle_binary_format_int64();

    test_diffrle_binary_format_uint8();
    test_diffrle_binary_format_uint16();
    test_diffrle_binary_format_uint32();
    test_diffrle_binary_format_uint64();
}

/***********************************************************************
 * Check that the RLE routines are able to properly decompress a
 * (long) sequence of data. */

#define DEFINE_DIFFRLE_DECOMPRESS_TEST(fn_name, pcomp_compr_fn,        \
                                       pcomp_decompr_fn, datatype_t)   \
    void fn_name(void)                                                 \
    {                                                                  \
        const size_t input_size = 10000;                               \
        datatype_t* input_buf                                          \
            = malloc(input_size * sizeof(datatype_t));                 \
        size_t compr_size = pcomp_diffrle_bufsize(input_size)          \
                            * sizeof(datatype_t);                      \
        datatype_t* compr_buf = malloc(compr_size);                    \
        datatype_t* decompr_buf                                        \
            = malloc(input_size * sizeof(datatype_t));                 \
        size_t decompr_size = input_size;                              \
        size_t idx;                                                    \
                                                                       \
        for (idx = 0; idx < input_size; ++idx) {                       \
            /* Pick a number between 0 and 9 */                        \
            input_buf[idx] = random() % 10;                            \
        }                                                              \
                                                                       \
        assert(pcomp_compr_fn(compr_buf, &compr_size, input_buf,       \
                              input_size) == PCOMP_STAT_SUCCESS);      \
        assert(pcomp_decompr_fn(decompr_buf, decompr_size, compr_buf,  \
                                compr_size) == PCOMP_STAT_SUCCESS);    \
                                                                       \
        assert(decompr_size == input_size);                            \
        for (idx = 0; idx < input_size; ++idx) {                       \
            assert(decompr_buf[idx] == input_buf[idx]);                \
        }                                                              \
                                                                       \
        free(input_buf);                                               \
        free(compr_buf);                                               \
        free(decompr_buf);                                             \
    }

DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_int8,
                               pcomp_compress_diffrle_int8,
                               pcomp_decompress_diffrle_int8, int8_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_int16,
                               pcomp_compress_diffrle_int16,
                               pcomp_decompress_diffrle_int16, int16_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_int32,
                               pcomp_compress_diffrle_int32,
                               pcomp_decompress_diffrle_int32, int32_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_int64,
                               pcomp_compress_diffrle_int64,
                               pcomp_decompress_diffrle_int64, int64_t)

DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_uint8,
                               pcomp_compress_diffrle_uint8,
                               pcomp_decompress_diffrle_uint8, uint8_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_uint16,
                               pcomp_compress_diffrle_uint16,
                               pcomp_decompress_diffrle_uint16,
                               uint16_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_uint32,
                               pcomp_compress_diffrle_uint32,
                               pcomp_decompress_diffrle_uint32,
                               uint32_t)
DEFINE_DIFFRLE_DECOMPRESS_TEST(test_diffrle_decompress_uint64,
                               pcomp_compress_diffrle_uint64,
                               pcomp_decompress_diffrle_uint64,
                               uint64_t)

void test_diffrle_decompression(void)
{
    test_diffrle_decompress_int8();
    test_diffrle_decompress_int16();
    test_diffrle_decompress_int32();
    test_diffrle_decompress_int64();

    test_diffrle_decompress_uint8();
    test_diffrle_decompress_uint16();
    test_diffrle_decompress_uint32();
    test_diffrle_decompress_uint64();
}

int main(void)
{
    test_diffrle_binary_format();
    test_diffrle_decompression();

    return 0;
}
