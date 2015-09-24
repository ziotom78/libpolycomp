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
#include <math.h>
#include <stdlib.h>

/***********************************************************************
 * Check that the quantization routines encode the output sequence in
 * the required format (see the "polycomp" paper, Fig. 3).
 */

#define DEFINE_QUANT_TEST(fn_name, pcomp_fn, datatype_t)               \
    void fn_name(void)                                                 \
    {                                                                  \
        datatype_t input_buf[] = { 3.06, 5.31, 2.25, 7.92, 4.86 };     \
        size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);  \
        size_t output_size = 4; /* Known by manual calculation */      \
        uint8_t* output_buf = malloc(output_size);                     \
                                                                       \
        pcomp_quant_params_t params;                                   \
        params.element_size = sizeof(float);                           \
        params.bits_per_sample = 5;                                    \
                                                                       \
        assert(pcomp_fn(output_buf, &output_size, input_buf,           \
                        input_size, &params) == PCOMP_STAT_SUCCESS);   \
                                                                       \
        assert(output_buf[0] == 36); /* Count */                       \
        assert(output_buf[1] == 65); /* Value */                       \
        assert(output_buf[2] == 247); /* Count */                      \
        assert(output_buf[3] == 0);                                    \
                                                                       \
        free(output_buf);                                              \
    }

DEFINE_QUANT_TEST(test_quant_compression_float,
                  pcomp_compress_quant_float, float)
DEFINE_QUANT_TEST(test_quant_compression_double,
                  pcomp_compress_quant_double, double)

/***********************************************************************
 * Check that the decompression of quantized data works as expected
 * (see the "polycomp" paper, Fig. 3).
 */

#define DEFINE_QUANT_DECOMPR_TEST(fn_name, compr_fn, decompr_fn,       \
                                  datatype_t)                          \
    void fn_name(void)                                                 \
    {                                                                  \
        datatype_t input_buf[] = { 3.06, 5.31, 2.25, 7.92, 4.86 };     \
        size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);  \
        size_t compr_size = 4; /* Known by manual calculation */       \
        uint8_t* compr_buf = malloc(compr_size);                       \
        datatype_t decompr_buf[5];                                     \
        size_t idx;                                                    \
                                                                       \
        pcomp_quant_params_t params;                                   \
        params.element_size = sizeof(float);                           \
        params.bits_per_sample = 5;                                    \
                                                                       \
        assert(compr_fn(compr_buf, &compr_size, input_buf, input_size, \
                        &params) == PCOMP_STAT_SUCCESS);               \
                                                                       \
        assert(decompr_fn(&decompr_buf[0], input_size, compr_buf,      \
                          compr_size, &params) == PCOMP_STAT_SUCCESS); \
                                                                       \
        for (idx = 0; idx < input_size; ++idx) {                       \
            assert(fabs(input_buf[idx] - decompr_buf[idx]) < 0.186);   \
        }                                                              \
                                                                       \
        free(compr_buf);                                               \
    }

DEFINE_QUANT_DECOMPR_TEST(test_quant_decompression_float,
                          pcomp_compress_quant_float,
                          pcomp_decompress_quant_float, float)
DEFINE_QUANT_DECOMPR_TEST(test_quant_decompression_double,
                          pcomp_compress_quant_double,
                          pcomp_decompress_quant_double, double)

int main(void)
{
    test_quant_compression_float();
    test_quant_compression_double();

    test_quant_decompression_float();
    test_quant_decompression_double();

    return 0;
}
