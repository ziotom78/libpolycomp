/* test_polycomp_high_level.c - Tests for high-level polynomial
 *                              compression functions
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

#include <stdio.h>

#define MAX_ERROR 0.1

void test_compression(void)
{
    double input[]
        = { 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 2.0, 6.0, 7.0, 9.0 };
    size_t input_size = sizeof(input) / sizeof(input[0]);
    double* decompr;
    size_t decompr_size;
    pcomp_polycomp_chunk_t** chunks;
    size_t num_of_chunks;
    pcomp_polycomp_t* params
        = pcomp_init_polycomp(4, 2, MAX_ERROR, PCOMP_ALG_USE_CHEBYSHEV);
    size_t idx;

    pcomp_compress_polycomp(&chunks, &num_of_chunks, input, input_size,
                            params);
    assert(num_of_chunks == 3);

    decompr_size = pcomp_total_num_of_samples(chunks, num_of_chunks);
    assert(decompr_size == input_size);

    decompr = malloc(decompr_size * sizeof(double));
    pcomp_decompress_polycomp(decompr, chunks, num_of_chunks);

    for (idx = 0; idx < decompr_size; ++idx) {
        assert(fabs(input[idx] - decompr[idx]) <= MAX_ERROR);
    }

    free(decompr);
    pcomp_free_polycomp(params);
    pcomp_free_chunks(chunks, num_of_chunks);
}

void test_encoding(void)
{
    double input[]
        = { 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 2.0, 6.0, 7.0, 9.0 };
    size_t input_size = sizeof(input) / sizeof(input[0]);
    double* decompr;
    size_t decompr_size;
    pcomp_polycomp_chunk_t** chunks;
    void* buf;
    size_t buf_size;
    size_t num_of_chunks;
    pcomp_polycomp_t* params
        = pcomp_init_polycomp(4, 2, MAX_ERROR, PCOMP_ALG_USE_CHEBYSHEV);
    size_t idx;

    pcomp_compress_polycomp(&chunks, &num_of_chunks, input, input_size,
                            params);

    buf = malloc(pcomp_chunks_num_of_bytes(chunks, num_of_chunks));
    assert(pcomp_encode_chunks(buf, &buf_size, chunks, num_of_chunks)
           == PCOMP_STAT_SUCCESS);
    pcomp_free_chunks(chunks, num_of_chunks);

    assert(pcomp_decode_chunks(&chunks, &num_of_chunks, buf)
           == PCOMP_STAT_SUCCESS);
    free(buf);

    decompr_size = pcomp_total_num_of_samples(chunks, num_of_chunks);
    assert(decompr_size == input_size);

    decompr = malloc(decompr_size * sizeof(double));
    pcomp_decompress_polycomp(decompr, chunks, num_of_chunks);

    for (idx = 0; idx < decompr_size; ++idx) {
        assert(fabs(input[idx] - decompr[idx]) <= MAX_ERROR);
    }

    free(decompr);
    pcomp_free_polycomp(params);
    pcomp_free_chunks(chunks, num_of_chunks);
}

int main(void)
{
    test_compression();
    test_encoding();
    return 0;
}
