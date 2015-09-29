/* test_polycomp_low_level.c - Tests for low-level polynomial
 *                             compression functions
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

#define EPSILON 1.0e-7
#define MAX_ALLOWABLE_ERROR 0.3

int test_no_compression(void)
{
    /* It is impossible to compress these data using the polynomial
     * compression algorithm, as they contain too many jumps */
    double input[]
        = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0 };
    size_t input_size = sizeof(input) / sizeof(input[0]);
    double* decompr = malloc(sizeof(double) * input_size);
    pcomp_polycomp_chunk_t* chunk = pcomp_init_chunk(input_size);
    double max_error;
    pcomp_chebyshev_t* inv_chebyshev
        = pcomp_init_chebyshev(input_size, PCOMP_TD_INVERSE);
    pcomp_polycomp_t* polycomp = pcomp_init_polycomp(
        input_size, 2, MAX_ALLOWABLE_ERROR, PCOMP_ALG_USE_CHEBYSHEV);
    size_t idx;

    pcomp_run_polycomp_on_chunk(polycomp, input, input_size, chunk,
                                &max_error);
    assert(pcomp_chunk_num_of_elements(chunk) == 10);
    assert(!pcomp_chunk_is_compressed(chunk));
    assert(pcomp_chunk_poly_coeffs(chunk) == NULL);
    assert(pcomp_chunk_cheby_coeffs(chunk) == NULL);

    pcomp_decompress_polycomp_chunk(decompr, chunk, inv_chebyshev);
    for (idx = 0; idx < input_size; ++idx) {
        assert(input[idx] == decompr[idx]);
    }

    free(decompr);
    pcomp_free_chunk(chunk);
    pcomp_free_polycomp(polycomp);

    return 0;
}

int test_no_chebyshev(void)
{
    /* These data are fitted perfectly by a straight line */
    double input[]
        = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    size_t input_size = sizeof(input) / sizeof(input[0]);
    double* decompr = malloc(sizeof(double) * input_size);
    pcomp_polycomp_chunk_t* chunk = pcomp_init_chunk(input_size);
    double max_error;
    pcomp_chebyshev_t* inv_chebyshev
        = pcomp_init_chebyshev(input_size, PCOMP_TD_INVERSE);
    pcomp_polycomp_t* polycomp = pcomp_init_polycomp(
        input_size, 2, MAX_ALLOWABLE_ERROR, PCOMP_ALG_USE_CHEBYSHEV);
    size_t idx;

    pcomp_run_polycomp_on_chunk(polycomp, input, input_size, chunk,
                                &max_error);
    assert(pcomp_chunk_num_of_elements(chunk) == 10);
    assert(pcomp_chunk_is_compressed(chunk));
    assert(pcomp_chunk_num_of_poly_coeffs(chunk) == 2);
    assert(pcomp_chunk_poly_coeffs(chunk) != NULL);
    assert(fabs(pcomp_chunk_poly_coeffs(chunk)[0] - 0.0) < EPSILON);
    assert(fabs(pcomp_chunk_poly_coeffs(chunk)[1] - 1.0) < EPSILON);
    assert(pcomp_chunk_cheby_coeffs(chunk) == NULL);

    pcomp_decompress_polycomp_chunk(decompr, chunk, inv_chebyshev);
    for (idx = 0; idx < input_size; ++idx) {
        assert(fabs(input[idx] - decompr[idx]) < EPSILON);
    }

    free(decompr);
    pcomp_free_chunk(chunk);
    pcomp_free_polycomp(polycomp);

    return 0;
}

int test_complete_compression_and_decompression(void)
{
    double input[]
        = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0 };
    size_t input_size = sizeof(input) / sizeof(input[0]);
    double* decompr = malloc(sizeof(double) * input_size);
    pcomp_polycomp_chunk_t* chunk = pcomp_init_chunk(input_size);
    double max_error;
    pcomp_chebyshev_t* inv_chebyshev
        = pcomp_init_chebyshev(input_size, PCOMP_TD_INVERSE);
    pcomp_polycomp_t* polycomp = pcomp_init_polycomp(
        input_size, 2, MAX_ALLOWABLE_ERROR, PCOMP_ALG_USE_CHEBYSHEV);
    size_t idx;

    pcomp_run_polycomp_on_chunk(polycomp, input, input_size, chunk,
                                &max_error);
    assert(pcomp_chunk_num_of_elements(chunk) == 10);
    assert(pcomp_chunk_is_compressed(chunk));
    assert(pcomp_chunk_num_of_poly_coeffs(chunk) == 2);
    assert(pcomp_chunk_poly_coeffs(chunk) != NULL);
    assert(pcomp_chunk_num_of_cheby_coeffs(chunk) < 10);

    pcomp_decompress_polycomp_chunk(decompr, chunk, inv_chebyshev);
    for (idx = 0; idx < input_size; ++idx) {
        assert(fabs(input[idx] - decompr[idx]) < MAX_ALLOWABLE_ERROR);
    }

    free(decompr);
    pcomp_free_chunk(chunk);
    pcomp_free_polycomp(polycomp);

    return 0;
}

int main(void)
{
    int result;

    result = test_no_compression();
    if (result != 0)
        return result;

    result = test_no_chebyshev();
    if (result != 0)
        return result;

    result = test_complete_compression_and_decompression();
    if (result != 0)
        return result;

    return 0;
}
