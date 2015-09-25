/* test_chebyshev.c - Tests for Chebyshev transform functions
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

#define EPSILON 1.0e-7

int main(void)
{
    double points[] = { 0.0, 1.0, 3.0 };
    double coeffs[3];
    double inv_coeffs[3];

    pcomp_chebyshev_t* chebyshev = pcomp_init_chebyshev(
        sizeof(points) / sizeof(points[0]), PCOMP_TD_DIRECT);
    assert(pcomp_run_chebyshev(chebyshev, PCOMP_TD_DIRECT, coeffs,
                               points) == PCOMP_STAT_SUCCESS);

    /* These values have been calculated by hand */
    assert(fabs(coeffs[0] - (+5.0 / 2.0)) < EPSILON);
    assert(fabs(coeffs[1] - (-3.0 / 2.0)) < EPSILON);
    assert(fabs(coeffs[2] - (+1.0 / 2.0)) < EPSILON);

    /* Now check that the inverse transform reconstructs the points
     * correctly */
    assert(pcomp_run_chebyshev(chebyshev, PCOMP_TD_INVERSE, inv_coeffs,
                               coeffs) == PCOMP_STAT_SUCCESS);

    assert(fabs(points[0] - inv_coeffs[0]) < EPSILON);
    assert(fabs(points[1] - inv_coeffs[1]) < EPSILON);
    assert(fabs(points[2] - inv_coeffs[2]) < EPSILON);

    pcomp_free_chebyshev(chebyshev);
    return 0;
}
