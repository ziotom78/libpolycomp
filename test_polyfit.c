/* test_polyfit.c - Tests for least-squares polynomial fits
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

int main(void)
{
    double coeffs[2];
    double points[] = { 1.0, 3.0, 5.0 };

    pcomp_poly_fit_data_t* poly_fit = pcomp_init_poly_fit(3, 2);
    assert(pcomp_run_poly_fit(poly_fit, coeffs, points)
           == PCOMP_STAT_SUCCESS);

    assert(fabs(coeffs[0] - (-1.0)) < 1.0e-7);
    assert(fabs(coeffs[1] - (+2.0)) < 1.0e-7);

    pcomp_free_poly_fit(poly_fit);
    return 0;
}
