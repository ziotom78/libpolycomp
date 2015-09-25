/* poly.c - Polynomial (de)compression routines
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

#include <gsl/gsl_block_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit.h>

#include <fftw3.h>

static double integer_power(int x, int y)
{
    double dbl_x = (double)x;

    if (y < 0)
        abort();

    if (y == 0)
        return 1.0;
    else if (y == 1)
        return dbl_x;
    else {
        double result = dbl_x * dbl_x;
        int cur_power = 2;
        while (2 * cur_power < y) {
            result *= result;
            cur_power *= 2;
        }

        if (y > cur_power)
            result *= integer_power(x, y - cur_power);

        return result;
    }
}

/***********************************************************************
 * Types and functions used for polynomial least-square fitting
 */

struct __pcomp_poly_fit_data_t {
    size_t num_of_points;
    size_t num_of_coeffs;
    gsl_multifit_linear_workspace* workspace;
    gsl_matrix* matrix;
    gsl_vector* y;
    gsl_vector* c;
    gsl_matrix* cov_matrix;
};

pcomp_poly_fit_data_t* pcomp_init_poly_fit(size_t num_of_points,
                                           size_t num_of_coeffs)
{
    size_t i, j;

    pcomp_poly_fit_data_t* poly_fit
        = malloc(sizeof(pcomp_poly_fit_data_t));
    poly_fit->num_of_points = num_of_points;
    poly_fit->num_of_coeffs = num_of_coeffs;
    poly_fit->workspace
        = gsl_multifit_linear_alloc(num_of_points, num_of_coeffs);
    poly_fit->matrix = gsl_matrix_alloc(num_of_points, num_of_coeffs);
    poly_fit->y = gsl_vector_alloc(num_of_points);
    poly_fit->c = gsl_vector_alloc(num_of_coeffs);
    poly_fit->cov_matrix
        = gsl_matrix_alloc(num_of_coeffs, num_of_coeffs);

    for (i = 0; i < num_of_points; ++i) {
        for (j = 0; j < num_of_coeffs; ++j) {
            gsl_matrix_set(poly_fit->matrix, i, j,
                           integer_power(i + 1, j));
        }
    }

    return poly_fit;
}

void pcomp_free_poly_fit(pcomp_poly_fit_data_t* poly_fit)
{
    if (poly_fit == NULL)
        return;

    gsl_matrix_free(poly_fit->matrix);
    gsl_vector_free(poly_fit->y);
    gsl_vector_free(poly_fit->c);
    gsl_matrix_free(poly_fit->cov_matrix);
    gsl_multifit_linear_free(poly_fit->workspace);

    free(poly_fit);
}

int pcomp_run_poly_fit(pcomp_poly_fit_data_t* poly_fit, double* coeffs,
                       const double* points)
{
    size_t idx;
    double chisq;

    if (poly_fit == NULL || coeffs == NULL || points == NULL)
        abort();

    for (idx = 0; idx < poly_fit->num_of_points; ++idx) {
        gsl_vector_set(poly_fit->y, idx, points[idx]);
    }

    if (gsl_multifit_linear(poly_fit->matrix, poly_fit->y, poly_fit->c,
                            poly_fit->cov_matrix, &chisq,
                            poly_fit->workspace) != 0) {
        return PCOMP_STAT_INVALID_FIT;
    }

    for (idx = 0; idx < poly_fit->num_of_coeffs; ++idx) {
        coeffs[idx] = gsl_vector_get(poly_fit->c, idx);
    }

    return PCOMP_STAT_SUCCESS;
}

/***********************************************************************
 * Types and functions used for computing Chebyshev transforms
 */

struct __pcomp_chebyshev_plan_t {
    double* input;
    double* output;
    fftw_plan fftw_plan_ptr;
    pcomp_transform_direction_t dir;
};

pcomp_chebyshev_plan_t*
pcomp_init_chebyshev_plan(size_t num_of_elements,
                          pcomp_transform_direction_t dir)
{
    pcomp_chebyshev_plan_t* plan
        = malloc(sizeof(pcomp_chebyshev_plan_t));

    plan->input = malloc(num_of_elements * sizeof(double));
    plan->output = malloc(num_of_elements * sizeof(double));
    plan->fftw_plan_ptr
        = fftw_plan_r2r_1d(num_of_elements, plan->input, plan->output,
                           FFTW_REDFT00, FFTW_ESTIMATE);
    plan->dir = dir;

    return plan;
}

void pcomp_free_chebyshev_plan(pcomp_chebyshev_plan_t* plan)
{
    if (plan == NULL)
        return;

    if (plan->fftw_plan_ptr != NULL)
        fftw_destroy_plan(plan->fftw_plan_ptr);

    free(plan);
}

/***********************************************************************
 * Types and functions used for applying the combined
 * fitting/Chebyshev transforms
 */

struct __pcomp_polycomp_data_t {
    pcomp_poly_fit_data_t* poly_fit;
    pcomp_chebyshev_plan_t* chebyshev_plan;
    pcomp_chebyshev_plan_t* inv_plan;
};

pcomp_polycomp_data_t* pcomp_init_polycomp_data(size_t num_of_samples,
                                                size_t num_of_coeffs)
{
    pcomp_polycomp_data_t* params
        = malloc(sizeof(pcomp_polycomp_data_t));

    params->poly_fit
        = pcomp_init_poly_fit(num_of_samples, num_of_coeffs);
    params->chebyshev_plan
        = pcomp_init_chebyshev_plan(num_of_samples, PCOMP_TD_DIRECT);
    params->inv_plan
        = pcomp_init_chebyshev_plan(num_of_samples, PCOMP_TD_INVERSE);

    return params;
}

void pcomp_free_polycomp_data(pcomp_polycomp_data_t* params)
{
    if (params == NULL)
        return;

    pcomp_free_poly_fit(params->poly_fit);
    pcomp_free_chebyshev_plan(params->chebyshev_plan);
    pcomp_free_chebyshev_plan(params->inv_plan);

    free(params);
}

/***********************************************************************
 * Evaluate the value of a polynomial at a point using Horner's formula
 */

static double eval_poly(double* coeffs, size_t num_of_coeffs, double x)
{
    if (coeffs == NULL)
        abort();

    if (num_of_coeffs >= 1) {
        size_t idx = num_of_coeffs - 1;
        double result = coeffs[idx];

        if (num_of_coeffs == 1)
            return result;

        --idx;
        while (1) {
            result = result * x + coeffs[idx];

            if (idx > 0)
                --idx;
            else
                break;
        }

        return result;
    }
    else
        return 0.0;
}

/***********************************************************************
 * Polynomial compression routines
 */

int pcomp_compress_poly_float(pcomp_poly_chunk_t** output_buf,
                              size_t* num_of_chunks,
                              const float* input_buf, size_t input_size,
                              const pcomp_poly_parameters* params)
{
    if (output_buf == NULL || num_of_chunks == NULL || input_buf == NULL
        || params == NULL) {
        abort();
    }

    return PCOMP_STAT_SUCCESS;
}
