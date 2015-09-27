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
#include <string.h>

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

struct __pcomp_chebyshev_t {
    double* input;
    double* output;
    size_t num_of_elements;
    fftw_plan fftw_plan_ptr;
    pcomp_transform_direction_t dir;
};

pcomp_chebyshev_t* pcomp_init_chebyshev(size_t num_of_elements,
                                        pcomp_transform_direction_t dir)
{
    pcomp_chebyshev_t* chebyshev = malloc(sizeof(pcomp_chebyshev_t));

    chebyshev->input = malloc(num_of_elements * sizeof(double));
    chebyshev->output = malloc(num_of_elements * sizeof(double));
    chebyshev->num_of_elements = num_of_elements;
    chebyshev->fftw_plan_ptr = fftw_plan_r2r_1d(
        num_of_elements, chebyshev->input, chebyshev->output,
        FFTW_REDFT00, FFTW_ESTIMATE);
    chebyshev->dir = dir;

    return chebyshev;
}

void pcomp_free_chebyshev(pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        return;

    if (plan->fftw_plan_ptr != NULL)
        fftw_destroy_plan(plan->fftw_plan_ptr);

    free(plan);
}

static double chebyshev_normalization(pcomp_transform_direction_t dir,
                                      size_t num_of_elements)
{
    if (dir == PCOMP_TD_DIRECT)
        return 1.0 / (((double)num_of_elements) - 1.0);
    else
        return 0.5;
}

int pcomp_run_chebyshev(pcomp_chebyshev_t* plan,
                        pcomp_transform_direction_t dir, double* output,
                        const double* input)
{
    double norm;
    size_t idx;

    if (plan == NULL)
        abort();

    if (input != NULL) {
        for (idx = 0; idx < plan->num_of_elements; ++idx) {
            plan->input[idx] = input[idx];
        }
    }

    fftw_execute(plan->fftw_plan_ptr);
    norm = chebyshev_normalization(dir, plan->num_of_elements);

    for (idx = 0; idx < plan->num_of_elements; ++idx) {
        plan->output[idx] *= norm;
    }

    if (output != NULL && output != plan->output) {
        for (idx = 0; idx < plan->num_of_elements; ++idx) {
            output[idx] = plan->output[idx];
        }
    }

    return PCOMP_STAT_SUCCESS;
}

/***********************************************************************
 * Types and functions used for applying the combined
 * fitting/Chebyshev transforms
 */

struct __pcomp_polycomp_t {
    pcomp_poly_fit_data_t* poly_fit;
    pcomp_chebyshev_t* chebyshev;
    pcomp_chebyshev_t* inv_chebyshev;
    double max_allowable_error;
    pcomp_polycomp_algorithm_t algorithm;
};

pcomp_polycomp_t*
pcomp_init_polycomp(size_t num_of_samples, size_t num_of_coeffs,
                    double max_allowable_error,
                    pcomp_polycomp_algorithm_t algorithm)
{
    pcomp_polycomp_t* params = malloc(sizeof(pcomp_polycomp_t));

    params->poly_fit
        = pcomp_init_poly_fit(num_of_samples, num_of_coeffs);
    params->chebyshev
        = pcomp_init_chebyshev(num_of_samples, PCOMP_TD_DIRECT);
    params->inv_chebyshev
        = pcomp_init_chebyshev(num_of_samples, PCOMP_TD_INVERSE);
    params->max_allowable_error = max_allowable_error;
    params->algorithm = algorithm;

    return params;
}

void pcomp_free_polycomp(pcomp_polycomp_t* params)
{
    if (params == NULL)
        return;

    pcomp_free_poly_fit(params->poly_fit);
    pcomp_free_chebyshev(params->chebyshev);
    pcomp_free_chebyshev(params->inv_chebyshev);

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
 * Remove sudden jumps in the data by applying an offset equal to +/-
 * "period".
 */

void pcomp_straighten(double* output, const double* input,
                      size_t num_of_elements, double period)
{
    size_t idx;
    double half_period = period * 0.5;
    double offset = 0.0;

    if (input == NULL || output == NULL)
        abort();

    if (period > 0) {
        for (idx = 1; idx < num_of_elements; ++idx) {
            double diff_with_previous = input[idx] - input[idx - 1];
            if (diff_with_previous > half_period)
                offset -= period;
            else
                offset += period;

            output[idx] = input[idx] + offset;
        }
    }
    else {
        for (idx = 0; idx < num_of_elements; ++idx)
            output[idx] = input[idx];
    }
}

/***********************************************************************
 * Chunk initialization/destruction
 */

pcomp_polycomp_chunk_t* pcomp_init_chunk(size_t num_of_elements)
{
    pcomp_polycomp_chunk_t* chunk
        = malloc(sizeof(pcomp_polycomp_chunk_t));

    chunk->num_of_elements = num_of_elements;

    chunk->is_compressed = 0;
    chunk->uncompressed
        = malloc(sizeof(double) * sizeof(chunk->num_of_elements));

    chunk->num_of_poly_coeffs = 0;
    chunk->poly_coeffs = NULL;

    chunk->num_of_cheby_coeffs = 0;
    chunk->cheby_coeffs = NULL;

    return chunk;
}

void pcomp_free_chunk(pcomp_polycomp_chunk_t* chunk)
{
    if (chunk->uncompressed != NULL)
        free(chunk->uncompressed);

    if (chunk->poly_coeffs != NULL)
        free(chunk->poly_coeffs);

    if (chunk->cheby_coeffs != NULL)
        free(chunk->cheby_coeffs);
}

/***********************************************************************
 * This routine implements the core of the "polynomial compression": a
 * polynomial fitting plus a Chebyshev transform on the residuals. It
 * is a static function because we're going to wrap a nicer interface
 * around it.
 */

static int polyfit_and_chebyshev(pcomp_polycomp_t* params,
                                 double* coeffs, const double* input,
                                 double* max_residual)
{
    size_t idx;
    int status;
    double running_max;

    status = pcomp_run_poly_fit(params->poly_fit, coeffs, input);
    if (status != PCOMP_STAT_SUCCESS)
        return status;

    for (idx = 0; idx < params->poly_fit->num_of_points; ++idx) {
        double abs_residual;

        params->chebyshev->input[idx]
            = input[idx]
              - eval_poly(coeffs, params->poly_fit->num_of_coeffs,
                          idx + 1.0);

        abs_residual = fabs(params->chebyshev->input[idx]);
        if (abs_residual > running_max)
            running_max = abs_residual;
    }

    if (max_residual != NULL)
        *max_residual = running_max;

    if (params->algorithm != PCOMP_ALG_NO_CHEBYSHEV) {
        status = pcomp_run_chebyshev(
            params->chebyshev, params->chebyshev->dir, NULL, NULL);
    }

    return status;
}

static size_t trunc_chebyshev(pcomp_chebyshev_t* chebyshev,
                              pcomp_chebyshev_t* inv_chebyshev,
                              double max_allowable_error,
                              double* max_error)
{
    size_t idx;
    size_t result = 0;
    double err;

    if (chebyshev == NULL || inv_chebyshev == NULL)
        abort();

    if (max_allowable_error <= 0.0)
        return 0;

    /* Start by setting all the coefficients to zero */
    for (idx = 0; idx < chebyshev->num_of_elements; ++idx) {
        inv_chebyshev->input[idx] = 0.0;
    }

    /* Add the coefficients one by one until the error is below
     * "max_allowable_error" */
    while (result < chebyshev->num_of_elements) {
        err = 0.0;

        pcomp_run_chebyshev(inv_chebyshev, inv_chebyshev->dir, NULL,
                            NULL);

        for (idx = 0; idx < chebyshev->num_of_elements; ++idx) {
            double cur_err = fabs(chebyshev->input[idx]
                                  - inv_chebyshev->output[idx]);
            if (cur_err > err)
                err = cur_err;
        }

        ++result;
        inv_chebyshev->input[result] = chebyshev->output[result];

        if (err < max_allowable_error)
            break;
    }

    if (max_error != NULL)
        *max_error = err;

    if (result >= chebyshev->num_of_elements)
        result = chebyshev->num_of_elements;

    return result;
}

/* This function is used internally by "pcomp_run_polycomp_on_chunk"
 * and "pcomp_decompress_poly_chunk" to make sure there is no memory
 * leak on the chunk passed as argument. */
static void clear_chunk(pcomp_polycomp_chunk_t* chunk)
{
    /* Leave chunk->uncompressed as it is, as it never changes */

    chunk->num_of_poly_coeffs = 0;
    if (chunk->poly_coeffs != NULL)
        free(chunk->poly_coeffs);

    chunk->num_of_cheby_coeffs = 0;
    if (chunk->cheby_coeffs != NULL)
        free(chunk->cheby_coeffs);
}

int pcomp_run_polycomp_on_chunk(pcomp_polycomp_t* params,
                                const double* input,
                                size_t num_of_elements,
                                pcomp_polycomp_chunk_t* chunk,
                                double* max_error)
{
    double* coeffs;
    size_t cheby_coeffs_to_retain = 0;
    double max_residual;
    int apply_chebyshev = 1;

    if (chunk == NULL || input == NULL || params == NULL
        || params->poly_fit == NULL || params->chebyshev == NULL
        || params->inv_chebyshev == NULL)
        abort();

    clear_chunk(chunk);

    if (num_of_elements != params->poly_fit->num_of_points)
        return PCOMP_STAT_INVALID_BUFFER;

    coeffs = malloc(sizeof(double) * params->poly_fit->num_of_coeffs);
    polyfit_and_chebyshev(params, coeffs, input, &max_residual);
    apply_chebyshev = (max_residual >= params->max_allowable_error)
                      && (params->algorithm != PCOMP_ALG_NO_CHEBYSHEV);

    if (apply_chebyshev) {
        cheby_coeffs_to_retain
            = trunc_chebyshev(params->chebyshev, params->inv_chebyshev,
                              params->max_allowable_error, max_error);

        chunk->is_compressed
            = (cheby_coeffs_to_retain + params->poly_fit->num_of_coeffs)
              < num_of_elements;
    }
    else {
        /* Assume that num_of_elements > deg(p) + 1 */
        chunk->is_compressed = 1;
    }

    chunk->num_of_elements = num_of_elements;
    if (chunk->is_compressed) {
        size_t idx;

        chunk->num_of_poly_coeffs = params->poly_fit->num_of_coeffs;
        chunk->poly_coeffs = coeffs;
        if (apply_chebyshev) {
            chunk->num_of_cheby_coeffs = cheby_coeffs_to_retain;
            chunk->cheby_coeffs
                = malloc(sizeof(double) * cheby_coeffs_to_retain);
            for (idx = 0; idx < cheby_coeffs_to_retain; ++idx) {
                chunk->cheby_coeffs[idx]
                    = params->chebyshev->output[idx];
            }
        }
        else {
            chunk->num_of_cheby_coeffs = 0;
            chunk->cheby_coeffs = NULL;
        }
    }
    else {
        size_t idx;

        free(coeffs);

        chunk->uncompressed = malloc(sizeof(double) * num_of_elements);
        for (idx = 0; idx < num_of_elements; ++idx)
            chunk->uncompressed[idx] = input[idx];

        if (max_error != NULL)
            *max_error = 0.0;
    }

    return PCOMP_STAT_SUCCESS;
}

int pcomp_decompress_polycomp_chunk(double* output,
                                    const pcomp_polycomp_chunk_t* chunk,
                                    pcomp_chebyshev_t* inv_chebyshev)
{
    if (output == NULL || chunk == NULL || inv_chebyshev == NULL)
        abort();

    if (chunk->is_compressed) {
        size_t idx;

        /* Compute the values of the polynomial at the points 1, 2, ...
         */
        for (idx = 0; idx < chunk->num_of_elements; ++idx) {
            output[idx] = eval_poly(chunk->poly_coeffs,
                                    chunk->num_of_poly_coeffs, idx + 1);
        }

        /* If present, add the contribution of the Chebyshev transform
         */
        if (chunk->num_of_cheby_coeffs > 0) {
            if (chunk->cheby_coeffs == NULL)
                abort();

            for (idx = 0; idx < chunk->num_of_cheby_coeffs; ++idx) {
                inv_chebyshev->input[idx] = chunk->cheby_coeffs[idx];
            }
            for (idx = chunk->num_of_cheby_coeffs;
                 idx < chunk->num_of_elements; ++idx) {
                inv_chebyshev->input[idx] = 0.0;
            }

            pcomp_run_chebyshev(inv_chebyshev, inv_chebyshev->dir, NULL,
                                NULL);

            for (idx = 0; idx < chunk->num_of_elements; ++idx) {
                output[idx] += inv_chebyshev->output[idx];
            }
        }
    }
    else {
        memcpy(output, chunk->uncompressed,
               sizeof(chunk->uncompressed[0]) * chunk->num_of_elements);
    }

    return PCOMP_STAT_SUCCESS;
}

/***********************************************************************
 * Polynomial compression routines
 */

int pcomp_compress_poly_float(pcomp_polycomp_chunk_t** output_buf,
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
