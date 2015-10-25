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
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_block_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit.h>

#include <fftw3.h>

/**********************************************************************/

/** \defgroup poly Polynomial compression
 *
 * ### The algorithm and its applicability
 *
 * Polynomial compression relies on a simple idea, that is to divide
 * the input data stream into subsets of consecutive samples (called
 * "chunks"), and to approximate each chunk by means of a polynomial.
 * Such compression is inherently lossy, as the residuals of the
 * fitting procedure are usually discarded.
 *
 * This idea has been widely applied in the literature. Libpolycomp
 * implements a significant improvement over it, because in case the
 * polynomial used for the fitting is not good enough (i.e., the
 * magnitude of the residuals is larger than the desired error
 * threshold) it saves selected components of the Chebyshev transform
 * of the residuals together with the polynomial.
 *
 * It is possible to avoid the usage of Chebyshev transforms. In this
 * case, if no polynomial of the desired degree are able to fit the
 * data with the given error threshold, the data for that chunk is
 * saved uncompressed.
 *
 * This works quite well for smooth data series, where changes between
 * consecutive samples are well described by slowly varying continuous
 * functions. It is not suitable if the signal contains noise, unless
 * this noise is significantly smaller than the signal and than the
 * error threshold.
 *
 * ### Implementation details
 *
 * The Libpolycomp library implements two families of functions which
 * implement polynomial compression:
 *
 * - Low-level functions are useful to work with each chunk
 *   separately, or when one needs to optimize the compression
 *   parameters.
 *
 * - High-level functions are used for compressing/decompressing long
 *   streams when no particular needs are required.
 *
 * Both the low-level and high-level functions use the \ref
 * pcomp_polycomp_t structure to determine which parameters to use for
 * the compression. The functions that allow to allocate/free/manage
 * this structure are the following:
 *
 * - \ref pcomp_init_polycomp
 * - \ref pcomp_free_polycomp
 * - \ref pcomp_polycomp_samples_per_chunk
 * - \ref pcomp_polycomp_num_of_poly_coeffs
 * - \ref pcomp_polycomp_max_error
 * - \ref pcomp_polycomp_algorithm
 * - \ref pcomp_polycomp_period
 * - \ref pcomp_polycomp_set_period
 *
 * Moreover, the function implements a number of facilities to compute
 * polynomial fits and Chebyshev transforms of discrete data.
 * Libpolycomp uses the GNU Scientific Library (GSL) to implement the
 * former, and the FFTW library for the latter.
 *
 * ### Low-level functions
 *
 * ### High-level functions
 *
 * ### Mathematical routines
 *
 * The following set of Libpolycomp routines allows to compute the
 * least squares best fit between a set of floating-point numbers
 * \f$x_i\f$ (with \f$i=1\ldots N\f$) and a polynomial, through the
 * points \f$(i, x_i)\f$:
 *
 * - \ref pcomp_init_poly_fit and \ref pcomp_free_poly_fit are used to
 *   initialize and free a \ref pcomp_poly_fit_data_t structure. Such
 *   structure can be used to apply repeatedly a fitting process to
 *   different sets of data, provided that the number of points
 *   \f$N\f$ and the degree of the fitting polynomial stay constant.
 *
 * - \ref pcomp_poly_fit_num_of_samples and \ref
 *   pcomp_poly_fit_num_of_coeffs are used to retrieve the fields of
 *   an instance of the opaque structure \ref pcomp_poly_fit_data_t.
 *
 * - \ref pcomp_run_poly_fit calculates the least-squares best fit
 *   between a set of points and a polynomial.
 *
 * The following set of routines compute the Chebyshev transform of a
 * set of floating-point numbers. They are a tiny wrapper around
 * analogous functions of the FFTW library, with the main purpose of
 * using the correct normalization constants in the forward and
 * inverse transforms:
 *
 * - \ref pcomp_init_chebyshev and \ref pcomp_free_chebyshev are used
 *   to allocate and free a \ref pcomp_chebyshev_t structure, which
 *   describes how to compute a Chebyshev transform for a set of
 *   \f$N\f$ numbers.
 *
 * - \ref pcomp_chebyshev_num_of_samples and \ref
 *   pcomp_chebyshev_direction are used to query the parameters of an
 *   instance of the opaque structure \ref pcomp_chebyshev_t.
 *
 * - \ref pcomp_run_chebyshev applies the Chebyshev transform (either in
 *   the forward or inverse direction) to a dataset.
 */

/**********************************************************************/

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
    size_t num_of_samples;
    size_t num_of_coeffs;
    gsl_multifit_linear_workspace* workspace;
    gsl_matrix* matrix;
    gsl_vector* y;
    gsl_vector* c;
    gsl_matrix* cov_matrix;
};

/** \ingroup poly
 *
 * \brief Allocate a new instance of the \ref pcomp_poly_fit_data_t
 * structure on the heap
 *
 * \param[in] num_of_samples Number of floating-point numbers that
 * must fit the polynomial
 *
 * \param[in] num_of_coeffs Number of coefficients of the
 * least-squares fitting polynomial \f$p(x)\f$. This is equal to
 * \f$\deg p(x) + 1\f$, where \f$\deg p(x)\f$ is the degree of the
 * polynomial. Thus, for a parabolic polynomial of the form \f$p(x) =
 * a x^2 + b x + c\f$, \a num_of_coeffs = 3.
 *
 * \returns A newly created instance of \ref pcomp_poly_fit_data_t
 * structure. This must be freed using \ref pcomp_free_poly_fit, once
 * it is no longer used.
 */
pcomp_poly_fit_data_t* pcomp_init_poly_fit(size_t num_of_samples,
                                           size_t num_of_coeffs)
{
    size_t i, j;

    pcomp_poly_fit_data_t* poly_fit
        = malloc(sizeof(pcomp_poly_fit_data_t));
    if (poly_fit == NULL)
        abort();

    poly_fit->num_of_samples = num_of_samples;
    poly_fit->num_of_coeffs = num_of_coeffs;
    poly_fit->workspace
        = gsl_multifit_linear_alloc(num_of_samples, num_of_coeffs);
    poly_fit->matrix = gsl_matrix_alloc(num_of_samples, num_of_coeffs);
    poly_fit->y = gsl_vector_alloc(num_of_samples);
    poly_fit->c = gsl_vector_alloc(num_of_coeffs);
    poly_fit->cov_matrix
        = gsl_matrix_alloc(num_of_coeffs, num_of_coeffs);

    for (i = 0; i < num_of_samples; ++i) {
        for (j = 0; j < num_of_coeffs; ++j) {
            gsl_matrix_set(poly_fit->matrix, i, j,
                           integer_power(i + 1, j));
        }
    }

    return poly_fit;
}

/** \ingroup poly
 *
 * \brief Free an instance of the \ref pcomp_poly_fit_data_t that has
 * been allocated via a call to \ref pcomp_init_poly_fit.
 *
 * \param[in] poly_fit Pointer to the structure to be freed
 */
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

/** \ingroup poly
 *
 * \brief Return the number of samples to be used in a polynomial fit
 *
 * \param[in] poly_fit Pointer to the structure detailing the fit
 *
 * \returns The number of samples that should be passed to a call to
 * \ref pcomp_run_poly_fit.
 */
size_t
pcomp_poly_fit_num_of_samples(const pcomp_poly_fit_data_t* poly_fit)
{
    if (poly_fit == NULL)
        abort();

    return poly_fit->num_of_samples;
}

/** \ingroup poly
 *
 * \brief Return the number of coefficients of the least-squares
 * fitting polynomial
 *
 * \param[in] poly_fit Pointer to the structure detailing the fit
 *
 * \returns The number of coefficients for the fitting polynomial (one
 * plus the polynomial degree)
 */
size_t
pcomp_poly_fit_num_of_coeffs(const pcomp_poly_fit_data_t* poly_fit)
{
    if (poly_fit == NULL)
        abort();

    return poly_fit->num_of_coeffs;
}

/** \ingroup poly
 *
 * \brief Calculates a polynomial least-squares fit.
 *
 * Compute a least-squares fit between the numbers \f$x_i\f$ (with
 * \f$i = 1 \ldots N\f$) and the polynomial \f$p(x)\f$ through the
 * points \f$(i, x_i)_{i=1}^N\f$. The coefficients of \f$p(x)\f$ are
 * saved in \a coeffs, from the least to the greatest degree.
 *
 * Here is an example of the usage of this function:
 *
 * \code{.c}
 * double points[] = { 1.0, 3.0, 5.0 };
 * double coeffs[2];
 * const size_t num_of_points = sizeof(points) / sizeof(points[0]);
 * const size_t num_of_coeffs = sizeof(coeffs) / sizeof(coeffs[0]);
 * pcomp_poly_fit_data_t* poly_fit;
 *
 * poly_fit = pcomp_init_poly_fit(num_of_points, num_of_coeffs);
 * pcomp_run_poly_fit(poly_fit, coeffs, points);
 * printf("The data are fitted by the polynomial y = %f + %f x\n",
 *        coeffs[0], coeffs[1]);
 * \endcode
 *
 * \param[in] poly_fit Pointer to a \ref pcomp_poly_fit_data_t
 * structure, created using the \ref pcomp_init_poly_fit function.
 *
 * \param[out] coeffs Pointer to an array where the coefficients of
 * the polynomial will be stored on exit. The array must have room for
 * a number of elements greater or equal than the value returned by
 * \ref pcomp_poly_fit_num_of_coeffs.
 *
 * \param[in] points Array of numbers \f$x_i\f$ to use in the fit. The
 * number of elements considered in the fit is equal to the return
 * value of \ref pcomp_poly_fit_num_of_samples.
 *
 * \returns \ref PCOMP_STAT_SUCCESS if the fit was computed
 * successfully, \ref PCOMP_STAT_INVALID_FIT if the data are incorrect
 * (e.g., there are fewer samples than unknowns).
 */
int pcomp_run_poly_fit(pcomp_poly_fit_data_t* poly_fit, double* coeffs,
                       const double* points)
{
    size_t idx;
    double chisq;

    if (poly_fit == NULL || coeffs == NULL || points == NULL)
        abort();

    for (idx = 0; idx < poly_fit->num_of_samples; ++idx) {
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
    size_t num_of_samples;
    fftw_plan fftw_plan_ptr;
    pcomp_transform_direction_t dir;
};

/** \ingroup poly
 *
 * \brief Allocate a new instance of the \ref pcomp_chebyshev_t
 * structure on the heap
 *
 * Despite the fact that this function takes the parameter \a dir, the
 * function which actually computes the Chebyshev transform (\ref
 * pcomp_run_chebyshev) allow to specify the desired direction. The
 * purpose of having \a dir encoded in \ref pcomp_chebyshev_t is that
 * sometimes it is useful to keep it memorized in the structure
 * itself.
 *
 * \param[in] num_of_samples Number of floating-point numbers that
 * will be transformed
 *
 * \param[in] dir Direction of the transform (either forward or
 * backward). This is used to determine the normalization constant of
 * the transform:
 * - If computing a forward transform, the normalization is \f$1 / (N
 *   - 1)\f$, with \f$N\f$ the number of samples.
 * - If computing a backward transform, the normalization is 1.
 *
 * \returns A newly created instance of \ref pcomp_poly_fit_data_t
 * structure. This must be freed using \ref pcomp_free_poly_fit, once
 * it is no longer used.
 */
pcomp_chebyshev_t* pcomp_init_chebyshev(size_t num_of_samples,
                                        pcomp_transform_direction_t dir)
{
    pcomp_chebyshev_t* chebyshev = malloc(sizeof(pcomp_chebyshev_t));
    if (chebyshev == NULL)
        abort();

    chebyshev->input = fftw_alloc_real(num_of_samples);
    chebyshev->output = fftw_alloc_real(num_of_samples);
    chebyshev->num_of_samples = num_of_samples;
    chebyshev->fftw_plan_ptr = fftw_plan_r2r_1d(
        num_of_samples, chebyshev->input, chebyshev->output,
        FFTW_REDFT00, FFTW_ESTIMATE);
    chebyshev->dir = dir;

    return chebyshev;
}

/** \ingroup poly
 *
 * \brief Free the memory allocated by a previous call to \ref
 * pcomp_init_chebyshev.
 *
 * \param[in] plan Pointer to the structure to be freed.
 */
void pcomp_free_chebyshev(pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        return;

    if (plan->input != NULL)
        fftw_free(plan->input);

    if (plan->output != NULL)
        fftw_free(plan->output);

    if (plan->fftw_plan_ptr != NULL)
        fftw_destroy_plan(plan->fftw_plan_ptr);

    free(plan);
}

/** \ingroup poly
 *
 * \brief Return the number of samples in a Chebyshev transform
 *
 * \param[in] plan Pointer to the Chebyshev plan.
 *
 * \returns The number of elements that are used in the Chebyshev
 * transform specified by \a plan.
 */
size_t pcomp_chebyshev_num_of_samples(const pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        abort();

    return plan->num_of_samples;
}

/** \ingroup poly
 *
 * \brief Return the direction of a Chebyshev transform
 *
 * \param[in] plan Pointer to the Chebyshev plan.
 *
 * \returns A \ref pcomp_transform_direction_t value specifying the
 * normalization used for the Chebyshev transform specified by \a
 * plan.
 */
pcomp_transform_direction_t
pcomp_chebyshev_direction(const pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        abort();

    return plan->dir;
}

static double chebyshev_normalization(pcomp_transform_direction_t dir,
                                      size_t num_of_samples)
{
    if (dir == PCOMP_TD_DIRECT)
        return 1.0 / (((double)num_of_samples) - 1.0);
    else
        return 0.5;
}

/** \ingroup poly
 *
 * \brief Compute a forward/backward Chebyshev discrete transform
 *
 * \param[in] plan Pointer to a Chebyshev plan created by \ref
 * pcomp_init_chebyshev
 *
 * \param[in] dir Direction of the transform. This parameter overrides
 * the internal direction of \a plan (returned by \ref
 * pcomp_chebyshev_direction).
 *
 * \param[out] output Pointer to an array of \c double values that will
 * contain the Chebyshev transform of \a input. It must have room for
 * a number of elements at least equal to the return value of \ref
 * pcomp_num_of_samples.
 *
 * \param[in] input Array of \c double values to be transformed. The
 * function will use the first N elements, where N is the return value
 * of \ref pcomp_num_of_samples.
 *
 * \returns \ref PCOMP_STAT_SUCCESS when successful.
 */
int pcomp_run_chebyshev(pcomp_chebyshev_t* plan,
                        pcomp_transform_direction_t dir, double* output,
                        const double* input)
{
    double norm;
    size_t idx;

    if (plan == NULL)
        abort();

    if (input != NULL) {
        for (idx = 0; idx < plan->num_of_samples; ++idx) {
            plan->input[idx] = input[idx];
        }
    }

    fftw_execute(plan->fftw_plan_ptr);
    norm = chebyshev_normalization(dir, plan->num_of_samples);

    for (idx = 0; idx < plan->num_of_samples; ++idx) {
        plan->output[idx] *= norm;
    }

    if (output != NULL && output != plan->output) {
        for (idx = 0; idx < plan->num_of_samples; ++idx) {
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
    size_t samples_per_chunk;
    pcomp_poly_fit_data_t* poly_fit;
    pcomp_chebyshev_t* chebyshev;
    pcomp_chebyshev_t* inv_chebyshev;
    double max_allowable_error;
    pcomp_polycomp_algorithm_t algorithm;
    double period;
};

/** \ingroup poly
 *
 * \brief Allocate space for a \ref pcomp_polycomp_t structure
 *
 * \param[in] samples_per_chunk Number of samples in each chunk
 *
 * \param[in] num_of_coeffs Number of polynomial coefficients to use
 *
 * \param[in] max_allowable_error Upper bound for the compression
 * error (positive value)
 *
 * \param[in] algorithm Kind of compression algorithm to use
 *
 * \returns A pointer to the newly allocate \ref pcomp_polycomp_t
 * structure. This must be freed using \ref pcomp_free_polycomp, once
 * it is no longer used.
 */
pcomp_polycomp_t*
pcomp_init_polycomp(pcomp_chunk_size_t samples_per_chunk,
                    pcomp_poly_size_t num_of_coeffs,
                    double max_allowable_error,
                    pcomp_polycomp_algorithm_t algorithm)
{
    pcomp_polycomp_t* params = malloc(sizeof(pcomp_polycomp_t));
    if (params == NULL)
	abort();

    params->samples_per_chunk = samples_per_chunk;
    params->poly_fit
        = pcomp_init_poly_fit(samples_per_chunk, num_of_coeffs);
    params->chebyshev
        = pcomp_init_chebyshev(samples_per_chunk, PCOMP_TD_DIRECT);
    params->inv_chebyshev
        = pcomp_init_chebyshev(samples_per_chunk, PCOMP_TD_INVERSE);
    params->max_allowable_error = max_allowable_error;
    params->algorithm = algorithm;
    params->period = 0.0;

    return params;
}

/** \ingroup poly
 *
 * \brief Free the memory allocated by \ref pcomp_init_polycomp for a
 * \ref pcomp_polycomp_t structure.
 *
 * \param[in] params Pointer to the structure to be freed
 */
void pcomp_free_polycomp(pcomp_polycomp_t* params)
{
    if (params == NULL)
        return;

    pcomp_free_poly_fit(params->poly_fit);
    pcomp_free_chebyshev(params->chebyshev);
    pcomp_free_chebyshev(params->inv_chebyshev);

    free(params);
}

/** \ingroup poly
 *
 * \brief Return the number of samples per chunk
 *
 * This function returns the size of each chunk but the last one in
 * the input data for a polynomial compression. Such chunks contain a
 * set of consecutive values in the input array passed to routines as
 * \ref pcomp_compress_polycomp.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The number of samples in each chunk.
 */
pcomp_chunk_size_t
pcomp_polycomp_samples_per_chunk(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->samples_per_chunk;
}

/** \ingroup poly
 *
 * \brief Return the number of coefficients for the fitting polynomial
 * used in the polynomial compression.
 *
 * The return value has the same meaning as the value returned by the
 * \ref pcomp_poly_fit_num_of_coeffs.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The number of coefficients of the fitting polynomial.
 */
pcomp_poly_size_t
pcomp_polycomp_num_of_poly_coeffs(const pcomp_polycomp_t* params)
{
    if (params == NULL || params->poly_fit == NULL)
        abort();

    return params->poly_fit->num_of_coeffs;
}

/** \ingroup poly
 *
 * \brief Return the upper bound on the error of the polynomial
 *compression.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The maximum allowable error for the polynomial compression.
 */
double pcomp_polycomp_max_error(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->max_allowable_error;
}

/** \ingroup poly
 *
 * \brief Return the kind of algorithm used for a polynomial
 *compression.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The algorithm to be used by the compressor.
 */
pcomp_polycomp_algorithm_t
pcomp_polycomp_algorithm(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->algorithm;
}

/** \ingroup poly
 *
 * \brief Return the period of the input data, or a number
 * less than or equal to 0 if the data have no periodicity.
 *
 * See also \ref pcomp_polycomp_set_period.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The periodicity. If zero or negative, no periodicity is
 * assumed in the data to be compressed.
 */
double pcomp_polycomp_period(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->period;
}

/** \ingroup poly
 *
 * \brief Set the periodicity of the data to be compressed
 *
 * If \a period is a value greater than zero, this is assumed to be
 * the periodicity of the input data: the value \a x is therefore
 * assumed equivalent to \a x + \a period and to \a x - \a period. It
 * is typically a multiple of Pi = 3.14159...
 *
 * The polynomial compressor can improve the compression ratio for
 * data if they have some form of periodicity.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \param[in] period The periodicity of the data, or a zero/negative
 * value if no periodicity should be assumed by the compressor.
 */
void pcomp_polycomp_set_period(pcomp_polycomp_t* params, double period)
{
    if (params == NULL)
        abort();

    params->period = period;
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

/***********************************************************************/

/** \ingroup poly
 *
 * \brief Remove sudden jumps from \a input
 *
 * Assuming that the data in the array \a input have a periodicity
 * equal to \a period, the function copies them to \a output while
 * applying a positive/negative offset equal to a multiple of \a
 * period.
 *
 * It is ok for \a input and \a output to point to the same memory
 * location.
 *
 * \param[out] output Pointer to the array that will contain the
 * result. It must have room for at least \a num_of_samples values.
 *
 * \param[in] input Array of \a num_of_samples values to process.
 *
 * \param[in] num_of_samples Number of samples to process in \a input
 *
 * \param[in] period Periodicity of the data. If less or equal to
 * zero, \a input is copied verbatim to \a output.
 */

void pcomp_straighten(double* output, const double* input,
                      size_t num_of_samples, double period)
{
    size_t idx;

    if (input == NULL || output == NULL)
        abort();

    if (period > 0) {
        double half_period = period * 0.5;
        double offset = 0.0;

        output[0] = input[0];

        for (idx = 1; idx < num_of_samples; ++idx) {
            double diff_with_previous = input[idx] - input[idx - 1];
            if (diff_with_previous > half_period)
                offset -= period;
            else if (diff_with_previous < -half_period)
                offset += period;

            output[idx] = input[idx] + offset;
        }
    }
    else {
        for (idx = 0; idx < num_of_samples; ++idx)
            output[idx] = input[idx];
    }
}

/***********************************************************************
 * Chunk initialization/destruction
 */

/* Information about a chunk of data compressed using the polynomial
 * compression */
struct __pcomp_polycomp_chunk_t {
    /* Number of samples in this chunk */
    size_t num_of_samples;

    /* Is this chunk compressed using polynomial/Chebyshev
     * coefficients? */
    int is_compressed;
    /* If the chunk is not compressed (is_compressed == 0), this
     * points to a buffer which holds "num_of_samples" uncompressed
     * samples */
    double* uncompressed;

    /* Polynomial coefficients, from the lowest-order to the
     * highest-order */
    size_t num_of_poly_coeffs;
    double* poly_coeffs;

    /* Chebyshev coefficients */
    uint8_t* cheby_mask;
    size_t num_of_cheby_coeffs; /* This is always less than
                                 * num_of_samples, as the Chebyshev
                                 * series is truncated. */
    double* cheby_coeffs;
};

pcomp_polycomp_chunk_t*
pcomp_init_chunk(pcomp_chunk_size_t num_of_samples)
{
    pcomp_polycomp_chunk_t* chunk
        = malloc(sizeof(pcomp_polycomp_chunk_t));
    if (chunk == NULL)
	abort();

    chunk->num_of_samples = num_of_samples;

    chunk->is_compressed = 0;
    chunk->uncompressed
        = malloc(sizeof(double) * sizeof(chunk->num_of_samples));

    chunk->num_of_poly_coeffs = 0;
    chunk->poly_coeffs = NULL;

    chunk->num_of_cheby_coeffs = 0;
    chunk->cheby_coeffs = NULL;
    chunk->cheby_mask = NULL;

    return chunk;
}

pcomp_polycomp_chunk_t*
pcomp_init_uncompressed_chunk(pcomp_chunk_size_t num_of_samples,
                              const double* samples)
{
    pcomp_polycomp_chunk_t* chunk
        = malloc(sizeof(pcomp_polycomp_chunk_t));
    const size_t num_of_bytes = sizeof(double) * num_of_samples;
    if (chunk == NULL)
	abort();

    chunk->num_of_samples = num_of_samples;

    chunk->is_compressed = 0;
    chunk->uncompressed = malloc(num_of_bytes);
    if (chunk->uncompressed == NULL)
	abort();
    memcpy(chunk->uncompressed, samples, num_of_bytes);

    chunk->num_of_poly_coeffs = 0;
    chunk->poly_coeffs = NULL;

    chunk->num_of_cheby_coeffs = 0;
    chunk->cheby_coeffs = NULL;
    chunk->cheby_mask = NULL;

    return chunk;
}

pcomp_polycomp_chunk_t* pcomp_init_compressed_chunk(
    pcomp_chunk_size_t num_of_samples,
    pcomp_poly_size_t num_of_poly_coeffs, const double* poly_coeffs,
    pcomp_chunk_size_t num_of_cheby_coeffs, const uint8_t* cheby_mask,
    const double* cheby_coeffs)
{
    size_t size;
    pcomp_polycomp_chunk_t* chunk;

    if (num_of_samples == 0 || poly_coeffs == NULL)
        abort();

    chunk = malloc(sizeof(pcomp_polycomp_chunk_t));
    if (chunk == NULL)
	abort();

    chunk->num_of_samples = num_of_samples;
    chunk->is_compressed = 1;
    chunk->uncompressed = NULL;

    chunk->num_of_poly_coeffs = num_of_poly_coeffs;
    size = num_of_poly_coeffs * sizeof(double);
    chunk->poly_coeffs = malloc(size);
    if (chunk->poly_coeffs == NULL)
	abort();
    memcpy(chunk->poly_coeffs, poly_coeffs, size);

    chunk->num_of_cheby_coeffs = num_of_cheby_coeffs;
    if (num_of_cheby_coeffs > 0) {
        size = num_of_cheby_coeffs * sizeof(double);
        chunk->cheby_coeffs = malloc(size);
	if (chunk->cheby_coeffs == NULL)
	    abort();
        memcpy(chunk->cheby_coeffs, cheby_coeffs, size);

        size = pcomp_chunk_cheby_mask_size(num_of_samples)
               * sizeof(uint8_t);
        chunk->cheby_mask = malloc(size);
	if (chunk->cheby_mask == NULL)
	    abort();
        memcpy(chunk->cheby_mask, cheby_mask, size);
    }
    else {
        chunk->cheby_mask = NULL;
        chunk->cheby_coeffs = NULL;
    }

    return chunk;
}

void pcomp_free_chunk(pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        return;

    if (chunk->uncompressed != NULL)
        free(chunk->uncompressed);

    if (chunk->poly_coeffs != NULL)
        free(chunk->poly_coeffs);

    if (chunk->cheby_coeffs != NULL)
        free(chunk->cheby_coeffs);

    if (chunk->cheby_mask != NULL)
        free(chunk->cheby_mask);

    free(chunk);
}

pcomp_chunk_size_t
pcomp_chunk_num_of_samples(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->num_of_samples;
}

size_t pcomp_chunk_num_of_bytes(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (chunk->is_compressed) {
        /* The size is calculated as follows:
         * - the "compressed" flag (int8_t)
         * - the number of samples (pcomp_chunk_size_t)
         * - the number N of polynomial coefficients (pcomp_poly_size_t)
         * - the size of the Chebyshev mask
         * - the number M of Chebyshev coefficients (pcomp_chunk_size_t)
         * - Nx8 bytes for the polynomial
         * - Mx8 bytes for the Chebyshev coefficients
         */
        return sizeof(int8_t) + sizeof(pcomp_chunk_size_t)
               + sizeof(pcomp_poly_size_t)
               + pcomp_chunk_cheby_mask_size(chunk->num_of_samples)
               + sizeof(pcomp_chunk_size_t)
               + (chunk->num_of_poly_coeffs
                  + chunk->num_of_cheby_coeffs) * sizeof(double);
    }
    else {
        /* The size is calculated as follows:
         * - 1 byte for the "uncompressed" flag
         * - 4 bytes for the number of samples
         * - Nx8 bytes for the samples
         */
        return sizeof(int8_t) + sizeof(pcomp_chunk_size_t)
               + chunk->num_of_samples * sizeof(double);
    }
}

int pcomp_chunk_is_compressed(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->is_compressed;
}

const double*
pcomp_chunk_uncompressed_data(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (chunk->is_compressed)
        return NULL;

    return chunk->uncompressed;
}

pcomp_poly_size_t
pcomp_chunk_num_of_poly_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->num_of_poly_coeffs;
}

const double*
pcomp_chunk_poly_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (chunk->num_of_poly_coeffs == 0)
        return NULL;

    return chunk->poly_coeffs;
}

pcomp_chunk_size_t
pcomp_chunk_num_of_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->num_of_cheby_coeffs;
}

const double*
pcomp_chunk_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (chunk->num_of_cheby_coeffs == 0)
        return NULL;

    return chunk->cheby_coeffs;
}

/* Return the number of bytes needed to store the mask of Chebyshev
 * coefficients */
size_t pcomp_chunk_cheby_mask_size(pcomp_chunk_size_t chunk_size)
{
    return chunk_size / CHAR_BIT
           + ((chunk_size % CHAR_BIT) > 0 ? 1 : 0);
}

const uint8_t*
pcomp_chunk_cheby_mask(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->cheby_mask;
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
    double running_max = -1.0; /* Negative stands for "uninitialized" */

    status = pcomp_run_poly_fit(params->poly_fit, coeffs, input);
    if (status != PCOMP_STAT_SUCCESS)
        return status;

    for (idx = 0; idx < params->samples_per_chunk; ++idx) {
        double abs_residual;

        params->chebyshev->input[idx]
            = input[idx]
              - eval_poly(coeffs, params->poly_fit->num_of_coeffs,
                          idx + 1.0);

        abs_residual = fabs(params->chebyshev->input[idx]);
        if (abs_residual > running_max || running_max < 0.0)
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

/***********************************************************************
 * Sort the array "positions" according to the absolute values of
 * "coeffs", in *descending* order. The function uses the merge sort
 * algorithm. */

static void sort_positions(pcomp_chunk_size_t positions[],
                           const double coeffs[], size_t num)
{
    size_t front, back;
    double pivot;
    pcomp_chunk_size_t temp;

    if (num < 2)
        return;

    pivot = fabs(coeffs[positions[num / 2]]);
    for (front = 0, back = num - 1;; front++, back--) {
        while (fabs(coeffs[positions[front]]) > pivot) {
            front++;
        }

        while (pivot > fabs(coeffs[positions[back]])) {
            if (back == 0)
                break;
            back--;
        }

        if (front >= back)
            break;

        temp = positions[front];
        positions[front] = positions[back];
        positions[back] = temp;
    }
    sort_positions(positions, coeffs, front);
    sort_positions(positions + front, coeffs, num - front);
}

static int get_bit(uint8_t* mask, size_t pos)
{
    return (mask[pos / CHAR_BIT] & (1 << (pos % CHAR_BIT))) != 0;
}

static void set_bit(uint8_t* mask, size_t pos)
{
    mask[pos / CHAR_BIT] |= (1 << (pos % CHAR_BIT));
}

static double compute_discrepancy(double a[], double b[], size_t num)
{
    size_t idx;
    double err = 0.0;
    for (idx = 0; idx < num; ++idx) {
        double cur_err = fabs(a[idx] - b[idx]);
        if (cur_err > err)
            err = cur_err;
    }

    return err;
}

/* On exit, the bits in "bitmask" will be set to 1 in correspondence of
 * every Chebyshev coefficient that must be retained. The function
 * returns the number of Chebyshev coefficients to retain (i.e., the
 * number of bits in "mask" that have been set to 1).*/
static size_t trunc_chebyshev(pcomp_chebyshev_t* chebyshev,
                              pcomp_chebyshev_t* inv_chebyshev,
                              double max_allowable_error, uint8_t* mask,
                              double* max_error)
{
    size_t idx;
    size_t cur_coeff = 0;
    pcomp_chunk_size_t* positions;
    double err;

    if (chebyshev == NULL || inv_chebyshev == NULL || mask == NULL
        || chebyshev->num_of_samples != inv_chebyshev->num_of_samples)
        abort();

    if (max_allowable_error <= 0.0)
        return 0;

    /* The "positions" array contains the indexes to the
     * chebyshev->output array, i.e., the list of Chebyshev
     * coefficients to sort in decreasing order. At the beginning each
     * entry is set to its own index:
     *
     *            +---+---+---+-----+
     * positions: | 0 | 1 | 2 | ... |
     *            +---+---+---+-----+
     *
     * After the call to "sort_positions", the "positions" array
     * contains the indexes to "chebyshev->output" ordered according
     * to the decreasing absolute value of the latters, e.g.:
     *
     *            +---+---+---+-----+
     * positions: | 7 | 2 | 5 | ... |
     *            +---+---+---+-----+
     */
    positions = malloc(chebyshev->num_of_samples
                       * sizeof(pcomp_chunk_size_t));
    if (positions == NULL)
	abort();
    for (idx = 0; idx < chebyshev->num_of_samples; ++idx) {
        positions[idx] = idx;
    }
    sort_positions(positions, chebyshev->output,
                   chebyshev->num_of_samples);

    /* Start by setting all the coefficients to zero */
    memset(&inv_chebyshev->input[0], 0,
           chebyshev->num_of_samples * sizeof(inv_chebyshev->input[0]));
    memset(&mask[0], 0,
           pcomp_chunk_cheby_mask_size(chebyshev->num_of_samples));

    /* Add the coefficients one by one until the error is below
     * "max_allowable_error" */
    cur_coeff = 0;
    while (cur_coeff < chebyshev->num_of_samples) {

        inv_chebyshev->input[positions[cur_coeff]]
            = chebyshev->output[positions[cur_coeff]];
        set_bit(mask, positions[cur_coeff]);
        ++cur_coeff;

        pcomp_run_chebyshev(inv_chebyshev, inv_chebyshev->dir, NULL,
                            NULL);
        err = compute_discrepancy(chebyshev->input,
                                  inv_chebyshev->output,
                                  chebyshev->num_of_samples);

        if (err < max_allowable_error)
            break;
    }

    if (max_error != NULL)
        *max_error = err;

    return cur_coeff;
}

/* This function is used internally by
 * "pcomp_run_polycomp_on_chunk"
 * and "pcomp_decompress_poly_chunk" to make sure there is no
 * memory
 * leak on the chunk passed as argument. */
static void clear_chunk(pcomp_polycomp_chunk_t* chunk)
{
    /* Leave chunk->uncompressed as it is, as it never changes
     */

    chunk->num_of_poly_coeffs = 0;
    if (chunk->poly_coeffs != NULL)
        free(chunk->poly_coeffs);

    chunk->num_of_cheby_coeffs = 0;
    if (chunk->cheby_coeffs != NULL)
        free(chunk->cheby_coeffs);
}

int pcomp_run_polycomp_on_chunk(pcomp_polycomp_t* params,
                                const double* input,
                                pcomp_chunk_size_t num_of_samples,
                                pcomp_polycomp_chunk_t* chunk,
                                double* max_error)
{
    uint8_t* mask = NULL;
    double* coeffs = NULL;
    size_t cheby_coeffs_to_retain = 0;
    int apply_chebyshev = 1;
    double* buf = NULL;
    const double* straightened_input;

    if (chunk == NULL || input == NULL || params == NULL
        || params->poly_fit == NULL || params->chebyshev == NULL
        || params->inv_chebyshev == NULL)
        abort();

    clear_chunk(chunk);

    if (num_of_samples != params->samples_per_chunk)
        return PCOMP_STAT_INVALID_BUFFER;

    if (params->period > 0.0) {
        buf = malloc(num_of_samples * sizeof(input[0]));
	if (buf == NULL)
	    abort();
        pcomp_straighten(buf, input, num_of_samples, params->period);
        straightened_input = buf; /* This preserve const-correctness */
    }
    else {
        straightened_input = input;
    }

    if (num_of_samples <= params->poly_fit->num_of_coeffs) {
        /* The number of element is so small that is better to
         * store
         * them uncompressed */
        chunk->is_compressed = 0;
    }
    else {
        double max_residual;

        /* Compute the polynomial fit and the full Chebyshev
         * transform */
        coeffs
            = malloc(sizeof(double) * params->poly_fit->num_of_coeffs);
	if (coeffs == NULL)
	    abort();
        polyfit_and_chebyshev(params, coeffs, straightened_input,
                              &max_residual);
        apply_chebyshev
            = (max_residual >= params->max_allowable_error)
              && (params->algorithm != PCOMP_ALG_NO_CHEBYSHEV);

        /* If the Chebyshev transform is needed, truncate it as much
         * as possible */
        if (apply_chebyshev) {
            mask = malloc(pcomp_chunk_cheby_mask_size(num_of_samples));
	    if (mask == NULL)
		abort();
            cheby_coeffs_to_retain = trunc_chebyshev(
                params->chebyshev, params->inv_chebyshev,
                params->max_allowable_error, mask, max_error);

            chunk->is_compressed = (cheby_coeffs_to_retain
                                    + params->poly_fit->num_of_coeffs)
                                   < num_of_samples;
            if (!chunk->is_compressed) {
                free(mask);
                mask = NULL;
            }
        }
        else {
            /* Assume that num_of_samples > deg(p) + 1 */
            chunk->is_compressed
                = (max_residual <= params->max_allowable_error);
        }
    }

    chunk->num_of_samples = num_of_samples;
    if (chunk->is_compressed) {
        size_t idx;

        chunk->num_of_poly_coeffs = params->poly_fit->num_of_coeffs;
        chunk->poly_coeffs = coeffs;
        if (apply_chebyshev) {
            size_t cheby_idx;

            chunk->num_of_cheby_coeffs = cheby_coeffs_to_retain;
            chunk->cheby_mask = mask;
            chunk->cheby_coeffs
                = malloc(sizeof(double) * cheby_coeffs_to_retain);
	    if (chunk->cheby_coeffs == NULL)
		abort();
            cheby_idx = 0;
            for (idx = 0; idx < params->chebyshev->num_of_samples;
                 ++idx) {
                if (get_bit(mask, idx)) {
                    chunk->cheby_coeffs[cheby_idx++]
                        = params->chebyshev->output[idx];
                }
            }
        }
        else {
            chunk->num_of_cheby_coeffs = 0;
            chunk->cheby_mask = NULL;
            chunk->cheby_coeffs = NULL;
        }
    }
    else {
        size_t idx;

        if (coeffs != NULL)
            free(coeffs);

        chunk->uncompressed = malloc(sizeof(double) * num_of_samples);
	if (chunk->uncompressed == NULL)
	    abort();
        for (idx = 0; idx < num_of_samples; ++idx)
            chunk->uncompressed[idx] = input[idx];

        if (max_error != NULL)
            *max_error = 0.0;
    }

    if (buf != NULL) {
        free(buf);
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

        /* Compute the values of the polynomial at the points 1,
         * 2, ...
         */
        for (idx = 0; idx < chunk->num_of_samples; ++idx) {
            output[idx] = eval_poly(chunk->poly_coeffs,
                                    chunk->num_of_poly_coeffs, idx + 1);
        }

        /* If present, add the contribution of the Chebyshev
         * transform
         */
        if (chunk->num_of_cheby_coeffs > 0) {
            size_t cur_cheby_idx = 0;
            if (chunk->cheby_coeffs == NULL)
                abort();

            for (idx = 0; idx < chunk->num_of_samples; ++idx) {
                if (get_bit(chunk->cheby_mask, idx)) {
                    if (cur_cheby_idx >= chunk->num_of_cheby_coeffs) {
                        abort();
                    }

                    inv_chebyshev->input[idx]
                        = chunk->cheby_coeffs[cur_cheby_idx++];
                }
                else {
                    inv_chebyshev->input[idx] = 0.0;
                }
            }
            pcomp_run_chebyshev(inv_chebyshev, inv_chebyshev->dir, NULL,
                                NULL);

            for (idx = 0; idx < chunk->num_of_samples; ++idx) {
                output[idx] += inv_chebyshev->output[idx];
            }
        }
    }
    else {
        memcpy(output, chunk->uncompressed,
               sizeof(chunk->uncompressed[0]) * chunk->num_of_samples);
    }

    return PCOMP_STAT_SUCCESS;
}

/***********************************************************************
 *
 */

int pcomp_compress_polycomp(pcomp_polycomp_chunk_t** output_buf[],
                            size_t* num_of_chunks,
                            const double* input_buf, size_t input_size,
                            const pcomp_polycomp_t* params)
{
    size_t idx;
    const double* cur_input = input_buf;
    pcomp_polycomp_t* chunk_params;

    if (output_buf == NULL || num_of_chunks == NULL || input_buf == NULL
        || params == NULL || params->poly_fit == NULL)
        abort();

    /* Calculate how many chunks we'll create */
    *num_of_chunks = input_size / params->samples_per_chunk;
    if (input_size % params->samples_per_chunk != 0)
        ++(*num_of_chunks);

    *output_buf
        = malloc(sizeof(pcomp_polycomp_chunk_t*) * (*num_of_chunks));
    if (*output_buf == NULL)
	abort();

    chunk_params = pcomp_init_polycomp(
        params->samples_per_chunk, params->poly_fit->num_of_coeffs,
        params->max_allowable_error, params->algorithm);

    /* Loop over the chunks and call
     * "pcomp_run_polycomp_on_chunk" for
     * each of them */
    for (idx = 0; idx < *num_of_chunks; ++idx) {
        size_t cur_chunk_size = params->samples_per_chunk;

        if ((cur_input - input_buf) + cur_chunk_size > input_size)
            cur_chunk_size
                = (size_t)(input_size - (cur_input - input_buf));

        if (cur_chunk_size != chunk_params->samples_per_chunk) {
            pcomp_free_polycomp(chunk_params);

            chunk_params = pcomp_init_polycomp(
                cur_chunk_size, params->poly_fit->num_of_coeffs,
                params->max_allowable_error, params->algorithm);
        }

        (*output_buf)[idx] = pcomp_init_chunk(cur_chunk_size);
        pcomp_run_polycomp_on_chunk(chunk_params, cur_input,
                                    cur_chunk_size, (*output_buf)[idx],
                                    NULL);

        cur_input += cur_chunk_size;
    }

    return PCOMP_STAT_SUCCESS;
}

size_t
pcomp_total_num_of_samples(pcomp_polycomp_chunk_t* const chunk_array[],
                           size_t num_of_chunks)
{
    size_t total = 0;
    size_t idx;
    if (chunk_array == NULL)
        abort();

    for (idx = 0; idx < num_of_chunks; ++idx) {
        if (chunk_array[idx] == NULL)
            abort();

        total += chunk_array[idx]->num_of_samples;
    }

    return total;
}

int pcomp_decompress_polycomp(
    double* output_buf, pcomp_polycomp_chunk_t* const chunk_array[],
    size_t num_of_chunks)
{
    size_t idx;
    double* cur_output_pos = output_buf;
    pcomp_chebyshev_t* inv_chebyshev = NULL;

    if (output_buf == NULL || chunk_array == NULL)
        abort();

    for (idx = 0; idx < num_of_chunks; ++idx) {
        if (chunk_array[idx] == NULL)
            abort();

        if (inv_chebyshev == NULL
            || inv_chebyshev->num_of_samples
                   != chunk_array[idx]->num_of_samples) {

            /* This does no harm if inv_chebyshev==NULL */
            pcomp_free_chebyshev(inv_chebyshev);

            inv_chebyshev = pcomp_init_chebyshev(
                chunk_array[idx]->num_of_samples, PCOMP_TD_INVERSE);
        }

        pcomp_decompress_polycomp_chunk(
            cur_output_pos, chunk_array[idx], inv_chebyshev);

        cur_output_pos += chunk_array[idx]->num_of_samples;
    }

    pcomp_free_chebyshev(inv_chebyshev);

    return PCOMP_STAT_SUCCESS;
}

void pcomp_free_chunks(pcomp_polycomp_chunk_t* chunk_array[],
                       size_t num_of_chunks)
{
    size_t idx;

    if (chunk_array == NULL)
        return;

    for (idx = 0; idx < num_of_chunks; ++idx) {
        pcomp_free_chunk(chunk_array[idx]);
    }

    free(chunk_array);
}

size_t pcomp_chunks_num_of_bytes(pcomp_polycomp_chunk_t* const chunks[],
                                 size_t num_of_chunks)
{
    size_t result = sizeof(size_t); /* Room for the number of chunks */
    size_t idx;

    for (idx = 0; idx < num_of_chunks; ++idx) {
        result += pcomp_chunk_num_of_bytes(chunks[idx]);
    }

    return result;
}

/***********************************************************************
 * Encode/decode a list of chunks into a raw stream of bytes (suitable
 * for I/O).
 */

#define SAVE_TO_PTR_AND_INCREMENT(buf, value, type)                    \
    {                                                                  \
        *((type*)buf) = value;                                         \
        buf = ((type*)buf) + 1;                                        \
    }

int pcomp_encode_chunks(void* buf, size_t* buf_size,
                        pcomp_polycomp_chunk_t* const chunk_array[],
                        size_t num_of_chunks)
{
    void* buf_ptr = buf;
    size_t chunk_idx;

    if (chunk_array == NULL || num_of_chunks == 0)
        abort();

    SAVE_TO_PTR_AND_INCREMENT(buf_ptr, num_of_chunks, size_t);

    for (chunk_idx = 0; chunk_idx < num_of_chunks; ++chunk_idx) {
        const pcomp_polycomp_chunk_t* cur_chunk
            = chunk_array[chunk_idx];
        size_t idx; /* Used for inner loops */

        SAVE_TO_PTR_AND_INCREMENT(buf_ptr, cur_chunk->is_compressed,
                                  uint8_t);
        SAVE_TO_PTR_AND_INCREMENT(buf_ptr, cur_chunk->num_of_samples,
                                  pcomp_chunk_size_t);

        if (cur_chunk->is_compressed) {
            size_t cheby_mask_size = pcomp_chunk_cheby_mask_size(
                cur_chunk->num_of_samples);

            SAVE_TO_PTR_AND_INCREMENT(buf_ptr,
                                      cur_chunk->num_of_poly_coeffs,
                                      pcomp_poly_size_t);
            for (idx = 0; idx < cur_chunk->num_of_poly_coeffs; idx++) {
                SAVE_TO_PTR_AND_INCREMENT(
                    buf_ptr, cur_chunk->poly_coeffs[idx], double);
            }

            SAVE_TO_PTR_AND_INCREMENT(buf_ptr,
                                      cur_chunk->num_of_cheby_coeffs,
                                      pcomp_chunk_size_t);
            if (cur_chunk->num_of_cheby_coeffs > 0) {
                /* Mask */
                for (idx = 0; idx < cheby_mask_size; ++idx) {
                    SAVE_TO_PTR_AND_INCREMENT(
                        buf_ptr, cur_chunk->cheby_mask[idx], uint8_t);
                }
                /* Chebyshev coefficients */
                for (idx = 0; idx < cur_chunk->num_of_cheby_coeffs;
                     idx++) {
                    SAVE_TO_PTR_AND_INCREMENT(
                        buf_ptr, cur_chunk->cheby_coeffs[idx], double);
                }
            }
        }
        else {
            for (idx = 0; idx < cur_chunk->num_of_samples; idx++) {
                SAVE_TO_PTR_AND_INCREMENT(
                    buf_ptr, cur_chunk->uncompressed[idx], double);
            }
        }
    }

    *buf_size = ((uint8_t*)buf_ptr) - ((uint8_t*)buf);
    return PCOMP_STAT_SUCCESS;
}

#define READ_FROM_PTR_AND_INCREMENT(var, pointer, type)                \
    {                                                                  \
        var = *((type*)(pointer));                                     \
        pointer = ((type*)(pointer)) + 1;                              \
    }

int pcomp_decode_chunks(pcomp_polycomp_chunk_t** chunk_array[],
                        size_t* num_of_chunks, const void* buf)
{
    const void* cur_ptr = buf;
    size_t chunk_idx;
    double* poly_buf = NULL;
    size_t poly_buf_size = 0;
    double* cheby_buf = NULL;
    size_t cheby_buf_size = 0;
    uint8_t* cheby_mask_buf = NULL;
    size_t cheby_mask_buf_size = 0;
    double* uncompr_buf = NULL;
    size_t uncompr_buf_size = 0;

    if (buf == NULL || chunk_array == NULL || num_of_chunks == NULL)
        abort();

    READ_FROM_PTR_AND_INCREMENT(*num_of_chunks, cur_ptr, size_t);
    *chunk_array
        = malloc(sizeof(pcomp_polycomp_chunk_t*) * (*num_of_chunks));
    if (*chunk_array == NULL)
	abort();

    for (chunk_idx = 0; chunk_idx < *num_of_chunks; ++chunk_idx) {
        uint8_t is_compressed;
        size_t num_of_samples;
        size_t idx; /* Used for inner loops */

        READ_FROM_PTR_AND_INCREMENT(is_compressed, cur_ptr, uint8_t);
        READ_FROM_PTR_AND_INCREMENT(num_of_samples, cur_ptr,
                                    pcomp_chunk_size_t);

        if (is_compressed) {
            size_t num_of_poly_coeffs;
            size_t num_of_cheby_coeffs;
            size_t cheby_mask_size
                = pcomp_chunk_cheby_mask_size(num_of_samples);

            READ_FROM_PTR_AND_INCREMENT(num_of_poly_coeffs, cur_ptr,
                                        pcomp_poly_size_t);
            if (num_of_poly_coeffs > poly_buf_size) {
                poly_buf_size = num_of_poly_coeffs;
                poly_buf
                    = realloc(poly_buf, poly_buf_size * sizeof(double));
            }
            for (idx = 0; idx < num_of_poly_coeffs; ++idx) {
                READ_FROM_PTR_AND_INCREMENT(poly_buf[idx], cur_ptr,
                                            double);
            }

            READ_FROM_PTR_AND_INCREMENT(num_of_cheby_coeffs, cur_ptr,
                                        pcomp_chunk_size_t);
            if (num_of_cheby_coeffs > 0) {
                if (cheby_mask_size > cheby_mask_buf_size) {
                    cheby_mask_buf_size = cheby_mask_size;
                    cheby_mask_buf = realloc(cheby_mask_buf,
                                             cheby_mask_buf_size
                                                 * sizeof(uint8_t));
                }
                for (idx = 0; idx < cheby_mask_size; ++idx) {
                    READ_FROM_PTR_AND_INCREMENT(cheby_mask_buf[idx],
                                                cur_ptr, uint8_t);
                }

                if (num_of_cheby_coeffs > cheby_buf_size) {
                    cheby_buf_size = num_of_cheby_coeffs;
                    cheby_buf = realloc(
                        cheby_buf, cheby_buf_size * sizeof(double));
                }
                for (idx = 0; idx < num_of_cheby_coeffs; ++idx) {
                    READ_FROM_PTR_AND_INCREMENT(cheby_buf[idx], cur_ptr,
                                                double);
                }
            }

            (*chunk_array)[chunk_idx] = pcomp_init_compressed_chunk(
                num_of_samples, num_of_poly_coeffs, poly_buf,
                num_of_cheby_coeffs, cheby_mask_buf, cheby_buf);
        }
        else {
            if (num_of_samples > uncompr_buf_size) {
                uncompr_buf_size = num_of_samples;
                uncompr_buf = realloc(
                    uncompr_buf, uncompr_buf_size * sizeof(double));
            }
            for (idx = 0; idx < num_of_samples; ++idx) {
                READ_FROM_PTR_AND_INCREMENT(uncompr_buf[idx], cur_ptr,
                                            double);
            }

            (*chunk_array)[chunk_idx] = pcomp_init_uncompressed_chunk(
                num_of_samples, uncompr_buf);
        }
    }

    if (poly_buf != NULL)
        free(poly_buf);

    if (cheby_buf != NULL)
        free(cheby_buf);

    if (uncompr_buf != NULL)
        free(uncompr_buf);

    return PCOMP_STAT_SUCCESS;
}
