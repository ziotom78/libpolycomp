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

#ifdef WITH_OPENMP

#include <omp.h>

#else

typedef int omp_int_t;
static omp_int_t omp_get_thread_num(void) { return 0; }
static omp_int_t omp_get_max_threads(void) { return 1; }

#endif

/**********************************************************************/

/** \defgroup poly Polynomial compression functions
 *
 * Polynomial compression relies on a simple idea, that is to divide
 * the input data stream into subsets of consecutive samples (called
 * "chunks"), and to approximate each chunk by means of a polynomial.
 * Such compression is inherently lossy, as the residuals of the
 * fitting procedure are usually discarded. If the polynomial used for
 * the fitting produces residuals that are too large, usually the
 * samples in the chunk are saved in uncompressed form.
 *
 * This idea has been widely applied in the literature. Libpolycomp
 * implements an improvement over it, because if the fit residuals are
 * too large, the library saves a chopped sequence of the Chebyshev
 * transform of the residuals. This allows to achieve better
 * compression ratios in those cases where polynomial fitting is not
 * always enough to keep compression errors below the desired
 * threshold. This kind of compression works quite well for smooth
 * data series, where changes between consecutive samples are well
 * described by slowly varying continuous functions. It is not
 * suitable if the signal contains noise, unless this noise is
 * significantly smaller than the signal and than the error threshold.
 *
 * Libpolycomp allows to avoid the usage of Chebyshev transforms. In
 * this case, if no polynomial of the desired degree are able to fit
 * the data with the given error threshold, the data for that chunk is
 * saved uncompressed.
 *
 * The typical workflow for applying polynomial compression is the
 * following:
 *
 * 1. Allocate a new \ref pcomp_polycomp_t object via a call to \ref
 *    pcomp_init_polycomp. Such object contains the parameters to be
 *    used for the compression, e.g., the size of each chunk, the
 *    degree of the fitting polynomial, whether to apply or not the
 *    Chebyshev transform to the residuals, etc.
 *
 * 2. Split the data into chunks and compress each of them using the
 *    function \ref pcomp_compress_polycomp.
 *
 * 3. Convert the list of chunks into a byte sequence using \ref
 *    pcomp_encode_chunks, typically with the purpose of saving it
 *    into a file or sending it through a pipe/socket/etc.
 *
 * The decompression workflow is specular:
 *
 * 1. Process the byte sequence containing the compressed data using
 *    \ref pcomp_decode_chunks. This will produce a list of chunks
 *    that are still compressed.
 *
 * 2. Decompress the chunks using the function \ref
 *    pcomp_decompress_polycomp.
 *
 * The compression functions described in this page use the \ref
 * pcomp_polycomp_t structure to determine which parameters to use for
 * the compression. The functions that allow to allocate/free/manage
 * this structure are the following:
 *
 * - \ref pcomp_init_polycomp and \ref pcomp_free_polycomp
 * - \ref pcomp_polycomp_samples_per_chunk
 * - \ref pcomp_polycomp_num_of_poly_coeffs
 * - \ref pcomp_polycomp_max_error
 * - \ref pcomp_polycomp_algorithm
 * - \ref pcomp_polycomp_period and \ref pcomp_polycomp_set_period
 *
 * It is possible to use a set of more low-level functions to use
 * polynomial compression. Refer to \ref poly_lowlevel for further
 * information.
 */

/** \defgroup poly_lowlevel Polynomial compression (low-level functions)
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

/** \ingroup polyfit
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

/** \ingroup polyfit
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

/** \ingroup polyfit
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

/** \ingroup polyfit
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

/** \ingroup polyfit
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

/** \ingroup cheby
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

/** \ingroup cheby
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

/** \ingroup cheby
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

/** \ingroup cheby
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

/** \ingroup cheby
 *
 * \brief Compute a forward/backward Chebyshev discrete transform
 *
 * \code{.c}
 * #define NUM_OF_POINTS 3
 * double points[NUM_OF_POINTS] = { 0.0, 1.0, 3.0 };
 * double transform[NUM_OF_POINTS];
 * pcomp_chebyshev_t* chebyshev;
 * size_t idx;
 *
 * chebyshev = pcomp_init_chebyshev(NUM_OF_POINTS, PCOMP_TD_DIRECT);
 * pcomp_run_chebyshev(chebyshev, PCOMP_TD_DIRECT, transform, points);
 *
 * puts("Transform:");
 * for (idx = 0; idx < NUM_OF_POINTS; ++idx) {
 *     printf("%f\t", transform[idx]);
 * }
 * puts("");
 * \endcode
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

/** \ingroup cheby
 *
 * \brief Return the input data used in the last call to \ref
 *pcomp_run_chebyshev
 *
 * If \ref pcomp_run_chebyshev was never called, the array returned by
 * this function contains garbage.
 *
 * \param[in] plan Pointer to a Chebyshev plan created by \ref
 * pcomp_init_chebyshev
 *
 * \return A pointer to the first element of the array of elements
 * used as input by the last call to \ref pcomp_run_chebyshev.
 */
const double* pcomp_chebyshev_input(const pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        abort();
    return plan->input;
}

/** \ingroup cheby
 *
 * \brief Return the output (Chebyshev transform) of the last call to
 *\ref pcomp_run_chebyshev
 *
 * If \ref pcomp_run_chebyshev was never called, the array returned by
 * this function contains garbage.
 *
 * \param[in] plan Pointer to a Chebyshev plan created by \ref
 * pcomp_init_chebyshev
 *
 * \return A pointer to the first element of the array of elements
 * containing the output of the last call to \ref pcomp_run_chebyshev.
 */
const double* pcomp_chebyshev_output(const pcomp_chebyshev_t* plan)
{
    if (plan == NULL)
        abort();
    return plan->output;
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

/** \ingroup poly_lowlevel
 *
 * \brief Return a pointer to a \ref pcomp_chebyshev_t structure
 * representing the forward Chebyshev transform.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The algorithm to be used by the compressor.
 */
pcomp_chebyshev_t*
pcomp_polycomp_forward_cheby(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->chebyshev;
}

/** \ingroup poly_lowlevel
 *
 * \brief Return a pointer to a \ref pcomp_chebyshev_t structure
 * representing the forward Chebyshev transform.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * containing the compression parameters
 *
 * \returns The algorithm to be used by the compressor.
 */
pcomp_chebyshev_t*
pcomp_polycomp_backward_cheby(const pcomp_polycomp_t* params)
{
    if (params == NULL)
        abort();

    return params->inv_chebyshev;
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
        int idx = num_of_coeffs - 1;
        double result = coeffs[idx];

        if (num_of_coeffs == 1)
            return result;

        for (idx = num_of_coeffs - 2; idx >= 0; --idx)
            result = result * x + coeffs[idx];

        return result;
    }
    else
        return 0.0;
}

/***********************************************************************/

/** \ingroup poly_lowlevel
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
                                 * series is chopped. */
    double* cheby_coeffs;
};

/** \ingroup poly_lowlevel
 *
 * \brief Allocate memory for a \ref pcomp_polycomp_chunk_t object
 *
 * \param[in] num_of_samples Number of samples that the chunk will be
 * capable to hold.
 *
 * \return A pointer to the newly allocated object. Use \ref
 * pcomp_free_chunk to free the memory once is no longer needed.
 */
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

/** \ingroup poly_lowlevel
 *
 * \brief Allocate memory for a \ref pcomp_polycomp_chunk_t object and
 * fill it with data in uncompressed form.
 *
 * \param[in] num_of_samples Number of samples that the chunk will be
 * capable to hold.
 *
 * \param[in] samples The (uncompressed) samples to copy into the
 * chunk. After the call, \a input is no longer needed and can be
 * freed without invalidating the pointer returned by the function.
 *
 * \return A pointer to the newly allocated object. Use \ref
 * pcomp_free_chunk to free the memory once is no longer needed.
 */
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

/** \ingroup poly_lowlevel
 *
 * \brief Allocate memory for a \ref pcomp_polycomp_chunk_t object and
 * fill it with data compressed using the polynomial compression
 * algorithm.
 *
 * \param[in] num_of_samples Number of samples that the chunk will be
 * capable to hold.
 *
 * \param[in] num_of_poly_coeffs Number of coefficients of the
 * interpolating polynomial.
 *
 * \param[in] poly_coeffs Pointer to the coefficients of the
 * interpolating polynomial. Their number must be equal to the
 * parameter \a num_of_poly_coeffs.
 *
 * \param[in] num_of_cheby_coeffs Number of nonzero Chebyshev
 * coefficients associated with the polynomial fit. This number is
 * always less than \a num_of_samples. Zero is allowed.
 *
 * \param[in] cheby_mask Bitmask representing the position of the
 * nonzero coefficients in \a cheby_coeffs within the full sequence.
 * (Use \ref pcomp_mask_get_bit and \ref pcomp_mask_set_bit to
 * read/write bits in the sequence.)
 *
 * \param[in] cheby_coeffs Array of nonzero Chebyshev coefficients.
 * Their number must be equal to \a num_of_cheby_coeffs.
 *
 * \return A pointer to the newly allocated object. Use \ref
 * pcomp_free_chunk to free the memory once is no longer needed.
 */
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

/** \ingroup poly_lowlevel
 *
 * \brief Free memory associated with a \ref pcomp_poly_chunk_t
 *
 * This function releases the memory allocated by one of the following
 * functions:
 * - \ref pcomp_init_chunk
 * - \ref pcomp_init_uncompressed_chunk
 * - \ref pcomp_init_compressed_chunk
 *
 * \param[in] chunk Pointer to the object to be freed.
 */
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

/** \ingroup poly_lowlevel
 *
 * \brief Return the number of samples in a chunk
 *
 * \param[in] chunk Pointer to the chunk data
 *
 * \returns The number of samples
 */
pcomp_chunk_size_t
pcomp_chunk_num_of_samples(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->num_of_samples;
}

/** \ingroup poly_lowlevel
 *
 * \brief Return the number of bytes necessary to encode a chunk
 *
 * Refer to \ref pcomp_encode_chunks and \ref pcomp_decode_chunks for
 * further details.
 *
 * \param[in] chunk Pointer to the chunk data
 *
 * \returns The number of bytes
 */
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

/** \ingroup poly_lowlevel
 *
 * \brief Return nonzero if the chunk holds data in uncompressed form.
 */
int pcomp_chunk_is_compressed(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->is_compressed;
}

/** \ingroup poly_lowlevel
 *
 * \brief If the chunks contain uncompressed data, returns a pointer
 * to the first element. Otherwise, return \c NULL.
 */
const double*
pcomp_chunk_uncompressed_data(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (chunk->is_compressed)
        return NULL;

    return chunk->uncompressed;
}

/** \ingroup poly_lowlevel
 *
 * \brief If the chunks contain compressed data, returns the number of
 * polynomial coefficients used in the compression. Otherwise, return
 * zero.
 */
pcomp_poly_size_t
pcomp_chunk_num_of_poly_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (!chunk->is_compressed)
        return 0;

    return chunk->num_of_poly_coeffs;
}

/** \ingroup poly_lowlevel
 *
 * \brief If the chunks contain compressed data, returns a pointer to
 * the first element of the array of coefficients of the interpolating
 * polynomial. Otherwise, return \c NULL.
 */
const double*
pcomp_chunk_poly_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (!chunk->is_compressed || chunk->num_of_poly_coeffs == 0)
        return NULL;

    return chunk->poly_coeffs;
}

/** \ingroup poly_lowlevel
 *
 * \brief If the chunks contain compressed data, returns the number of
 * nonzero Chebyshev coefficients held in the chunk. Otherwise, return
 * zero.
 */
pcomp_chunk_size_t
pcomp_chunk_num_of_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (!chunk->is_compressed)
        return 0;

    return chunk->num_of_cheby_coeffs;
}

/** \ingroup poly_lowlevel
 *
 * \brief If the chunks contain compressed data, returns a pointer to
 * the first element of the Chebyshev transform of the fit residuals.
 * Otherwise, return \c NULL.
 */
const double*
pcomp_chunk_cheby_coeffs(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    if (!chunk->is_compressed || chunk->num_of_cheby_coeffs == 0)
        return NULL;

    return chunk->cheby_coeffs;
}

/** \ingroup poly_lowlevel
 *
 * \brief Return the number of bytes required for the bitmask of
 * nonzero Chebyshev coefficients.
 *
 * The polynomial compression compresses Chebyshev transforms by
 * saving only those coefficients that are significantly different
 * from zero. In order to keep track of the position of such
 * coefficients in the full array, a bit mask is used. This function
 * determines how many bytes are required for such mask, which is
 * internally represented by Libpolycomp as an array of \c uint8_t
 * values.
 *
 * \param[in] chunk_size Number of samples in the chunk
 *
 * \returns The number of bytes (\c uint8_t values) required for the
 * mask.
 */
size_t pcomp_chunk_cheby_mask_size(pcomp_chunk_size_t chunk_size)
{
    return chunk_size / CHAR_BIT
           + ((chunk_size % CHAR_BIT) > 0 ? 1 : 0);
}

/** \ingroup poly_lowlevel
 *
 * \brief Return a pointer to the bitmask of nonzero Chebyshev
 * coefficients for a chunk
 *
 * \param[in] chunk Pointer to the chunk
 *
 * \returns A pointer to the array of bytes which make up the mask.
 * Use \ref pcomp_mask_get_bit to access the values of each bit.
 */
const uint8_t*
pcomp_chunk_cheby_mask(const pcomp_polycomp_chunk_t* chunk)
{
    if (chunk == NULL)
        abort();

    return chunk->cheby_mask;
}

/**********************************************************************/

/** \ingroup poly_lowlevel
 *
 * \brief Compute a polynomial fit of the data in \a input and a
 * Chebyshev transform of the residuals
 *
 * Note that this function *always* computes the Chebyshev transform
 * of the data, even if there is a perfect fit between the polynomial
 * and the input data.
 *
 * \param[in] params Pointer to a \ref pcomp_polycomp_t structure
 * initialized by \ref pcomp_init_polycomp.
 *
 * \param[out] coeffs Pointer to the array that on exit will hold the
 * coefficients of the best-fit polynomial. It must have enough room
 * for a number of elements equal to the return value of \ref
 * pcomp_polycomp_num_of_poly_coeffs.
 *
 * \param[out] cheby_residuals Pointer to an array that on exit will
 * contain the Chebyshev transform of the residuals of the fit. It can
 * be \c NULL; in any case, these numbers can be obtained by the use
 * of a call to \ref pcomp_polycomp_forward_cheby and \ref
 * pcomp_chebyshev_output.
 *
 * \param[in] input Pointer to the array of values to be transformed.
 * The number of values used is equal to the return value of the
 * function \ref pcomp_polycomp_samples_per_chunk.
 *
 * \param[out] max_residual Pointer to a variable that will hold the
 * maximum absolute value of the discrepancy between each sample in \a
 * input and the polynomial fit. It can be \c NULL.
 *
 * \returns If no errors occurred, \ref PCOMP_STAT_SUCCESS. Otherwise,
 * the function returns the code of the error.
 */

int pcomp_polyfit_and_chebyshev(pcomp_polycomp_t* params,
                                double* coeffs, double* cheby_residuals,
                                const double* input,
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

        if (cheby_residuals != NULL
            && cheby_residuals != params->chebyshev->output) {
            memcpy(cheby_residuals, params->chebyshev->output,
                   sizeof(params->chebyshev->output[0])
                       * params->chebyshev->num_of_samples);
        }
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

/** \ingroup poly_lowlevel
 *
 * \brief Return the value of the bit at the position \a pos in the
 * bitmask \a mask.
 *
 * \param[in] mask Pointer to the first byte of the mask
 *
 * \param[in] pos Zero-based index of the bit in the mask
 *
 * \returns Either 0 or 1, depending on the value of the bit.
 */
int pcomp_mask_get_bit(const uint8_t* mask, size_t pos)
{
    return (mask[pos / CHAR_BIT] & (1 << (pos % CHAR_BIT))) != 0;
}

/** \ingroup poly_lowlevel
 *
 * \brief Set the value of the bit at position \a pos in the bitmask \a
 *
 * \param[inout] mask The bitmask to modify
 *
 * \param[in] pos Zero-based index of the byte to set
 *
 * \param[in] value Value of the bit (either 0 or 1; any value
 * different from zero is treated as equal to 1)
 */
void pcomp_mask_set_bit(uint8_t* mask, size_t pos, int value)
{
    if (value != 0) {
        mask[pos / CHAR_BIT] |= (1 << (pos % CHAR_BIT));
    }
    else {
        mask[pos / CHAR_BIT] &= ~(((uint8_t)1) << (pos % CHAR_BIT));
    }
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

/** \ingroup poly_lowlevel
 *
 * \brief Find the smallest subset of Chebyshev coefficients that can
 * approximate a Chebyshev transform with an error less than \a
 * max_allowable_error.
 *
 * On exit, the bits in \a bitmask will be set to 1 in correspondence
 * of every Chebyshev coefficient that must be retained. The function
 * returns the number of Chebyshev coefficients to retain (i.e., the
 * number of bits in \a mask that have been set to 1).
 *
 * \param[in] chebyshev Pointer to a \ref pcomp_chebyshev_t structure
 * used to compute the forward Chebyshev transform (from the space of
 * the fit residuals to the Chebyshev space)
 *
 * \param[in] inv_chebyshev Pointer to a \ref pcomp_chebyshev_t
 * structure representing the inverse transform of \a chebyshev.
 *
 * \param[in] max_allowable_error The maximum allowed discrepancy for
 * the chopped Chebyshev, as measured in the space of the fit
 * residuals.
 *
 * \param[out] mask Bitmask that will contain the position of the
 * unchopped Chebyshev terms of the transform of the fit residuals.
 * Use \ref pcomp_mask_get_bit to access each element.
 *
 * \param[out] max_error If not \c NULL, it will be set to the maximum
 * error due to the chopping of the Chebyshev transform represented by
 * \a mask.
 *
 * \returns The number of bits equal to one in \a mask, i.e., the
 * number of unchopped Chebyshev coefficients.
 */
size_t pcomp_find_chebyshev_mask(pcomp_chebyshev_t* chebyshev,
                                 pcomp_chebyshev_t* inv_chebyshev,
                                 double max_allowable_error,
                                 uint8_t* mask, double* max_error)
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
        pcomp_mask_set_bit(mask, positions[cur_coeff], 1);
        ++cur_coeff;

        pcomp_run_chebyshev(inv_chebyshev, inv_chebyshev->dir, NULL,
                            NULL);
        err = compute_discrepancy(chebyshev->input,
                                  inv_chebyshev->output,
                                  chebyshev->num_of_samples);

        if (err < max_allowable_error)
            break;
    }

    free(positions);
    if (max_error != NULL)
        *max_error = err;

    return cur_coeff;
}

/* This function is used internally by "pcomp_run_polycomp_on_chunk"
 * and "pcomp_decompress_poly_chunk" to make sure there is no memory
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

/** \ingroup poly_lowlevel
 *
 * \brief Compress the first \a num_of_samples elements in \a input
 * and store them in \a chunk.
 *
 * The function determines if the data in \a input can be efficiently
 * compressed using polynomial compression with the parameters
 * described by \a params. If it is so, it stores the compressed data
 * in \a chunk. If the compression ratio is small or equal to one, or
 * if the compression error is too large, the function copies the data
 * in \a input into \a chunk in uncompressed format.
 *
 * The following example shows how to use this function together with
 * \ref pcomp_init_chunk:
 *
 * \code{.c}
 * double input[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
 *                    7.0, 8.0, 9.0, 11.0 };
 * size_t input_size = sizeof(input) / sizeof(input[0]);
 * pcomp_polycomp_chunk_t* chunk;
 * pcomp_polycomp_t* polycomp;
 * double max_error;
 * size_t idx;
 *
 * polycomp = pcomp_init_polycomp(input_size, 2, 1.0e-5,
 *                                PCOMP_ALG_USE_CHEBYSHEV);
 * chunk = pcomp_init_chunk(input_size);
 *
 * pcomp_run_polycomp_on_chunk(polycomp, input, input_size, chunk,
 *                             &max_error);
 * \endcode
 *
 * \param[in] params Pointer to a \a pcomp_polycomp_t structure
 * (created using \ref pcomp_init_polycomp) which provides the
 * parameters of the compression.
 *
 * \param[in] input The sequence of numbers to compress. Their number
 * is equal to the parameter \a num_of_samples
 *
 * \param[in] num_of_samples Number of values in \a input to compress
 *
 * \param[inout] chunk The chunk that will contain the data, either in
 * compressed or uncompressed format. It must have already been
 * initialized via a call to \ref pcomp_init_chunk.
 *
 * \param[out] max_error On exit, the function writes the compression
 * error here. It can be \c NULL.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS (if no errors occurred), or
 * the error code.
 */
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
        /* The number of element is so small that is better to store
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
        pcomp_polyfit_and_chebyshev(params, coeffs, NULL,
                                    straightened_input, &max_residual);
        apply_chebyshev
            = (max_residual >= params->max_allowable_error)
              && (params->algorithm != PCOMP_ALG_NO_CHEBYSHEV);

        /* If the Chebyshev transform is needed, chop it as much as
         * possible */
        if (apply_chebyshev) {
            mask = malloc(pcomp_chunk_cheby_mask_size(num_of_samples));
            if (mask == NULL)
                abort();
            cheby_coeffs_to_retain = pcomp_find_chebyshev_mask(
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
                if (pcomp_mask_get_bit(mask, idx)) {
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

/** \ingroup poly_lowlevel
 *
 * \brief Decompress the data in a chunk
 *
 * This function performs the decompression of a chunk, and it is the
 * counterpart of \ref pcomp_run_polycomp_on_chunk. Here is an
 * example:
 *
 * \code{.c}
 * double* decompr;
 * pcomp_chebyshev_t* inv_chebyshev;
 *
 * // We assume that "chunk" has already been initialized somewhere
 * decompr = malloc(sizeof(double) *
 *                  pcomp_chunk_num_of_samples(chunk));
 *
 * inv_chebyshev = pcomp_init_chebyshev(input_size,
 *                                      PCOMP_TD_INVERSE);
 * pcomp_decompress_polycomp_chunk(decompr, chunk, inv_chebyshev);
 * \endcode
 *
 * \param[out] output Pointer to the array that will contain the
 * uncompressed data
 *
 * \param[in] chunk The chunk to decompress
 *
 * \param[in] inv_chebyshev Pointer to a \ref pcomp_chebyshev_t object
 * that performs the inverse Chebyshev transform. The function does
 * not allocate an object of this kind because in this way such
 * objects can be reused on subsequent calls to \ref
 * pcomp_decompress_polycomp_chunk.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS if no error occurred, or
 * the error code.
 */
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
                if (pcomp_mask_get_bit(chunk->cheby_mask, idx)) {
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

/** \ingroup poly
 *
 * \brief Compress the array \a input_buf using polynomial compression
 *
 * This function compresses the first \a input_size elements of the
 * array \a input_buf using the polynomial compression scheme. The
 * output is an array of chunks saved in \a output_buf (the number of
 * elements of this array is saved in \a num_of_chunks). The \a params
 * variable specifies the parameters used by the compression
 * algorithm.
 *
 * Here is an example showing how to compress a sequence of numbers in
 * the variable \a input:
 *
 * \code{.c}
 * double input[] = { 1.0, 2.0, 3.0, 4.0, 3.0, 2.0,
 *                    1.0, 2.0, 6.0, 7.0, 9.0 };
 * size_t input_size = sizeof(input) / sizeof(input[0]);
 * double* decompr;
 * size_t decompr_size;
 * pcomp_polycomp_chunk_t** chunks;
 * size_t num_of_chunks;
 * pcomp_polycomp_t* params;
 * size_t idx;
 *
 * params = pcomp_init_polycomp(4, 2, 1.0e-5, PCOMP_ALG_USE_CHEBYSHEV);
 * pcomp_compress_polycomp(&chunks, &num_of_chunks, input, input_size,
 *                         params);
 *
 * // Print some information for each chunk
 * for(idx = 0; idx < num_of_chunks; ++idx) {
 *     printf("Chunk %lu of %lu: %s\n", idx + 1, num_of_chunks,
 *            pcomp_chunk_is_compressed(chunks[idx]) ?
 *                "compressed" : "uncompressed");
 * }
 * \endcode
 *
 * Once the sequence \a input_buf is compressed, the array of chunks
 * can either be analyzed (e.g., using a \c for loop as in the example
 * above) or encoded using the \ref pcomp_encode_chunks. Once the
 * variable \a output_buf is no longer used, it should be freed via a
 * call to \ref pcomp_free_chunks.
 *
 * \param[out] output_buf Pointer to a variable that will receive the
 * address of an array of \ref pcomp_polycomp_chunk_t variables
 * created by the function. Such array contains the whole set of data
 * in \a input in compressed format. The array can be freed via a call
 * to \ref pcomp_free_chunks.
 *
 * \param[out] num_of_chunks On output, the variable will contain the
 * number of chunks saved in \a output_buf.
 *
 * \param[in] input_buf Pointer to the array of numbers to compress.
 *
 * \param[in] input_size Number of elements in \a input_buf to
 * compress.
 *
 * \param[in] params Parameters used for the compression. The variable
 * must have been created via a call to \ref pcomp_init_polycomp.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS if no error occurred, or
 * the error code.
 */
int pcomp_compress_polycomp(pcomp_polycomp_chunk_t** output_buf[],
                            size_t* num_of_chunks,
                            const double* input_buf, size_t input_size,
                            const pcomp_polycomp_t* params)
{
    size_t idx;
    pcomp_polycomp_t** chunk_params;
    pcomp_polycomp_t* last_chunk_params;
    size_t samples_in_last_chunk = 0;

    if (output_buf == NULL || num_of_chunks == NULL || input_buf == NULL
        || params == NULL || params->poly_fit == NULL)
        abort();

    /* Calculate how many chunks we'll create */
    *num_of_chunks = input_size / params->samples_per_chunk;
    samples_in_last_chunk = input_size % params->samples_per_chunk;
    if (samples_in_last_chunk != 0)
        ++(*num_of_chunks);

    *output_buf
        = malloc(sizeof(pcomp_polycomp_chunk_t*) * (*num_of_chunks));
    if (*output_buf == NULL)
        abort();

    /* Allocate a pcomp_polycomp_t structure for each of the OpenMP
     * threads */
    chunk_params
        = malloc(sizeof(pcomp_polycomp_t*) * omp_get_max_threads());
    for (idx = 0; idx < omp_get_max_threads(); ++idx) {
        chunk_params[idx] = pcomp_init_polycomp(
            params->samples_per_chunk, params->poly_fit->num_of_coeffs,
            params->max_allowable_error, params->algorithm);
    }

    if (samples_in_last_chunk != 0) {
        /* This is going to be used by just *one* OpenMP process */
        last_chunk_params = pcomp_init_polycomp(
            samples_in_last_chunk, params->poly_fit->num_of_coeffs,
            params->max_allowable_error, params->algorithm);
    }
    else {
        last_chunk_params = NULL;
    }

#pragma omp parallel for
    for (idx = 0; idx < *num_of_chunks; ++idx) {
        const double* cur_input = input_buf
                                  + params->samples_per_chunk * idx;
        pcomp_polycomp_t* cur_params;
        size_t cur_chunk_size;

        if (idx + 1 < *num_of_chunks || last_chunk_params == NULL) {
            cur_params = chunk_params[omp_get_thread_num()];
            cur_chunk_size = params->samples_per_chunk;
        }
        else {
            cur_params = last_chunk_params;
            cur_chunk_size = samples_in_last_chunk;
        }

        if (cur_params == NULL)
            abort();

        (*output_buf)[idx] = pcomp_init_chunk(cur_chunk_size);
        pcomp_run_polycomp_on_chunk(cur_params, cur_input,
                                    cur_chunk_size, (*output_buf)[idx],
                                    NULL);
    }

    for (idx = 0; idx < omp_get_max_threads(); ++idx)
        pcomp_free_polycomp(chunk_params[idx]);
    free(chunk_params);

    pcomp_free_polycomp(last_chunk_params);

    return PCOMP_STAT_SUCCESS;
}

/** \ingroup poly
 *
 * \brief Compute the sum of the number of samples encoded in \a
 *chunk_array
 *
 * \param[in] chunk_array Array of \ref pcomp_polycomp_chunk_t
 * variables. Typically, such array is created via a call to \ref
 * pcomp_compress_polycomp.
 *
 * \param[in] num_of_chunks Number of elements in \a chunk_array
 *
 * \returns The overall number of samples encoded in the sequence of
 * chunks
 */
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

/** \ingroup poly
 *
 * \brief Decompress a sequence of chunks
 *
 * This function is the counterpart for \ref pcomp_compress_polycomp.
 *
 * \param[out] output_buf Pointer to the variable that will hold the
 * uncompressed data. It must have room for a number of elements at
 * least equal to the return value of \ref pcomp_total_num_of_samples.
 *
 * \param[in] chunk_array Array of chunks holding the data in
 * compressed format.
 *
 * \param[in] num_of_chunks Number of elements in the array \a
 * chunk_array.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS if no error occurred, or
 * the error code.
 */
int pcomp_decompress_polycomp(
    double* output_buf, pcomp_polycomp_chunk_t* const chunk_array[],
    size_t num_of_chunks)
{
    size_t idx;
    pcomp_chebyshev_t** inv_chebyshev_array = NULL;
    size_t* start_pos;

    if (output_buf == NULL || chunk_array == NULL || num_of_chunks == 0)
        abort();

    /* Compute where the function should save the data after having
     * decompressed each chunk. Here we do *not* assume that all the
     * chunks but the last have the same size */
    start_pos = malloc(sizeof(size_t) * num_of_chunks);
    if (start_pos == NULL)
        abort();

    start_pos[0] = 0;
    for (idx = 0; idx < num_of_chunks - 1; ++idx) {
        start_pos[idx + 1] = start_pos[idx]
                             + chunk_array[idx]->num_of_samples;
    }

    /* "inv_chebyshev_array" is an array of pointers to
     * pcomp_chebyshev_t structures. We create them in advance so that
     * OpenMP threads can work without interferences among them. NULL
     * values are placed for those chunks whose size is equal to the
     * size of the first chunk: in this case, we can just reuse the
     * first pcomp_chebyshev_t array (the FFTW documentation says it's
     * safe to call fftw_execute at the same time from many threads,
     * see http://www.fftw.org/doc/Thread-safety.html). */
    inv_chebyshev_array
        = malloc(sizeof(pcomp_chebyshev_t) * num_of_chunks);
    if (inv_chebyshev_array == NULL)
        abort();
    for (idx = 0; idx < num_of_chunks; ++idx) {
        if (idx == 0 || (chunk_array[idx]->num_of_samples
                         != chunk_array[0]->num_of_samples)) {
            inv_chebyshev_array[idx] = pcomp_init_chebyshev(
                chunk_array[idx]->num_of_samples, PCOMP_TD_INVERSE);
        }
        else {
            inv_chebyshev_array[idx] = NULL;
        }
    }

/* Now decompress the chunks. The code is optimized for the case
 * where most of the chunks have the same size as the first, but
 * it works in the general case too */
#pragma omp parallel for
    for (idx = 0; idx < num_of_chunks; ++idx) {
        pcomp_chebyshev_t* cur_cheby;

        if (chunk_array[idx] == NULL)
            abort();

        if (inv_chebyshev_array[idx] != NULL)
            cur_cheby = inv_chebyshev_array[idx];
        else
            cur_cheby = inv_chebyshev_array[0];

        pcomp_decompress_polycomp_chunk(output_buf + start_pos[idx],
                                        chunk_array[idx], cur_cheby);
    }

    free(start_pos);
    for (idx = 0; idx < num_of_chunks; ++idx) {
        pcomp_free_chebyshev(inv_chebyshev_array[idx]);
    }
    free(inv_chebyshev_array);

    return PCOMP_STAT_SUCCESS;
}

/** \ingroup poly
 *
 * \brief Free an array of chunks
 *
 * \param[in] chunk_array An array of chunks. This variable must have
 * been allocated by a call to \ref pcomp_compress_polycomp.
 *
 * \param[in] num_of_chunks Number of elements in the array \a
 * chunk_array
 */
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

/** \ingroup poly
 *
 * \brief Number of bytes required by \ref pcomp_encode_chunks
 *
 * This function computes the number of bytes required to encode the
 * array of chunks in the variable \a chunks. Unlike functions like
 * \ref pcomp_rle_bufsize, this function provides an exact estimate,
 * not an upper bound.
 *
 * \param[in] chunks Array of chunks to encode. This should have been
 * initialized via a call to \ref pcomp_compress_polycomp.
 *
 * \param[in] num_of_chunks Number of elements in \a chunks.
 *
 * \returns The number of bytes required for the output buffer used by
 * \ref pcomp_encode_chunks.
 */
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

/** \ingroup poly
 *
 * \brief Encode a list of chunks into a sequence of raw bytes
 *
 * This function transforms an array of instances to \ref
 * pcomp_polycomp_chunk_t variables into a sequence of raw bytes,
 * suitable for I/O. It can be used together with \ref
 * pcomp_compress_polycomp to compress a dataset and save it into a
 * binary file.
 *
 * To decode byte sequences produced by this function, use \ref
 * pcomp_decode_chunks.
 *
 * \param[out] buf Pointer to a memory buffer that will receive the
 * result of the encoding. It must have room for a number of bytes (\c
 * uint8_t) at least equal to the return value of \ref
 * pcomp_chunks_num_of_bytes.
 *
 * \param[out] buf_size On exit, it will contain the number of bytes
 * actually written in \a buf. The latter number is equal to the value
 * returned by \ref pcomp_chunks_num_of_bytes.
 *
 * \param[in] chunk_array Array of chunks to encode
 *
 * \param[in] num_of_chunks Number of elements in the array \a
 * chunk_array.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS if no error occurred, or
 * the error code.
 */
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

/** \ingroup poly
 *
 * \brief Decode a byte sequence created by \ref pcomp_encode_chunks
 * into an array of chunks.
 *
 * This function can be used to read from a binary file or a socket a
 * sequence of chunks encoded by \ref pcomp_encode_chunks. The
 * function allocates memory for an array of \ref
 * pcomp_polycomp_chunk_t structures and returns it in the variable \a
 * chunk_array. The latter variable must be freed using \ref
 * pcomp_free_chunks once it is no longer needed.
 *
 * This function is the counterpart for \ref pcomp_encode_chunks.
 *
 * \param[out] chunk_array Pointer to an array that will contain the
 * chunks decoded from \a buf.
 *
 * \param[out] num_of_chunks On exit, this variable will hold the
 * number of chunks saved in \a chunk_array.
 *
 * \param[in] Pointer to the byte sequence to decode.
 *
 * \returns Either \ref PCOMP_STAT_SUCCESS if no error occurred, or
 * the error code.
 */
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
