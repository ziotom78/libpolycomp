/* This program has been designed to test the speed of the polynomial
 * compression. It is *not* included in the test suite, as it does not
 * check the results. If you want to hack libpolycomp's source code to
 * make it faster, use this program as a benchmark. */

#include "libpolycomp.h"
#include <stdio.h>
#include <sys/time.h>
#include <fitsio.h>
#include <stdlib.h>

static struct timeval tstart, tstop;

void start_timer() { gettimeofday(&tstart, NULL); }

double stop_timer()
{
    gettimeofday(&tstop, NULL);
    return (tstop.tv_sec + tstop.tv_usec / 1000000.0)
           - (tstart.tv_sec + tstart.tv_usec / 1000000.0);
}

double* read_data_from_fits(const char* file_name,
                            size_t* num_of_elements)
{
    fitsfile* file = NULL;
    LONGLONG num_of_rows;
    int status = 0;
    double* data;

    if (fits_open_table(&file, file_name, READONLY, &status) != 0)
        goto error;

    if (fits_get_num_rowsll(file, &num_of_rows, &status) != 0)
        goto error;

    data = malloc(sizeof(double) * num_of_rows);
    if (fits_read_col(file, TDOUBLE, 1, 1, 1, num_of_rows, NULL, data,
                      NULL, &status) != 0)
        goto error;

    if (fits_close_file(file, &status) != 0)
        fits_report_error(stderr, status);

    *num_of_elements = (size_t)num_of_rows;
    return data;

error:
    fits_report_error(stderr, status);
    if (file != NULL)
        fits_close_file(file, &status);
    exit(1);
}

void compress_data(pcomp_chunk_size_t samples_per_chunk,
                   pcomp_poly_size_t num_of_coeffs, double max_error,
                   const double* input, size_t num_of_elements)
{
    pcomp_polycomp_chunk_t** chunks;
    size_t num_of_chunks;

    pcomp_polycomp_t* params
        = pcomp_init_polycomp(samples_per_chunk, num_of_coeffs,
                              max_error, PCOMP_ALG_USE_CHEBYSHEV);

    pcomp_compress_polycomp(&chunks, &num_of_chunks, input,
                            num_of_elements, params);

    pcomp_free_polycomp(params);
    pcomp_free_chunks(chunks, num_of_chunks);
}

int main(int argc, const char* argv[])
{
    double* input;
    size_t num_of_elements;
    pcomp_chunk_size_t samples_per_chunk;
    pcomp_poly_size_t num_of_coeffs;
    double max_error;
    double duration;

    if (argc != 5) {
        fputs("usage: test_speed FITS_FILE SAMPLES_PER_CHUNK "
              "NUM_COEFFS MAX_ERROR\n\n",
              stderr);
        fputs("error: you must supply the name of a FITS file\n"
              "containing a tabular HDU where to read the data\n"
              "to be used in the test\n",
              stderr);
        return 1;
    }

    input = read_data_from_fits(argv[4], &num_of_elements);
    samples_per_chunk = atoi(argv[1]);
    num_of_coeffs = atoi(argv[2]);
    max_error = atof(argv[3]);

    printf("input file name: %s\n", argv[4]);
    printf("number of samples: %lu\n", num_of_elements);
    printf("samples per chunk: %u\n", samples_per_chunk);
    printf("number of coefficients: %u\n", num_of_coeffs);
    printf("maximum error: %e\n", max_error);

    start_timer();
    compress_data(samples_per_chunk, num_of_coeffs, max_error, input,
                  num_of_elements);
    duration = stop_timer();
    printf("elapsed time (s): %lf\n", duration);

    free(input);

    return 0;
}
