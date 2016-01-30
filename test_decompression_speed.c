/* This program has been designed to test the speed of the polynomial
 * compression. It is *not* included in the test suite, as it does not
 * check the results. If you want to hack libpolycomp's source code to
 * make it faster, use this program as a benchmark. */

#include "libpolycomp.h"
#include <stdio.h>
#include <sys/time.h>
#include <fitsio.h>
#include <stdlib.h>

#ifdef WITH_OPENMP

#include <omp.h>

#else

typedef int omp_int_t;
static omp_int_t omp_get_max_threads(void) { return 1; }

#endif

static struct timeval tstart, tstop;

void start_timer() { gettimeofday(&tstart, NULL); }

double stop_timer()
{
    gettimeofday(&tstop, NULL);
    return (tstop.tv_sec + tstop.tv_usec / 1000000.0)
           - (tstart.tv_sec + tstart.tv_usec / 1000000.0);
}

void* read_file_data(const char* file_name)
{
    FILE* f;
    long file_size;
    void* buf;

    f = fopen(file_name, "rb");
    if (f == NULL) {
        fprintf(stderr, "error opening file \"%s\": ", file_name);
        perror(NULL);
        exit(1);
    }

    fseek(f, 0, SEEK_END);
    file_size = ftell(f);
    fseek(f, 0, SEEK_SET);

    buf = malloc(file_size);
    if (fread(buf, 1, file_size, f) < file_size) {
        fprintf(stderr, "error reading file \"%s\": ", file_name);
        perror(NULL);
        exit(1);
    }

    fclose(f);

    return buf;
}

void save_data_in_fits_file(const char* file_name, double* data,
                            size_t num_of_elements)
{
    fitsfile* file;
    int status = 0;
    char* ttype[] = { "SAMPLE" };
    char* tform[] = { "1D" };
    char* tunit[] = { "" };

    if (fits_create_file(&file, file_name, &status) != 0)
        goto error;

    if (fits_create_tbl(file, BINARY_TBL, 0, 1, ttype, tform, tunit,
                        "DECOMPR", &status) != 0)
        goto error;

    if (fits_write_col(file, TDOUBLE, 1, 1, 1, num_of_elements, data,
                       &status) != 0)
        goto error;

    fits_close_file(file, &status);
    return;

error:
    fits_report_error(stderr, status);
    if (file != NULL)
        fits_close_file(file, &status);
    exit(1);
}

int main(int argc, const char* argv[])
{
    void* buf;
    double* decompr;
    size_t num_of_elements;
    size_t num_of_chunks;
    pcomp_polycomp_chunk_t** chunks;
    double duration;

    if (argc != 3) {
        fputs("usage: test_decompression_speed BINARY_FILE "
              "OUTPUT_FITS_FILE\n",
              stderr);
        return 1;
    }

#if WITH_OPENMP
    printf("OpenMP enabled: 1\n");
    printf("number of OpenMP threads: %d\n", omp_get_max_threads());
#else
    printf("OpenMP enabled: 0\n");
#endif
    printf("input file name: %s\n", argv[1]);
    printf("output file name: %s\n", argv[2]);

    buf = read_file_data(argv[1]);
    start_timer();
    pcomp_decode_chunks(&chunks, &num_of_chunks, buf);
    duration = stop_timer();
    printf("decoding elapsed time (s): %lf\n", duration);

    num_of_elements = pcomp_total_num_of_samples(chunks, num_of_chunks);
    decompr = malloc(sizeof(double) * num_of_elements);

    start_timer();
    pcomp_decompress_polycomp(decompr, chunks, num_of_chunks);
    duration = stop_timer();
    printf("decompression elapsed time (s): %lf\n", duration);

    save_data_in_fits_file(argv[2], decompr, num_of_elements);
    return 0;
}
