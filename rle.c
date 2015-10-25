/* rle.c - Run-Length encoding
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

/** \defgroup RLE Run-Length Encoding
 *
 * ### The algorithm and its applicability
 *
 * Libpolycomp implements routines for compressing and decompressing
 * streams of data using the Run-Length Encoding (RLE) scheme. This
 * kind of compression is perfect for data streams which contain long
 * sequences of repeated data, e.g.,
 * \verbatim 1041 1041 1041 1041 1280 1280 1041 1041 1041 1041
 \endverbatim
 *
 * The RLE scheme works by encoding each sample together with the
 * number of consecutive repeats found. Therefore, for the previous
 * example the encoding would be
 * \verbatim 1041 4  1280 2  1041 4 \endverbatim
 * In certain cases, RLE can outperform other well-known compression
 * schemes like gzip and bzip2.
 *
 * The RLE scheme can be applied reliably only on sequences of integer
 * data. The algorithm involves the comparison between consecutive
 * values, and this cannot be done reliably with floating-point
 * numbers because of round-off errors. The Libpolycomp library
 * implements many similar functions (e.g., \ref
 * pcomp_compress_rle_int8) for signed and unsigned integers, with
 * sizes of 1, 2, 4, and 8 bytes.
 *
 * Libpolycomp correctly takes into account possible overflows.
 * Suppose for instance that the input datastream contains 1000
 * repetitions of the 8-bit value "152". In this case, Libpolycomp
 * produces the following output:
 * \verbatim 152 256  152 256  152 256  152 232 \endverbatim
 * as 1000 = 256 * 3 + 232.
 *
 * ### Implementation details
 *
 * The compression and decompression routines require the caller to
 * pre-allocate the memory which will contain the output. For
 * compression routines, one typically uses the function \ref
 * pcomp_rle_bufsize to get an upper estimate to the size of the
 * output buffer, and after the compression resizes the buffer (is
 * needed) in order to free the unused memory. See the documentation
 * for \ref pcomp_compress_rle_int8 and \ref pcomp_decompress_rle_int8
 * for examples.
 */

#include "libpolycomp.h"
#include <stdlib.h>

/** \ingroup RLE
 *
 *\brief Calculate an upper limit for the size of a buffer holding
 * RLE-encoded streams.
 *
 * Return the number of elements required for a buffer used to hold
 * the RLE-compressed version of a datastream. It is typically used
 * together with functions like \ref pcomp_compress_rle_int8 to
 * pre-allocate the buffer that will contain the result.
 *
 * The number returned by this function might be severely
 * overestimated. Therefore, when using this function to allocate a
 * buffer for functions like \ref pcomp_compress_rle_int8, it is
 * recommended to claim any unused space at the end of the buffer.
 * (Functions like \ref pcomp_compress_rle_int8 inform the caller
 * about the number of elements actually written in the buffer.)
 *
 * \param[in] input_size Number of elements of the data stream to
 *compress
 * \returns The number of *elements* (not bytes) of the buffer
 */
size_t pcomp_rle_bufsize(size_t input_size) { return 2 * input_size; }

/***********************************************************************
 * Run-length compression routines
 */

#define IMPLEMENT_RLE_COMPR_FN(name, datatype_t, max_value)            \
    int name(datatype_t* output_buf, size_t* output_size,              \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t true_output_size = 0;                                   \
        size_t input_idx = 0;                                          \
                                                                       \
        if (output_buf == NULL || output_size == NULL                  \
            || input_buf == NULL)                                      \
            abort();                                                   \
                                                                       \
        if (input_size == 0) {                                         \
            *output_size = 0;                                          \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        if (*output_size < 2) {                                        \
            return PCOMP_STAT_INVALID_BUFFER;                          \
        }                                                              \
                                                                       \
        while (input_idx < input_size) {                               \
            datatype_t first_datum_in_the_run = input_buf[input_idx];  \
            datatype_t count = 0;                                      \
            while (count < max_value && input_idx < input_size         \
                   && input_buf[input_idx]                             \
                          == first_datum_in_the_run) {                 \
                input_idx++;                                           \
                count++;                                               \
            }                                                          \
                                                                       \
            if (true_output_size >= *output_size - 2) {                \
                return PCOMP_STAT_INVALID_BUFFER;                      \
            }                                                          \
                                                                       \
            output_buf[true_output_size++] = count;                    \
            output_buf[true_output_size++] = first_datum_in_the_run;   \
        }                                                              \
                                                                       \
        *output_size = true_output_size;                               \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int8, int8_t, INT8_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int16, int16_t, INT16_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int32, int32_t, INT32_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_int64, int64_t, INT64_MAX)

IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint8, uint8_t, UINT8_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint16, uint16_t, UINT16_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint32, uint32_t, UINT32_MAX)
IMPLEMENT_RLE_COMPR_FN(pcomp_compress_rle_uint64, uint64_t, UINT64_MAX)

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_int8(int8_t* output_buf,
 *                                 size_t* output_size,
 *                                 const int8_t* input_buf,
 *                                 size_t input_size)
 *
 * \brief Compress an array of int8_t values using the RLE
 * compression
 *
 * The size of the output buffer is typically guessed using \ref
 * pcomp_rle_bufsize. After the call to the function, the buffer
 * should be resized to claim unused space at the end. Here is an
 * example:
 *
 * \code{.c}
 * int8_t input_buf[] = { 10, 10, 20, 30, 30, 30 };
 * size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);
 * size_t output_size = pcomp_rle_bufsize(input_size) *
 *                      sizeof(int8_t);
 * int8_t* output_buf = malloc(output_size);
 *
 * pcomp_compress_rle_int8(output_buf, &output_size, input_buf,
 *                         input_size);
 *
 * output_buf = realloc(output_buf, output_size * sizeof(int8_t));
 * \endcode
 *
 * \param[out] output_buf The buffer that will hold the compressed
 * stream. It must have space for a number of elements at least equal
 * to the value returned by pcomp_rle_bufsize.
 *
 * \param[inout] output_size On entering the function, this must
 * specify the number of *elements* (not bytes) that can be written in
 * \a output_buf. On exit, it will contain the actual number of
 * elements written.
 *
 * \param[in] input_buf The array of values to compress
 *
 * \param[in] input_size Number of *elements* (not bytes) in the array
 * \a input_buf
 *
 * \returns \ref PCOMP_STAT_SUCCESS if the encoding completed
 * successfully. Otherwise, the error code specifies the kind of error
 * occurred during the call.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_int16(int16_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int16_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int16_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_int32(int32_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int32_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int32_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_int64(int64_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int64_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int64_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_uint8(uint8_t* output_buf,
 *                                 size_t* output_size,
 *                                 const uint8_t* input_buf,
 *                                 size_t input_size)
 *
 * \brief Compress an array of uint8_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_uint16(uint16_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint16_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint16_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_uint32(uint32_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint32_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint32_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for more
 * information.
 */

/** \ingroup RLE
 *
 * \fn int pcomp_compress_rle_uint64(uint64_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint64_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint64_t values using the RLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_rle_int8 for
 * more information.
 */

/***********************************************************************
 * Run-Length decompression routines
 */

#define IMPLEMENT_RLE_DECOMPR_FN(name, datatype_t)                     \
    int name(datatype_t* output_buf, size_t output_size,               \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t input_idx = 0;                                          \
        size_t output_idx = 0;                                         \
                                                                       \
        if (output_buf == NULL || input_buf == NULL)                   \
            abort();                                                   \
                                                                       \
        if (input_size == 0) {                                         \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        if (input_size % 2 != 0) {                                     \
            return PCOMP_STAT_INVALID_ENCODING;                        \
        }                                                              \
                                                                       \
        while (output_idx < output_size                                \
               && input_idx < input_size - 1) {                        \
            datatype_t count = input_buf[input_idx];                   \
            datatype_t value = input_buf[input_idx + 1];               \
            datatype_t idx;                                            \
                                                                       \
            for (idx = 0; idx < count; ++idx) {                        \
                output_buf[output_idx++] = value;                      \
            }                                                          \
                                                                       \
            input_idx += 2;                                            \
        }                                                              \
                                                                       \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int8, int8_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int16, int16_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int32, int32_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_int64, int64_t)

IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint8, uint8_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint16, uint16_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint32, uint32_t)
IMPLEMENT_RLE_DECOMPR_FN(pcomp_decompress_rle_uint64, uint64_t)
