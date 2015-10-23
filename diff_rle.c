/* diff_rle.c - Differenced run-length encoding
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

/** \defgroup diffRLE Differenced Run-Length Encoding
 *
 * ### The algorithm and its applicability
 *
 * The differenced RLE compression is like RLE encoding, but it is
 * applied to the stream of consecutive differences of the original
 * stream. The first value of the input stream is prepended to the
 * sequence of run-length encoded values, in order to allow proper
 * decompression.
 *
 * For instance, given the following sequence:
 * \verbatim 10 11 12 13 15 17 19 20 21 \endverbatim
 * the compressed stream is
 * \verbatim 10  1 3  2 3  1 2 \endverbatim
 *
 * The main domain of applicability of this variant of RLE is for
 * integer quantities that increase regularly, like the number of
 * clock ticks in a digital chronometer.
 *
 * ### Implementation details
 *
 * As for RLE routines like \ref pcomp_compress_rle_int8, the caller
 * is expected to pre-allocate the memory which will contain the
 * output. The routine \ref pcomp_diffrle_bufsize is meant to provide
 * an upper estimate for the number of *elements* that will be saved
 * in the compressed stream.
 */

/** \ingroup diffRLE
 *
 *\brief Calculate an upper limit for the size of a buffer holding
 * streams encoded using differenced RLE.
 *
 * Return the number of elements required for a buffer used to hold
 * the RLE-compressed version of a datastream. It is typically used
 * together with functions like \ref pcomp_compress_diffrle_int8 to
 * pre-allocate the buffer that will contain the result.
 *
 * The number returned by this function might be severely
 * overestimated. Therefore, when using this function to allocate a
 * buffer for functions like \ref pcomp_compress_diffrle_int8, it is
 * recommended to claim any unused space at the end of the buffer.
 * (Functions like \ref pcomp_compress_diffrle_int8 inform the caller
 * about the number of elements actually written in the buffer.)
 *
 * \param[in] input_size Number of elements of the data stream to
 *compress
 * \returns The number of *elements* (not bytes) of the buffer
 */
size_t pcomp_diffrle_bufsize(size_t input_size)
{
    return 1 + 2 * input_size;
}

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_int8(int8_t* output_buf,
 *                                     size_t* output_size,
 *                                     const int8_t* input_buf,
 *                                     size_t input_size)
 *
 * \brief Compress an array of int8_t values using the diffRLE
 * compression
 *
 * The size of the output buffer is typically guessed using \ref
 * pcomp_diffrle_bufsize. After the call to the function, the buffer
 * should be resized to claim unused space at the end. Here is an
 * example:
 *
 * \code{.c}
 * int8_t input_buf[] = { 10, 10, 20, 30, 30, 30 };
 * size_t input_size = sizeof(input_buf) / sizeof(input_buf[0]);
 * size_t output_size = pcomp_diffrle_bufsize(input_size) *
 *                      sizeof(int8_t);
 * int8_t* output_buf = malloc(output_size);
 *
 * pcomp_compress_diffrle_int8(output_buf, &output_size, input_buf,
 *                             input_size);
 *
 * output_buf = realloc(output_buf, output_size * sizeof(int8_t));
 * \endcode
 *
 * \param[out] output_buf The buffer that will hold the compressed
 * stream. It must have space for a number of elements at least equal
 * to the value returned by pcomp_diffrle_bufsize.
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
 * \returns PCOMP_STAT_SUCCESS if the encoding completed successfully.
 * Otherwise, the error code specifies the kind of error occurred
 * during the call.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_int16(int16_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int16_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int16_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_int32(int32_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int32_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int32_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_int64(int64_t* output_buf,
 *                                  size_t* output_size,
 *                                  const int64_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of int64_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_uint8(uint8_t* output_buf,
 *                                 size_t* output_size,
 *                                 const uint8_t* input_buf,
 *                                 size_t input_size)
 *
 * \brief Compress an array of uint8_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_uint16(uint16_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint16_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint16_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_uint32(uint32_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint32_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint32_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 *more
 * information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_compress_diffrle_uint64(uint64_t* output_buf,
 *                                  size_t* output_size,
 *                                  const uint64_t* input_buf,
 *                                  size_t input_size)
 *
 * \brief Compress an array of uint64_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_compress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_int8(int8_t* output_buf,
 *                                     size_t output_size,
 *                                     const int8_t* input_buf,
 *                                     size_t input_size)
 *
 * \brief Decompress an array of int8_t values encoded using the
 * diffRLE compression
 *
 * The size of the output buffer must have been saved somewhere
 * instead with the compressed data. It is currently not possible to
 * retrieve it directly from \a input_buf.
 *
 * \param[out] output_buf The buffer that will hold the decompressed
 * stream. The space required for the decompressed chunks must have
 * been already allocated.
 *
 * \param[in] output_size Number of *elements* (not bytes) in the
 * decompressed stream \a output_buf.
 *
 * \param[in] input_buf The array of values to compress
 *
 * \param[in] input_size Number of *elements* (not bytes) in the array
 * \a input_buf
 *
 * \returns PCOMP_STAT_SUCCESS if the decoding completed successfully.
 * Otherwise, the error code specifies the kind of error occurred
 * during the call.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_int16(int16_t* output_buf,
 *                                        size_t output_size,
 *                                        const int16_t* input_buf,
 *                                        size_t input_size)
 *
 * \brief Compress an array of int16_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_int32(int32_t* output_buf,
 *                                        size_t output_size,
 *                                        const int32_t* input_buf,
 *                                        size_t input_size)
 *
 * \brief Compress an array of int32_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_int64(int64_t* output_buf,
 *                                        size_t output_size,
 *                                        const int64_t* input_buf,
 *                                        size_t input_size)
 *
 * \brief Compress an array of int64_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_uint8(uint8_t* output_buf,
 *                                        size_t output_size,
 *                                        const uint8_t* input_buf,
 *                                        size_t input_size)
 *
 * \brief Compress an array of uint8_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_uint16(uint16_t* output_buf,
 *                                         size_t output_size,
 *                                         const uint16_t* input_buf,
 *                                         size_t input_size)
 *
 * \brief Compress an array of uint16_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_uint32(uint32_t* output_buf,
 *                                         size_t output_size,
 *                                         const uint32_t* input_buf,
 *                                         size_t input_size)
 *
 * \brief Compress an array of uint32_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 *more
 * information.
 */

/** \ingroup diffRLE
 *
 * \fn int pcomp_decompress_diffrle_uint64(uint64_t* output_buf,
 *                                         size_t output_size,
 *                                         const uint64_t* input_buf,
 *                                         size_t input_size)
 *
 * \brief Compress an array of uint64_t values using the diffRLE
 * compression
 *
 * Refer to the documentation for \ref pcomp_decompress_diffrle_int8 for
 * more information.
 */

/***********************************************************************
 * Differenced run-length compression routines
 */

#define IMPLEMENT_DIFFRLE_COMPR_FN(name, rle_compr_fn, datatype_t)     \
    int name(datatype_t* output_buf, size_t* output_size,              \
             const datatype_t* input_buf, size_t input_size)           \
    {                                                                  \
        size_t rle_output_size;                                        \
        datatype_t* diff_buf = NULL;                                   \
        int rle_compr_result;                                          \
        size_t idx;                                                    \
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
        output_buf[0] = input_buf[0];                                  \
        if (input_size == 1) {                                         \
            *output_size = 1;                                          \
            return PCOMP_STAT_SUCCESS;                                 \
        }                                                              \
                                                                       \
        diff_buf = malloc((input_size - 1) * sizeof(datatype_t));      \
        for (idx = 0; idx < input_size - 1; ++idx) {                   \
            diff_buf[idx] = input_buf[idx + 1] - input_buf[idx];       \
        }                                                              \
                                                                       \
        rle_output_size = *output_size - 1;                            \
        rle_compr_result                                               \
            = rle_compr_fn(output_buf + 1, &rle_output_size, diff_buf, \
                           input_size - 1);                            \
        free(diff_buf);                                                \
        if (rle_compr_result != PCOMP_STAT_SUCCESS)                    \
            return rle_compr_result;                                   \
                                                                       \
        *output_size = rle_output_size + 1;                            \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int8,
                           pcomp_compress_rle_int8, int8_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int16,
                           pcomp_compress_rle_int16, int16_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int32,
                           pcomp_compress_rle_int32, int32_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_int64,
                           pcomp_compress_rle_int64, int64_t)

IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint8,
                           pcomp_compress_rle_uint8, uint8_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint16,
                           pcomp_compress_rle_uint16, uint16_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint32,
                           pcomp_compress_rle_uint32, uint32_t)
IMPLEMENT_DIFFRLE_COMPR_FN(pcomp_compress_diffrle_uint64,
                           pcomp_compress_rle_uint64, uint64_t)

/***********************************************************************
 * Differenced run-length decompression routines
 */

#define IMPLEMENT_DIFFRLE_DECOMPR_FN(name, datatype_t)                 \
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
        if (input_size % 2 != 1) {                                     \
            return PCOMP_STAT_INVALID_ENCODING;                        \
        }                                                              \
                                                                       \
        output_buf[output_idx++] = input_buf[input_idx++];             \
                                                                       \
        while (output_idx < output_size                                \
               && input_idx < input_size - 1) {                        \
            datatype_t count = input_buf[input_idx];                   \
            datatype_t incr = input_buf[input_idx + 1];                \
            datatype_t idx;                                            \
                                                                       \
            for (idx = 0; idx < count; ++idx) {                        \
                output_buf[output_idx] = output_buf[output_idx - 1]    \
                                         + incr;                       \
                output_idx++;                                          \
            }                                                          \
                                                                       \
            input_idx += 2;                                            \
        }                                                              \
                                                                       \
        return PCOMP_STAT_SUCCESS;                                     \
    }

IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int8, int8_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int16, int16_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int32, int32_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_int64, int64_t)

IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint8, uint8_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint16, uint16_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint32, uint32_t)
IMPLEMENT_DIFFRLE_DECOMPR_FN(pcomp_decompress_diffrle_uint64, uint64_t)
