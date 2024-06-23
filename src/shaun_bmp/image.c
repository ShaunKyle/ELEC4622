#include <stdlib.h> // malloc
#include <stdint.h>

#include "image.h"
#include "bmp_io.h"

//! \brief Reads image data from an open bitmap file.
//!
//! Only unaccessed rows of the bmp will be read. The bmp file handle will be 
//! closed by this function.
//!
//! Pixels in the boundary extension border are currently uninitialized.
//!
//! @param image_info   Image struct to store pixel data and image info
//! @param bmp_info     Bitmap to read
//! @param border       Width of boundary extension in pixels
//!
//! @return Bitmap IO_ERR
int read_image_from_bmp(image *image_info, bmp *bmp_info, int border) {
    const int width = bmp_info->cols;
    const int height = bmp_info->rows;
    const int planes = bmp_info->num_components;
    const int stride = width + 2*border;

    // Original image info
    image_info->cols = width;
    image_info->rows = height;
    image_info->num_components = planes;

    // Boundary extension info
    image_info->border = border;
    image_info->stride = stride;

    // Allocate memory to hold bmp image data
    uint8_t *handle = malloc(stride * (height + 2*border) * planes);
    image_info->handle = handle;
    image_info->buf = handle + ((border*stride) + border) * planes;

    // Read every available line into image_info
    // Note: Rows that have already been accessed won't be read.
    uint8_t *current_line_ptr = image_info->buf; // start of actual image
    while(bmp_info->num_unaccessed_rows > 0) {
        int fileReadErr = read_bmp_line(bmp_info, current_line_ptr);
        if (fileReadErr != 0) {
            // print_bmp_file_error(fileReadErr);
            return fileReadErr;
        }
        current_line_ptr += stride * planes; // next line of image
    }

    // TODO: Boundary extension (zero padding by default?)

    // TODO: a stripe is added to the bottom of 8-bit images. Why?

    // Close bmp file
    close_bmp(bmp_info);

    return 0;
}

//! \brief Writes image data to a new bitmap file
//!
//! @param image_info   Image struct with pixel data to create bmp from
//! @param fname        Name of bitmap file to output
//!
//! @return Bitmap IO_ERR
int export_image_as_bmp(image *image_info, const char *fname) {
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int planes = image_info->num_components;
    const int stride = image_info->stride;

    // Create new bitmap file
    bmp output_bmp;
    create_bmp(&output_bmp, fname, width, height, planes);

    // Write pixel data to bitmap file until file is full
    uint8_t *current_line_ptr = image_info->buf;
    while (output_bmp.num_unaccessed_rows > 0) {
        int fileWriteErr = write_bmp_line(&output_bmp, current_line_ptr);
        if (fileWriteErr != 0) {
            // print_bmp_file_error(fileWriteErr);
            return fileWriteErr;
        }
        current_line_ptr += stride * planes;
    }

    // Close bmp file
    close_bmp(&output_bmp);

    return 0;
}
