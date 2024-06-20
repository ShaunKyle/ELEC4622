#include <stdlib.h> // malloc
#include <stdint.h>

#include "image.h"
#include "bmp_io.h"

//! \brief Reads image data from an open bitmap file.
//!
//! Only unaccessed rows of the bmp will be read. The bmp file handle will be 
//! closed by this function.
//!
//! @param image_info   Image struct to store pixel data and image info
//! @param bmp_info     Bitmap to read
//!
//! @return Bitmap IO_ERR
int read_image_from_bmp(image *image_info, bmp *bmp_info) {
    const int width = bmp_info->cols;
    const int height = bmp_info->rows;
    const int planes = bmp_info->num_components;

    // Populate image struct
    image_info->cols = width;
    image_info->rows = height;
    image_info->num_components = planes;

    // Allocate memory to hold bmp image data
    image_info->data = malloc(width * height * planes);

    // Read every available line into image_info
    // Note: Rows that have already been accessed won't be read.
    uint8_t *current_line_ptr = image_info->data; // start of pixel data
    while(bmp_info->num_unaccessed_rows > 0) {
        int fileReadErr = read_bmp_line(bmp_info, current_line_ptr);
        if (fileReadErr != 0) {
            // print_bmp_file_error(fileReadErr);
            return fileReadErr;
        }
        current_line_ptr += width * planes; // next line of image_info
    }

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

    // Create new bitmap file
    bmp output_bmp;
    create_bmp(&output_bmp, fname, width, height, planes);

    // Write pixel data to bitmap file until file is full
    uint8_t *current_line_ptr = image_info->data;
    while (output_bmp.num_unaccessed_rows > 0) {
        int fileWriteErr = write_bmp_line(&output_bmp, current_line_ptr);
        if (fileWriteErr != 0) {
            // print_bmp_file_error(fileWriteErr);
            return fileWriteErr;
        }
        current_line_ptr += width * planes;
    }

    // Close bmp file
    close_bmp(&output_bmp);

    return 0;
}
