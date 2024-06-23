#include <stdlib.h> // malloc
#include <stdint.h>
#include <string.h> // memset

#include "image.h"
#include "bmp_io.h"

static uint8_t convert_from_pixel_to_byte(pixel_t pixel_value) {
    // TODO: Change this when pixel_t type changes to float or fixed point.
    return pixel_value;
}

static pixel_t convert_from_byte_to_pixel(uint8_t pixel_byte) {
    return pixel_byte;
}


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
    // Note that we need to multiply by num of components...

    // Original image info
    image_info->cols = width;
    image_info->rows = height;
    image_info->num_components = planes;

    // Boundary extension info
    image_info->border = border;
    image_info->stride = stride;

    // Allocate memory to hold bmp image data
    pixel_t *handle = malloc(
        stride * (height + 2*border) * planes * sizeof(pixel_t));
    image_info->handle = handle;
    image_info->buf = handle + (((border*stride) + border) * planes);

    // Read every available line into image_info
    // Note: Rows that have already been accessed won't be read.
    pixel_t *current_line_ptr = image_info->buf; // start of actual image
    uint8_t *line_from_bmp = malloc(width * planes);
    while(bmp_info->num_unaccessed_rows > 0) {
        // Read line
        memset(line_from_bmp, 0, width);
        int fileReadErr = read_bmp_line(bmp_info, line_from_bmp);
        if (fileReadErr != 0) {
            return fileReadErr;
        }
        
        // Convert pixels in line from byte (uint8_t) to pixel_t
        for (int col = 0; col < width; col++) {
            if (planes == 3) {
                current_line_ptr[col*planes] = 
                    convert_from_byte_to_pixel(line_from_bmp[col*planes]);
                current_line_ptr[col*planes+1] =
                    convert_from_byte_to_pixel(line_from_bmp[col*planes+1]);
                current_line_ptr[col*planes+2] = 
                    convert_from_byte_to_pixel(line_from_bmp[col*planes+2]);
            } else {
                current_line_ptr[col] = 
                    convert_from_byte_to_pixel(line_from_bmp[col]);
            }
        }

        // Next line of image
        current_line_ptr += stride * planes;
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
    pixel_t *current_line_ptr = image_info->buf;
    uint8_t *line_to_bmp = malloc(width * planes);
    while (output_bmp.num_unaccessed_rows > 0) {
        // Convert pixels in current line from pixel_t to byte (uint8_t)
        for (int col = 0; col < width; col++) {
            if (planes == 3) {
                line_to_bmp[col*planes] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes]);
                line_to_bmp[col*planes+1] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes+1]);
                line_to_bmp[col*planes+2] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes+2]);
            } else {
                line_to_bmp[col] = 
                    convert_from_pixel_to_byte(current_line_ptr[col]);
            }
        }

        // Write line
        int fileWriteErr = write_bmp_line(&output_bmp, line_to_bmp);
        if (fileWriteErr != 0) {
            return fileWriteErr;
        }

        // Next line of image
        current_line_ptr += stride * planes;
    }

    // Close bmp file
    close_bmp(&output_bmp);

    return 0;
}

int export_image_and_border_as_bmp(image *image_info, const char *fname) {
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int planes = image_info->num_components;
    const int stride = image_info->stride;
    const int border = image_info->border;

    // Create new bitmap file
    bmp output_bmp;
    create_bmp(&output_bmp, fname, stride, height+2*border, planes);  //change width and height

    // Write pixel data to bitmap file until file is full
    pixel_t *current_line_ptr = image_info->handle; // change from buf to handle
    uint8_t *line_to_bmp = malloc(stride * planes);  // change from width to stride
    while (output_bmp.num_unaccessed_rows > 0) {
        // Convert pixels in current line from pixel_t to byte (uint8_t)
        for (int col = 0; col < width; col++) {
            if (planes == 3) {
                line_to_bmp[col*planes] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes]);
                line_to_bmp[col*planes+1] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes+1]);
                line_to_bmp[col*planes+2] = 
                    convert_from_pixel_to_byte(current_line_ptr[col*planes+2]);
            } else {
                line_to_bmp[col] = 
                    convert_from_pixel_to_byte(current_line_ptr[col]);
            }
        }

        // Write line
        int fileWriteErr = write_bmp_line(&output_bmp, line_to_bmp);
        if (fileWriteErr != 0) {
            return fileWriteErr;
        }

        // Next line of image
        current_line_ptr += stride * planes;
    }

    // Close bmp file
    close_bmp(&output_bmp);

    return 0;
}


//////////////////////
// Image processing //
//////////////////////

void perform_boundary_extension(image *image_info) {
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int stride = image_info->stride;
    const int border = image_info->border;
    const int planes = image_info->num_components;

    // Extend upwards
    pixel_t *first_line = image_info->buf;
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int pixel_index = (-row * stride + col) * planes;
            if (planes == 3) {
                first_line[pixel_index] = 0;
                first_line[pixel_index+1] = 0;
                first_line[pixel_index+2] = 0;
            } else {
                first_line[pixel_index] = 0;
            }
            
        }
    }

    // Extend downwards
    pixel_t *last_line = image_info->buf + ((height-1)*stride) * planes;
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int pixel_index = (row * stride + col) * planes;
            if (planes == 3) {
                last_line[pixel_index] = 0;
                last_line[pixel_index+1] = 0;
                last_line[pixel_index+2] = 0;
            } else {
                last_line[pixel_index] = 0;
            }
        }
    }

    // Extend rows to the left and right (not working as intended)
    pixel_t *left_edge = image_info->buf - ((border*stride) * planes);
    pixel_t *right_edge = left_edge + ((width-1) * planes);
    for (int row = height + 2*border; row > 0; row--) {
        for (int col = 1; col <= border; col++) {
            int left_index = -col*planes;
            int right_index = col*planes;
            if (planes == 3) {
                left_edge[left_index] = 0;
                left_edge[left_index+1] = 0;
                left_edge[left_index+2] = 0;

                right_edge[right_index] = 0;
                right_edge[right_index+1] = 0;
                right_edge[right_index+2] = 0;
            } else {
                left_edge[left_index] = 0;
                right_edge[right_index] = 0;
            }
        }

        left_edge += stride*planes;
        right_edge += stride*planes;
    }

}
