#include <stdlib.h> // malloc
#include <stdint.h>
#include <string.h> // memset

#include "image.h"
#include "bmp_io.h"

static uint8_t convert_from_pixel_to_byte(pixel_t pixel_value) {
    const int max_byte_value = 255;
    uint8_t byte_value = (uint8_t) (pixel_value * max_byte_value);
    return byte_value;
}

static pixel_t convert_from_byte_to_pixel(uint8_t pixel_byte) {
    const int max_byte_value = 255;
    pixel_t float_value = (float) pixel_byte / max_byte_value;
    return float_value;
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

//! \brief Writes image data, including border extensions, to a new bitmap file
//!
//! Useful for debugging and illustrating border extension methods.
//!
//! @param image_info   Image struct with pixel data to create bmp from
//! @param fname        Name of bitmap file to output
//!
//! @return Bitmap IO_ERR
int export_image_and_border_as_bmp(image *image_info, const char *fname) {
    // We use stride instead of width, to include the border.
    // const int width = image_info->cols;
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
        for (int col = 0; col < stride; col++) {
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

//! \brief Symmetric odd extension of image boundaries
//!
//! The thickness of border is determined by image_info->border.
void perform_boundary_extension(image *image_info) {
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int stride = image_info->stride;
    const int border = image_info->border;
    const int planes = image_info->num_components;

    // Precondition: Border thickness should be less than image dimensions
    if ((border >= width) || (border >= height)) {
        fprintf(stderr,
            "Image too small/border too thick for symmetric extension.\n");
        // TODO: return error?
        return;
    }

    // Extend upward
    // x[-n1, -n2] = x[n1, n2]
    pixel_t *first_line = image_info->buf;  // Upper boundary (row 0)
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int border_index = (-row * stride + col) * planes;  // row -n
            int image_index = (row * stride + col) * planes;    // row +n
            if (planes == 3) {
                first_line[border_index] = first_line[image_index];
                first_line[border_index+1] = first_line[image_index+1];
                first_line[border_index+2] = first_line[image_index+2];
            } else {
                first_line[border_index] = first_line[image_index];
            }
        }
    }

    // Extend downward
    // x[(N-1)+n1, n2] = x[(N-1)-n1, n2]
    // last_line is row (N-1)
    pixel_t *last_line = image_info->buf + ((height-1)*stride) * planes;
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int border_index = (row * stride + col) * planes;
            int image_index = (-row * stride + col) * planes;
            if (planes == 3) {
                last_line[border_index] = last_line[image_index];
                last_line[border_index+1] = last_line[image_index+1];
                last_line[border_index+2] = last_line[image_index+2];
            } else {
                last_line[border_index] = last_line[image_index];
            }
        }
    }

    // Extend rows left and right
    pixel_t *left_edge = image_info->buf - ((border*stride) * planes);
    pixel_t *right_edge = left_edge + ((width-1) * planes);
    for (int row = 0; row < height + 2*border; row++) {
        for (int col = 1; col <= border; col++) {
            int left_index = (-col + row*stride) * planes;
            int right_index = (col + row*stride) * planes;
            if (planes == 3) {
                left_edge[left_index] = left_edge[right_index];
                left_edge[left_index+1] = left_edge[right_index+1];
                left_edge[left_index+2] = left_edge[right_index+2];
                right_edge[right_index] = right_edge[left_index];
                right_edge[right_index+1] = right_edge[left_index+1];
                right_edge[right_index+2] = right_edge[left_index+2];
            } else {
                left_edge[left_index] = left_edge[right_index];
                right_edge[right_index] = right_edge[left_index];
            }
        }
    }
}

void perform_boundary_extension_zero_padding(image *image_info) {
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int stride = image_info->stride;
    const int border = image_info->border;
    const int planes = image_info->num_components;

    const int pad = 0;  // Constant value to pad border with
                        // TODO: symmetric extension instead

    // Extend upwards
    pixel_t *first_line = image_info->buf;
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int pixel_index = (-row * stride + col) * planes;
            if (planes == 3) {
                first_line[pixel_index] = pad;
                first_line[pixel_index+1] = pad;
                first_line[pixel_index+2] = pad;
            } else {
                first_line[pixel_index] = pad;
            }
        }
    }

    // Extend downwards
    pixel_t *last_line = image_info->buf + ((height-1)*stride) * planes;
    for (int row = 1; row <= border; row++) {
        for (int col = 0; col < width; col++) {
            int pixel_index = (row * stride + col) * planes;
            if (planes == 3) {
                last_line[pixel_index] = pad;
                last_line[pixel_index+1] = pad;
                last_line[pixel_index+2] = pad;
            } else {
                last_line[pixel_index] = pad;
            }
        }
    }

    // Extend rows to the left and right
    pixel_t *left_edge = image_info->buf - ((border*stride) * planes);
    pixel_t *right_edge = left_edge + ((width-1) * planes);
    for (int row = height + 2*border; row > 0; row--) {
        for (int col = 1; col <= border; col++) {
            int left_index = -col*planes;
            int right_index = col*planes;
            if (planes == 3) {
                left_edge[left_index] = pad;
                left_edge[left_index+1] = pad;
                left_edge[left_index+2] = pad;

                right_edge[right_index] = pad;
                right_edge[right_index+1] = pad;
                right_edge[right_index+2] = pad;
            } else {
                left_edge[left_index] = pad;
                right_edge[right_index] = pad;
            }
        }

        left_edge += stride*planes;
        right_edge += stride*planes;
    }

}

//! \brief Convolves an image with a point spread function (PSF)
//!
//! "Extent" refers to the number of elements extending from the center point 
//! in a direction. (e.g. a 3x3 PSF has an extent of 1 in both horizontal and 
//! vertical directions).
//!
//! @param image_in             Input image to convolve
//! @param image_out            Result of convolution
//! @param psf_values           Array of PSF values
//! @param extent_horizontal    Horizontal extent of the PSF
//! @param extent_vertical      Vertical extent of the PSF
void apply_filter(image *image_in, image *image_out, pixel_t *psf_values, 
int extent_horizontal, int extent_vertical) {
    // Image constants
    const int height = image_in->rows;
    const int width = image_in->cols;
    const int planes = image_in->num_components;
    const int stride = image_in->stride;

    // Filter constants
    const int extent_cols = extent_horizontal;
    const int extent_rows = extent_vertical;
    const int dim_cols = 2*extent_horizontal + 1;   // horizontal dimension
    // const int dim_rows = 2*extent_vertical + 1;     // vertical dimension
    // const int taps = dim_cols * dim_rows;           // number of elements

    // TODO: Precondition: Check if border is large enough for PSF overhang

    // TODO: Precondition: Check if dc gain is 1
    // sum of taps / num of taps

    // Initialize output image
    image_out->rows = height;
    image_out->cols = width;
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->stride = width;

    pixel_t *out_buf = malloc(height * width * planes * sizeof(pixel_t));
    image_out->handle = out_buf;
    image_out->buf = out_buf;

    // Find center of PSF
    pixel_t *h = psf_values + dim_cols*extent_vertical + extent_horizontal;

    // Perform convolution
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_in->buf + (row*stride + col)*planes;
            pixel_t *out_p = out_buf + (row*width + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            pixel_t sum[3] = {0,0,0};
            for (int y = -extent_rows; y < extent_rows+1; y++) {
                for (int x = -extent_cols; x < extent_cols+1; x++) {
                    int index_in = (y*stride + x) * planes;
                    int index_h = (-y*dim_cols - x); // mirrored
                    // printf("h[%d, %d] = %f\n", x, y, h[index_h]);
                    if (planes == 3) {
                        sum[0] += in_p[index_in] * h[index_h];
                        sum[1] += in_p[index_in+1] * h[index_h];
                        sum[2] += in_p[index_in+2] * h[index_h];
                    } else {
                        sum[0] += in_p[index_in] * h[index_h];
                    }
                }
            }

            // Store result of inner product
            if (planes == 3) {
                out_p[0] = sum[0];
                out_p[1] = sum[1];
                out_p[2] = sum[2];
            } else {
                *out_p = sum[0];
            }
        }
    }

    // TODO: What about "unused" portions of border?
    // e.g. a 1D horizontal PSF will leave the top and bottom borders untouched
    // Keeping these values in the output image could be useful for separable
    // filter implementation.
}
