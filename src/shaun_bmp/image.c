#include <stdlib.h> // malloc
#include <stdint.h>
#include <string.h> // memset
#include <assert.h> // hacky preconditions for the jank function at end of file

#include "image.h"
#include "bmp_io.h"

static uint8_t convert_from_pixel_to_byte(pixel_t pixel_value) {
    const int max_byte_value = 255;

    // Clip max/min pixel values to avoid numerical wrap-around.
    if (pixel_value > 1.0) {
        pixel_value = 1.0;
    }
    else if (pixel_value < 0.0) {
        pixel_value = 0.0;
    }

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

//! \brief Copy image data to a new image struct
//!
//! Useful when an image requires boundary extension for a filtering step
//!
//! @param image_info   Image to copy pixel data from
//! @param image_copy   Image to copy pixel data to
//! @param border       Thickness of border in pixels
void copy_image(image *image_info, image *image_copy, int border) {
    // Constants between the source and destination images
    const int width = image_info->cols;
    const int height = image_info->rows;
    const int planes = image_info->num_components;

    // Copied image info
    image_copy->cols = width;
    image_copy->rows = height;
    image_copy->num_components = planes;
    image_copy->border = border;
    image_copy->stride = width + 2*border;

    // Note: there are 2 different strides... don't mix them up.

    // Allocate memory to hold boundary extended image
    pixel_t *handle = malloc(
        image_copy->stride * (height + 2*border) * planes * sizeof(pixel_t));
    image_copy->handle = handle;
    image_copy->buf = 
        handle + (((border * image_copy->stride) + border) * planes);

    // Copy image data from image_in to image_copy, one row at a time.
    // Note: Can't copy more than a single row, because border gets in the way.
    for (int row = 0; row < height; row++) {
        pixel_t *source = image_info->buf + row * image_info->stride * planes;
        pixel_t *dest = image_copy->buf + row * image_copy->stride * planes;
        // memcpy(source, dest, width * planes * sizeof(pixel_t));
        for (int col = 0; col < width; col++) {
            if (planes == 3) {
                const int n = col*planes;
                dest[n] = source[n];
                dest[n+1] = source[n+1];
                dest[n+2] = source[n+2];
            } else {
                dest[col] = source[col];
            }
        }
    }

    // Remember to perform boundary extension on image_copy...
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
            "%dx%d image is too small for %d thick odd symmetric border.\n", 
            width, height, border);
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

//! \brief Shift intensity of all pixel values in image by a constant value.
//!
//! y[n] = x[n] + shift
void perform_level_shift(image *image_info, pixel_t shift) {
    const int height = image_info->rows;
    const int stride = image_info->stride;
    const int true_height = height + 2*image_info->border;
    const int planes = image_info->num_components;

    pixel_t *handle = image_info->handle;

    for (int row = 0; row < true_height; row++) {
        for (int col = 0; col < stride; col++) {
            const int n = (row * stride * planes) + (col * planes);
            for (int i = 0; i < planes; i++) {
                handle[n+i] = handle[n+i] + shift;
            }
        }
    }
}

//! \brief Scale intensity of all pixel values in image by a constant value.
//!
//! y[n] = scale * x[n]
//!
//! Note that this operation can also be achieved by scaling filter tap values.
void perform_scaling(image *image_info, float scale) {
    const int height = image_info->rows;
    const int stride = image_info->stride;
    const int true_height = height + 2*image_info->border;
    const int planes = image_info->num_components;

    pixel_t *handle = image_info->handle;

    for (int row = 0; row < true_height; row++) {
        for (int col = 0; col < stride; col++) {
            const int n = (row * stride * planes) + (col * planes);
            for (int i = 0; i < planes; i++) {
                handle[n+i] = scale * handle[n+i];
            }
        }
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

//! \brief Cascaded convolution of two 1D point spread functions (PSF).
//!
//! "Extent" refers to the number of elements extending from the center point 
//! in a direction. (e.g. a 3x3 PSF has an extent of 1 in both horizontal and 
//! vertical directions).
//!
//! @param image_in             Input image to convolve
//! @param image_out            Result of convolution
//! @param h1_values            Values of 1D PSF in horizontal direction
//! @param h2_values            Values of 1D PSF in vertical direction
//! @param extent1              Extent of h1  i.e. region of support [H1,1]
//! @param extent1              Extent of h2  i.e. region of support [1,H2]
void apply_separable_filters(image *image_in, image *image_out, 
pixel_t *h1_values, pixel_t *h2_values, int extent1, int extent2) {
    // Image constants
    const int height = image_in->rows;
    const int width = image_in->cols;
    const int planes = image_in->num_components;
    const int stride = image_in->stride;
    const int border = image_in->border;

    // Filter h1 constants
    const int EXTENT1 = extent1;
    // const int DIM1 = 2*EXTENT1 + 1;

    // Filter h2 constants
    const int EXTENT2 = extent2;
    // const int DIM2 = 2*EXTENT2 + 1;

    // Find center of PSF h1 and h2
    pixel_t *h1 = h1_values + EXTENT1;
    pixel_t *h2 = h2_values + EXTENT2;

    // TODO: Precondition: Check if border is large enough for PSF overhang

    // Initialize intermediate image  i.e. image_h1 = ( image_in * h1 )[n]
    image image_h1;
    image_h1.rows = height;
    image_h1.cols = width;
    image_h1.num_components = planes;
    image_h1.border = border;   // We keep the border in this image
    image_h1.stride = stride;   // 
    
    image_h1.handle = malloc(
        (height+2*border) * stride * planes * sizeof(pixel_t));
    image_h1.buf = image_h1.handle + (stride*border + border) * planes;

    // Initialize output image
    image_out->rows = height;
    image_out->cols = width;
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->stride = width;

    pixel_t *out_buf = malloc(height * width * planes * sizeof(pixel_t));
    image_out->handle = out_buf;
    image_out->buf = out_buf;

    // Perform convolution with h1
    // TODO: Keep border from image_in for next convolution?
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_in->buf + (row*stride + col)*planes;
            pixel_t *out_p = image_h1.buf + (row*stride + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h1 is horizontal vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT1; n < EXTENT1+1; n++) {
                int index_in = (n) * planes;
                int index_h1 = (-n); // mirrored
                if (planes == 3) {
                    sum[0] += in_p[index_in] * h1[index_h1];
                    sum[1] += in_p[index_in+1] * h1[index_h1];
                    sum[2] += in_p[index_in+2] * h1[index_h1];
                } else {
                    sum[0] += in_p[index_in] * h1[index_h1];
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

    // Hack: Just extend border again, rather than reusing border from image_in
    perform_boundary_extension(&image_h1);

    // Perform convolution with h2
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_h1.buf + (row*stride + col)*planes;
            pixel_t *out_p = image_out->buf + (row*width + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h2 is vertical vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT2; n < EXTENT2+1; n++) {
                int index_in = (n*stride) * planes;
                int index_h2 = (-n); // mirrored
                if (planes == 3) {
                    sum[0] += in_p[index_in] * h2[index_h2];
                    sum[1] += in_p[index_in+1] * h2[index_h2];
                    sum[2] += in_p[index_in+2] * h2[index_h2];
                } else {
                    sum[0] += in_p[index_in] * h2[index_h2];
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
}

// Same as apply_separable_filters, but only applies filter to every 2nd pixel.
// Kind of like decimation.
void apply_separable_filters_2n(image *image_in, image *image_out, 
pixel_t *h1_values, pixel_t *h2_values, int extent1, int extent2) {
    // Image constants
    const int height = image_in->rows;
    const int width = image_in->cols;
    const int planes = image_in->num_components;
    const int stride = image_in->stride;
    const int border = image_in->border;

    // Filter h1 constants
    const int EXTENT1 = extent1;
    // const int DIM1 = 2*EXTENT1 + 1;

    // Filter h2 constants
    const int EXTENT2 = extent2;
    // const int DIM2 = 2*EXTENT2 + 1;

    // Find center of PSF h1 and h2
    pixel_t *h1 = h1_values + EXTENT1;
    pixel_t *h2 = h2_values + EXTENT2;

    // TODO: Precondition: Check if border is large enough for PSF overhang

    // Initialize intermediate image  i.e. image_h1 = ( image_in * h1 )[n]
    image image_h1;
    image_h1.rows = height;
    image_h1.cols = width / 2;      // Decimate in horizontal direction
    image_h1.num_components = planes;
    image_h1.border = border;   // We keep the border in this image
    image_h1.stride = image_h1.cols + 2*border;
    
    image_h1.handle = malloc(
        (image_h1.rows+2*border) * image_h1.stride * planes * sizeof(pixel_t));
    image_h1.buf = image_h1.handle + (image_h1.stride*border + border) * planes;

    // Initialize output image
    image_out->rows = height / 2;   // Decimate in vertical direction
    image_out->cols = width / 2;
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->stride = image_out->cols;

    pixel_t *out_buf = malloc(height/2 * width/2 * planes * sizeof(pixel_t));
    image_out->handle = out_buf;
    image_out->buf = out_buf;

    // Perform convolution with h1
    // TODO: Keep border from image_in for next convolution?
    for (int row = 0; row < (height); row++) {
        for (int col = 0; col < (width/2); col++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_in->buf + (row*stride + col*2)*planes;
            pixel_t *out_p = image_h1.buf + (row*image_h1.stride + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h1 is horizontal vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT1; n < EXTENT1+1; n++) {
                int index_in = (n) * planes;
                int index_h1 = (-n); // mirrored
                if (planes == 3) {
                    sum[0] += in_p[index_in] * h1[index_h1];
                    sum[1] += in_p[index_in+1] * h1[index_h1];
                    sum[2] += in_p[index_in+2] * h1[index_h1];
                } else {
                    sum[0] += in_p[index_in] * h1[index_h1];
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

    // Hack: Just extend border again, rather than reusing border from image_in
    perform_boundary_extension(&image_h1);

    // Perform convolution with h2
    for (int row = 0; row < (height/2); row++) {
        for (int col = 0; col < (width/2); col++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_h1.buf + (row*2*image_h1.stride + col)*planes;
            pixel_t *out_p = image_out->buf + (row*image_out->stride + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h2 is vertical vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT2; n < EXTENT2+1; n++) {
                int index_in = (n*image_h1.stride) * planes;
                int index_h2 = (-n); // mirrored
                if (planes == 3) {
                    sum[0] += in_p[index_in] * h2[index_h2];
                    sum[1] += in_p[index_in+1] * h2[index_h2];
                    sum[2] += in_p[index_in+2] * h2[index_h2];
                } else {
                    sum[0] += in_p[index_in] * h2[index_h2];
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
}


// Same as apply_separable_filters(), but only applie filters to a single 
// component of the image
void apply_separable_filters_to_comp(image *image_in, image *image_out, 
pixel_t *h1_values, pixel_t *h2_values, int extent1, int extent2, 
int component) {
    // Copy paste from apply_separable_filters
    // Yeah, it's jank. But this is due in 2 hours.

    // Image constants
    const int height = image_in->rows;
    const int width = image_in->cols;
    const int planes = image_in->num_components;
    const int stride = image_in->stride;
    const int border = image_in->border;

    // Filter h1 constants
    const int EXTENT1 = extent1;
    // const int DIM1 = 2*EXTENT1 + 1;

    // Filter h2 constants
    const int EXTENT2 = extent2;
    // const int DIM2 = 2*EXTENT2 + 1;

    // Find center of PSF h1 and h2
    pixel_t *h1 = h1_values + EXTENT1;
    pixel_t *h2 = h2_values + EXTENT2;

    // TODO: Precondition: Check if border is large enough for PSF overhang

    // Initialize intermediate image  i.e. image_h1 = ( image_in * h1 )[n]
    image image_h1;
    image_h1.rows = height;
    image_h1.cols = width;
    image_h1.num_components = planes;
    image_h1.border = border;   // We keep the border in this image
    image_h1.stride = stride;   // 
    
    image_h1.handle = malloc(
        (height+2*border) * stride * planes * sizeof(pixel_t));
    image_h1.buf = image_h1.handle + (stride*border + border) * planes;

    // Initialize output image
    image_out->rows = height;
    image_out->cols = width;
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->stride = width;

    pixel_t *out_buf = malloc(height * width * planes * sizeof(pixel_t));
    image_out->handle = out_buf;
    image_out->buf = out_buf;

    // Perform convolution with h1
    // TODO: Keep border from image_in for next convolution?
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_in->buf + (row*stride + col)*planes;
            pixel_t *out_p = image_h1.buf + (row*stride + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h1 is horizontal vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT1; n < EXTENT1+1; n++) {
                int index_in = (n) * planes;
                int index_h1 = (-n); // mirrored
                sum[component] += in_p[index_in+component] * h1[index_h1];
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

    // Hack: Just extend border again, rather than reusing border from image_in
    perform_boundary_extension(&image_h1);

    // Perform convolution with h2
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_h1.buf + (row*stride + col)*planes;
            pixel_t *out_p = image_out->buf + (row*width + col)*planes;

            // Inner product of image with mirrored PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // Reminder that h2 is vertical vector
            pixel_t sum[3] = {0,0,0};
            for (int n = -EXTENT2; n < EXTENT2+1; n++) {
                int index_in = (n*stride) * planes;
                int index_h2 = (-n); // mirrored
                sum[component] += in_p[index_in+component] * h2[index_h2];
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
}


void apply_optimized_moving_average_filter(image *image_in, image *image_out, 
int extent1, int extent2) {
    // Image constants
    const int height = image_in->rows;
    const int width = image_in->cols;
    const int planes = image_in->num_components;
    const int stride = image_in->stride;
    const int border = image_in->border;

    // Filter h1 constants
    const int EXTENT1 = extent1;
    const int DIM1 = 2*EXTENT1 + 1;

    // Filter h2 constants
    const int EXTENT2 = extent2;
    const int DIM2 = 2*EXTENT2 + 1;

    // Find value for PSF taps
    pixel_t h = 1.0 / (DIM1 * DIM2);

    // TODO: Precondition: Check if border is large enough for PSF overhang

    // Initialize intermediate image  i.e. image_h1 = ( image_in * h1 )[n]
    image image_h1;
    image_h1.rows = height;
    image_h1.cols = width;
    image_h1.num_components = planes;
    image_h1.border = border;   // We keep the border in this image
    image_h1.stride = stride;   // 
    
    image_h1.handle = malloc(
        (height+2*border) * stride * planes * sizeof(pixel_t));
    image_h1.buf = image_h1.handle + (stride*border + border) * planes;

    // Initialize output image
    image_out->rows = height;
    image_out->cols = width;
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->stride = width;

    pixel_t *out_buf = malloc(height * width * planes * sizeof(pixel_t));
    image_out->handle = out_buf;
    image_out->buf = out_buf;

    // Perform convolution with h1
    // TODO: Keep border from image_in for next convolution?
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_in->buf + (row*stride + col)*planes;
            pixel_t *out_p = image_h1.buf + (row*stride + col)*planes;

            // Optimization: Moving average inner products
            //
            // Step 1: Take inner product for pixel c
            // +---------+
            // |a b c d e|f g
            // +---------+
            //
            // The inner product for step 1 requires summing the product of 
            // pixels a to e with a PSF value five times.
            //
            // sum1 = a*PSF[-2] + b*PSF[-1] + c*PSF[0] + d*PSF[1] + e*PSF[2]
            //
            // Step 2: Take inner product for pixel d
            //   +---------+
            //  a|b c d e f|g
            //   +---------+
            //
            // For a moving average filter, the PSF values are all the same. 
            // Therefore, we can reuse part of the calculation from step 1.
            //
            // sum2 = sum1 - a*PSF + f*PSF
            //

            // Inner product of image with mirrored moving average PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // horizontal
            pixel_t sum[3] = {0,0,0};
            pixel_t prev_sum[3] = {0,0,0};

            if (col == 0) {
                // For first pixel, sum all pixel values
                for (int n = -EXTENT1; n < EXTENT1+1; n++) {
                    int index_in = (n) * planes;
                    if (planes == 3) {
                        sum[0] += in_p[index_in] * h;
                        sum[1] += in_p[index_in+1] * h;
                        sum[2] += in_p[index_in+2] * h;
                    } else {
                        sum[0] += in_p[index_in] * h;
                    }
                }
            } else {
                // For rest of row, use optimization
                int n_sub = (-EXTENT1 - 1) * planes;
                int n_add = EXTENT1 * planes;
                if (planes == 3) {
                    sum[0] = prev_sum[0] + h * (in_p[n_add] - in_p[n_sub]);
                    sum[1] = prev_sum[1] + h * (in_p[n_add+1] - in_p[n_sub+1]);
                    sum[2] = prev_sum[2] + h * (in_p[n_add+2] - in_p[n_sub+2]);
                } else {
                    sum[0] = prev_sum[0] + h * (in_p[n_add] - in_p[n_sub]);
                }
            }

            // Store result of inner product
            if (planes == 3) {
                out_p[0] = prev_sum[0] = sum[0];
                out_p[1] = prev_sum[1] = sum[1];
                out_p[2] = prev_sum[2] = sum[2];
            } else {
                *out_p = prev_sum[0] = sum[0];
            }

        }
    }

    // Hack: Just extend border again, rather than reusing border from image_in
    perform_boundary_extension(&image_h1);

    // Perform convolution with h2
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col ++) {
            // Pixel coordinate on input and output images
            pixel_t *in_p = image_h1.buf + (row*stride + col)*planes;
            pixel_t *out_p = image_out->buf + (row*width + col)*planes;

            // Inner product of image with mirrored moving average PSF
            // `h_mirror[n1, n2] = h[-n1, -n2]`
            // vertical
            pixel_t sum[3] = {0,0,0};
            pixel_t prev_sum[3] = {0,0,0};

            if (col == 0) {
                // For first pixel, sum all pixel values
                for (int n = -EXTENT2; n < EXTENT2+1; n++) {
                    int index_in = (n*stride) * planes;
                    if (planes == 3) {
                        sum[0] += in_p[index_in] * h;
                        sum[1] += in_p[index_in+1] * h;
                        sum[2] += in_p[index_in+2] * h;
                    } else {
                        sum[0] += in_p[index_in] * h;
                    }
                }
            } else {
                // For rest of row, use optimization
                int n_sub = (-EXTENT2 - 1) * planes;
                int n_add = EXTENT2 * planes;
                if (planes == 3) {
                    sum[0] = prev_sum[0] + h * (in_p[n_add] - in_p[n_sub]);
                    sum[1] = prev_sum[1] + h * (in_p[n_add+1] - in_p[n_sub+1]);
                    sum[2] = prev_sum[2] + h * (in_p[n_add+2] - in_p[n_sub+2]);
                } else {
                    sum[0] = prev_sum[0] + h * (in_p[n_add] - in_p[n_sub]);
                }
            }

            // Store result of inner product
            if (planes == 3) {
                out_p[0] = prev_sum[0] = sum[0];
                out_p[1] = prev_sum[1] = sum[1];
                out_p[2] = prev_sum[2] = sum[2];
            } else {
                *out_p = prev_sum[0] = sum[0];
            }
        }
    }
}





// // For project 1 task3b
// // Takes a color plane from each of the input images, and splices them all 
// // together in image_out.
// void hacky_RGB_image_splice(image *image_out, 
// image *image_R, image *image_G, image *image_B) {
//     const int width = image_R->cols;
//     const int height = image_R->rows;
//     const int planes = 3;

//     // Precondition: Check if images are the same size
//     assert(image_B->rows == height);
//     assert(image_G->rows == height);
//     assert(image_B->cols == width);
//     assert(image_G->cols == width);

//     // Precondtion: Check if images have BGR planes
//     assert(image_R->num_components == planes);
//     assert(image_G->num_components == planes);
//     assert(image_B->num_components == planes);

//     // Initialize image_out
//     image_out->num_components = planes;
//     image_out->border = 0;
//     image_out->cols = width;
//     image_out->rows = height;
//     image_out->handle = malloc(width * height * planes);
//     image_out->buf = image_out->handle;

//     // Extract color planes
//     for (int row = 0; row < height; row++) {
//         for (int col = 0; col < width; col++) {
//             const int nB = (row*image_B->stride + col) * planes;
//             const int nG = (row*image_G->stride + col) * planes + 1;
//             const int nR = (row*image_R->stride + col) * planes + 2;

//             // TODO: image_out index is wrong. stride = width.
            
//             image_out->buf[nB] = image_B->buf[nB];
//             image_out->buf[nG] = image_G->buf[nG];
//             image_out->buf[nR] = image_R->buf[nR];
//         }
//     }
// }

void extract_component(image *image_out, image *image_BGR, int component) {
    const int width = image_BGR->cols;
    const int height = image_BGR->rows;
    const int planes = image_BGR->num_components;
    const int stride = image_BGR->stride;

    // Precondition: RGB image
    assert(image_BGR->num_components == planes);

    // Initialize image_out
    image_out->num_components = 1;
    image_out->border = 0;
    image_out->cols = width;
    image_out->rows = height;
    image_out->handle = malloc(width * height * 1);
    image_out->buf = image_out->handle;

    // Extract color planes
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            pixel_t *in_p = image_BGR->buf + (row*stride + col)*planes + component;
            pixel_t *out_p = image_out->buf + (row*width + col)*1;
            *out_p = *in_p;

            // const int nBGR = (row*image_BGR->stride+col) * planes + component;
            // const int nOut = (row*width + col) * 1;
            // image_out->buf[nOut] = image_BGR->buf[nBGR];
            
            
        }
    }
}

void hacky_combine_planes_into_RGB(image *image_out, 
image *image_R, image *image_G, image *image_B) {
    const int width = image_B->cols;
    const int height = image_B->rows;
    const int planes = 3;

    assert(image_R->num_components == 1);
    assert(image_G->num_components == 1);
    assert(image_B->num_components == 1);

    // Initialize image_out
    image_out->num_components = planes;
    image_out->border = 0;
    image_out->cols = width;
    image_out->rows = height;
    image_out->handle = malloc(width * height * planes);
    image_out->buf = image_out->handle;

    // Combine color planes
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            pixel_t *inR_p = image_R->buf + (row*image_R->stride + col)*1;
            pixel_t *inG_p = image_G->buf + (row*image_G->stride + col)*1;
            pixel_t *inB_p = image_B->buf + (row*image_B->stride + col)*1;
            pixel_t *out_p = image_out->buf + (row*width + col)*planes;
            
            out_p[0] = *inB_p;
            out_p[1] = *inG_p;
            out_p[2] = *inR_p;
        }
    }
}
