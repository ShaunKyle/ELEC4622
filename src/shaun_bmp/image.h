#ifndef SHAUN_BMP_IMAGE_H
#define SHAUN_BMP_IMAGE_H

#include "bmp_io.h"

////////////////////////////
// Data type declarations //
////////////////////////////

typedef struct image image;
typedef float pixel_t;      // Normalized value (0 to 255) -> (0.0 to 1.0)
// TODO: try fixed point.
// see convert functions in image.c
// When processing image, values may go outside 0-255 range. float is easy.
// Penalties associated with floating-point are covered in Ch2, pg 15.

////////////////
// Structures //
////////////////

struct image {
    // Image info
    int num_components, rows, cols;

    // Boundary extension border info
    int border, stride;

    // Memory
    pixel_t *handle;    // Points to start of pixel data (including border)
    pixel_t *buf;       // Points to start of actual image data
};

////////////////
// Public API //
////////////////

// Image data I/O
int read_image_from_bmp(image *image_info, bmp *bmp_info, int border);
int export_image_as_bmp(image *image_info, const char *fname);
int export_image_and_border_as_bmp(image *image_info, const char *fname);
void copy_image(image *image_info, image *image_copy, int border);
// TODO: delete_image?

// Image processing
// TODO: Is this pre/post processing?
void perform_boundary_extension(image *image_info);
void perform_level_shift(image *image_info, pixel_t shift);
void perform_scaling(image *image_info, float scale);

// TODO: Move filtering operations to a separate file? filter.c
void apply_filter(image *image_in, image *image_out, pixel_t *psf_values, 
int extent_horizontal, int extent_vertical);
void apply_separable_filters(image *image_in, image *image_out, 
pixel_t *h1_values, pixel_t *h2_values, int extent1, int extent2);

// TODO: rethink filter API. For now, here's separate function for applying a 
// filter to a single component of an image
void apply_separable_filters_to_comp(image *image_in, image *image_out, 
pixel_t *h1_values, pixel_t *h2_values, int extent1, int extent2, 
int component);
void apply_optimized_moving_average_filter(image *image_in, image *image_out, 
int extent1, int extent2);
// void hacky_RGB_image_splice(image *image_out, 
// image *image_R, image *image_G, image *image_B);

// I misread the task. more hacky. Also very broken.
void hacky_combine_planes_into_RGB(image *image_out, 
image *image_R, image *image_G, image *image_B);
void extract_component(image *image_out, image *image_BGR, int component);

#endif // SHAUN_BMP_IMAGE_H
