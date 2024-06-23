#ifndef SHAUN_BMP_IMAGE_H
#define SHAUN_BMP_IMAGE_H

#include "bmp_io.h"

////////////////////////////
// Data type declarations //
////////////////////////////

typedef struct image image;

////////////////
// Structures //
////////////////

struct image {
    // Image info
    int num_components, rows, cols;

    // Boundary extension border info
    int border, stride;

    // Memory
    uint8_t *handle;    // Points to start of pixel data (including border)
    uint8_t *buf;       // Points to start of actual image data
    // TODO: Change data type. When processing image, values may go outside 0-255 range. float is easy.
};

////////////////
// Public API //
////////////////

// Image data I/O
int read_image_from_bmp(image *image_info, bmp *bmp_info, int border);
int export_image_as_bmp(image *image_info, const char *fname);

// Image processing
// void apply_boundary_extension(image *image_info, int border_size);
// void apply_filter(image *image_in, image *image_out);

#endif // SHAUN_BMP_IMAGE_H
