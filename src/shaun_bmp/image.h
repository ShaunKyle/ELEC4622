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

    // TODO: Boundary extension info?

    // Memory
    uint8_t *data;  // Points to start of pixel data
    // TODO: Change data type. When processing image, values may go outside 0-255 range.
};

////////////////
// Public API //
////////////////

// Image data I/O
int read_image_from_bmp(image *image_info, bmp *bmp_info);
int export_image_as_bmp(image *image_info, const char *fname);

// Image processing
void apply_boundary_extension(image *image_info, int border_size);

#endif // SHAUN_BMP_IMAGE_H
