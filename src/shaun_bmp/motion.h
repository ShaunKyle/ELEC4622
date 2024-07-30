#ifndef SHAUN_BMP_MOTION_H
#define SHAUN_BMP_MOTION_H

#include "image.h"

struct motion_vector {
    int x;
    int y;
};

enum norm_type {
    MAD,    // Mean absolute difference
    MSE,    // Mean squared error
};

typedef struct motion_vector mvector_t;
typedef enum norm_type norm_e;

mvector_t estimate_motion_block(
    image *source, image *target, int start_row, int start_col, 
    int block_size, int search_bounds, norm_e norm);
void compensate_for_motion(image *source, image *target, mvector_t vec);

#endif // SHAUN_BMP_MOTION_H
