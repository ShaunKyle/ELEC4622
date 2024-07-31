#ifndef SHAUN_BMP_MOTION_H
#define SHAUN_BMP_MOTION_H

#include "image.h"

struct motion_vector {
    int x;
    int y;
};

struct motion_vector_fixed_precision {
    int x, y;
    int precision;
};

enum norm_type {
    MAD,    // Mean absolute difference
    MSE,    // Mean squared error
};

typedef struct motion_vector mvector_t;
typedef struct motion_vector_fixed_precision mvector2_t;
typedef enum norm_type norm_e;

mvector_t estimate_motion_block(
    image *source, image *target, int start_row, int start_col, 
    int block_size, int search_bounds, norm_e norm,
    double *mse_out // hack
);
mvector2_t estimate_motion_block_bilinear(
    image *source, image *target, int start_row, int start_col, 
    int block_size, int search_bounds, norm_e norm, int g_inv,
    double *mse_out // hack
);
mvector2_t estimate_motion_block_sinc(
    image *source, image *target, int start_row, int start_col, 
    int block_size, int search_bounds, norm_e norm
);
mvector2_t estimate_motion_block_telescopic(
    image *source, image *target, int start_row, int start_col, 
    int block_size, int search_bounds, norm_e norm,
    double *mse_out // hack
);
void compensate_motion_block(
    image *source, image *target, mvector_t vec, 
    int start_row, int start_col, int block_size
);
void compensate_motion_block_fractional(
    image *source, image *target, mvector2_t vec, 
    int start_row, int start_col, int block_size
);

#endif // SHAUN_BMP_MOTION_H
