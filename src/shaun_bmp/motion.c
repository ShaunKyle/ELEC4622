#include "motion.h"
#include <math.h>

//! \brief Finds the motion vector which best describes the movement of a block 
//! in the target frame.
//!
//! @param source       Reference image (frame 1)
//! @param target       Target image (frame 2)
//! @param start_row    Origin of block in target frame
//! @param start_col    Origin of block in target frame
//! @param block_size   Size of motion block (pels)
//! @param search_bound Search range (pels) for motion vector in each direction
//! @param norm_choice  Norm used for minimization criteria
mvector_t estimate_motion_block(image *source, image *target, int start_row, 
int start_col, int block_size, int search_bound, norm_e norm_choice) {
    const int B = block_size;
    const int S = search_bound;
    const int planes = source->num_components;

    mvector_t vec, best_vec;
    pixel_t norm, best_norm, mse, best_mse;

    best_norm = 256*B*B;  //TODO: Why 256?

    // Search for best vector (x,y) using some optimization criterion
    // The search range is [-S,+S] pels in each direction
    for (vec.y = -S; vec.y <= S; vec.y++) {
        for (vec.x = -S; vec.x <= S; vec.x++) {
            // Check if block translated by motion vector is contained within 
            // the source frame. Skip this motion vector if outside.
            // TODO: Leave a "border" so that all motion vectors are possible?
            int src_row = start_row - vec.y;
            int src_col = start_col - vec.x;
            if ((src_row < 0) || (src_col < 0) || 
                ((src_row+B) > source->rows) || 
                ((src_col+B) > source->cols)) {
                // printf("(%d, %d) translated out of source frame\n", vec.x, vec.y);
                continue;
            }
            // printf("(%d, %d) within source frame\n", vec.x, vec.y);

            // Calculate the norm of the displaced frame difference (DFD) 
            // between the target block and displaced source block.
            pixel_t *source_p = source->buf + 
                (src_row * source->stride + src_col) * planes;
            pixel_t *target_p = target->buf + 
                (start_row * target->stride + start_col) * planes;
            norm = 0;
            mse = 0;
            for (int r = 0; r < B; r++) {
                for (int c = 0; c < B; c++) {
                    const int stride = source->stride;  //assume target is same
                    const int index = (c + (r * stride)) * planes;

                    // Calculate DFD
                    pixel_t diff = target_p[index] - source_p[index];
                    // Side note: Don't bother using more than one color plane
                    // We care about the shapes, not the color.

                    // Calculate norm
                    // Note: In both cases, use sum rather than mean.
                    if (norm_choice == MSE) {
                        // puts("Chosen norm: Mean squared error");
                        norm += (diff*diff);
                    }
                    else {
                        // puts("Chosen norm: Mean absolute difference");
                        norm += (diff < 0) ? (-diff) : diff;
                    }
                    mse += (diff*diff);
                }
            }

            // Check if norm has improved for chosen motion vector
            if (norm < best_norm) {
                best_norm = norm;
                best_vec = vec;
                // printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
                best_mse = mse;
            }
        }
    }

    // printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
    printf("MSE (tgt to displaced ref): %f\n", best_mse);
    return best_vec;
}


//! \brief Same as estimate_motion_block, but uses bilinear interpolation to 
// increase accuracy of motion vector search.
//
// `g` is grid spacing in pels (e.g. 1/2 pel grid). 
// `g_inv` is 1/g (e.g. 1/2 pel grid corresponds to g_inv=2).
mvector2_t estimate_motion_block_bilinear(image *source, image *target, int start_row, 
int start_col, int block_size, int search_bound, norm_e norm_choice, int g_inv) {
    // const int g_inv = 2;    // g = pixel spacing (half-pel grid)

    const int B = block_size;
    const int S = search_bound * g_inv;
    const int planes = source->num_components;

    mvector2_t vec, best_vec;
    vec.precision = g_inv;  // Half-pel precision

    pixel_t norm, best_norm, mse, best_mse;

    best_norm = 256*B*B*2;  //TODO: Why 256?

    // Search for best vector (x,y) using some optimization criterion
    // The search range is [-S,+S] pels in each direction
    for (vec.y = -S; vec.y <= S; vec.y++) {
        for (vec.x = -S; vec.x <= S; vec.x++) {
            // Check if block translated by motion vector is contained within 
            // the source frame. Skip this motion vector if outside.
            // TODO: Leave a "border" so that all motion vectors are possible?
            double src_row = start_row - vec.y/(double)g_inv;
            double src_col = start_col - vec.x/(double)g_inv;
            if ((src_row < 0) || (src_col < 0) || 
                ((src_row+B) > source->rows) || 
                ((src_col+B) > source->cols)) {
                // printf("(%f, %f) translated out of source frame\n", vec.x/(double)g_inv, vec.y/(double)g_inv);
                continue;
            }
            // printf("(%f, %f) within source frame\n", vec.x/(double)g_inv, vec.y/(double)g_inv);

            // Calculate the norm of the displaced frame difference (DFD) 
            // between the target block and displaced source block.
            // pixel_t *source_p;
            pixel_t *target_p = target->buf + 
                (start_row * target->stride + start_col) * planes;
            double sigma1 = fmod(src_col, 1.0);
            double sigma2 = fmod(src_row, 1.0);
            pixel_t *source_p = source->buf + 
                    ((int)src_row * source->stride + (int)src_col) * planes;
            
            norm = 0;
            mse = 0;
            for (int r = 0; r < B; r++) {
                for (int c = 0; c < B; c++) {
                    const int stride = source->stride;  //assume target is same
                    const int index = (c + (r * stride)) * planes;

                    // Interpolate pixel on source frame
                    pixel_t interp = 
                        source_p[index] * (1.0-sigma2)*(1.0-sigma1) +
                        source_p[index+planes] * (1-sigma2)*sigma1 +
                        source_p[index+stride] * (1-sigma1)*sigma2 +
                        source_p[index+stride+planes] * sigma2*sigma1;


                    // Calculate DFD
                    pixel_t diff = target_p[index] - interp;
                    // Side note: Don't bother using more than one color plane
                    // We care about the shapes, not the color.

                    // Calculate norm
                    // Note: In both cases, use sum rather than mean.
                    if (norm_choice == MSE) {
                        // puts("Chosen norm: Mean squared error");
                        norm += (diff*diff);
                    }
                    else {
                        // puts("Chosen norm: Mean absolute difference");
                        norm += (diff < 0) ? (-diff) : diff;
                    }
                    mse += (diff*diff);
                }
            }

            // Check if norm has improved for chosen motion vector
            if (norm < best_norm) {
                best_norm = norm;
                best_vec = vec;
                // printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
                best_mse = mse;
            }
        }
    }

    // printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
    printf("Bilinear%d MSE (tgt to displaced ref): %f\n", g_inv, best_mse);
    return best_vec;
}


//! \brief Motion compensation for target frame by copying displaced blocks 
// from reference frame.
void compensate_motion_block(image *source, image *target, mvector_t vec, 
int start_row, int start_col, int block_size) {
    // Coordinate for start of displaced block in reference frame
    const int src_row = start_row - vec.y;
    const int src_col = start_col - vec.x;

    // Pointer to start of target block, and start of displaced reference block
    const int planes = source->num_components;
    pixel_t *source_p = 
        source->buf + (src_row * source->stride + src_col) * planes;
    pixel_t *target_p = 
        target->buf + (start_row * target->stride + start_col) * planes;
    
    // Overwrite target block with displaced reference block
    for (int r = 0; r < block_size; r++) {
        for (int c = 0; c < block_size; c++) {
            const int index = (r*source->stride + c) * planes;
            for (int p = 0; p < planes; p++) {
                target_p[index+p] = source_p[index+p];
            }
        }
    }
}


//! \brief Same as compensate_motion_block, but for fractional motion vectors.
void compensate_motion_block_fractional(image *source, image *target, 
mvector2_t vec, int start_row, int start_col, int block_size) {
    const int g_inv = vec.precision;

    // Coordinate for start of displaced block in reference frame
    double src_row = start_row - vec.y/(double)g_inv;
    double src_col = start_col - vec.x/(double)g_inv;
    double sigma1 = fmod(src_col, 1.0);
    double sigma2 = fmod(src_row, 1.0);

    // Pointer to start of target block, and start of displaced reference block
    const int planes = source->num_components;
    pixel_t *source_p = source->buf + 
        ((int)src_row * source->stride + (int)src_col) * planes;
    pixel_t *target_p = 
        target->buf + (start_row * target->stride + start_col) * planes;
    
    // Overwrite target block with displaced reference block
    for (int r = 0; r < block_size; r++) {
        for (int c = 0; c < block_size; c++) {
            const int index = (r*source->stride + c) * planes;
            for (int p = 0; p < planes; p++) {
                // Interpolate pixel on source frame
                const int stride = source->stride;
                pixel_t interp = 
                    source_p[index+p] * (1.0-sigma2)*(1.0-sigma1) +
                    source_p[index+planes+p] * (1-sigma2)*sigma1 +
                    source_p[index+stride+p] * (1-sigma1)*sigma2 +
                    source_p[index+stride+planes+p] * sigma2*sigma1;
                
                // Write
                target_p[index+p] = interp;
            }
        }
    }
}
