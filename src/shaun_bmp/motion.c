#include "motion.h"
#include <math.h>
#include <stdlib.h> // abs


// For windowed sinc design...
static float hanning_window(int n, int extent) {
    // If extent is 0, we don't need a window.
    if (extent == 0) {
        return 1.0;
    }

    // Hanning window
    if (abs(n) < extent) {
        // printf("%f\n", cos(M_PI * n / extent));
        return (1 + cos(M_PI * n / extent)) / 2.0;
    }
    else {
        // puts("wall");
        return 0.0;
    }
}



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

    best_norm = 256*B*B*2*1000;  //TODO: Why 256?

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


//! \brief Same as estimate_motion_block, but uses windowed sinc interpolation 
// to increase accuracy of motion vector search.
//
// `g` is grid spacing in pels (e.g. 1/2 pel grid). 
// `g_inv` is 1/g (e.g. 1/2 pel grid corresponds to g_inv=2).
mvector2_t estimate_motion_block_sinc(image *source, image *target, int start_row, 
int start_col, int block_size, int search_bound, norm_e norm_choice) {
    const int g_inv = 4;    // Only for 1/4 pel grids.

    const int B = block_size;
    const int S = search_bound * g_inv;
    const int planes = source->num_components;

    mvector2_t vec, best_vec;
    vec.precision = g_inv;  // fractional pel precision

    pixel_t norm, best_norm, mse, best_mse;

    best_norm = 256*B*B*2*1000;  //TODO: Why 256?


    // Design shifted windowed sincs for interpolation
    // 1/4, 2/4 and 3/4 shifts
    const int EXTENT = 3;
    const int TAPS = 2*EXTENT + 1;
    pixel_t h1[EXTENT];
    pixel_t h2[EXTENT];
    pixel_t h3[EXTENT];
    for (int t = 0; t < TAPS; t++) {
        const int n = t - EXTENT;
        const double nshift1 = n - 0.25;
        const double nshift2 = n - 0.5;
        const double nshift3 = n - 0.75;
        h1[TAPS-t] = sin(M_PI * nshift1) / (M_PI * nshift1);
        h1[TAPS-t] *= hanning_window(n, EXTENT);
        h2[TAPS-t] = sin(M_PI * nshift2) / (M_PI * nshift2);
        h2[TAPS-t] *= hanning_window(n, EXTENT);
        h3[TAPS-t] = sin(M_PI * nshift3) / (M_PI * nshift3);
        h3[TAPS-t] *= hanning_window(n, EXTENT);
        // Using TAPS-t instead of t for indexing
        // To implement mirror PSF.
    }

    // Check DC gain of windowed sinc
    double dc_gain1 = 0;
    double dc_gain2 = 0;
    double dc_gain3 = 0;
    for (int row = 0; row < TAPS; row++) {
        for (int col = 0; col < TAPS; col++) {
            dc_gain1 += (h1[row] * h1[col]);
            dc_gain2 += (h2[row] * h2[col]);
            dc_gain3 += (h3[row] * h3[col]);
        }
    }
    
    // Normalize dc_gain to 1
    // Note that this design is for separable filters. Sqrt.
    for (int tap = 0; tap < TAPS; tap++) {
        h1[tap] /= sqrt(dc_gain1);
        h2[tap] /= sqrt(dc_gain2);
        h3[tap] /= sqrt(dc_gain3);

        // printf("h1[%d] = %f\n", tap, h1[tap]);
        // printf("h2[%d] = %f\n", tap, h1[tap]);
        // printf("h3[%d] = %f\n", tap, h1[tap]);
    }


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
                    // pixel_t interp = 
                    //     source_p[index] * (1.0-sigma2)*(1.0-sigma1) +
                    //     source_p[index+planes] * (1-sigma2)*sigma1 +
                    //     source_p[index+stride] * (1-sigma1)*sigma2 +
                    //     source_p[index+stride+planes] * sigma2*sigma1;
                    // TODO: bilinear interp stride needs to be multiplied by planes...
                    // Hasn't been an issue, because planes = 1.

                    pixel_t interp_col[TAPS];
                    for (int row = 0; row < TAPS; row++) {
                        interp_col[row] = 0;
                        for (int t = 0; t < TAPS; t++) {
                            const int n1 = t - EXTENT;
                            
                            // Decide which shifted sinc to use
                            int tap1;
                            if (sigma1 < 0.5) {
                                tap1 = h1[t];
                            } else if (sigma1 == 0.5) {
                                tap1 = h2[t];
                            } else {
                                tap1 = h3[t];
                            }

                            interp_col[row] += source_p[index + (row*stride + n1) * planes] * tap1;
                        }
                    }
                    

                    pixel_t interp = 0;
                    for (int t = 0; t < TAPS; t++) {
                        // const int n2 = t - EXTENT;

                        // Decide which shifted sinc to use
                        int tap2;
                        if (sigma2 < 0.5) {
                            tap2 = h1[t];
                        } else if (sigma2 == 0.5) {
                            tap2 = h2[t];
                        } else {
                            tap2 = h3[t];
                        }

                        interp += interp_col[t] * tap2;
                    }
                    // Further optimisation would be to save intermediate sums.

                    // if (interp > 0.0) {
                    //     printf("interp %f\n", interp);
                    // }
                    


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
    printf("Sinc%d MSE (tgt to displaced ref): %f\n", g_inv, best_mse);
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
