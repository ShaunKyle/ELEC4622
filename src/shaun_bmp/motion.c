#include "motion.h"

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
    pixel_t norm, best_norm;

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
                }
            }

            // Check if norm has improved for chosen motion vector
            if (norm < best_norm) {
                best_norm = norm;
                best_vec = vec;
                // printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
            }
        }
    }

    printf("DEBUG: (%d, %d) norm is %f\n", best_vec.x, best_vec.y, best_norm);
    return best_vec;
}
