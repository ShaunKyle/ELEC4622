#include "motion.h"

//! \brief Finds the motion vector which best describes the motion between the 
//! reference and target frames.
mvector estimate_motion(image *source, image *target, int search_bounds) {
    // Assume that input images have already been sliced into equal blocks
    const int block_width = source->cols;
    const int block_height = source->rows;

    mvector vec, best_vec;
    int sad, best_sad;

    best_sad = 999999999;

    // Search for best vector (x,y) using some optimization criterion
    // The search range is [-S,+S] pels in each direction
    for (vec.y = -search_bounds; vec.y <= search_bounds; vec.y++) {
        for (vec.x = -search_bounds; vec.x <= search_bounds; vec.x++) {

            // TODO: modification to ptrs for given vec?

            // Calculate sum-of-abs-differences (SAD) for given vec
            pixel_t *source_p = source->buf;
            pixel_t *target_p = target->buf;
            sad = 0;
            for (int r = 0; r < block_height; r++) {
                for (int c = 0; c < block_width; c++) {
                    const int index = 
                        (c + (r * source->stride)) * source->num_components;
                    // Side note: Don't bother using more than one color plane
                    // We care about the shapes, not the color.
                    int diff = target_p[index] - source_p[index];
                    sad += (diff < 0)?(-diff):diff;
                }
            }

            // Check if SAD improved
            if (sad < best_sad) {
                best_sad = sad;
                best_vec = vec;
            }
        }
    }

    return best_vec;
}