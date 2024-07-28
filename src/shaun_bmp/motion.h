#ifndef SHAUN_BMP_MOTION_H
#define SHAUN_BMP_MOTION_H

#include "image.h"

typedef struct motion_vector mvector;

struct motion_vector {
    int x;
    int y;
};

mvector estimate_motion(image *source, image *target, int search_bounds);
void compensate_for_motion(image *source, image *target, mvector vec);

#endif // SHAUN_BMP_MOTION_H
