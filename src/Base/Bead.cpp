//
// Created by Robert Rambo on 12/01/2017.
//

#include "Bead.h"

Bead::Bead(){};
Bead::Bead(float x, float y, float z, float contrast) : x(x), y(y), z(z), contrast(contrast){
    this->vec.x = x;
    this->vec.y = y;
    this->vec.z = z;
//    this->x = x;
//    this->y = y;
//    this->z = z;
//    this->contrast = contrast;

    /**
     * (x+y+z)
     * (x-y-z)
     * (x+y-z)
     * (x-y+z)
     * octant 1 (q1) (x+y+z)
     * octant 2 (q2) (x-y-z)
     * octant 3 (q3) (x+y-z)
     * octant 4 (q4) (x-y+z)
     */

    q1 = x+y+z;
    q2 = x-y-z;
    q3 = x+y-z;
    q4 = x-y+z;
}

void Bead::setContrast(float value) { this->contrast = value; }
