//
// Created by Robert Rambo on 12/01/2017.
//

#include "Bead.h"

Bead::Bead(float x, float y, float z, float contrast) {
    this->vec.x = x;
    this->vec.y = y;
    this->vec.z = z;
    this->x = x;
    this->y = y;
    this->z = z;
    this->contrast = contrast;
}

void Bead::setContrast(float value) { this->contrast = value; }
