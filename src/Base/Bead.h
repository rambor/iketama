//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_BEAD_H
#define IKETAMA_BEAD_H

#include "vector3.h"

class Bead {

    vector3 vec;
    float x;
    float y;
    float z;
    float contrast;

public:
    Bead(float x, float y, float z, float contrast);
    float const getContrast() const {return contrast;} // likely relative to protein ?
    void setContrast(float contrast);

    float const & getX() const {return vec.x;}
    float const & getY() const {return vec.y;}
    float const & getZ() const {return vec.z;}
    vector3 const & getVec() const {return vec;}
};

#endif //IKETAMA_BEAD_H
