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

    /**
      * octant 1 (q1) (x+y+z)
      * octant 2 (q2) (x-y-z)
      * octant 3 (q3) (x+y-z)
      * octant 4 (q4) (x-y+z)
      *
      */

    float q1, q2, q3, q4;

public:
    Bead();
    Bead(float x, float y, float z, float contrast);
    float const getContrast() const {return contrast;} // likely relative to protein ?
    void setContrast(float contrast);

    float const & getX() const {return vec.x;}
    float const & getY() const {return vec.y;}
    float const & getZ() const {return vec.z;}

    /**
     * x+y+z
     */
    float const & getQ1() const {return q1;}

    /**
     * x-y-z
     */
    float const & getQ2() const {return q2;}

    /**
    * x+y-z
    */
    float const & getQ3() const {return q3;}

    /**
     * x-y+z
     */
    float const & getQ4() const {return q4;}

    vector3 const & getVec() const {return vec;}
};

#endif //IKETAMA_BEAD_H
