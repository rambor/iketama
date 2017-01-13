//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_COORDS_H
#define IKETAMA_COORDS_H

#include <iostream>
#include <cstdio>
#include <string>


class Coords {

public:
    Coords();
    Coords(float x, float y, float z, std::string atomType, float occ);

    float x;
    float y;
    float z;
    std::string type;
    float occ;

};


#endif //IKETAMA_COORDS_H
