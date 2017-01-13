//
// Created by Robert Rambo on 12/01/2017.
//

#include "Coords.h"

using namespace std;

Coords::Coords(){
    this->x = 0;
    this->y = 0;
    this->z = 0;
    this->type = "";
    this->occ = 0;
}

Coords::Coords(float x, float y, float z, std::string atomType, float occ) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->type = atomType;
    this->occ = occ;
}