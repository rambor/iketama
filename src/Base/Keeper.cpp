//
// Created by xos81802 on 01/12/2017.
//

#include "Keeper.h"

Keeper::Keeper(int upperLimit) : upperLimit(upperLimit) {

    indices.resize(upperLimit);
    for(int i=0; i<upperLimit; i++){
        indices[i]=0;
    }
}

/**
 * sorted vector
 * @param upperLimit
 * @param vertexIndices
 */
void Keeper::updateIndices(int upperLimit, std::vector<int> &vertexIndices, double d_kl, double contactsE, double totalE) {
    this->upperLimit = upperLimit;
    int * ptr = (this->upperLimit != 0) ? &indices.front() : NULL;
    for(int i = 0; i < this->upperLimit; i++) {
        ptr[i] = vertexIndices[i];
    }

    this->d_kl = d_kl;
    this->contactsPotential = contactsE;
    this->totalEnergy = totalE;
}

void Keeper::recalculateEnergy(double constant){
    this->totalEnergy = d_kl + constant*contactsPotential;
}