//
// Created by Robert Rambo on 12/01/2017.
//


#include "Phase.h"


using namespace std;

Phase::Phase(int volume, float sigma, float contrast){
    this->volume = volume;
    this->contrast = contrast;
    this->sigma = sigma;
}


Phase::Phase(int volume, float sigma, float contrast, string id) : Phase(volume, sigma, contrast) {
    this->id = id;
}

void Phase::setNumberOfBeads(float bead_radius) {
    float value = 4.0/3.0*M_PI*bead_radius*bead_radius*bead_radius;

    this->bead_radius = bead_radius;
    this->number_of_beads = (int)(volume/value);
    this->sigma_number_of_beads = sigma*this->number_of_beads;
    //cout << "Number of Beads " << this->number_of_beads << " => " << volume << endl;
    this->electrons_per_bead = (value * (this->contrast));
}

/**
 * startIndex is inclusive
 * workingLimit is exclusive
 */
void  Phase::setLimits(int startIndex, int workingLimit) {
    this->startIndex = startIndex;
    this->workingLimit = workingLimit;
}