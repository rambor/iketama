//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_PHASE_H
#define IKETAMA_PHASE_H


#include <string>
#include <vector>
#include "math.h"
#include <iostream>
#include <cstdio>
#include <fstream>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <chrono>
#include <random>

class Phase {

private:
    int volume;

    float sigma;
    float contrast;
    std::string id; // ID is a string identifier set in the input text file for multi-component analysis
    //   std::vector<float> statusVector;
    int workingLimit;
    int number_of_beads;
    float bead_radius;
    float electrons_per_bead;
    float sigma_number_of_beads;
    int startIndex;

public:
    Phase(int volume, float sigma, float contrast);
    Phase(int volume, float sigma, float contrast, std::string id);

    void setSigma(float sigma);
    float getSigma() const {return sigma;}

//    void setStatusVector(int totalBeads){ statusVector.resize(totalBeads,0); }

    std::string getID() const {return id;}
    float getContrast() const {return contrast;}

    void setWorkingLimit(int value){ workingLimit = value;}

    void setLimits(int startIndex, int workingLimit);

    int getWorkingLimit() const { return workingLimit;}
    int getStartIndex() const { return startIndex;}

    void setNumberOfBeads(float bead_radius);

    int getVolume() const {return volume;}

    float getElectronsPerBead() const { return electrons_per_bead; }

    int getNumberOfBeads() const { return number_of_beads;}
    float getSigmaNumberOfBeads() const {return sigma_number_of_beads;}
};

#endif //IKETAMA_PHASE_H
