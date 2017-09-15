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
#include <set>
#include "Segment.h"

class Model;

class Phase {

private:
    int volume;

    float contrast, sigma;
    std::string id; // ID is a string identifier set in the input text file for multi-component analysis
    //   std::vector<float> statusVector;
    int workingLimit;
    int number_of_beads;
    float bead_radius;
    float electrons_per_bead;
    float sigma_number_of_beads;
    int startIndex;
    int multiplicity=1;
    float bead_resolution=100.0f;
    bool contiguous=true;

    int indexOfBaseModel=0;

    /*
     * Phases manages all the Universes of a particular phase
     *
     * PRDATA => file1.dat NAME => A
     * PRDATA => file2.dat NAME => AB
     * PRDATA => file3.dat NAME => AC
     *
     * We should expect three phases to be specified
     * In this case phases A, B and C
     * A will be assigned to 3 universes (Model classes)
     *
     */



    // Phase class should be used to make the components it needs

public:

    Phase(int volume, float sigma, float contrast);

    Phase(int volume, float sigma, float contrast, int multiplicity, std::string id, bool contiguous);

    // create component based on volume and multiplicity

    // if volume is negative and small, we fix the number of indices for the component


    // void setStatusVector(int totalBeads){ statusVector.resize(totalBeads,0); }

    std::string getID() const {return id;}

    float getContrast() const {return contrast;}

    void setNumberOfBeads(float bead_radius);

    int getVolume() const {return volume;}

    float getElectronsPerBead() const { return electrons_per_bead; }

    // methods are add/remove, move or swap
    // move is to find a new position that is not selected
    // swap is find another bead that is occupied by a different phase and swap


};

#endif //IKETAMA_PHASE_H
