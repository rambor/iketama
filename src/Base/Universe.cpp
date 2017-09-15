//
// Created by Robert Rambo on 15/09/2017.
//

#include "Universe.h"

// create one universe per PRDATA file
Universe::Universe(std::string belongs_to_name, Data * pData, int searchSpace, std::string mode) {

    this->name = belongs_to_name;
    float bead_radius = 0.499999999999995*(pData->getBinWidth());
    interconnectivityCutOff = bead_radius*2.001;

    this->model = Model(searchSpace, bead_radius, true, mode);
}


// for each phase, need to make instance of Segment
void Universe::createSegment(std::string phaseName, int volumeOfSegment, bool contiguous){

    // how to handle multiplicity ?
    // Do I want to try to match them? For every bead I add/remove, do same to its cousins?
    segments.insert(std::pair<std::string, Segment> (phaseName, Segment(phaseName, volumeOfSegment, &(this->model), contiguous)));

}