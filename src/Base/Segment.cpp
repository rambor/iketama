//
// Created by Robert Rambo on 14/09/2017.
//

#include "Segment.h"
#include "../Model.h"

using namespace std;

Segment::Segment(std::string id, int volume, Model *pModel, bool contiguous) {

    this->id = id;
    this->volume = volume;
    this->contiguous = contiguous;
    this->pModel = pModel;

    //
    if (volume > 0){
        this->upper_number_of_beads = std::ceil((double)volume/(double)pModel->getBeadVolume());
        this->lower_number_of_beads = std::ceil(this->upper_number_of_beads - 0.15*this->upper_number_of_beads);
    } else {
        this->fixed = true;
        this->upper_number_of_beads = -1*volume;
    }

    // is segment a single beads_in_use,
    // if multiplicity > 1, adding or removing a lattice point means simultaneous add/remove from all segments that define multiplicity

    // how are different phases in same Universe managed collectively?

    // adding a bead means selecting one that has not been selected in another phase
}