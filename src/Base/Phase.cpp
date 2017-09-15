//
// Created by Robert Rambo on 12/01/2017.
//


#include "Phase.h"
//#include "../Model.h"

using namespace std;

Phase::Phase(int volume, float sigma, float contrast){
    this->volume = volume;
    this->contrast = contrast;
    this->sigma = sigma;
}


Phase::Phase(int volume, float sigma, float contrast, int multiplicity, std::string id, bool contiguous) : Phase(volume, sigma, contrast) {
    this->multiplicity = multiplicity;
    this->contiguous = contiguous;
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
     * PRDATA => file1.dat NAME => A
     * PRDATA => file2.dat NAME => AB
     * PRDATA => file3.dat NAME => AC
     *
     * We should expect three phases to be specified
     * In this case phases A, B and C
     * A will be assigned to 3 universes (Model classes)
     *
     * PHASE => A VOL => X00000 CONTRAST => 1 MULTIPLITICY => 1 CONTINGUOUS => true BELONGS_TO => A BELONGS_TO => AB BELONGS_TO => AC
     * PHASE => B VOL => X0000 CONTRAST => 1  MULTIPLITICY => 1 CONTINGUOUS => true  BELONGS_TO => AB
     * PHASE => C VOL => X0000 CONTRAST => 1  MULTIPLITICY => 1 CONTINGUOUS => true  BELONGS_TO => AC
     *
     *
 */
//void Phase::addSegment(std::string phaseName, Model * pBelongs_To) {
//
//    // volume divided by the datasets resolution determines upper number of beads
//
//    /*
//     * if multiplicity is > 1, volume should represent the subunit which will be multiplied
//     * Component(std::string id, float volume, Model *pModel, bool contiguous);
//     *
//     * can not have multiplicity > 1 and contiguous
//     */
//    segments.push_back(Segment(phaseName, volume, multiplicity, pBelongs_To, contiguous));
//
//    // need to set the highest resolution Model
//    if (pBelongs_To->getBeadRadius() < bead_resolution ){
//        bead_resolution = pBelongs_To->getBeadRadius();
//        indexOfBaseModel = segments.size()-1;
//    }
//}