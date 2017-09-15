//
// Created by Robert Rambo on 07/08/2017.
//

#include "MultiDataComponent.h"


// have to set the data with the highest resolution as the base universe
// each dataset gets its own universe




MultiDataComponent::MultiDataComponent(std::string setup_filename) {

    this->setupFile = setup_filename;

    // parse the file and add create Universes (based on PRDATA)

    /*
     * read setup file
     * Steps:
     * 1. parse the file for PRDATA and create Universe for each PRDATA
     * 2. parse the file and make a phase for each
     */

}


void MultiDataComponent::addDataset(std::string filename, std::string name) {

    // open file
    // create IofQ object
    datasets.insert(pair<std::string, Data> (name, Data()));
    datasets[name].addPofRData(filename);

    //
}





void MultiDataComponent::createInitialModel() {



    // for each round
    //    randomize phases list
    //    for through each phase
    //         1. add or remove
    //         2. positional
    //         3. swap

    // for each phase in each universe, make randomize selection limited by upperlimit
    // * phase is not aware of the other phases
    // * universe has list of beads_in_use


}