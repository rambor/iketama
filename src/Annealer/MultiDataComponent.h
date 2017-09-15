//
// Created by Robert Rambo on 07/08/2017.
//

#ifndef IKETAMA_MULTIDATACOMPONENT_H
#define IKETAMA_MULTIDATACOMPONENT_H

#include <string>
#include <vector>
#include <bead.h>
#include <functions.h>
#include <Data.h>
#include <math.h>
#include <regex>
#include <iostream>

// create single universe using smallest possible lattice dimensions
// each model belongs to one or more distributions

class MultiDataComponent {

    std::vector<Component> components;

    Objective datasets;
    std::map<std::string, Universe> universes;
    std::map<std::string, Phase> phases;

    std::string setupFile;

public:

    MultiDataComponent(std::string setupfilename);

    void addDataset(std::string filename, std::string name);

    bool createInitialModel();

};


#endif //IKETAMA_MULTIDATACOMPONENT_H
