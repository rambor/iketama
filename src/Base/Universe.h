//
// Created by Robert Rambo on 15/09/2017.
//

#ifndef IKETAMA_UNIVERSE_H
#define IKETAMA_UNIVERSE_H
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include "Phase.h"
#include "../Model.h"


class Universe {

    Model model;
    std::string name;
    float interconnectivityCutOff;
    std::set<int> total_beads_in_use;

    std::map <std::string, Segment> segments;


public:

    /**
     * name is derived from NAME that matches PRDATA
     */
    Universe(std::string belongs_to_name, Data * pData, int searchSpace, std::string mode);

    void createSegment(std::string phaseName, int volumeOfSegment, bool contiguous);

};


#endif //IKETAMA_UNIVERSE_H
