//
// Created by Robert Rambo on 14/09/2017.
//

#ifndef IKETAMA_SEGMENT_H
#define IKETAMA_SEGMENT_H

#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

class Model;

class Segment {

    std::string id;
    int volume;
    int upper_number_of_beads;
    int lower_number_of_beads;
    bool contiguous;
    bool fixed = false;
    Model * pModel;  // for the euler tour
    std::vector< std::set<int>> beads_in_use;

public:

    Segment(std::string id, int volume, Model *pModel, bool contiguous);

    void createInitialModel(std::set<int> total_beads_in_use);

};


#endif //IKETAMA_SEGMENT_H
