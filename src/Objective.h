//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_OBJECTIVE_H
#define IKETAMA_OBJECTIVE_H


#include <string>
#include <vector>
#include <map>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <iostream>
//#include "Base/Data.h"
#include "Model.h"
#include "Data.h"

class Objective {

    // HASH KEY => DATA
    //std::unordered_map <std::string, Data> datasets;
    std::map <std::string, Data> datasets;
    std::vector<std::string> keys;
    std::vector<float> weights;


public:
    Objective();
    ~Objective(){
        // remove all
//        std::map<std::string, Data>::iterator it = datasets.begin();
//        while(it != datasets.end()){
//            it = datasets.erase(it);
//        }

        // clear pointers to nodes
//        for (std::map<std::string, Data>::iterator it=datasets.begin(); it!=datasets.end(); ++it){
//            // delete points in tour
//            delete(*it);
//            *it = NULL;
//        }

        datasets.clear();


    }

    void addDataObject(std::string iofqfile, std::string pofrfile);

    void addModel(std::string iofqfileName, Model * model);

    //void createWorkingSets(Model & model);

    Data * getDataObject (std::string key) {return &datasets[key];}
    Data * getDataObjectByIndex (int index) {return &datasets[keys[index]];}

    Data * getLastDataset();
    Data * getMainDataset();

    int getTotalDatasets(){ return datasets.size();}

};



#endif //IKETAMA_OBJECTIVE_H
