//
// Created by Robert Rambo on 15/02/2017.
//

#ifndef IKETAMA_COMPONENT_H
#define IKETAMA_COMPONENT_H
#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <iostream>
#include <set>
#include "../EulerTour/EulerTour.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <qhull_a.h>
#include <libqhull.h>
#ifdef __cplusplus
}
#endif

class Model;

class Component {

private:

    std::string id;
    std::set<int> anchors;            // can be empty
    std::vector<int> resids;          // can be empty
    std::vector<std::string> chains;  // can be empty
    std::set<int> beads_in_use;       // never empty
    std::set<int> best;       // never empty
    std::vector<int> cvxPoints;       // never empty
    EulerTour tour;
    Model * pModel;  // for the euler tour
    int totalComponents;
    int targetNumberOfBeads=0;
    float volume;
    bool empty = true;

    float totalContactSum;

public:

    Component(std::string id, float volume, Model *pModel);

    ~Component(){
    }

    bool checkID(std::string str);
    float numberOfLatticePointsPotential(float value);
    int removeLatticePoint(int index);
    int addLatticePoint(int index); // add to beads_in_use and

    bool indexInUse(int index);

    void addResid(int index, std::string chain);

    std::string const getID() const { return id;}
    int getTotalResids() { return resids.size();}
    int getResidByIndex(int index) { return resids[index];}
    std::string getChainsByIndex(int index) { return chains[index];}
    void addAnchor(int index);

    // desired number of beads but exclude the anchors
    void setTargetNumberOfLatticePoints(float value);
    int getTargetNumberOfLatticePoints();

    bool anchorsEmpty(){ return empty;}

    std::set<int> * getAnchors(){ return &anchors; }
    std::set<int> * getBeadsInUse(){ return &beads_in_use; }
    int getTotalNumberOfComponents(){ return tour.getNumberOfComponents();}

    bool inUse(int index){ return !(beads_in_use.find(index) == beads_in_use.end()); }
    bool isAnchor(int index);
    float potential();
    void printAnchors();
    void printSet();
    void printBest();
    void writeToFile(std::string nameOf);

    void copyBestToInUse();
    float calculateCVXVolume();
    void setBest();

    void setTotalContactSum(float value){ totalContactSum = value;}
    float getTotalContactSum(){ return totalContactSum;}
    void printConstraints();
};


#endif //IKETAMA_COMPONENT_H
