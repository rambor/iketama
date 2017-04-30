//
// Created by Robert Rambo on 17/04/2017.
//

#ifndef IKETAMA_EULERTOURSMART_H
#define IKETAMA_EULERTOURSMART_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "SmartNode.h"
#include "Tour.h"
#include <algorithm>
#include <ctime>

class Model;
class SmartNode;


struct setcmp {
    bool operator()( Tour &tour1, Tour &tour2) {
        return (tour1.getRootNode() < tour2.getRootNode());
    }
};


class EulerTourSmart {

    int totalComponents;

    std::map<int, std::shared_ptr<SmartNode> > nodes; // maintain ownership of Node

    std::set< std::shared_ptr<Tour>, setcmp > tours;    // list of weak pointers, key is root of tour
    // collection of points, N, can have 1 to N tours

    //void createSubTour(SmartNode * pNode, Tour * subTourToLoad);

    // use a set of tours and custom comparator
    void createInitialTour(int workingLimit, Model *pModel, std::vector<int>::iterator beginIt);
    bool addToTour(int nodeToAdd);

    struct by_root_node_of_tour { // functor
        by_root_node_of_tour(int const & keytofind) : key(keytofind) { }
        bool operator () (Tour & p) { return p.getRootNode() == key; }
    private:
        int key;
    };

public:
    EulerTourSmart(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel);
    EulerTourSmart();

    int addNode(int latticePoint, Model *pModel);
    int removeNode(int indexOfNode);




};


#endif //IKETAMA_EULERTOURSMART_H
