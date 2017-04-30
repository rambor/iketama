//
// Created by Robert Rambo on 17/04/2017.
//

#ifndef IKETAMA_TOURS_H
#define IKETAMA_TOURS_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "SmartNode.h"
#include <ctime>

class SmartNode;


class Tour {

    struct by_root_node_of_tour : std::unary_function<std::weak_ptr<SmartNode>, bool>{
        by_root_node_of_tour(int keyToFind) : key(keyToFind){}
        bool operator () (std::weak_ptr<SmartNode> p) { return p.lock().get()->getKey() == key; }
    private:
        int key;
    };

    int rootNode;

    std::list< std::weak_ptr<SmartNode> > eulerTour;    // list of weak pointers, key is root of tour

    void createSubTour(std::weak_ptr<SmartNode> nodeToAdd);
    void updateRootNodes(int indexOfRoot);

    void popOff();

public:
    Tour(int key); // root of tour should be key
    Tour(std::weak_ptr<SmartNode> nodeToAdd);

    const int getRootNode();
    const int getKeyOfRootNode();

    void addToTour(std::shared_ptr<SmartNode> nodeToAdd);

    std::list< std::weak_ptr<SmartNode> > * getTour(){ return &eulerTour;}

    std::list< std::weak_ptr<SmartNode> >::iterator getStartOfTour(){ return eulerTour.begin();}
    std::list< std::weak_ptr<SmartNode> >::iterator getEndOfTour(){ return eulerTour.end();}

    int getSizeOfTour(){return eulerTour.size();}
    void reRootTour(std::shared_ptr<SmartNode> newRootNode);

    void mergeTours(std::weak_ptr<SmartNode> insertionPoint, Tour & subTour);
};



#endif //IKETAMA_TOURS_H
