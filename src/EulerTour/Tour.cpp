//
// Created by Robert Rambo on 17/04/2017.
//

#include "Tour.h"

Tour::Tour(int key){
    this->rootNode = key;
}


Tour::Tour(std::weak_ptr<SmartNode> nodeToAdd){
    createSubTour(nodeToAdd);
    //this->rootNode = nodeToAdd.lock().get()->getKey();
}

void Tour::addToTour(std::shared_ptr<SmartNode> nodeToAdd){
    eulerTour.push_back(nodeToAdd); // casts nodeToAdd to weak_ptr
}


void Tour::reRootTour(std::shared_ptr<SmartNode> newRootNode){
    // find root in
    int newNode = newRootNode.get()->getKey();

    // search tour for
    std::list< std::weak_ptr<SmartNode> >::iterator inTour = std::find_if(eulerTour.begin(), eulerTour.end(), by_root_node_of_tour(newNode));
    std::list< std::weak_ptr<SmartNode> > tempList;

    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    if ( eulerTour.front().lock().get()->getRootNodeOfTour() != newNode && (eulerTour.size() > 2)){

        eulerTour.pop_front();

        tempList.splice(tempList.begin(), eulerTour, eulerTour.begin(), inTour);
        tempList.push_back(newRootNode);
        // transfer components from tempList to subTour
        eulerTour.splice(eulerTour.end(), tempList);
    }

    //update nodes pointer to tour
    for(std::list< std::weak_ptr<SmartNode> >::iterator it = eulerTour.begin(); it != eulerTour.end(); ++it){
//        (*it).lock().get()->setPointerToTour(this);
    }

    // clear tempList
    tempList.clear();
}


/**
 * Creates proper Euler Tour where first and last Node are the same
 */
void Tour::createSubTour(std::weak_ptr<SmartNode> nodeToAdd){

    eulerTour.push_back(nodeToAdd); // casts nodeToAdd to weak_ptr

    int totalNeighbors = nodeToAdd.lock().get()->getTotalNeighbors();
    int totalIterations = 2*totalNeighbors + 1;
    int index=0;

    // go through each neighbor of pNode and add to new Tour
    for(int i=1; i<totalIterations; i++){
        if (i%2 == 0){
            eulerTour.push_back(nodeToAdd);
        } else {
            // add from neighborhood and not Node
            eulerTour.push_back(nodeToAdd.lock().get()->getPointerToNeighborByIndex(index));
            index++;
        }
    }
}

void Tour::popOff(){
    this->eulerTour.pop_back();
}


void Tour::mergeTours(std::weak_ptr<SmartNode> insertionPoint, Tour & subTour){

    subTour.reRootTour(insertionPoint.lock().get());
    subTour.popOff(); // popoff

    //update nodes in subTour to point to tour
    int root = this->getKeyOfRootNode();
    for(std::list< std::weak_ptr<SmartNode> >::iterator it = subTour.getStartOfTour(); it != subTour.getEndOfTour(); ++it){
        if ( (*it).lock().get()->getRootNodeOfTour() != root){
            (*it).lock().get()->setPointerToTour(this);
        }
    }

    // find insertion point in this tour
    std::list< std::weak_ptr<SmartNode> >::iterator inTour = std::find_if (eulerTour.begin(), eulerTour.end(), by_root_node_of_tour(insertionPoint.lock().get()->getKey()));

    eulerTour.splice(inTour, *subTour.getTour());
}


const int Tour::getRootNode() {return eulerTour.front().lock().get()->getKey(); }
const int Tour::getKeyOfRootNode(){return eulerTour.front().lock().get()->getKey(); }