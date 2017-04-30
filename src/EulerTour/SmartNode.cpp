//
// Created by Robert Rambo on 17/04/2017.
//

#include "SmartNode.h"
#include "Tour.h"
// Constructor
//SmartNode::SmartNode(){};

// Constructor
SmartNode::SmartNode(int key) {
    this->key = key;
    totalNeighbors = 0;
    //adjacencyList.reserve(40); // size of possible neighbors in neighborhood
    adjacencyList.resize(40);
}

void SmartNode::setPointerToTour(std::shared_ptr<Tour> pointer){
    pointerToTour = pointer;
}

int SmartNode::getRootNodeOfTour(){
    return pointerToTour.lock().get()->getRootNode();
}


void SmartNode::addNeighbor(const std::shared_ptr<SmartNode> & pNode) {

    adjacencyList[totalNeighbors] = pNode; // makes weak pointer from pNode since adjacencyList is a vector of weak_ptr
    adjacencyListIndex.insert(pNode->getKey()); // set can never have duplicated values
    totalNeighbors++;

    // add this->Node if not present in Neighbor
    if (!(pNode->isNeighborPresent(this->getKey()))){ // need a weak ptr from this?
        pNode->addNeighbor(shared_from_this());
    }
}

//


bool SmartNode::isNeighborPresent(int index){
    if (adjacencyListIndex.find(index) != adjacencyListIndex.end()) {
        return true;
    }
    return false;
}



void SmartNode::removeNeighbor(const std::shared_ptr<SmartNode> & pNode) {

    for(int i=0;i < totalNeighbors; i++){

        if (adjacencyList[i].lock() == pNode){ // do they have the same address (comparing pointers)
            // move to end and decrement
            // if we need to remove it, swap to the last (in use position), and decrement count
            std::iter_swap(adjacencyList.begin()+i, adjacencyList.begin()+totalNeighbors-1);
            adjacencyListIndex.erase(pNode->getKey()); //vector
            totalNeighbors--;
            break;
        }
    }

    // remove itself from pNode's neighborhood list
    if (pNode->isNeighborPresent(this->getKey())){
        pNode->removeNeighbor(shared_from_this());
    }
}