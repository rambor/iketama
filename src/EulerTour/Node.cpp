//
// Created by Robert Rambo on 11/01/2017.
//

#include "Node.h"

// Constructor
Node::Node(){};

// Constructor
Node::Node(int key) {
    this->key = key;
    totalNeighbors = 0;
    adjacencyList.reserve(40); // size of possible neighbors in neighborhood
    //adjacencyList.resize(12);
}


void Node::addNeighbor(Node * pNode) {
    adjacencyList.push_back(pNode);
    //adjacencyList[totalNeighbors] = pNode;
    adjacencyListIndex.insert(pNode->getKey()); // set can never have duplicated values
    totalNeighbors++;

    if (!(pNode->isNeighborPresent(this->getKey()))){
        pNode->addNeighbor(this);
    }
}


bool Node::isNeighborPresent(int index){
    if (adjacencyListIndex.find(index) != adjacencyListIndex.end()) {
        return true;
    }
    return false;
}


void Node::removeNeighbor(Node * pNode) {

    for(int i=0;i < totalNeighbors; i++){
        //std::cout << this->getKey() << " - " << pNode->getKey() << " " << adjacencyList[i]->getKey() << " " << (adjacencyList[i] == pNode) << std::endl;
        //if (adjacencyList[i]->getKey() == pNode->getKey()){
        if (adjacencyList[i] == pNode){
            adjacencyList.erase(adjacencyList.begin()+i);
            adjacencyListIndex.erase(pNode->getKey());
            totalNeighbors--;
            break;
        }
    }

    // remove itself from pNode's neighborhood list
    if (pNode->isNeighborPresent(this->getKey())){
        pNode->removeNeighbor(this);
    }
}

void Node::setFirst(std::list<Node>::iterator *pFirst) {
    this->pFirst = pFirst;
}

void Node::setLast(std::list<Node>::iterator *pLast) {
    this->pLast = pLast;
}

bool Node::validate(){

    if (adjacencyList.size() != adjacencyListIndex.size()){
        std::cout << "Possible repeats "<< std::endl;
        std::cout << "     ADJACENCYLIST : " << adjacencyList.size() << " " << totalNeighbors << std::endl;
        std::cout << "ADJACENCYLISTINDEX : " << adjacencyListIndex.size() << totalNeighbors << std::endl;

        for(int i=0;i < adjacencyList.size(); i++){
            std::cout << i << " " << adjacencyList[i]->getKey() << std::endl;
        }

        return false;
    }

    return true;
}
