//
// Created by Robert Rambo on 11/01/2017.
//

#include "Node.h"

// Constructor
/**
 * CAN NOT HAVE A NODE WITH NO KEY!
 */

Node::Node(int key) : key(key) {
    totalNeighbors = 0;
    adjacencyList.resize(40);
}


void Node::addNeighbor(Node * pNode) {

        adjacencyList[totalNeighbors] = pNode->getKey();
        adjacencyListIndex.insert(pNode->getKey()); // set can never have duplicated values
        totalNeighbors++;

    if (!(pNode->isNeighborPresent(this->getKey()))){
        pNode->addNeighbor(this);
    }
}



void Node::removeNeighbor(Node * pNode) {

    for(int i=0;i < totalNeighbors; i++){
        //std::cout << this->getKey() << " - " << pNode->getKey() << " " << adjacencyList[i]->getKey() << " " << (adjacencyList[i] == pNode) << std::endl;
        if (adjacencyList[i] == pNode->getKey()){ // do they have the same address (comparing pointers)
            // adjacencyList.erase(adjacencyList.begin()+i);
            // move to end and decrement
            // if we need to remove it, swap to the last (in use position), and decrement count
            //if (i < (totalNeighbors - 1)){
                std::iter_swap(adjacencyList.begin()+i, adjacencyList.begin()+totalNeighbors-1);
            //}

            adjacencyListIndex.erase(pNode->getKey()); // remove from set
            totalNeighbors--;
            break;
        }
    }

    // remove itself from pNode's neighborhood list
    if (pNode->isNeighborPresent(this->getKey())){
        pNode->removeNeighbor(this);
    }
}

//void Node::setFirst(std::list<Node>::iterator *pFirst) {
//    this->pFirst = pFirst;
//}
//
//void Node::setLast(std::list<Node>::iterator *pLast) {
//    this->pLast = pLast;
//}

bool Node::validate(){

    if (totalNeighbors != adjacencyListIndex.size()){
        std::cout << "Possible repeats "<< std::endl;
        std::cout << "           ADJACENCYLIST : " << adjacencyList.size() << " >= " << totalNeighbors << std::endl;
        std::cout << "ADJACENCYLISTINDEX (SET) : " << adjacencyListIndex.size() << " != " << totalNeighbors << std::endl;

        for(int i=0; i < totalNeighbors; i++){
            std::cout << i << " " << adjacencyList[i] << std::endl;
        }

        return false;
    }
    //this->printNeighbors();

    return true;
}


void Node::printNeighbors(){
    int count=1;
    std::cout << " SIZE OF ADJACENCY LIST " << adjacencyList.size() << " " << totalNeighbors << " "  << std::endl;

    for(int i=0; i<totalNeighbors; i++){
        std::cout << "      " << count << " NEIGHBOR OF " << key << " " << adjacencyList[i] << " " << totalNeighbors << std::endl;
        count++;
    }

    // print adjacency list
    this->printAdjacencyList("FROM PRINTNEIGHBORS");
}


void Node::printAdjacencyList(std::string text){

    std::cout << "      ADJACENCY LIST " << std::endl;
    int count=1;
    for(std::set<int>::iterator it= adjacencyListIndex.begin(); it!=adjacencyListIndex.end(); ++it){
        std::cout << "      " << count << " ADJACENCY LIST : " << key << " " << *it << std::endl;
        count++;
    }
}


// setting pointer to tour should override previous root node
// void Node::setPointerToTour(std::list < Node * > * pointer){ // point to a location in map
//    pointerToTour = pointer;
//    rootNode = (*pointerToTour).front()->getKey();
////    std::cout << key << "          ROOT : SETTING POINT TO TOUR ( SIZE : " << pointer->size() << " ) "<< rootNode << std::endl;  // no size means list is empty
//}

bool Node::isNeighborPresent(int index) {

    if (adjacencyListIndex.find(index) != adjacencyListIndex.end()){
        return true;
    }
    return false;
}
