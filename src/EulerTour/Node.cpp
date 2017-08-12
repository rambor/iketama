//
// Created by Robert Rambo on 11/01/2017.
//

#include "Node.h"

// Constructor
/**
 * CAN NOT HAVE A NODE WITH NO KEY!
 */
Node::Node(int key) : key(key) {
    //this->key = key;
    totalNeighbors = 0;
    adjacencyList.resize(40);
}


void Node::addNeighbor(Node * pNode) {

    //adjacencyList.push_back(pNode);
    adjacencyList[totalNeighbors] = pNode;
    adjacencyListIndex.insert(pNode->getKey()); // set can never have duplicated values
    totalNeighbors++;

//    if (adjacencyList.size() != totalNeighbors){
//        std::cout << "addNeighbor NODE ERROR => DO NOT EQUAL " << adjacencyList.size() <<  " != " << totalNeighbors << std::endl;
//        exit(0);
//    }

    if (!(pNode->isNeighborPresent(this->getKey()))){
        pNode->addNeighbor(this);
        //std::cout << this->getKey() << " ADDING NEIGHBOR => " << pNode->getKey() << std::endl;
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

        if (adjacencyList[i] == pNode){ // do they have the same address (comparing pointers)
            //adjacencyList.erase(adjacencyList.begin()+i);
            // move to end and decrement
//            if (adjacencyListIndex.size() != totalNeighbors){
//                std::cout << "SIZES DO NOT MATCH NODE::removeNeighbor " << adjacencyListIndex.size() << " != " << totalNeighbors << std::endl;
//                exit(0);
//            }
            // if we need to remove it, swap to the last (in use position), and decrement count
            std::iter_swap(adjacencyList.begin()+i, adjacencyList.begin()+totalNeighbors-1);

//            printAdjacencyList("FROM REMOVE NEIGHBOR ");
            adjacencyListIndex.erase(pNode->getKey()); //vector
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
        std::cout << "     ADJACENCYLIST : " << adjacencyList.size() << " " << totalNeighbors << std::endl;
        std::cout << "ADJACENCYLISTINDEX : " << adjacencyListIndex.size() << totalNeighbors << std::endl;

        for(int i=0;i < totalNeighbors; i++){
            std::cout << i << " " << adjacencyList[i]->getKey() << std::endl;
        }

        return false;
    }

    return true;
}


void Node::printNeighbors(){
    int count=1;
    std::cout << " SIZE OF ADJACENCY LIST " << adjacencyList.size() << " " << totalNeighbors << std::endl;
    for(int i=0; i<totalNeighbors; i++){
        std::cout << "      " << count << " NEIGHBOR OF " << key << " " << adjacencyList[i]->getKey() << std::endl;
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
void Node::setPointerToTour(std::list < Node * > * pointer){
    pointerToTour = pointer;
    //std::cout << key << " ROOT : SETTING POINT TO TOUR ( SIZE : " << pointer->size() << " ) "<<  std::endl;
    rootNode = (*pointerToTour).front()->getKey();
    //rootNode = (*pointerToTour).front()->key;
    //std::cout << "        key : " << this->getKey() << "  ROOT NODE SET => " << rootNode << std::endl;
}
