//
// Created by Robert Rambo on 11/01/2017.
//

#ifndef IKETAMA_NODE_H
#define IKETAMA_NODE_H

#include <string>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <iostream>

class Node {

    int key;
    std::vector< Node * > adjacencyList;
    std::set< int > adjacencyListIndex;
    int totalNeighbors = 0;
    //int tour_list_index;
    int rootNode;
    std::list < Node * > * pointerToTour;

//    std::list < Node >::iterator * pFirst;
//    std::list < Node >::iterator * pLast;
    bool accessed = false;

    void printAdjacencyList(std::string text);

public:
    //Node();
    Node(int key);

    ~Node(){
        //std::cout << key << " ADJACENCY LIST SIZE " << totalNeighbors << std::endl;
        for (auto it = adjacencyList.begin() ; it != adjacencyList.end(); ++it) {
            *it = NULL;
            it = adjacencyList.erase(it);
        }
        adjacencyListIndex.clear();
        adjacencyList.clear();
        pointerToTour = NULL;
    }

    int getKey(){return key;}
    void addNeighbor(Node * pNode);
    void removeNeighbor(Node * pNode);
    void printNeighbors();

//    void setFirst(std::list < Node >::iterator * pFirst);
//    void setLast(std::list < Node >::iterator * pLast);
    void setAccessed(bool value){this->accessed = value;}

    void setPointerToTour(std::list < Node * > * pointer);//{ pointerToTour = pointer;}
    std::list < Node * > * getPointerToTour(){ return pointerToTour;}

    //void setRootNodeOfTour(int index){rootNode = index;}
    int getRootNodeOfTour(){return rootNode;}

    bool getAccessed(){return this->accessed;}

    int getTotalNeighbors(){return totalNeighbors;}
    Node * getPointerToNeighborByIndex(int index){return adjacencyList[index];}
    bool isNeighborPresent(int index);

    std::set< int > * getIteratorToIndices(){return &adjacencyListIndex;}
//    std::list < Node >::iterator * getFirst(){return pFirst;}
//    std::list < Node >::iterator * getLast(){return pLast;}

    bool validate();
};

#endif //IKETAMA_NODE_H
