//
// Created by Robert Rambo on 17/04/2017.
//

#ifndef IKETAMA_SMARTNODE_H
#define IKETAMA_SMARTNODE_H

#include <string>
#include <vector>
#include <list>
#include <string>
#include <memory>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <iostream>

class Tour;


class SmartNode : public std::enable_shared_from_this<SmartNode> {
    int key;
    std::vector< std::weak_ptr<SmartNode> > adjacencyList;

    std::set< int > adjacencyListIndex;
    int totalNeighbors = 0;
    //int tour_list_index;
    //std::list < SmartNode * > * pointerToTour;
    std::weak_ptr<Tour> pointerToTour;
    bool accessed = false;

public:

    //SmartNode();
    SmartNode(int key);
    ~SmartNode(){
        //std::cout << key << " ADJACENCY LIST SIZE " << totalNeighbors << std::endl;
        for (auto it = adjacencyList.begin() ; it != adjacencyList.end(); ++it) {
            it = adjacencyList.erase(it);
        }
        adjacencyList.clear();
    }

    int getKey(){return key;}
    int getTotalNeighbors(){return totalNeighbors;}

    void addNeighbor(const std::shared_ptr<SmartNode> & pNode);
    void removeNeighbor(const std::shared_ptr<SmartNode> & pNode);

    std::weak_ptr<SmartNode> getPointerToNeighborByIndex(int index){return adjacencyList[index];}

    bool isNeighborPresent(int index);

    // if a Node is changed to a new tour, these have to be updated
    void setPointerToTour(std::shared_ptr<Tour> pointer);
    std::weak_ptr<Tour> getPointerToTour(){ return pointerToTour;}

    int getRootNodeOfTour();
};


#endif //IKETAMA_SMARTNODE_H
