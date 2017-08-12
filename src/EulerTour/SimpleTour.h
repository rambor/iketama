//
// Created by Robert Rambo on 10/05/2017.
//

#ifndef IKETAMA_SIMPLETOUR_H
#define IKETAMA_SIMPLETOUR_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "Node.h"
#include <ctime>

struct setcmpNodes {
    bool operator()( Node * node1, Node * node2) {
        return (node1->getKey() < node2->getKey());
    }
};

class SimpleTour {

    std::list< Node *> tour;
    std::set< Node *, setcmpNodes > nodesInTour;

public:
    SimpleTour(Node * pNode); // root of tour should be key
    SimpleTour(std::list< Node *> pTempTour); // root of tour should be key

    std::list< Node *> * getPointerToTour() { return &tour; }
    int getRootNodeOfTour(){ return tour.front()->getKey(); }

    std::set< Node *> * getPointerToNodesInTour(){ return (std::set<Node *> *) &nodesInTour;}

    void removeFromSet(Node * pNode);
    void updateRootNode();

    void mergeTour(Node * pNeighbor, SimpleTour * pTour);
    void insertSubTour(Node * pNeighbor, std::list< Node * > * subTour);

    void rerootTour(Node * pRootNode);
    void mergeNodeIndices(SimpleTour * pTour);
    // merging tours would require merging sets and lists
};


#endif //IKETAMA_SIMPLETOUR_H
