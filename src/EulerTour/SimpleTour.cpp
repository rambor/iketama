//
// Created by Robert Rambo on 10/05/2017.
//

#include "SimpleTour.h"

SimpleTour::SimpleTour(Node * pNode) {
    tour.push_back(pNode);
    nodesInTour.insert(pNode);
    this->updateRootNode();
}


SimpleTour::SimpleTour(std::list<Node *> pTempTour) : tour(pTempTour) {
    // build set
    for(std::list<Node *>::iterator it = tour.begin(); it!= tour.end(); ++it){
        nodesInTour.insert(*it);
    }
    this->updateRootNode();
}


void SimpleTour::removeFromSet(Node * pNode){
   nodesInTour.erase(pNode);
}


void SimpleTour::updateRootNode() {
    std::set<Node *>::iterator it;
    for (it = nodesInTour.begin(); it != nodesInTour.end(); ++it) {
        //std::cout << "         SIZE of tour " << tour.size() << std::endl;
        (*it)->setPointerToTour(&tour);
        //std::cout << "   UPDATING ROOT NODE " << (*it)->getKey() << "  " << (*it)->getRootNodeOfTour() << std::endl;
    }
}


// subTour must be rooted to pNeighbor
void SimpleTour::insertSubTour(Node * pNeighbor, std::list<Node *> * subTour) {
    // subTour is rooted to pNode
    // pNeighbor => 6
    //               subTour => 58565
    //
    // pExistingNeighborTour => 124313426797621   (root node is 1)
    //
    // reroot        subTour => 65856
    //
    std::list< Node * >::iterator inTour = std::find (tour.begin(), tour.end(), pNeighbor);

    if(inTour == tour.end()){
        std::cout << "SPLING ERROR" << std::endl;
        exit(0);
    }

    //
    // locate 6 at ...267...
    //
    subTour->pop_back();

    // add to set
    for(std::list<Node *>::iterator it = subTour->begin(); it!= subTour->end(); ++it){
        nodesInTour.insert(*it);
    }
    //
    // pop_back      subTour => 6585
    // merge tour
    // merge subTour of newNode with existing tour
    tour.splice(inTour, *subTour); // Add additional tours to base
//    count = 0;
//    for(std::list<Node *>::iterator it=tour.begin(); it!=tour.end(); ++it){
//        std::cout << "      AFTER  TOUR " << count << " " << (*it)->getKey() << std::endl;
//        count++;
//    }
    // subTour should be empty now
    // nothing prevents inTour from pointing to end() and then having a phony splice, need to check that inTour never points to end
}


void SimpleTour::mergeTour(Node * pNeighbor, SimpleTour * pTour){

    // reroot pTour to pNeighbor
    pTour->rerootTour(pNeighbor);
    std::list< Node * > * pList = pTour->getPointerToTour();
    // find pNeighbor in this->tour
    std::list< Node * >::iterator inTour = std::find (tour.begin(), tour.end(), pNeighbor);
    // better to iterate over node list and check for tour membership
    // for all nodes in this tour, reset tourIndex
    //std::cout << "\n  ========     POINTERS " << pExistingNeighborTour << " " << &tours[rootToBaseTour] << " ++++\n"<< std::endl;
    // reassign before popping, in case pList size is 1
    pList->pop_back();
    // merge pList into existingNeighborTour
    tour.splice(inTour, *pList);

    //update nodes in pTour before erasing Tour
    for (std::set<Node *>::iterator it = pTour->getPointerToNodesInTour()->begin(); it != pTour->getPointerToNodesInTour()->end(); ++it) {
        (*it)->setPointerToTour(this->getPointerToTour());
        this->nodesInTour.insert(*it);
    }
}



/**
 * set all Nodes in the merged tour to point to this->tour
 */
void SimpleTour::mergeNodeIndices(SimpleTour * pTour){
    std::set<Node *>::iterator it;
    for (it = pTour->getPointerToNodesInTour()->begin(); it != pTour->getPointerToNodesInTour()->end(); ++it) {
        (*it)->setPointerToTour(&tour);
    }
}


void SimpleTour::rerootTour(Node * pRootNode){
    // 124313426797621 => reroot to 6
    std::list< Node * >::iterator inTour = std::find (tour.begin(), tour.end(), pRootNode);
    // inTour points to 6 at ..4267...

    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    if (tour.front()->getKey() != pRootNode->getKey() && (tour.size() > 2)){
        tour.pop_front();
        // subTourToLoad => 24313426797621
        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), tour, tour.begin(), inTour);
        // tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(pRootNode);
        //      tempList => 24313426
        // transfer components from tempList to subTour
        tour.splice(tour.end(), tempList);
        // subTourToLoad => 6797621-24313426
        // clear tempList
        tempList.clear();
    }
}