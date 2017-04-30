//
// Created by Robert Rambo on 17/04/2017.
//
#include "EulerTourSmart.h"
#include "../Model.h"

// Constructor
EulerTourSmart::EulerTourSmart(){
}

EulerTourSmart::EulerTourSmart(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel){
    //pSelectedLattice = &beginIt;
    totalComponents=0;
    this->createInitialTour(workingLimit, pModel, beginIt);
}

/**
 * critical, new node must already be added to nodes!
 */
void EulerTourSmart::createInitialTour(int workingLimit, Model *pModel, std::vector<int>::iterator beginIt) {

    // for each node, make adjacency list
    std::vector<int>::iterator it;
    int neighbor, nodeToInsert;
    SmartNode * pNode;


    for(int i=0; i < workingLimit; i++){ // iterate over selected beads and make Node and add to nodes
        // make new Node from lattice/bead
        nodeToInsert = *(beginIt + i);
        //nodes.insert ( std::pair<int,Node>(nodeToInsert, Node(nodeToInsert) ) );

        this->addNode(nodeToInsert, pModel);

        // won't need this anymore, moved to addNode
//        nodes.insert ( std::pair<int, std::shared_ptr<SmartNode> > (nodeToInsert, std::make_shared<SmartNode>(nodeToInsert) )); // make object on HEAP
//        pNode = nodes[nodeToInsert].get();
//
//        // create the neighborhood for nodes in nodes list
//        // as I add each node, create neighborhood
//
//        it = pModel->getPointerToNeighborhood(pNode->getKey());
//        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
//            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
//            neighbor = *(it+j);
//
//            // if not found, itIndex will report last
//            if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){
//                nodes[neighbor].get()->addNeighbor(nodes[nodeToInsert]);
//                // call neighbor and add pNode as neighbor
//            } else if (neighbor == -1) {
//                break;
//            }
//        }
//
//        addToTour(nodeToInsert);
        // create subtour
    } // end of adding beads

    // make Nodes
    // then for each Node, create neighborhood?

    std::cout << " INITIAL TOURS SIZE " << tours.size() << std::endl;
    //totalComponents = tours.size();
}



/**
 * add a node to existing set of nodes
 * first, create neighborhood from existing nodes
 * find which tour to tours to add to and merge any tours that are bridged by new node
 * returns number of nodes in tour
 */
int EulerTourSmart::addNode(int newNode, Model *pModel) {

    //std::vector<int>::iterator itIndex;
    int neighbor;

    nodes.insert ( std::pair<int, std::shared_ptr<SmartNode> > (newNode, std::make_shared<SmartNode>(newNode) )); // make object on HEAP
    SmartNode * pNode = nodes[newNode].get();

    // build neighborhood (check from existing nodes or beads in use)
    auto it = pModel->getPointerToNeighborhood(newNode);
    for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
        // if neighbor is inside workinglimit, don't add
        neighbor = *(it+j); // index of bead used as key
        // if not found, itIndex will report last
        if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){ // look in nodes list
            //nodes[neighbor].get()->addNeighbor(nodes[nodeToInsert]);
            pNode->addNeighbor(nodes[neighbor]); // add neighbor to both parent and child
        } else if (neighbor == -1){
            break;
        }
    }

    // add to Tour
    //addToTour(newNode);

    totalComponents = tours.size();
    return totalComponents;
}


/**
 * Find node in nodes list
 * Remove node from neighbors and break into subTours if necessary
 * for each of its neighbor in neighborhood
 * if removing node creates a new tour, add to tours
 */
int EulerTourSmart::removeNode(int indexOfNode){

    SmartNode * pNodeToRemove = nodes[indexOfNode].get(), * pNeighbor;
    std::vector<int> checkedNodes;

    std::list< SmartNode * > * pTour;// = &tours[pNodeToRemove->getRootNodeOfTour()];  // tour with node to remove

    // delete pointer in vector
    // std::string result;

    // if node has neighbors, need to remove itself from each neighbor and recalculate tours

//    if (pNodeToRemove->getTotalNeighbors() > 0){
//        // remove node
//        std::unordered_set <int> remainderNodes;
//
//        while (pNodeToRemove->getTotalNeighbors() > 0) {  // for each neighbor, remove edge with node_to_remove
//
//            pNeighbor = pNodeToRemove->getPointerToNeighborByIndex(0);
//            pTour = &tours[pNodeToRemove->getRootNodeOfTour()];  // get tour that contains node to remove
//            this->rerootAndReassignSubTour(pNeighbor, pTour);    // reroot and reassign subTour in Tours which means pointer changes
//            // pTour may no longer exists
//            pTour = &tours[pNeighbor->getRootNodeOfTour()];      // reassign pointer
//
//            std::list< Node * >::iterator pFirst = getFirstOccurrenceInTour(pNodeToRemove, pTour);
//            std::list< Node * >::reverse_iterator pLast = getLastOccurrenceInTour(pNodeToRemove, pTour);
//            // search from start iterate through
//            bool terminateNow = false;
//            while ( (pFirst = std::find(pFirst, pTour->end(), pNodeToRemove)) != pTour->end()) {
//                // Do something with iter
//                auto pPreviousFirst = std::next(pFirst, -1);
//                auto pNextFirst = std::next(pFirst, 1);
//
//                if (((*pPreviousFirst)->getKey() == pNeighbor->getKey()) || ((*pNextFirst)->getKey() == pNeighbor->getKey())){ //
//                    // search for Neighbor in reverse
//                    while ( (pLast = std::find(pLast, pTour->rend(), pNodeToRemove)) != pTour->rend()) { // search for Neighbor in reverse
//                        // Do something with iter
//                        auto pPreviousLast = std::next(pLast, 1);
//                        auto pNextLast = std::next(pLast, -1);
//                        // from start
//                        if (((*pPreviousLast)->getKey() == pNeighbor->getKey()) || ((*pNextLast)->getKey() == pNeighbor->getKey())){ // if true, implies edge found
//                            //std::cout << (*pPreviousFirst)->getKey() << " pPreviousLast SEQUENCE :" << (*pPreviousLast)->getKey() << std::endl;
//                            //std::cout << (*pFirst)->getKey() << "         pLast SEQUENCE :" << (*pLast)->getKey() << std::endl;
//                            //std::cout << (*pNextFirst)->getKey() << "     pNextLast SEQUENCE :" << (*pNextLast)->getKey()<< std::endl;
//                            //std::cout << " LAST SEQUENCE BASE :" << (*pNextLast.base())->getKey()<< std::endl;
//                            // have two cases to distinguish
//                            // case 1:   Y ... Y-X ... X-Y ... Y
//                            // case 2:   Y ... X-Y ... Y-X ... Y
//                            if (((*pPreviousFirst)->getKey() == pNeighbor->getKey()) && ((*pNextLast)->getKey() == pNeighbor->getKey())) { // case 1:   Y ... Y-X ... X-Y ... Y
//                                //             Last -> NodeToRemove
//                                //        pNextLast -> Node +1 to Last
//                                // pNextLast.base() -> Node +2 to Last
//                                // pLast is a reverse iterator, decrementing moves the location up, so pod points to pNextLast
//                                auto pod = pPreviousLast.base(); // will potential go past the list, equivalent to list.end()
//
//                                if (pFirst == pPreviousLast.base()){ // do they point to same place?
//                                    // XYX = > X where Y is neighbor of X, meaning do they point to Y
//                                    pTour->erase(pPreviousFirst, pLast.base());
//                                } else {
//                                    remainderNodes.clear();
//                                    // populate remainderNodes map
//                                    for(std::list< Node * >::iterator it = pTour->begin(); it != pFirst; ++it) {
//                                        remainderNodes.insert((*it)->getKey());
//                                    }
//
//                                    for(std::list< Node * >::reverse_iterator it = pTour->rbegin(); it != pLast; ++it) {
//                                        remainderNodes.insert((*it)->getKey());
//                                    }
//
//                                    bool found = false;
//                                    std::unordered_set<int>::iterator got;
//
//                                    auto pSearch = pFirst; // make copy of pFirst, pFirst should be pointing to Node_to_Remove
//                                    // pod points to
//                                    for(; pSearch != pod; ++pSearch) {
//                                        // search to see if a bounded node can be found in remainderNodes
//                                        // if true, break out and substitute
//                                        if (!((*pSearch)->getAccessed())){
//                                            got = remainderNodes.find ( (*pSearch)->getKey() );
//                                            (*pSearch)->setAccessed(true);
//                                            checkedNodes.push_back((*pSearch)->getKey());
//
//                                            if (got != remainderNodes.end()){
//                                                found = true;
//                                                break;
//                                            }
//                                        }
//                                    }
//
//                                    std::list< Node * > tempList;
//                                    // find first occurrence of newRoot
//                                    tempList.splice(tempList.begin(), *pTour, pFirst, --pNextLast.base());
//                                    pTour->erase(--pNextLast.base()); // removes duplicate nodes next two each other
//
//                                    if (found){ // make substitution
//                                        // reroot subtour bounded by neighbors and substitute
//                                        rerootSubTour(&nodes[*got], &tempList);
//                                        tempList.pop_back();
//                                        std::list< Node * >::iterator gotIt = std::find(pTour->begin(), pTour->end(), &nodes[*got]);
//                                        pTour->splice(gotIt, tempList);
//
//                                    } else { // create new list, split out and reset root nodes of items in new list
//                                        // reassign the root node in Nodes of tempList before adding to tours
//                                        resetRootNodesInSubTour(&tempList);
//                                        // add tempList to tours
//                                        tours.emplace ( std::pair< int, std::list<Node *> >((*tempList.begin())->getKey(), tempList) );  // key is the root of the tour
//                                    }
//
//                                    // clear tempList
//                                    for (auto it = tempList.begin() ; it != tempList.end(); ++it) {
//                                        it = tempList.erase(it);
//                                    }
//                                    tempList.clear();
//                                }
//
//                            } else if (((*pNextFirst)->getKey() == pNeighbor->getKey()) && ((*pPreviousLast)->getKey() == pNeighbor->getKey())) { // case 2:   Y ... X-Y ... Y-X ... Y
//                                // does remainderNodes need to be cleared for each neighbor
//                                // populate remainderNodes map
//                                // pFirst and pLast points to NodeToRemove
//
//                                remainderNodes.clear();
//                                for(std::list< Node * >::iterator it = pTour->begin(); it != pFirst; ++it) {
//                                    remainderNodes.insert((*it)->getKey());
//                                }
//
//                                for(std::list< Node * >::reverse_iterator it = pTour->rbegin(); it != pLast; ++it) {
//                                    remainderNodes.insert((*it)->getKey());
//                                }
//
//                                bool found = false;
//                                std::unordered_set<int>::iterator got;
//
//                                auto pSearch = pFirst; // make copy of pFirst, pFirst should be pointing to NeighborNode
//                                auto pod = pPreviousLast.base();
//                                // pod points to
//                                ++pSearch;
//                                for(; pSearch != pod; ++pSearch) {
//                                    // search to see if a bounded node can be found in remainderNodes
//                                    // if true, break out and substitute
//                                    if (!((*pSearch)->getAccessed()) ){
//                                        got = remainderNodes.find ( (*pSearch)->getKey() );
//                                        (*pSearch)->setAccessed(true);
//                                        checkedNodes.push_back((*pSearch)->getKey());
//
//                                        if (got != remainderNodes.end()){
//                                            found = true;
//                                            break;
//                                        }
//                                    }
//                                }
//
//                                std::list< Node * > tempList;
//                                // find first occurrence of newRoot
//                                tempList.splice(tempList.begin(), *pTour, pNextFirst, pPreviousLast.base()); // exclusive of pPreviousLast.base()
//                                pTour->erase(--pLast.base()); // removes duplicate nodes next two each other
//
//                                if (found){ // make substitution
//                                    // reroot subtour bounded by neighbors and substitute
//                                    if (tempList.size() > 1){
//                                        rerootSubTour(&nodes[*got], &tempList);
//                                        tempList.pop_back();
//                                        std::list< Node * >::iterator gotIt = std::find(pTour->begin(), pTour->end(), &nodes[*got]);
//                                        pTour->splice(gotIt, tempList);
//                                    }
//                                } else { // break out and return false
//                                    // pFirst points to node to remove
//                                    // pLast points to node to remove
//                                    resetRootNodesInSubTour(&tempList);
//                                    tours.emplace ( std::pair< int, std::list<Node *> >((*tempList.begin())->getKey(), tempList) );  // key is the root of the tour
//                                }
//
////                                std::cout << "tempList " << tempList.size() << std::endl;
//                                // clear tempList
//                                for (auto it = tempList.begin() ; it != tempList.end(); ++it) {
//                                    it = tempList.erase(it);
//                                }
//                                tempList.clear();
//
//                            } // end of if statement
//
//                            pNodeToRemove->removeNeighbor(pNeighbor);  // nodeToRemove exists until we have removed all neighbors
//                            resetAccessed(&checkedNodes);
//                            pLast = --(pTour->rend());
//                            terminateNow = true;
//                            break;
//                        } // end if statement, if above statement excuted, it means we made the substitution and move to next neighbor
//                        ++pLast;
//                    } // while NodeToRemove loop LAST
//                } // end of if statement that checks if neighbor node is next to node_to_remove
//                if (terminateNow){
//                    break;
//                }
//                pFirst++;
//            } // while NodeToRemove loop FIRST
//            //result = "remove node " + std::to_string(pNodeToRemove->getKey()) + " " + std::to_string(pNodeToRemove->getTotalNeighbors());
//        } // while loop for each neighbors
//    } else { // no neighbors
//        //
//        // tours is a map => std::map<int, std::list< Node *> > tours; // key is the root of the tour
//        // list must be cleared first
//        pTour = &tours[pNodeToRemove->getRootNodeOfTour()];
//
//        // std::map<int, std::list< Node *> > tours; // key is the root of the tour
////        for(std::list<Node *>::iterator lit = pTour->begin(); lit != pTour->end(); ++lit){
////            //it->second.erase(lit);
////            pTour->erase(lit);
////        }
//
//        for (auto it = pTour->begin() ; it != pTour->end(); ++it) {
//            it = pTour->erase(it); // erase invalidates existing iterators
//        }
//
//        pTour->clear();
//        tours.erase(pNodeToRemove->getRootNodeOfTour());
//    }

    // remove the node from nodes
//    if (pNodeToRemove->getTotalNeighbors() == 0){
//        //tours.erase(pNodeToRemove->getRootNodeOfTour());
//        if (nodes.find(pNodeToRemove->getKey()) == nodes.end()){
//            std::cout << "NOT FOUND in NODES LIST!  CANT REMOVE " << pNodeToRemove->getKey() << std::endl;
//        }
//        // nodes is a map => std::map<int, Node> nodes;
//        nodes.erase(pNodeToRemove->getKey()); // this should destroy the object
//    }

    nodes.erase(pNodeToRemove->getKey()); // this should destroy the object

    totalComponents = tours.size();
    return totalComponents;
}




/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 */
bool EulerTourSmart::addToTour(int nodeToAdd){

    SmartNode * pNode = nodes[nodeToAdd].get();
    std::shared_ptr<SmartNode> pNeighbor;

    // add to Tour
    //if (validateList("********* BEFORE ADDING NEW TOUR MEMBERS")){
    //    return false;
    //}

    if (pNode->getTotalNeighbors() == 0){
        // add node to tours as new tour
        // set node tour index mapping
        // in old Euler Tour, tours was a map whose key is the root node
        // std::map<int, std::list< Node *> > tours; // key is the root of the tour
        //tours.insert ( std::pair< int, std::list<Node *> >(pNode->getKey(), std::list< Node *>()) );  // key is the root of the tour
        tours.insert(std::make_shared<Tour>( nodes[nodeToAdd].get() ));
        // add node to tour
        //tours[pNode->getKey()].push_back(&nodes[nodeToAdd]);
        std::set< std::shared_ptr<Tour> >::iterator inTour = std::find_if(tours.begin(), tours.end(), by_root_node_of_tour(pNode->getKey()));
        pNode->setPointerToTour(*inTour); // sets pointer to list< Node *> in vector
        //pNode->setRootNodeOfTour(pNode->getKey()); // if no neighbors, root node is the key of only node

        // pNode has_many relation with Tours
    } else {
        // create subtour of new Node
        // then check if new node is found an previous tour and merge
        Tour subTour(nodes[nodeToAdd]); // make a Tour with Node using neighbors to make the tour

        int totalNeighbors = pNode->getTotalNeighbors();
        int keyOfNeighborNode = pNode->getKey();
        // new Node could bridge 1 or more prior nodes
        std::set< int > numberOfNeighborhoods;

        for (int j=0; j < totalNeighbors; j++){
            pNeighbor = pNode->getPointerToNeighborByIndex(j).lock().get();
            // get tour index to merge into
            numberOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour()); // each neighbor could already belong to another tour
            // get pointer to list of neighbor, is it the same?
            keyOfNeighborNode = pNeighbor->getKey();
        }

        // set tour membership of new Node, we use the last neighbor as the key
        pNeighbor = nodes[keyOfNeighborNode];
//        pNode->setPointerToTour(pNeighbor->getPointerToTour()); // sets pointer to list< Node *> in vector
        const int rootToBaseTour = pNeighbor->getRootNodeOfTour();
        // find Tour that has keyOfNeighborNode
        std::set< std::shared_ptr<Tour> >::iterator pExistingNeighborTour = std::find_if (tours.begin(), tours.end(), by_root_node_of_tour(rootToBaseTour));
        // merge subTour with tour that contains neighbor
        (*pExistingNeighborTour).get()->mergeTours(pNeighbor, subTour);
        // subTour should be empty now

        // any additional tours will be greater than minTourIndex
        if (numberOfNeighborhoods.size() > 1){ // if number of neighborhoods is greater than 1, it means the new node bridges at least two
            // find tour greater than minTourIndex, reroot tour and merge
            for (int j=0; j < (totalNeighbors-1); j++){

                //pNeighbor = pNode->getPointerToNeighborByIndex(j);
                pNeighbor = pNode->getPointerToNeighborByIndex(j).lock().get();
                int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();

                // if they have the same root, already in BaseTour
                // if different, new Node bridges two tours
                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again

                    std::set< std::shared_ptr<Tour> >::iterator pCurrentTour = std::find_if (tours.begin(), tours.end(), by_root_node_of_tour(rootOfCurrentTour));

                    (*pExistingNeighborTour).get()->mergeTours(pNeighbor, *(*pCurrentTour).get());

                    std::cout << " TOUR " << (*pCurrentTour).get()->getSizeOfTour() << std::endl;
                    tours.erase(pCurrentTour);

                    //std::cout << " TOURS SIZE addToTour : " << tours.size() << std::endl;
                    //tours[rootToBaseTour] =  *pExistingNeighborTour;
                }
            }
        } // pointer deallocation problem exiting block


    } // end of creating subtour

    return true;
}




