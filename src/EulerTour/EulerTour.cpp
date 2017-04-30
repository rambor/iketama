//
// Created by Robert Rambo on 11/01/2017.
//
#include "EulerTour.h"
#include "../Model.h"

// Constructor
EulerTour::EulerTour(){
}

EulerTour::EulerTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel){
    this->createInitialTour(workingLimit, pModel, beginIt);
}


/**
 * add a node to existing set of nodes
 * first, create neighborhood from existing nodes
 * find which tour to tours to add to and merge any tours that are bridged by new node
 * returns number of nodes in tour
 *
 * Possible undefined condition, using find on nodes can lead to nodes.end()
 *
 * the control of active nodes is maintained by beads_in_use and bead_indices
 */
int EulerTour::addNode(int newNode, Model *pModel) {

    //std::vector<int>::iterator itIndex;
    int neighbor;

    nodes.insert ( std::pair<int,Node>(newNode, Node(newNode) ) ); // map -> node_index -> pointer
    Node * pNode = &(nodes.find(newNode)->second); // retrieve the newly made node

   // Node * pNode = &(nodes[newNode]); // was causing issues, an empty constructor would be made
    //
    // can lead to nodes.end() which doesn't contain a node
    // somehow getting a nonsense node to return
    // build neighborhood (check from existing nodes or beads in use)
    auto it = pModel->getPointerToNeighborhood(newNode);
    for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
        // if neighbor is inside workinglimit, don't add
        neighbor = *(it+j); // index of bead used as key
        // if not found, itIndex will report last
        if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){ // look in nodes list
            // nodes.find(0) returns true but nodeKey is something not possible
            // how does this return true?
            pNode->addNeighbor(&(nodes.find(neighbor)->second)); // add neighbor to both parent and child
            // again nothing prevents me from finding nodes.end()
            //pNode->addNeighbor(&(nodes[neighbor])); // add neighbor to both parent and child
            // has no idea which list it is in
        } else if (neighbor == -1){
            break;
        }
    }

    // add to Tour
    addToTour(newNode);
    totalComponents = tours.size();
//    validateNodes("AFTER ADDING NODE : " + newNode);
    return totalComponents;
}

// reuse Tour and reset
int EulerTour::newTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel) {

    // clear pointers to nodes
    for (std::map<int, std::list< Node *>>::iterator it=tours.begin(); it!=tours.end(); ++it){
        // delete points in tour
       for(std::list<Node *>::iterator lit = it->second.begin(); lit != it->second.end(); ++lit){
           it->second.erase(lit);
       }
        it->second.clear();
    }

    // remove all nodes
    std::map<int,Node>::iterator it = nodes.begin();
    while(it != nodes.end()){
        it = nodes.erase(it);
    }

    tours.clear();
    nodes.clear();

    //pSelectedLattice = &beginIt;
    this->createInitialTour(workingLimit, pModel, beginIt);
    std::cout << "Tour size after initialization of new tour: " << tours.size() << std::endl;
    return totalComponents;
}



//void EulerTour::createBackUp(){
//
////    std::map<int, Node> backUpedNodes; // backup of Nodes
//    std::clock_t start = std::clock();
//
//    //cout << "REMOVE old : " << oldmeth/(double) CLOCKS_PER_SEC << " new : " << newmeth/(double) CLOCKS_PER_SEC << endl;
//    backedUpNodes.clear();
//    for(std::map<int, Node>::iterator it = nodes.begin(); it != nodes.end(); it++){
//        backedUpNodes.insert( std::pair<int,Node>(it->first, it->second )) ;
//    }
//
//    std::cout << "COPY MAP TIME : " << (std::clock() - start)/(double) CLOCKS_PER_SEC << std::endl;
//
////    std::list< Node * > backUpedTours; //
//
//}

/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 */
bool EulerTour::addToTour(int nodeToAdd){

    //Node *pNeighbor, * pNode = &(nodes[nodeToAdd]);
    Node *pNeighbor, * pNode = &(nodes.find(nodeToAdd)->second);
    // add to Tour
    //if (validateList("********* BEFORE ADDING NEW TOUR MEMBERS")){
    //    return false;
    //}

    if (pNode->getTotalNeighbors() == 0){
        // add node to tours as new tour
        // set node tour index mapping
        // create tour (empty)
        tours.insert ( std::pair< int, std::list<Node *> >(pNode->getKey(), std::list< Node *>()) );  // key is the root of the tour
        // add node
        tours[pNode->getKey()].push_back(pNode); // push pointer of node into tour
        //tours[pNode->getKey()].push_back(&nodes[nodeToAdd]); // push pointer of node into tour
        pNode->setPointerToTour(&tours[pNode->getKey()]); // sets pointer to list< Node *> in vector
    } else {
        // create subtour of new Node
        // then check if new node is found an previous tour and merge
        // std::list < std::unique_ptr<Node> > myList;
        std::list< Node * > subTour;
        createSubTour(pNode, &subTour); // create subtour using pNode's neighbors

        //printList("/nsubtour of : " + std::to_string(pNode->getKey()), &subTour);
        int totalNeighbors = pNode->getTotalNeighbors(); // neighborhood was created when adding Node, checks existing nodes for proximity
        int keyOfNeighborNode = pNode->getKey();
        // new Node could bridge 1 or more prior nodes
        std::set< int > rootOfNeighborhoods;

        for (int j=0; j < totalNeighbors; j++){
            pNeighbor = pNode->getPointerToNeighborByIndex(j);
            // get tour index to merge into
            rootOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour()); // how many neighborhoods does newNode belong to?
            // get pointer to list of neighbor, is it the same?
            keyOfNeighborNode = pNeighbor->getKey();
        }

        // what happens if node that is removed is root node for a tour?
        // set tour membership of new Node
        pNeighbor = &(nodes.find(keyOfNeighborNode)->second);
        //pNeighbor = &nodes[keyOfNeighborNode];

        rerootSubTour(pNeighbor, &subTour); //prepare for merging in pNeighbor tour (doesn't change tour membership of nodes, still points to old location)
        pNode->setPointerToTour(pNeighbor->getPointerToTour()); // sets pointer to list< Node *> in vector

        const int rootToBaseTour = pNode->getRootNodeOfTour();
        // re-root to node found incommon to a previous tour in tours
        // merge in tours[i]
        // get base tour of neighbor to merge into
        //std::list< Node * > * pExistingNeighborTour = &tours[rootToBaseTour];
        std::list< Node * > * pExistingNeighborTour = pNeighbor->getPointerToTour();
        // subTour is rooted to pNode
        //               subTour => 58565
        //
        // pExistingNeighborTour => 124313426797621   (root node is 1)
        //
        // reroot        subTour => 65856
        //
        std::list< Node * >::iterator inTour = std::find (pExistingNeighborTour->begin(), pExistingNeighborTour->end(), pNeighbor);
        //
        // locate 6 at ...267...
        //
        subTour.pop_back();
        //
        // pop_back      subTour => 6585
        // merge tour
        // merge subTour of newNode with existing tour
        pExistingNeighborTour->splice(inTour, subTour); // Add additional tours to base
        // subTour should be empty now
        //
        // pExistingNeighborTour => 1243134265856797621   (root node is 1)
        //
        // nodes in this tour will have mixed roots with pNeighbor
        // if new node bridges two or more tours, pExistingNeighborTour will be mixed
        if (rootOfNeighborhoods.size() > 1){ // if number of neighborhoods is greater than 1, it means the new node bridges at least two
            // find tour greater than minTourIndex, reroot tour and merge
            for (int j=0; j < (totalNeighbors-1); j++){
                // if two neighbors are in same neighborhood, and I add the first one
                // then i have a problem if rootNodeOfTour is not updated before next neighbor is checked
                pNeighbor = pNode->getPointerToNeighborByIndex(j);
                int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();

                // if two nodes belong to same tour, only need to add tour once
                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again
                    std::list< Node * > * pList = &(tours.find(rootOfCurrentTour)->second);

                    // pList should not be zero size
                    //std::cout << j << " rootOfCurrentTour " << rootOfCurrentTour << " pList SIZE => " << pList->size() << std::endl;
                    rerootSubTour(pNeighbor, pList);

                    inTour = std::find (pExistingNeighborTour->begin(), pExistingNeighborTour->end(), pNeighbor);

                    // better to iterate over node list and check for tour membership
                    // for all nodes in this tour, reset tourIndex
                    //std::cout << "\n  ========     POINTERS " << pExistingNeighborTour << " " << &tours[rootToBaseTour] << " ++++\n"<< std::endl;

                    for(std::list< Node * >::iterator pIt = pList->begin(); pIt != pList->end(); ++pIt){
                        if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
                            (*pIt)->setPointerToTour(pExistingNeighborTour);
                        }
                    }
                    // reassign before popping, in case pList size is 1
                    pList->pop_back();
                    // merge pList into existingNeighborTour
                    pExistingNeighborTour->splice(inTour, *pList);
                    // remove, erasing from vector changes address of the other elements,
                    // pList must be zero in size
                    tours.erase(rootOfCurrentTour); // calls destructor on list held by tours map
                }
            }
        }

        // after all neighbors have been added, we then update the nodes tour
        for(std::list< Node * >::iterator pIt = pExistingNeighborTour->begin(); pIt != pExistingNeighborTour->end(); ++pIt){
            if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
                (*pIt)->setPointerToTour(pExistingNeighborTour);
            }
        }

    } // end of creating subtour

    return true;
}




/**
 * Print a tour, diagnostic tool
 */
void EulerTour::printList(std::string text, std::list< Node * > * list){
    int tourLength = list->size();

    int count=1;
    std::cout << text << std::endl;
    std::cout << "      LIST START Length : " << tourLength << std::endl;
    for(std::list< Node * >::iterator iterator = list->begin(); iterator != list->end(); iterator++) {
        std::cout << "                   LIST : (" << count << ") " << (*iterator)->getKey() << std::endl;
        count++;
    }
    std::cout << " _________________________ " << std::endl;
}


/**
 * passing a pointer to a list of Nodes
 *
 */
void EulerTour::rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    if (subTourToLoad->front()->getKey() != newRoot->getKey()){
        // printList("================> REROOT_AND_REASSIGN_SUBTOUR REROOT Before", subTourToLoad);
        //
        // subTourToLoad => 124313426797621   (root node is 1)
        //       newRoot => 6
        std::list< Node * >::iterator inTour = std::find (subTourToLoad->begin(), subTourToLoad->end(), newRoot);
        subTourToLoad->pop_front();
        //
        //        inTour => 6 at ...4267...
        // subTourToLoad => 24313426797621
        //
        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour); // exclude inTour position, transfer is upto
        //      tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(newRoot);
        //      tempList => 24313426
        // copy to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
        // subTourToLoad => 679762124313426
        resetRootNodesInSubTourOfTour(subTourToLoad);
        // clear tempList
        tempList.clear();
    }
}

/**
 *               9
 *               |
 *   1 - 2 - 6 - 7
 *   |   |
 *   3 - 4
 *   124313426797621
 *  smallest tour size is 1, then 3, there is no two
 *  Only reroot the tour, there is no updating of the nodes root tour
 */
void EulerTour::rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    // 124313426797621 => reroot to 6
    std::list< Node * >::iterator inTour = std::find (subTourToLoad->begin(), subTourToLoad->end(), newRoot);
    // inTour points to 6 at ..4267...
    std::list< Node * > tempList;
    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    if (subTourToLoad->front()->getKey() != newRoot->getKey() && (subTourToLoad->size() > 2)){
        subTourToLoad->pop_front();
        // subTourToLoad => 24313426797621
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour);
        // tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(newRoot);
        //      tempList => 24313426
        // transfer components from tempList to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
        // subTourToLoad => 6797621-24313426
    }

    // clear tempList
    tempList.clear();
}

int EulerTour::removeNode(int indexOfNode){

    //Node * pNodeToRemove = &(nodes[indexOfNode]);
    Node * pNodeToRemove = &(nodes.find(indexOfNode)->second);


//if (nodes.find(indexOfNode) == nodes.end()){
//    std::cout << "LOOKING FOR A NODE NOT FOUND " << indexOfNode << std::endl;
//    exit(0);
//}

    // check neighbors
//    validateNodes();
    //validateNodesAndTours("FROM REMOVE ");

    std::list< Node * > * pTour;// = &tours[pNodeToRemove->getRootNodeOfTour()];  // tour with node to remove

    if (pNodeToRemove->getTotalNeighbors() > 0){

        if (pNodeToRemove->getTotalNeighbors() == 1){ // root to neighbor

            pTour = &(tours.find(pNodeToRemove->getRootNodeOfTour())->second);

            //pTour = &tours[pNodeToRemove->getRootNodeOfTour()];    // get tour that contains node to remove
            //std::cout << "(REMOVE SINGLE NODE) KEY " << pNodeToRemove->getKey() << "( " << pNodeToRemove->getRootNodeOfTour() << " ) " << " TOTAL NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << " SIZE " << pTour->size() << std::endl;
            // pointer to this tour
            // reroot it
            // tour still points to it after popping ends
            rerootSubTour(pNodeToRemove, pTour);

            pTour->pop_front();
            pTour->pop_back();

            pTour->front()->removeNeighbor(pNodeToRemove);

            resetRootNodesInSubTourOfTour(pTour); // reassigns tour and should clear old map entry
        } else {

            pTour = &tours[pNodeToRemove->getRootNodeOfTour()];    // get tour that contains node to remove
            // reroot tour to NodeToRemove and remove old tour
//            std::cout << "(REMOVE NODE) KEY " << pNodeToRemove->getKey() << "( " << pNodeToRemove->getRootNodeOfTour() << " ) " << " TOTAL NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << " SIZE " << pTour->size() << std::endl;

            rerootAndReassignSubTour(pNodeToRemove, pTour);  // remove old tour and set as newly rooted tour
            pTour = &tours[pNodeToRemove->getRootNodeOfTour()];    // get tour that contains node to remove

            std::list< Node * >::iterator  pNext;
            Node * pNeighbor;
            int keyOfNodeToRemove = pNodeToRemove->getKey();

            while (pNodeToRemove->getTotalNeighbors() > 0) {

                //rerootSubTour(pNodeToRemove, pTour);
                // printList("NODETOREMOVE ", pTour);
                //std::cout << " => " << pNodeToRemove->getKey() << std::endl;
                // who is the neighbor
                std::list< Node * >::iterator pRightNeighbor = std::next(pTour->begin(), 1);
                std::list< Node * >::iterator pStartHere = std::next(pTour->begin(), 2);

                pNeighbor = *pRightNeighbor; // get pointer to Neighbor

                // find next neighbor pair
                if ((*pStartHere)->getKey() == keyOfNodeToRemove){ // XYX where X is node to remove and Y is neighbor
                    // check if Y is in pTour, if it is delete first two nodes from Tour
                    pNext = pStartHere;
                    // if not, single node becomes a new tour
                } else {

                    std::list< Node * >::iterator pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor); // look for occurrences of pNeighbor
                    while(pLeftNeighbor != pTour->end()){ // find neighbor and check to see if next to it is NodeToRemove

                        pNext = std::next(pLeftNeighbor, 1); // could be .end() ?
                        if ((*pNext)->getKey() == pNodeToRemove->getKey()) {
                            break;
                        }
                        // find next
                        pStartHere = std::next(pLeftNeighbor, 1);
                        pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor);
                    }

//                    if (std::distance(pLeftNeighbor, pNext) == 0){
//                        printList("pLeft - pNext is zero", pTour);
//                        exit(0);
//                    }
                }

//                if (std::distance(pTour->begin(), pNext) > pTour->size()){
//                    printList("TOO FAR pNext", pTour);
//                        exit(0);
//                }

                // potential pLeftNeighbor will reach end?
                pTour->pop_front();

                std::list< Node * > subTour;
                subTour.splice(subTour.begin(), *pTour, pTour->begin(), pNext); // upTo but not including pNext, pNext should point to node to remove
//                if (subTour.size() == 0){
//                    printList("SUBTOUR IS ZERO ", pTour);
//                    std::cout << pNodeToRemove->getKey() << " pNext distance " << std::distance(pTour->begin(), pNext) << std::endl;
//                    std::cout << pNodeToRemove->getKey() << "      pNeighbor " << pNeighbor->getKey() << std::endl;
//                }
                //
                // if pNext is last node, then pTour will have only one element
                //
                // now I have two tours
                // check if I can sub back into remaining pTour
                std::set<int> nodesToCheck; // how could subTour be empty?
                for(std::list< Node * >::iterator sit = subTour.begin(); sit != subTour.end(); ++sit){
                    nodesToCheck.insert((*sit)->getKey());
                }

                // merge if possible, check each node in nodesToCheck
                std::list< Node * >::iterator locationOfCommonNode;
                for(std::set< int >::iterator sit = nodesToCheck.begin(); sit != nodesToCheck.end(); ++sit){
                    locationOfCommonNode = std::find_if(pTour->begin(), pTour->end(), find_by_key(*sit));
                    if (locationOfCommonNode != pTour->end() ){
                        // merge subTour with pTour
                        rerootSubTour(*locationOfCommonNode, &subTour);
                        // root pTour to common node
                        subTour.pop_back();
                        pTour->splice(locationOfCommonNode, subTour);
                        break;
                    }
                }

                // make a new Tour from subTour if it couldn't be merged(linked)
                if (subTour.size() > 0){
                    int newRootNodeForTour = subTour.front()->getKey();
                    //tours.insert ( std::pair< int, std::list<Node *> >(newRootNodeForTour, subTour) );  // key is the root of the tour
                    tours[newRootNodeForTour] =subTour;  // key is the root of the tour

                    auto pNewRootNode = &tours[newRootNodeForTour];
                    // should only iterate of the set of unique nodes in the tour, not all of them
                    // make a tour object that includes the list and set of nodes in tour?
                    for(std::list< Node * >::iterator it = pNewRootNode->begin(); it != pNewRootNode->end(); ++it) {
                        if ((*it)->getPointerToTour() != pNewRootNode){
                            (*it)->setPointerToTour(pNewRootNode);
                        }
                    }
                }
//                if (pTour->size() <= 1){ // move subTour into pTour, neighbor count should be 1
//                    std::cout << " SIZE of pTour => " << pTour->size() << " NODETOREMOVE => " << pNodeToRemove->getKey() << " NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << std::endl;
//                    printList("pTour size 1 : ROOT " + std::to_string(pNodeToRemove->getRootNodeOfTour()), pTour);
//                }


                // pTour should be rooted at NodeToRemove at this point
                // should stay that way until pNodeToRemove is gone
                pNodeToRemove->removeNeighbor(pNeighbor); // should be in balance after removing
//                validateNodes("INSIDE REMOVE NODES");
                // if no more pNodeToRemove in pTour, can't reroot to node to remove
            }

//            std::cout << "     NODE TO REMOVE : " << pNodeToRemove->getKey() << " " << pNodeToRemove->getTotalNeighbors() << std::endl;

            // remove pTour (should contain only a single node and be rooted to nodeToRemove)
//            if (pTour->size() > 1){
//                printList(" pTOUR size Greater than 1", pTour);
//                exit(0);
//            }

            resetRootNodesInSubTourOfTour(pTour); // pTour should always contain the node to remove as a single tour
            tours[pNodeToRemove->getRootNodeOfTour()].clear(); // clear List
            tours.erase(pNodeToRemove->getRootNodeOfTour()); // remove from tour
        }

    } else { // no neighbors
        //
        // tours is a map => std::map<int, std::list< Node *> > tours; // key is the root of the tour
        // list must be cleared first
        pTour = &tours[pNodeToRemove->getRootNodeOfTour()];
        pTour->erase(pTour->begin(), pTour->end()); // erase invalidates existing iterators
        pTour->clear();
        tours.erase(pNodeToRemove->getRootNodeOfTour());
    }

    // remove node
    // remove the node from nodes
    if (pNodeToRemove->getTotalNeighbors() == 0){
//        std::cout << " REMOVING NODE : " << pNodeToRemove->getKey() << " ROOTNODE :" << pNodeToRemove->getRootNodeOfTour() << std::endl;
//        if (nodes.find(pNodeToRemove->getKey()) == nodes.end()){
//            std::cout << "ERROR NOT FOUND in NODES LIST!  CANT REMOVE " << pNodeToRemove->getKey() << std::endl;
//            exit(0);
//        }
        // nodes is a map => std::map<int, Node> nodes;
        nodes.erase(pNodeToRemove->getKey());

//        validateNodesAndTours("END OF REMOVE NODE");

    } else {
        std::cout << " NODE NEIGHBOR LIST NOT EMPTY PROBLEM WITH CODE?" << std::endl;
        exit(0);
    }
//    validateNodes("AFTER REMOVE NODES");
    totalComponents = tours.size();
    return totalComponents;
}



bool EulerTour::validateTour(std::list<Node *> * tourtocheck){
    // each pair of numbers needs a reverse sequence
    // go through each element of list and validate it is a node
    //std::map<int, int> pairs;
    std::list<Node *>::iterator stopIt = std::next(tourtocheck->end(),-1);
    for(std::list<Node *>::iterator nodeIt=tourtocheck->begin(); nodeIt != stopIt; nodeIt++){

        std::list<Node *>::iterator next = std::next(nodeIt, 1);

        int firstKey = (*nodeIt)->getKey();
        int secondKey = (*next)->getKey();

        bool isTrue = true;

        for(std::list<Node *>::iterator nodeIt2=tourtocheck->begin(); nodeIt2 != stopIt; nodeIt2++){

            std::list<Node *>::iterator next2 = std::next(nodeIt2, 1);
            int firstKey2  = (*nodeIt2)->getKey();
            int secondKey2 = (*next2)->getKey();
            if (firstKey == secondKey2 && secondKey == firstKey2){
                isTrue = false;
                break;
            }
        }

        if (isTrue){
            std::cout << "InVALID TOur " << firstKey << " " << secondKey << std::endl;
            return false;
        }
    }

    return true;
}

bool EulerTour::validateNodes(std::string str){

    for (it_type it = nodes.begin(); it != nodes.end(); it++){
        // check for repeats
        Node tempNode = it->second;
        if (!(it->second).validate()){
            std::cout << str << std::endl;
            std::cout << " INVALID NODE " << it->first << " " << tempNode.getKey() << std::endl;
            exit(0);
            return false;
        }

        // for each node, check that his neighbors are present in node list
        int totalNeighbors = tempNode.getTotalNeighbors();
        std::set<int> * pSet = tempNode.getIteratorToIndices();
        int count = 0;

        for(std::set<int>::iterator it = pSet->begin(); it != pSet->end(); ++it){ // go through neighbors of selected node and check for reciprocity
            if (nodes.find(*it) == nodes.end()){
                std::cout << str << std::endl;
                std::cout << count << " ( " << totalNeighbors << " ) " << " NEIGHBOR NODE NOT FOUND IN NODES LIST " << tempNode.getKey() << " neigh => " << *it<< std::endl;
                std::cout << count << " BASE NODE NEIGHBORS : " << std::endl;
                tempNode.printNeighbors();
                exit(0);
                return false;
            }

            // check if neighbor has current node as neighbor
            if ( !(nodes.find(*it)->second.isNeighborPresent(tempNode.getKey())) ){
                std::cout << str << std::endl;
                std::cout << count << " NODE - NEIGHBOR relationship invalid " << std::endl;
                std::cout << " BASE NODE : " << tempNode.getKey() << " has neighbor <=> " << nodes.find(*it)->second.getKey() << " lack of reciprocity" << std::endl;
                // print neighborhood of node
                std::cout << count << " BASE NODE NEIGHBORS : " << std::endl;
                tempNode.printNeighbors();
                std::cout << count << "  NEIGHBOR NEIGHBORS : " << std::endl;
                nodes.find(*it)->second.printNeighbors();
                exit(0);
                return false;
            }
            count++;
        }
    }
    return true;
}

bool EulerTour::validateNodesAndTours(std::string text){


    int count=1;
    for(std::map<int, Node>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        // validate tour
        //std::cout << count  << " EXAMINING NODE " << (*it).second.getKey() << std::endl;
        if ((*it).first != (*it).second.getKey()){
            std::cout << " NODE keys do not match : " << text << std::endl;
            return true;
        }

        if ((*it).second.getPointerToTour()->size() == 0){
            std::cout << (*it).first << " ZERO SIZE TOUR " << text << std::endl;
            return true;
        }

       ++count;
    } // use shared pointer? make instance on heap

    return false;
}


bool EulerTour::validateList(std::string text){

    std::cout << " VALIDATING LISTS IN TOURS :::::: " << text << std::endl;
    int count = 1;
    for(std::map<int, std::list< Node *> >::iterator it = tours.begin(); it != tours.end(); it++){

        int tourIndex = it->first; // tour index should be an existing node
        std::list<Node *> nodelist = it->second;

        // check if the root of tour is in the nodes list
        if (nodes.find(tourIndex) == nodes.end()){
            std::cout << " @ " << text << " Node Does Not exist => " << tourIndex << " tour size : " << nodelist.size()<< std::endl;
            return true;
        }

        // check if nodelist size violates limit
        if (nodelist.size() == 2){
            std::cout << " @ " << text << " LIST SIZE == 2 ************************************************* tourindex => " << tourIndex << std::endl;
            return true;
        }

        // check if nodelist contains a legitimate tour (every edge contains reciprocal pair)
        if (nodelist.size() > 2 && !validateTour(&nodelist)){
            std::cout << " @ " << text << " INVALID LIST  " << tourIndex << std::endl;
            printList("INVALID LIST or TOUR ", &nodelist);
            return true;
        }

        // checks that first and last nodes of tour are the same
        if (*(nodelist.begin()) != nodelist.back()){
            printList(" DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR ", &nodelist);
            std::cout << " DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR " << std::endl;
            return false;
        }


        // go through each element of list and validate it is a node
        for(std::list<Node *>::iterator nodeIt=nodelist.begin(); nodeIt != nodelist.end(); nodeIt++){
            int key = (*nodeIt)->getKey();
            if (nodes.find(key) == nodes.end()){
                std::cout << count << " LOCALE : " << text << " ======> Analyzing tour : " << tourIndex << " of " << tours.size() << std::endl;
                std::cout << count << " LOCALE : " << text << " ======> NUMBER OF ELES : " << nodelist.size() << std::endl;
                std::cout << count << " LOCALE : " << text << " ======> NUMBER OF KEY  : " << key << std::endl;
                printList("MESSED UP TOUR ", &nodelist);
                return true;
            }
        }
        count++;
    }

    return false;
}

///**
// * subTour is a pointer to tour
// * newRootNode
// */
//void EulerTour::resetRootNodesInSubTour(std::list<Node *> * subTour){
//
//    // get first node in tour
//    int newRootNodeForTour = (*(subTour->begin()))->getKey();
//
//    for(std::list< Node * >::iterator it = subTour->begin(); it != subTour->end(); ++it) {
//        if ((*it)->getRootNodeOfTour() != newRootNodeForTour){
//            (*it)->setRootNodeOfTour(newRootNodeForTour);
//        }
//    }
//}


/**
 * subTour is a pointer to tour
 * removes old key from the tours map.
 * what happens if they oldNodeKey and newRootNodeForTour are the same?
 */
void EulerTour::resetRootNodesInSubTourOfTour(std::list<Node *> * subTour){

    int oldNodeKey = (*(subTour->begin()))->getRootNodeOfTour();
    int newRootNodeForTour = subTour->front()->getKey();

    //resetRootNodesInSubTour(subTour);
//    std::cout << newRootNodeForTour << " SUBTOUR " << subTour->size() << std::endl;

//    if (tours.find(newRootNodeForTour) != tours.end()){
//        std::cout << newRootNodeForTour << "     EXISTS " << std::endl;
//        std::cout << newRootNodeForTour << "       size " << tours[newRootNodeForTour].size() << std::endl;
//    }
    tours[newRootNodeForTour] = *subTour;
    // if oldNodeKey and newRootNodeForTour are not the same, then insert will not insert, must delete first and create new entry
    if (oldNodeKey != newRootNodeForTour){

        std::list< Node *> * ptempList = &tours[oldNodeKey];
        ptempList->erase(ptempList->begin(), ptempList->end());
        tours.erase(oldNodeKey); // remove tour from MAP based on key.
    }

    std::list< Node * > * pointerToTour = &tours[newRootNodeForTour];
    for(std::list< Node * >::iterator it = pointerToTour->begin(); it != pointerToTour->end(); ++it) {
        //if ((*it)->getRootNodeOfTour() != newRootNodeForTour){
            (*it)->setPointerToTour(pointerToTour);
        //}
    }

    // std::map<int, std::list< Node *> > tours; // key is the root of the tour
    // std::cout << " SIZE OF TOUR BEFORE ERASE " << tours[oldNodeKey].size() << std::endl;
    // if oldNodeKey == newRootKey cause the subTour is size 1, we got a problem
    // std::cout << " SIZE OF TOUR  AFTER ERASE " << tours[oldNodeKey].size() << " " << tours[newRootNodeForTour].size() << std::endl;

//    if (oldNodeKey != newRootNodeForTour){
//        std::cout << "        OLDNODEKEY " << oldNodeKey << " NEW : " << newRootNodeForTour << std::endl;
//        std::cout << "      SIZE OF NEW REROOT TOUR " << tours[newRootNodeForTour].size() << std::endl;
//        //tours.erase(oldNodeKey); // remove tour from MAP based on key.
//    }

}


/**
 *
 */
void EulerTour::resetAccessed(std::vector<int> * checkedNodes){
    int totalChecked = checkedNodes->size();
    for(int i=0; i<totalChecked; i++){
        nodes.find((*checkedNodes)[i])->second.setAccessed(false);
        //nodes[(*checkedNodes)[i]].setAccessed(false);
    }
}


/**
 * Creates proper Euler Tour where first and last Node are the same
 */
void EulerTour::createSubTour(Node * pNode, std::list< Node * > * subTourToLoad){

    subTourToLoad->push_back(pNode);
    int totalNeighbors = pNode->getTotalNeighbors();
    int totalIterations = 2*totalNeighbors + 1;
    int index=0;

    for(int i=1; i<totalIterations; i++){
        if (i%2 == 0){
            subTourToLoad->push_back(pNode);
        } else {
            // add from neighborhood and not Node
            subTourToLoad->push_back(pNode->getPointerToNeighborByIndex(index));
            index++;
        }
    }
}


/**
 * critical, new node must already be added to nodes!
 */
void EulerTour::createInitialTour(int workingLimit, Model *pModel, std::vector<int>::iterator beginIt) {

    // for each node, make adjacency list
    std::vector<int>::iterator it;
    int neighbor, nodeToInsert;
    Node * pNode;

    for(int i=0; i < workingLimit; i++){
        // make new Node from lattice/bead
        //nodeToInsert = *(*pSelectedLattice + i);
        nodeToInsert = *(beginIt + i);

        nodes.insert ( std::pair<int,Node>(nodeToInsert, Node(nodeToInsert) ) );
        //pNode = &(nodes[nodeToInsert]);
        pNode = &(nodes.find(nodeToInsert)->second);

        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        it = pModel->getPointerToNeighborhood(pNode->getKey());
        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            neighbor = *(it+j);
            //itIndex = std::find(*pSelectedLattice, *pSelectedLattice + currentNodesSize, neighbor);
            // if not found, itIndex will report last
            if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){
                //std::cout << i << " Adding distance : " << distance << " nodes : " << currentNodesSize << " workinglimit : " << workingLimit << std::endl;
                //std::cout << (nodes.find(neighbor)->second).getKey() << std::endl;
                pNode->addNeighbor(&(nodes.find(neighbor)->second)); // if node is found to have a neighbor, add to node's adjacency list
                //pNode->addNeighbor(&(nodes[neighbor])); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == -1) {
                break;
            }
        }

        addToTour(nodeToInsert);
        // create subtour
    } // end of adding beads

    std::cout << " INITIAL TOURS SIZE " << tours.size() << std::endl;
    totalComponents = tours.size();
}

/**
 * search tour for first occurrence and return iterator to the list
 */
std::list< Node * >::iterator EulerTour::getFirstOccurrenceInTour(Node * pNode, std::list<Node *> * pTour){
    return std::find(pTour->begin(), pTour->end(), pNode);
};

std::list< Node * >::reverse_iterator EulerTour::getLastOccurrenceInTour(Node * pNode, std::list<Node *> * pTour){
    return std::find(pTour->rbegin(), pTour->rend(), pNode);
};


/**
 *               9 - 10
 *               | / |
 *   1 - 2 - 6 - 7 - 8
 *   | 5 |
 *   3 - 4
 */
void EulerTour::testSetOne(){

//    nodes.clear();
//    nodes.insert ( std::pair<int,Node>(1, Node(1))  );
//    addToTour(1);
//    nodes.insert ( std::pair<int,Node>(2, Node(2))  );
//    nodes[2].addNeighbor(&(nodes[1]));
//    //nodes[1].addNeighbor(&(nodes[2])); // similarly, add the node to neighbor's list
//    addToTour(2);
//    nodes.insert ( std::pair<int,Node>(3, Node(3))  );
//    nodes[3].addNeighbor(&(nodes[1]));
//    //nodes[1].addNeighbor(&(nodes[3])); // similarly, add the node to neighbor's list
//    addToTour(3);
//    nodes.insert ( std::pair<int,Node>(4, Node(4))  );
//    nodes[4].addNeighbor(&(nodes[2]));
//    //nodes[2].addNeighbor(&(nodes[4]));
//    nodes[4].addNeighbor(&(nodes[3]));
//    //nodes[3].addNeighbor(&(nodes[4]));
//    addToTour(4);
//
//    nodes.insert ( std::pair<int,Node>(5, Node(5))  );
//    nodes[5].addNeighbor(&(nodes[1]));
//    nodes[1].addNeighbor(&(nodes[5]));
//    nodes[5].addNeighbor(&(nodes[2]));
//    nodes[2].addNeighbor(&(nodes[5]));
//    nodes[5].addNeighbor(&(nodes[3]));
//    nodes[3].addNeighbor(&(nodes[5]));
//    nodes[5].addNeighbor(&(nodes[4]));
//    nodes[4].addNeighbor(&(nodes[5]));
//    addToTour(5);
//
//
//    nodes.insert ( std::pair<int,Node>(7, Node(7))  );
//    //nodes[7].addNeighbor(&(nodes[6]));
////    nodes[6].addNeighbor(&(nodes[7]));
//    addToTour(7);
//
//    nodes.insert ( std::pair<int,Node>(8, Node(8))  );
//    nodes[8].addNeighbor(&(nodes[7]));
//    nodes[7].addNeighbor(&(nodes[8]));
//    addToTour(8);
//
//    nodes.insert ( std::pair<int,Node>(9, Node(9))  );
//    nodes[9].addNeighbor(&(nodes[7]));
//    nodes[7].addNeighbor(&(nodes[9]));
//    addToTour(9);
//
//    nodes.insert ( std::pair<int,Node>(10, Node(10))  );
//    nodes[10].addNeighbor(&(nodes[7]));
//    //nodes[7].addNeighbor(&(nodes[10]));
//    nodes[10].addNeighbor(&(nodes[8]));
//    //nodes[8].addNeighbor(&(nodes[10]));
//    nodes[10].addNeighbor(&(nodes[9]));
//    //nodes[9].addNeighbor(&(nodes[10]));
//    addToTour(10);
//
//    nodes.insert ( std::pair<int,Node>(6, Node(6))  );
//    nodes[6].addNeighbor(&(nodes[2]));
//    nodes[6].addNeighbor(&(nodes[7]));
//    nodes[6].addNeighbor(&(nodes[3]));
//    nodes[6].addNeighbor(&(nodes[10]));
//    addToTour(6);
//
//    std::cout << " BEFORE TOTAL NUMBER OF COMPONENTS " << tours.size() << std::endl;
//    this->removeNode(6);
//    this->removeNode(10);
//    std::cout << "  AFTER TOTAL NUMBER OF COMPONENTS " << tours.size() << std::endl;
//    for(std::map<int, Node>::iterator it = nodes.begin(); it!=nodes.end(); ++it){
//        std::cout << " NODE : " << (*it).first << std::endl;
//    }
//
//
//    nodes.insert ( std::pair<int,Node>(6, Node(6))  );;
//    nodes[6].addNeighbor(&(nodes[3]));
//    nodes[6].addNeighbor(&(nodes[7]));
//    addToTour(6);

    // should point same thing
//    std::list< Node * > * pExistingNeighborTour = &tours[7];
//    std::cout << "FRONT " << pExistingNeighborTour->front()->getKey() << std::endl;
//    std::cout << "BEFORE ADDRESS " << pExistingNeighborTour << " == " << &tours[7] << std::endl;
//    pExistingNeighborTour->pop_front();
//    std::cout << "FRONT " << pExistingNeighborTour->front()->getKey() << " == " << tours[7].front()->getKey() << std::endl;
//    std::cout << "AFTER  ADDRESS " << pExistingNeighborTour << " == " << &tours[7] << std::endl;

    std::cout << " ******* EULER TOUR ******* " << std::endl;
    std::cout << " TOTAL NUMBER OF COMPONENTS " << tours.size() << std::endl;
    std::cout << " ******* EULER TOUR ******* " << std::endl;

}

// check against nodes list, should match
bool EulerTour::checkNodesList(std::set<int> * beads_in_use){
    int size = beads_in_use->size();
    int nodeSize = nodes.size();

    if (nodeSize != size){
        std::cout << " SIZES DO NOT MATCH IN TOUR : beads_in_use => " << size << " != " << nodeSize << std::endl;
        return true;
    }

    for(std::set<int>::iterator it = beads_in_use->begin(); it != beads_in_use->end(); ++it){
        if (nodes.find(*it) == nodes.end()){
            std::cout << " NODE NOT FOUND : bead_index => " << *it << std::endl;
          return true;
        }
    }

    return false;
}

void EulerTour::checkTourSize(std::string note){

    for(std::map<int, std::list< Node *> >::iterator it=tours.begin(); it!=tours.end(); ++it){
        if ((*it).second.size() == 0){
            std::cout << "EMPTY TOUR " << (*it).first << std::endl;
            std::cout << note <<  std::endl;
            exit(0);
        }
    } // key is the root of the tour

}

