//
// Created by Robert Rambo on 11/01/2017.
//
#include "EulerTour.h"
#include "../Model.h"

// Constructor
EulerTour::EulerTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel){
    pSelectedLattice = &beginIt;
    this->createInitialTour(workingLimit, pModel);
}


/**
 * add a node to existing set of nodes
 * first, create neighborhood from existing nodes
 * find which tour to tours to add to and merge any tours that are bridged by new node
 * returns number of nodes in tour
 */
int EulerTour::addNode(int latticePoint, Model *pModel) {

    //std::vector<int>::iterator itIndex;
    int neighbor;
    int newNode = latticePoint;

    nodes.insert ( std::pair<int,Node>(newNode, Node(newNode) ) );
    Node * pNode = &(nodes[newNode]);

    // build neighborhood
    auto it = pModel->getPointerToNeighborhood(newNode);
    for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
        // if neighbor is inside workinglimit, don't add
        neighbor = *(it+j); // index of bead used as key
        // if not found, itIndex will report last
        if (nodes.find(neighbor) != nodes.end()){ // look in nodes list

            pNode->addNeighbor(&(nodes[neighbor])); // add neighbor to both parent and child
            // nodes[neighbor].addNeighbor(pNode);     // has no idea which list it is in

        } else if (neighbor == -1){
            break;
        }
    }

    // add to Tour
    addToTour(newNode);
    totalComponents = tours.size();
    return totalComponents;

}

void EulerTour::createBackUp(){

//    std::map<int, Node> backUpedNodes; // backup of Nodes
    std::clock_t start = std::clock();

    //cout << "REMOVE old : " << oldmeth/(double) CLOCKS_PER_SEC << " new : " << newmeth/(double) CLOCKS_PER_SEC << endl;
    backedUpNodes.clear();
    for(std::map<int, Node>::iterator it = nodes.begin(); it != nodes.end(); it++){
        backedUpNodes.insert( std::pair<int,Node>(it->first, it->second )) ;
    }

    std::cout << "COPY MAP TIME : " << (std::clock() - start)/(double) CLOCKS_PER_SEC << std::endl;

//    std::list< Node * > backUpedTours; //

}

/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 */
bool EulerTour::addToTour(int nodeToAdd){

    Node *pNeighbor, * pNode = &(nodes[nodeToAdd]);
    // add to Tour

    //if (validateList("********* BEFORE ADDING NEW TOUR MEMBERS")){
    //    return false;
    //}

    if (pNode->getTotalNeighbors() == 0){
        // add node to tours as new tour
        // set node tour index mapping
        tours.insert ( std::pair< int, std::list<Node *> >(pNode->getKey(), std::list< Node *>()) );  // key is the root of the tour
        tours[pNode->getKey()].push_back(&nodes[nodeToAdd]);

        pNode->setPointerToTour(&tours[pNode->getKey()]); // sets pointer to list< Node *> in vector
        pNode->setRootNodeOfTour(pNode->getKey());

    } else {
        // create subtour of new Node
        std::list< Node * > subTour;
        createSubTour(pNode, &subTour); // create subtour using pNode's neighbors

        int totalNeighbors = pNode->getTotalNeighbors();
        int keyOfNeighborNode = pNode->getKey();
        // new Node could bridge 1 or more prior nodes
        std::set< int > numberOfNeighborhoods;

        for (int j=0; j < totalNeighbors; j++){
            pNeighbor = pNode->getPointerToNeighborByIndex(j);
            // get tour index to merge into
            numberOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour());
            // get pointer to list of neighbor, is it the same?
            keyOfNeighborNode = pNeighbor->getKey();
        }

        // set tour membership of new Node
        pNeighbor = &nodes[keyOfNeighborNode];
        pNode->setPointerToTour(pNeighbor->getPointerToTour()); // sets pointer to list< Node *> in vector
        pNode->setRootNodeOfTour(pNeighbor->getRootNodeOfTour());

        const int rootToBaseTour = pNode->getRootNodeOfTour();
        // re-root to node found incommon to a previous tour in tours
        rerootSubTour(&nodes[keyOfNeighborNode], &subTour);
        // substitute in tours[i]
        std::list< Node * > * pExistingNeighborTour = &tours[rootToBaseTour];

        std::list< Node * >::iterator inTour = std::find (pExistingNeighborTour->begin(), pExistingNeighborTour->end(), pNeighbor);
        subTour.pop_back();
        // merge subTour of newNode with existing tour
        pExistingNeighborTour->splice(inTour, subTour); // Add additional tours to base
        // any additional tours will be greater than minTourIndex
        if (numberOfNeighborhoods.size() > 1){ // if number of neighborhoods is greater than 1, it means the new node bridges at least two
            // find tour greater than minTourIndex, reroot tour and merge
            for (int j=0; j < totalNeighbors; j++){

                pNeighbor = pNode->getPointerToNeighborByIndex(j);
                int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();

                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again
                    std::list< Node * > * pList = &tours[rootOfCurrentTour];

                    rerootSubTour(pNeighbor, pList);

                    inTour = std::find (pExistingNeighborTour->begin(), pExistingNeighborTour->end(), pNeighbor);
                    pList->pop_back();
                    // better to iterate over node list and check for tour membership
                    // for all nodes in this tour, reset tourIndex
                    for (std::map<int, Node>::iterator iterator = nodes.begin(); iterator != nodes.end(); ++iterator) {
                        if (((*iterator).second).getRootNodeOfTour() == rootOfCurrentTour){
                            (*iterator).second.setRootNodeOfTour(rootToBaseTour);
                        }
                    }
                    // merge pList into existingNeighborTour
                    pExistingNeighborTour->splice(inTour, *pList);
                    // if pExistingNeighborTour is single element, splicing pushes the node to the back
                    // remove, erasing from vector changes address of the other elements,
                    tours.erase(rootOfCurrentTour);
                    tours[rootToBaseTour] =  *pExistingNeighborTour;
                    //if (validateList("FROM ADDING NEW TOUR MEMBERS")){
                    //    return false;
                    //};
                }
            }
        } // pointer deallocation problem exiting block

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
 *
 *
 */
void EulerTour::rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    if (subTourToLoad->front()->getKey() != newRoot->getKey()){
        //printList("================> REROOT_AND_REASSIGN_SUBTOUR REROOT Before", subTourToLoad);
        std::list< Node * >::iterator inTour = std::find (subTourToLoad->begin(), subTourToLoad->end(), newRoot);
        subTourToLoad->pop_front();
        //printList("SUBTOUR REROOT ", subTourToLoad);

        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour);
        tempList.push_back(newRoot);
        // copy to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);

        resetRootNodesInSubTourOfTour(subTourToLoad);
    }
}

/**
 *
 *  smallest tour size is 1, then 3, there is no two
 */
void EulerTour::rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    std::list< Node * >::iterator inTour = std::find (subTourToLoad->begin(), subTourToLoad->end(), newRoot);
    std::list< Node * > tempList;

    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    if (subTourToLoad->front()->getKey() != newRoot->getKey() && (subTourToLoad->size() > 2)){

        subTourToLoad->pop_front();
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour);
        tempList.push_back(newRoot);
        // copy to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
    }
}


/**
 * Find node in nodes list
 * for each of its neighbor in neighborhood
 * if removing node creates a new tour, add to tours
 */
int EulerTour::removeNode(int indexOfNode){

    Node * pNodeToRemove = &(nodes[indexOfNode]), * pNeighbor;
    std::vector<int> checkedNodes;
    std::list< Node * > * pTour;// = &tours[pNodeToRemove->getRootNodeOfTour()];  // tour with node to remove
    // delete pointer in vector
    std::string result;

    if (pNodeToRemove->getTotalNeighbors() > 0){
        // remove node
        std::unordered_set <int> remainderNodes;

        while (pNodeToRemove->getTotalNeighbors() > 0) {  // for each neighbor, remove edge with node_to_remove

            pNeighbor = pNodeToRemove->getPointerToNeighborByIndex(0);
            pTour = &tours[pNodeToRemove->getRootNodeOfTour()];  // get tour that contains node to remove
            this->rerootAndReassignSubTour(pNeighbor, pTour);    // reroot and reassign subTour in Tours which means pointer changes
            // pTour may no longer exists
            pTour = &tours[pNeighbor->getRootNodeOfTour()];      // reassign pointer

            std::list< Node * >::iterator pFirst = getFirstOccurrenceInTour(pNodeToRemove, pTour);
            std::list< Node * >::reverse_iterator pLast = getLastOccurrenceInTour(pNodeToRemove, pTour);
            // search from start iterate through
            bool terminateNow = false;
            while ( (pFirst = std::find(pFirst, pTour->end(), pNodeToRemove)) != pTour->end()) {
                // Do something with iter
                auto pPreviousFirst = std::next(pFirst, -1);
                auto pNextFirst = std::next(pFirst, 1);

                if (((*pPreviousFirst)->getKey() == pNeighbor->getKey()) || ((*pNextFirst)->getKey() == pNeighbor->getKey())){ //
                    // search for Neighbor in reverse
                    while ( (pLast = std::find(pLast, pTour->rend(), pNodeToRemove)) != pTour->rend()) { // search for Neighbor in reverse
                        // Do something with iter
                        auto pPreviousLast = std::next(pLast, 1);
                        auto pNextLast = std::next(pLast, -1);
                        // from start
                        if (((*pPreviousLast)->getKey() == pNeighbor->getKey()) || ((*pNextLast)->getKey() == pNeighbor->getKey())){ // if true, implies edge found
                            //std::cout << (*pPreviousFirst)->getKey() << " pPreviousLast SEQUENCE :" << (*pPreviousLast)->getKey() << std::endl;
                            //std::cout << (*pFirst)->getKey() << "         pLast SEQUENCE :" << (*pLast)->getKey() << std::endl;
                            //std::cout << (*pNextFirst)->getKey() << "     pNextLast SEQUENCE :" << (*pNextLast)->getKey()<< std::endl;
                            //std::cout << " LAST SEQUENCE BASE :" << (*pNextLast.base())->getKey()<< std::endl;
                            // have two cases to distinguish
                            // case 1:   Y ... Y-X ... X-Y ... Y
                            // case 2:   Y ... X-Y ... Y-X ... Y
                            if (((*pPreviousFirst)->getKey() == pNeighbor->getKey()) && ((*pNextLast)->getKey() == pNeighbor->getKey())) { // case 1:   Y ... Y-X ... X-Y ... Y
                                //             Last -> NodeToRemove
                                //        pNextLast -> Node +1 to Last
                                // pNextLast.base() -> Node +2 to Last
                                // pLast is a reverse iterator, decrementing moves the location up, so pod points to pNextLast
                                auto pod = pPreviousLast.base(); // will potential go past the list, equivalent to list.end()

                                if (pFirst == pPreviousLast.base()){ // do they point to same place?
                                    pTour->erase(pPreviousFirst, pLast.base());
                                    // XYX = > X where Y is neighbor of X
                                } else {
                                    remainderNodes.clear();
                                    // populate remainderNodes map
                                    for(std::list< Node * >::iterator it = pTour->begin(); it != pFirst; ++it) {
                                        remainderNodes.insert((*it)->getKey());
                                    }

                                    for(std::list< Node * >::reverse_iterator it = pTour->rbegin(); it != pLast; ++it) {
                                        remainderNodes.insert((*it)->getKey());
                                    }

                                    bool found = false;
                                    std::unordered_set<int>::iterator got;

                                    auto pSearch = pFirst; // make copy of pFirst, pFirst should be pointing to Node_to_Remove
                                    // pod points to
                                    for(; pSearch != pod; ++pSearch) {
                                        // search to see if a bounded node can be found in remainderNodes
                                        // if true, break out and substitute
                                        if (!((*pSearch)->getAccessed())){
                                            got = remainderNodes.find ( (*pSearch)->getKey() );
                                            (*pSearch)->setAccessed(true);
                                            checkedNodes.push_back((*pSearch)->getKey());

                                            if (got != remainderNodes.end()){
                                                found = true;
                                                break;
                                            }
                                        }
                                    }

                                    std::list< Node * > tempList;
                                    // find first occurrence of newRoot
                                    tempList.splice(tempList.begin(), *pTour, pFirst, --pNextLast.base());
                                    pTour->erase(--pNextLast.base()); // removes duplicate nodes next two each other

                                    if (found){ // make substitution
                                        // reroot subtour bounded by neighbors and substitute
                                        rerootSubTour(&nodes[*got], &tempList);
                                        tempList.pop_back();
                                        std::list< Node * >::iterator gotIt = std::find(pTour->begin(), pTour->end(), &nodes[*got]);
                                        pTour->splice(gotIt, tempList);

                                    } else { // create new list, split out and reset root nodes of items in new list
                                        // reassign the root node in Nodes of tempList before adding to tours
                                        resetRootNodesInSubTour(&tempList);
                                        // add tempList to tours
                                        tours.emplace ( std::pair< int, std::list<Node *> >((*tempList.begin())->getKey(), tempList) );  // key is the root of the tour
                                    }
                                }

                            } else if (((*pNextFirst)->getKey() == pNeighbor->getKey()) && ((*pPreviousLast)->getKey() == pNeighbor->getKey())) { // case 2:   Y ... X-Y ... Y-X ... Y
                                // does remainderNodes need to be cleared for each neighbor
                                // populate remainderNodes map
                                // pFirst and pLast points to NodeToRemove

                                remainderNodes.clear();
                                for(std::list< Node * >::iterator it = pTour->begin(); it != pFirst; ++it) {
                                    remainderNodes.insert((*it)->getKey());
                                }

                                for(std::list< Node * >::reverse_iterator it = pTour->rbegin(); it != pLast; ++it) {
                                    remainderNodes.insert((*it)->getKey());
                                }

                                bool found = false;
                                std::unordered_set<int>::iterator got;

                                auto pSearch = pFirst; // make copy of pFirst, pFirst should be pointing to NeighborNode
                                auto pod = pPreviousLast.base();
                                // pod points to
                                ++pSearch;
                                for(; pSearch != pod; ++pSearch) {
                                    // search to see if a bounded node can be found in remainderNodes
                                    // if true, break out and substitute
                                    if (!((*pSearch)->getAccessed()) ){
                                        got = remainderNodes.find ( (*pSearch)->getKey() );
                                        (*pSearch)->setAccessed(true);
                                        checkedNodes.push_back((*pSearch)->getKey());

                                        if (got != remainderNodes.end()){
                                            found = true;
                                            break;
                                        }
                                    }
                                }

                                std::list< Node * > tempList;
                                // find first occurrence of newRoot
                                tempList.splice(tempList.begin(), *pTour, pNextFirst, pPreviousLast.base()); // exclusive of pPreviousLast.base()
                                pTour->erase(--pLast.base()); // removes duplicate nodes next two each other

                                if (found){ // make substitution
                                    // reroot subtour bounded by neighbors and substitute
                                    if (tempList.size() > 1){
                                        rerootSubTour(&nodes[*got], &tempList);
                                        tempList.pop_back();
                                        std::list< Node * >::iterator gotIt = std::find(pTour->begin(), pTour->end(), &nodes[*got]);
                                        pTour->splice(gotIt, tempList);
                                    }
                                } else { // break out and return false
                                    // pFirst points to node to remove
                                    // pLast points to node to remove
                                    resetRootNodesInSubTour(&tempList);
                                    tours.emplace ( std::pair< int, std::list<Node *> >((*tempList.begin())->getKey(), tempList) );  // key is the root of the tour
                                }

                            } // end of if statement

                            pNodeToRemove->removeNeighbor(pNeighbor);  // nodeToRemove exists until we have removed all neighbors
                            resetAccessed(&checkedNodes);
                            pLast = --(pTour->rend());
                            terminateNow = true;
                            break;
                        } // end if statement, if above statement excuted, it means we made the substitution and move to next neighbor
                        ++pLast;
                    } // while NodeToRemove loop LAST
                } // end of if statement that checks if neighbor node is next to node_to_remove
                if (terminateNow){
                    break;
                }
                pFirst++;
            } // while NodeToRemove loop FIRST
            //result = "remove node " + std::to_string(pNodeToRemove->getKey()) + " " + std::to_string(pNodeToRemove->getTotalNeighbors());
        } // while loop for each neighbors
    } else {
        tours.erase(pNodeToRemove->getRootNodeOfTour());
    }

    // remove the node from nodes
    if (pNodeToRemove->getTotalNeighbors() == 0){
        //tours.erase(pNodeToRemove->getRootNodeOfTour());
        nodes.erase(pNodeToRemove->getKey());
    }

    totalComponents = tours.size();
    return totalComponents;
}


bool EulerTour::validateTour(std::list<Node *> * tourtocheck){
    // each pair of numbers needs a reverse sequence
    // go through each element of list and validate it is a node
    std::map<int, int> pairs;
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

bool EulerTour::validateNodes(){
    for (it_type it = nodes.begin(); it != nodes.end(); it++){
        // check for repeats
        Node tempNode = it->second;
        if (!(it->second).validate()){
            std::cout << " INVALID NODE " << it->first << " " << tempNode.getKey() << std::endl;
            return false;
        }
    }
    return true;
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

/**
 * subTour is a pointer to tour
 * newRootNode
 */
void EulerTour::resetRootNodesInSubTour(std::list<Node *> * subTour){

    int newRootNodeForTour = (*(subTour->begin()))->getKey();

    for(std::list< Node * >::iterator it = subTour->begin(); it != subTour->end(); ++it) {
        (*it)->setRootNodeOfTour(newRootNodeForTour);
    }
}


/**
 * subTour is a pointer to tour
 * removes old key from the tours map.
 */
void EulerTour::resetRootNodesInSubTourOfTour(std::list<Node *> * subTour){

    int oldNodeKey = (*(subTour->begin()))->getRootNodeOfTour();
    int newRootNodeForTour = (*(subTour->begin()))->getKey();

    resetRootNodesInSubTour(subTour);

    std::list< Node *> tempList = tours[oldNodeKey];
    tours[newRootNodeForTour] = tempList;
    tours.erase(oldNodeKey);
    // reassign pointer
    // subTour = &tours[newRootNodeForTour];
}


/**
 *
 */
void EulerTour::resetAccessed(std::vector<int> * checkedNodes){
    int totalChecked = checkedNodes->size();
    for(int i=0; i<totalChecked; i++){
        nodes[(*checkedNodes)[i]].setAccessed(false);
    }
}


/**
 *
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
void EulerTour::createInitialTour(int workingLimit, Model *pModel) {

    // for each node, make adjacency list
    std::vector<int>::iterator it;
    int neighbor, nodeToInsert;
    Node * pNode;

    for(int i=0; i < workingLimit; i++){
        // make new Node from lattice/bead
        nodeToInsert = *(*pSelectedLattice + i);
        nodes.insert ( std::pair<int,Node>(nodeToInsert, Node(nodeToInsert) ) );
        pNode = &(nodes[nodeToInsert]);

        // create the neighborhood for nodes in nodes list
        it = pModel->getPointerToNeighborhood(pNode->getKey());
        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            neighbor = *(it+j);
            //itIndex = std::find(*pSelectedLattice, *pSelectedLattice + currentNodesSize, neighbor);
            // if not found, itIndex will report last
            //distance = std::distance(*pSelectedLattice, itIndex);
            //if ( (neighbor != -1) && (distance < currentNodesSize) ) {
            if (nodes.find(neighbor) != nodes.end()){
                //std::cout << i << " Adding distance : " << distance << " nodes : " << currentNodesSize << " workinglimit : " << workingLimit << std::endl;
                //std::cout << (nodes.find(neighbor)->second).getKey() << std::endl;
                pNode->addNeighbor(&(nodes[neighbor])); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == -1) {
                break;
            }
        }

        addToTour(nodeToInsert);
        // create subtour
    } // end of adding beads

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

    nodes.clear();
    nodes.insert ( std::pair<int,Node>(1, Node(1))  );
    addToTour(1);

    nodes.insert ( std::pair<int,Node>(2, Node(2))  );
    nodes[2].addNeighbor(&(nodes[1]));
    nodes[1].addNeighbor(&(nodes[2])); // similarly, add the node to neighbor's list
    addToTour(2);

    nodes.insert ( std::pair<int,Node>(3, Node(3))  );
    nodes[3].addNeighbor(&(nodes[1]));
    nodes[1].addNeighbor(&(nodes[3])); // similarly, add the node to neighbor's list
    addToTour(3);

    nodes.insert ( std::pair<int,Node>(4, Node(4))  );
    nodes[4].addNeighbor(&(nodes[2]));
    nodes[2].addNeighbor(&(nodes[4]));
    nodes[4].addNeighbor(&(nodes[3]));
    nodes[3].addNeighbor(&(nodes[4]));
    addToTour(4);

    nodes.insert ( std::pair<int,Node>(5, Node(5))  );
    nodes[5].addNeighbor(&(nodes[1]));
    nodes[1].addNeighbor(&(nodes[5]));
    nodes[5].addNeighbor(&(nodes[2]));
    nodes[2].addNeighbor(&(nodes[5]));
    nodes[5].addNeighbor(&(nodes[3]));
    nodes[3].addNeighbor(&(nodes[5]));
    nodes[5].addNeighbor(&(nodes[4]));
    nodes[4].addNeighbor(&(nodes[5]));
    addToTour(5);

    nodes.insert ( std::pair<int,Node>(6, Node(6))  );
    nodes[6].addNeighbor(&(nodes[2]));
    nodes[2].addNeighbor(&(nodes[6]));
    addToTour(6);

    nodes.insert ( std::pair<int,Node>(7, Node(7))  );
    nodes[7].addNeighbor(&(nodes[6]));
    nodes[6].addNeighbor(&(nodes[7]));
    addToTour(7);

    nodes.insert ( std::pair<int,Node>(8, Node(8))  );
    nodes[8].addNeighbor(&(nodes[7]));
    nodes[7].addNeighbor(&(nodes[8]));
    addToTour(8);

    nodes.insert ( std::pair<int,Node>(9, Node(9))  );
    nodes[9].addNeighbor(&(nodes[7]));
    nodes[7].addNeighbor(&(nodes[9]));
    addToTour(9);

    nodes.insert ( std::pair<int,Node>(10, Node(10))  );
    nodes[10].addNeighbor(&(nodes[7]));
    nodes[7].addNeighbor(&(nodes[10]));
    nodes[10].addNeighbor(&(nodes[8]));
    nodes[8].addNeighbor(&(nodes[10]));
    nodes[10].addNeighbor(&(nodes[9]));
    nodes[9].addNeighbor(&(nodes[10]));
    addToTour(10);

    std::cout << " ******* EULER TOUR ******* " << std::endl;
    std::cout << " TOTAL NUMBER OF COMPONENTS " << tours.size() << std::endl;
    std::cout << " ******* EULER TOUR ******* " << std::endl;

}


