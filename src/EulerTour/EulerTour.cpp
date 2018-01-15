//
// Created by Robert Rambo on 11/01/2017.
//
#include "EulerTour.h"
#include "../Model.h"

// Constructor
EulerTour::EulerTour(){};

EulerTour::EulerTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel){
    //nodes.clear();
    //tours.clear();
    this->createInitialTour(workingLimit, pModel, beginIt);
}

EulerTour::EulerTour(std::set<int> & indices, Model *pModel){
    //nodes.clear();
    //tours.clear();
    // for each node, make adjacency list
    //std::vector<int>::iterator vit;
    int neighbor, nodeToInsert;
    Node * pNode;

    for(std::set<int>::iterator sit = indices.begin(); sit != indices.end(); ++sit){
        // make new Node from lattice/bead
        nodeToInsert = *sit;

        nodes.insert ( std::pair<int, Node>(nodeToInsert, Node(nodeToInsert) ) );
        pNode = &(nodes.find(nodeToInsert)->second);

        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        auto vit = pModel->getPointerToNeighborhood(pNode->getKey());
        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            neighbor = *(vit+j);
            // if not found, itIndex will report last
            if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){
                pNode->addNeighbor(&(nodes.find(neighbor)->second)); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == -1) {
                break;
            }
        }

        addToTour(nodeToInsert);
    } // end of adding beads

    //std::cout << " => INITIAL TOURS SIZE " << tours.size() << std::endl;
    totalComponents = tours.size();
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
int EulerTour::addNode(int newNode, Model * pModel) {

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
            pNode->addNeighbor(&(nodes.find(neighbor)->second)); // add neighbor to both parent and child
            // again nothing prevents me from finding nodes.end()
            // has no idea which list it is in
        } else if (neighbor == -1){
            break;
        }
    }

    //printTourSizes();
    //std::cout << " NEW NODE " << newNode << " => " << pNode->getTotalNeighbors() <<  "  key of new node " << (nodes.find(newNode)->second).getKey() << std::endl;
    // add to Tour
    addToTour(newNode);
    totalComponents = tours.size();
    //totalComponents = simpleTours.size();
//    validateNodes("AFTER ADDING NODE : " + std::to_string(newNode));
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


/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 */
bool EulerTour::addToTour(int nodeToAdd){

    // get the new node from NODES map list
    Node *pNeighbor, * pNode = &(nodes.find(nodeToAdd)->second);
    // add to Tour
    //if (validateList("********* BEFORE ADDING NEW TOUR MEMBERS")){
    //    return false;
    //}

    if (pNode->getTotalNeighbors() == 0){ // if no existing neighbors, make a new tour
        // add node to tours as new tour
        // set node tour index mapping
        // create tour (empty)
//        std::cout << " new tour " << std::endl;
        tours.insert ( std::pair< int, std::list<Node *> >(pNode->getKey(), std::list< Node *>()) );  // key is the root of the tour
        pNode->setRootNodeOfTour(pNode->getKey());
        tours.find(pNode->getKey())->second.push_back(pNode);
//        printList("NEWLY ADDED SINGLE ", &tours[nodes.find(nodeToAdd)->second.getRootNodeOfTour()]);
    } else {
        // create subtour of new Node using Node's found neighbors
        std::list< Node * > subTour;  // subtour is created and has no membership to existing tour
        createSubTour(pNode, &subTour); // create subtour using pNode's neighbors

        int totalNeighbors = pNode->getTotalNeighbors(); // neighborhood was created when adding Node, checks existing nodes for proximity
        int keyOfNeighborNode = pNode->getKey();
        // new Node could bridge 1 or more prior nodes
        std::set< int > rootOfNeighborhoods;

        for (int j=0; j < totalNeighbors; j++){
            pNeighbor = &nodes.find(pNode->getPointerToNeighborByIndex(j))->second;
            // get tour index to merge into
            rootOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour()); // how many neighborhoods does newNode belong to?
            // get pointer to list of neighbor, is it the same?
            keyOfNeighborNode = pNeighbor->getKey(); // use last neighbor as start
        }
        // what happens if node that is removed is root node for a tour?
        // set tour membership of new Node
        pNeighbor = &(nodes.find(keyOfNeighborNode)->second);
        rerootSubTour(pNeighbor, &subTour); //prepare for merging in pNeighbor tour (doesn't change tour membership of nodes, still points to old location)
        // subTour nodes do not have root at this point
//        std::cout << " size of pNeighboor tour " << tours.find(pNeighbor->getRootNodeOfTour())->second.size() << std::endl;

        //pNode->setPointerToTour(pNeighbor->getPointerToTour()); // sets pointer to list< Node *> in vector
        //nodes.find(nodeToAdd)->second.setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
//        pNode->setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
        pNode->setRootNodeOfTour(pNeighbor->getRootNodeOfTour()); // new
        const int rootToBaseTour = pNode->getRootNodeOfTour();
        // Merge the new subTour with existing tour specified by Neighbor
        // re-root to node found incommon to a previous tour in tours
        // merge in tours[i]
        // get base tour of neighbor to merge into
        //std::list< Node * > * pExistingNeighborTour = pNeighbor->getPointerToTour();
        std::list< Node * > * pExistingNeighborTour = &tours.find(pNeighbor->getRootNodeOfTour())->second;
        int existingRootNode = tours.find(pNeighbor->getRootNodeOfTour())->first;

        // subTour is rooted to pNode
        //               subTour => 58565
        //
        // pExistingNeighborTour => 124313426797621   (root node is 1)
        //
        // reroot        subTour => 65856
        //
        //std::list< Node * >::iterator inTour = std::find (pExistingNeighborTour->begin(), pExistingNeighborTour->end(), pNeighbor);
//        std::cout << " EulerTour::addToTour ADD " << nodeToAdd << " pNeighbor " << pNeighbor->getKey() << " tour size : " << pExistingNeighborTour->size() << std::endl;

        std::list< Node * >::iterator inTour = std::find_if(pExistingNeighborTour->begin(), pExistingNeighborTour->end(), find_Node_by_key(pNeighbor->getKey()));
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
//                pNeighbor = pNode->getPointerToNeighborByIndex(j);
                pNeighbor = &nodes.find(pNode->getPointerToNeighborByIndex(j))->second;
                int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();
//                std::cout << j << " neighbor " << pNeighbor->getKey() << std::endl;
                // if two nodes belong to same tour, only need to add tour once
                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again
                    std::list< Node * > * pList = &(tours.find(rootOfCurrentTour)->second);
                    // pList should not be zero size
                    rerootSubTour(pNeighbor, pList);
                    inTour = std::find_if(pExistingNeighborTour->begin(), pExistingNeighborTour->end(),find_Node_by_key(pNeighbor->getKey()));
                    // better to iterate over node list and check for tour membership
                    // for all nodes in this tour, reset tourIndex
                    for(std::list< Node * >::iterator pIt = pList->begin(); pIt != pList->end(); ++pIt){
                        // SLOW STEP (CAN BE OPTIMIZED)
                        //if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
                        //  (*pIt)->setPointerToTour(pExistingNeighborTour);
                        //nodes.find((*pIt)->getKey())->second.setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
//                        nodes.find((*pIt)->getKey())->second.setPointerToTour(pExistingNeighborTour); //update node
                        nodes.find((*pIt)->getKey())->second.setRootNodeOfTour(existingRootNode);
                        //(*pIt)->setRootNodeOfTour(pExistingNeighborTour->front()->getKey()); // new
                        //}
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

//        printList("NEWLY ADDED", &tours[nodes.find(nodeToAdd)->second.getRootNodeOfTour()]);
        // after all neighbors have been added, we then update the nodes tour
        // only want to iterate over the unique nodes in the tour rather than all members of the tour!
        //
//        for(std::list< Node * >::iterator pIt = pExistingNeighborTour->begin(); pIt != pExistingNeighborTour->end(); ++pIt){
//            // SLOW STEP (CAN BE OPTIMIZED)
//            if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
//                (*pIt)->setPointerToTour(pExistingNeighborTour);
//                std::cout << " UPDATING " << std::endl;
//                exit(0);
//            }
//        }

    } // end of creating subtour

    return true;
}



/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 */
//bool EulerTour::addToSimpleTour(int nodeToAdd){
//
//    Node *pNeighbor, * pNode = &(nodes.find(nodeToAdd)->second);
//    // add to Tour
//    if (pNode->getTotalNeighbors() == 0){
//        // add node to tours as new tour
//        // set node tour index mapping
//        // create tour (empty)
//        simpleTours.insert( std::pair< int, SimpleTour>(pNode->getKey(), SimpleTour(pNode)));
//    } else {
//        // create subtour of new Node
//        // then check if new node is found an previous tour and merge
//        std::list< Node * > subTour;
//        //std::cout << " Creating subTour for " << nodeToAdd << " T_n : " << pNode->getTotalNeighbors() << std::endl;
//        createSubTour(pNode, &subTour); // create subtour using pNode's neighbors
//        //printList("/nsubtour of : " + std::to_string(pNode->getKey()), &subTour);
//        int totalNeighbors = pNode->getTotalNeighbors(); // neighborhood was created when adding Node, checks existing nodes for proximity
//        int keyOfNeighborNode = pNode->getKey();
//        // new Node could bridge 1 or more prior nodes
//        std::set< int > rootOfNeighborhoods;
//
//        for (int j=0; j < totalNeighbors; j++){
//            pNeighbor = pNode->getPointerToNeighborByIndex(j);
//            // get tour index to merge into
//            rootOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour()); // how many neighborhoods does newNode belong to?
//            // get pointer to list of neighbor, is it the same?
//            keyOfNeighborNode = pNeighbor->getKey();
//        }
//
//        // what happens if node that is removed is root node for a tour?
//        // set tour membership of new Node
//        pNeighbor = &(nodes.find(keyOfNeighborNode)->second);
//        // get SimpleTour of pNeighbor
//        SimpleTour * pExistingNeighborTour = &(simpleTours.find(pNeighbor->getRootNodeOfTour())->second);
//
//        rerootSubTour(pNeighbor, &subTour); //prepare for merging in pNeighbor tour (doesn't change tour membership of nodes, still points to old location)
//        // insert subTour in pExistingTour at pNeighbor position
//        pExistingNeighborTour->insertSubTour(pNeighbor, &subTour); // don't update the neighbors, neighbors carry the link to the other tours
//
//        pNode->setPointerToTour(pExistingNeighborTour->getPointerToTour()); // sets pointer to list< Node *> in vector
//        const int rootToBaseTour = pNode->getRootNodeOfTour();
//
//        // nodes in this tour will have mixed roots with pNeighbor
//        // if new node bridges two or more tours, pExistingNeighborTour will be mixed
//        if (rootOfNeighborhoods.size() > 1){ // if number of neighborhoods is greater than 1, it means the new node bridges at least two
//            // find tour greater than minTourIndex, reroot tour and merge
//            for (int j=0; j < (totalNeighbors-1); j++){
//                // if two neighbors are in same neighborhood, and I add the first one
//                // then i have a problem if rootNodeOfTour is not updated before next neighbor is checked
//                pNeighbor = pNode->getPointerToNeighborByIndex(j);
//                int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();
//                // if two nodes belong to same tour, only need to add tour once
//                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again
//                    SimpleTour * tempTour = &(simpleTours.find(rootOfCurrentTour)->second);
//                    if (simpleTours.find(rootOfCurrentTour) ==  simpleTours.end()){
//                        std::cout << " REACHED END OF FIND " << rootOfCurrentTour << " " <<  std::endl;
//                        exit(0);
//                    }
//                    pExistingNeighborTour->mergeTour(pNeighbor, tempTour);
//                    // pList must be zero in size
//                    simpleTours.erase(rootOfCurrentTour);
//                }
//            }
//        }
//
//        //check node test
////        int count=0;
////        for (std::list<Node *>::iterator it = pExistingNeighborTour->getPointerToTour()->begin(); it != pExistingNeighborTour->getPointerToTour()->end(); ++it){
////            std::cout << count << " " << (*it)->getKey() << " ROOT => " << (*it)->getRootNodeOfTour() << std::endl;
////            if (rootToBaseTour != (*it)->getRootNodeOfTour()){
////                exit(0);
////            }
////            count++;
////        }
//
//    } // end of creating subtour
//
//    return true;
//}




/**
 * Print a tour, diagnostic tool
 */
void EulerTour::printList(std::string text, std::list< Node * > * list){
    int tourLength = list->size();

    int count=1;
    std::cout << text << std::endl;
    std::cout << "      LIST START Length : " << tourLength << std::endl;
    std::cout << "              ROOT NODE : " << list->front()->getKey() << std::endl;
    for(std::list< Node * >::iterator iterator = list->begin(); iterator != list->end(); iterator++) {
        std::cout << "                   LIST : (" << count << ") " << (*iterator)->getKey() << " root -> " << (*iterator)->getRootNodeOfTour() << std::endl;
        count++;
    }
    std::cout << " _________________________ " << std::endl;
}


/**
 * passing a pointer to a list of Nodes
 *
 */
void EulerTour::rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    //std::cout << newRoot->getKey() << " subTour size : " << subTourToLoad->size() << std::endl;
    if (subTourToLoad->front()->getKey() != newRoot->getKey()){
        // printList("================> REROOT_AND_REASSIGN_SUBTOUR REROOT Before", subTourToLoad);
        //
        // subTourToLoad => 124313426797621   (root node is 1)
        //       newRoot => 6
        std::list< Node * >::iterator inTour = std::find_if(subTourToLoad->begin(), subTourToLoad->end(), find_Node_by_key(newRoot->getKey()));
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
//        printList("          => rerootAndReassignSubTour", subTourToLoad);
        // subTourToLoad => 679762124313426
        resetRootNodesInSubTourOfTour(subTourToLoad);
        // clear tempList
        tempList.clear();
    }
}

/**
 * reroot an existing tour that is in simpleTours map
 */
//void EulerTour::rerootActiveTour(Node * newRoot){
//
//    int oldNodeKey = newRoot->getRootNodeOfTour();
//    SimpleTour * pTempTour = &(simpleTours.find(oldNodeKey)->second);
//
//    if (pTempTour->getRootNodeOfTour() != newRoot->getKey()){ //
//        pTempTour->rerootTour(newRoot);
//        SimpleTour newTour = simpleTours.find(oldNodeKey)->second;
//        newTour.updateRootNode();
//
//        simpleTours.insert( std::pair< int, SimpleTour>(newTour.getRootNodeOfTour(), SimpleTour(*newTour.getPointerToTour())));
//        simpleTours.erase(oldNodeKey);
//        // if oldNodeKey and newRootNodeForTour are not the same, then insert will not insert, must delete first and create new entry
//    }
//}

/**
 *               9
 *               |
 *   1 - 2 - 6 - 7
 *   |   |
 *   3 - 4
 *   124313426797621
 *  smallest tour size is 1, then 3, there is no two
 *
 *  Only reroot the tour, there is no updating of the nodes root tour
 */
void EulerTour::rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    // 124313426797621 => reroot to 6
    std::list< Node * >::iterator inTour = std::find_if(subTourToLoad->begin(), subTourToLoad->end(), find_Node_by_key(newRoot->getKey()));
    // inTour points to 6 at ..4267...
    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    int sizeOfTour = subTourToLoad->size();
    if (subTourToLoad->front()->getKey() != newRoot->getKey() && (sizeOfTour > 2)){

        subTourToLoad->pop_front();
        // subTourToLoad => 24313426797621
        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour);
        // tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(newRoot);
        //      tempList => 24313426
        // transfer components from tempList to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
        //printList("rerooted ", subTourToLoad);
        // subTourToLoad => 6797621-24313426
        // clear tempList
        tempList.clear();
    }

//    else if (sizeOfTour == 1) {
//        std::cout << " TOUR SIZE is ONE !!!!!  " << std::endl;
//        std::cout << " TOUR SIZE is ONE !!!!!  " << subTourToLoad->front()->getKey() << std::endl;
//    }
}


int EulerTour::removeNode(int indexOfNode){

    Node * pNodeToRemove = &(nodes.find(indexOfNode)->second);
    // check neighbors
    std::list< Node * > * pTour;// = &tours[pNodeToRemove->getRootNodeOfTour()];  // tour with node to remove
    // break the euler tour and remove the pNodeToRemove
    if (pNodeToRemove->getTotalNeighbors() > 0){

//        std::cout << " NODE TO REMOVE " << pNodeToRemove->getKey() << " rooted => " << pNodeToRemove->getRootNodeOfTour() << std::endl;

            pTour = &tours.find(pNodeToRemove->getRootNodeOfTour())->second;    // get tour that contains node to remove

            // reroot tour to NodeToRemove and remove old tour
//            printList(" REMOVE BEFORE REROOT ", pTour);
//            std::cout << "(REMOVE NODE) KEY " << pNodeToRemove->getKey() << "( " << pNodeToRemove->getRootNodeOfTour() << " ) " << " TOTAL NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << " SIZE " << pTour->size() << std::endl;
            rerootAndReassignSubTour(pNodeToRemove, pTour);  // remove old tour and set as newly rooted tour

            //pTour = &tours[pNodeToRemove->getRootNodeOfTour()];    // get tour that contains node to remove
            pTour = &tours.find(pNodeToRemove->getRootNodeOfTour())->second;    // get tour that contains node to remove

//            printList("        AFTER  REROOT ", pTour);

            std::list< Node * >::iterator  pNext;
            Node * pNeighbor;
            int keyOfNodeToRemove = pNodeToRemove->getKey();

            while (pNodeToRemove->getTotalNeighbors() > 0) {
                //rerootSubTour(pNodeToRemove, pTour);
                //printList("NODETOREMOVE ", pTour);
                //std::cout << " => " << pNodeToRemove->getKey() << std::endl;
                // who is the neighbor
                std::list< Node * >::iterator pRightNeighbor = std::next(pTour->begin(), 1);
                std::list< Node * >::iterator pStartHere = std::next(pTour->begin(), 2);

//                std::cout <<  " right =>" << (*pRightNeighbor)->getKey() <<  " next => " << (*pStartHere)->getKey() << std::endl;

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
                }

                // potential pLeftNeighbor will reach end?
                pTour->pop_front();
                std::list< Node * > subTour;
                //std::cout << "     1b pTour popped " << pTour->size() << std::endl;
//                std::cout << "   1b pTour neighbor " << pTour->size() << " pNeighbor => " << (*pNeighbor).getKey() << std::endl;
                subTour.splice(subTour.begin(), *pTour, pTour->begin(), pNext); // upTo but not including pNext, pNext should point to node to remove
//                std::cout << "2 after splice pTour " << pTour->size() << " subTour " << subTour.size() << std::endl;
                //
                // if pNext is last node, then pTour will have only one element
                //
                // now I have two tours
                // 1. subTour
                // 2. pTour
                // check if I can sub back into remaining pTour
                std::set<int> nodesToCheck; // how could subTour be empty?
                for(std::list< Node * >::iterator sit = subTour.begin(); sit != subTour.end(); ++sit){
                    nodesToCheck.insert((*sit)->getKey());
                }

                // merge if possible, check each node in nodesToCheck
                std::list< Node * >::iterator locationOfCommonNode;
                for(std::set< int >::iterator sit = nodesToCheck.begin(); sit != nodesToCheck.end(); ++sit){
                    locationOfCommonNode = std::find_if(pTour->begin(), pTour->end(), find_Node_by_key(*sit));
                    if (locationOfCommonNode != pTour->end()){
                        // merge subTour with pTour
                        rerootSubTour(*locationOfCommonNode, &subTour);
                        // root pTour to common node
                        subTour.pop_back();
                        pTour->splice(locationOfCommonNode, subTour);
                        break;
                    }
                }

                // issue is if pTour size is 1, what happens to node subTour?
//                std::cout << "3 after reroot splice pTour " << pTour->size() << " subTour " << subTour.size() << std::endl;
                // make a new Tour from subTour if it couldn't be merged(linked)
                if (subTour.size() > 0){
                    int newRootNodeForTour = subTour.front()->getKey();
//                    std::cout << " -----> REMOVE SUBTOUR SIZE  " << subTour.size() << " " << newRootNodeForTour <<  std::endl;
                    //tours.insert ( std::pair< int, std::list<Node *> >(newRootNodeForTour, subTour) );  // key is the root of the tour
                    if (tours.find(newRootNodeForTour) == tours.end()){
                        tours.insert ( std::pair< int, std::list<Node *> >(newRootNodeForTour, std::list< Node *>()) );  // key is the root of the tour
                    } else {
                        std::cout << "pExisting Tour " << std::endl;
                        tours.find(newRootNodeForTour)->second.clear();
                        exit(0);
                    }

                    // get the pointer to the new Node List
                    std::list< Node *> * pNewRootNodeTour = &(tours.find(newRootNodeForTour)->second); // std::map<int, std::list< Node *> > tours;
                    // empty tour at this point

                    for(std::list< Node * >::iterator it = subTour.begin(); it != subTour.end(); ++it) {
                        auto pNodeTemp = &nodes.find((*it)->getKey())->second;
                        pNewRootNodeTour->push_back(pNodeTemp);
//                        pNodeTemp->setPointerToTour(&(tours.find(newRootNodeForTour)->second)); //update node
                        pNodeTemp->setRootNodeOfTour(newRootNodeForTour);

//                        pNewRootNodeTour->push_back(*it);
//                        nodes.find((*it)->getKey())->second.setPointerToTour(&(tours.find(newRootNodeForTour)->second)); //update node
                        // (*it)->setPointerToTour(&(tours.find(newRootNodeForTour)->second));
                    }
//                    printList(" AFETR REMOVE ", pNewRootNodeTour);
                }
                // pTour should be rooted at NodeToRemove at this point
                // should stay that way until pNodeToRemove is gone
                //nodes.find(pNeighbor->getKey())->second
                pNodeToRemove->removeNeighbor(&nodes.find(pNeighbor->getKey())->second); // should be in balance after removing
                // if no more pNodeToRemove in pTour, can't reroot to node to remove
            }

            // remove pTour (should contain only a single node and be rooted to nodeToRemove)
           // resetRootNodesInSubTourOfTour(pTour); // pTour should always contain the node to remove as a single tour
//                if (pTour->size() <= 1){ // move subTour into pTour, neighbor count should be 1
//                    std::cout << " SIZE of pTour => " << pTour->size() << " NODETOREMOVE => " << pNodeToRemove->getKey() << " NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << std::endl;
//                    printList("pTour size 1 : ROOT " + std::to_string(pNodeToRemove->getRootNodeOfTour()), pTour);
//                }


        if (tours.find(pNodeToRemove->getRootNodeOfTour()) != tours.end()){
            tours.find(pNodeToRemove->getRootNodeOfTour())->second.clear(); // clear List
            tours.erase(pNodeToRemove->getRootNodeOfTour()); // remove from tour
        }

//        }
    } else { // no neighbors
        //
        // tours is a map => std::map<int, std::list< Node *> > tours; // key is the root of the tour
        // list must be cleared first
        if (tours.find(pNodeToRemove->getRootNodeOfTour()) != tours.end()){
            pTour = &tours.find(pNodeToRemove->getRootNodeOfTour())->second; // list
            pTour->erase(pTour->begin(), pTour->end()); // erase invalidates existing iterators
            pTour->clear();
            tours.erase(pNodeToRemove->getRootNodeOfTour());
        }
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
    } else {
        std::cout << " NODE NEIGHBOR LIST NOT EMPTY PROBLEM WITH CODE?" << std::endl;
        exit(0);
    }

    //printTourSizes();
    //validateNodes("FROM REMOVE NODE AT END : " +  std::to_string(indexOfNode));

    totalComponents = tours.size();
    return totalComponents;
}



//int EulerTour::removeNodeFromSimpleTour(int indexOfNode){
//
//    Node * pNodeToRemove = &(nodes.find(indexOfNode)->second);
//
////if (nodes.find(indexOfNode) == nodes.end()){
////    std::cout << "LOOKING FOR A NODE NOT FOUND " << indexOfNode << std::endl;
////    exit(0);
////}
//    // check neighbors
//    if (pNodeToRemove->getTotalNeighbors() > 0){
//
////        if (pNodeToRemove->getTotalNeighbors() == 1){ // root to neighbor
////
////            SimpleTour * pTour = &(simpleTours.find(pNodeToRemove->getRootNodeOfTour())->second);
////            //pTour = &(tours.find(pNodeToRemove->getRootNodeOfTour())->second);
////
////            //pTour = &tours[pNodeToRemove->getRootNodeOfTour()];    // get tour that contains node to remove
////            //std::cout << "(REMOVE SINGLE NODE) KEY " << pNodeToRemove->getKey() << "( " << pNodeToRemove->getRootNodeOfTour() << " ) " << " TOTAL NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << " SIZE " << pTour->size() << std::endl;
////            // pointer to this tour
////            // reroot it
////            // tour still points to it after popping ends
////            rerootSubTour(pNodeToRemove, pTour);
////
////            pTour->pop_front();
////            pTour->pop_back();
////
////            pTour->front()->removeNeighbor(pNodeToRemove);
////
////            resetRootNodesInSubTourOfTour(pTour); // reassigns tour and should clear old map entry
////
////        } else {
//
//            rerootActiveTour(pNodeToRemove);
//            SimpleTour * pSimpleTour = &(simpleTours.find(pNodeToRemove->getRootNodeOfTour())->second);
//
//            std::list< Node * > * pTour = pSimpleTour->getPointerToTour();
//
//            std::list< Node * >::iterator pNext;
//            Node * pNeighbor;
//            int keyOfNodeToRemove = pNodeToRemove->getKey();
//
//            while (pNodeToRemove->getTotalNeighbors() > 0) {
//                // who is the neighbor
//                std::list< Node * >::iterator pRightNeighbor = std::next(pTour->begin(), 1);
//                std::list< Node * >::iterator pStartHere = std::next(pTour->begin(), 2);
//
//                pNeighbor = *pRightNeighbor; // get pointer to Neighbor
//                // find next neighbor pair
//                if ((*pStartHere)->getKey() == keyOfNodeToRemove){ // XYX where X is node to remove and Y is neighbor
//                    // check if Y is in pTour, if it is delete first two nodes from Tour
//                    pNext = pStartHere;
//                    // if not, single node becomes a new tour
//                } else {
//
//                    std::list< Node * >::iterator pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor); // look for occurrences of pNeighbor
//                    while(pLeftNeighbor != pTour->end()){ // find neighbor and check to see if next to it is NodeToRemove
//
//                        pNext = std::next(pLeftNeighbor, 1); // could be .end() ?
//                        if ((*pNext)->getKey() == pNodeToRemove->getKey()) {
//                            break;
//                        }
//                        // find next
//                        pStartHere = std::next(pLeftNeighbor, 1);
//                        pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor);
//                    }
//                }
//
//                // 12343545321
//                // remove node 3
//                // 34354532123
//                // potential pLeftNeighbor will reach end?
//                pTour->pop_front();
//                std::list< Node * > subTour;
//                subTour.splice(subTour.begin(), *pTour, pTour->begin(), pNext); // upTo but not including pNext, pNext should point to node to remove
//                //
//                // if pNext is last node, then pTour will have only one element
//                //
//                // now I have two tours
//                // check if I can sub back into remaining pTour
//                std::set<int> nodesToCheck; // how could subTour be empty?
//                for(std::list< Node * >::iterator sit = subTour.begin(); sit != subTour.end(); ++sit){
//                    nodesToCheck.insert((*sit)->getKey());
//                }
//
//                // merge if possible, check each node in nodesToCheck
//                // see if nodes in subTour can map back to parental pTour
//                std::list< Node * >::iterator locationOfCommonNode;
//                for(std::set< int >::iterator sit = nodesToCheck.begin(); sit != nodesToCheck.end(); ++sit){
//                    locationOfCommonNode = std::find_if(pTour->begin(), pTour->end(), find_by_key(*sit));
//                    if (locationOfCommonNode != pTour->end() ){
//                        if (subTour.size() == 1){
//                            subTour.pop_back();
//                        } else {
//                            // merge subTour with pTour
//                            rerootSubTour(*locationOfCommonNode, &subTour);
//                            // root pTour to common node
//                            subTour.pop_back();
//                            pTour->splice(locationOfCommonNode, subTour);
//                        }
//                        break;
//                    }
//                }
//
//                // make a new Tour from subTour if it couldn't be merged(linked)
//                if (subTour.size() > 0){
//                    int newRootNodeForTour = subTour.front()->getKey();
//                    simpleTours.insert( std::pair< int, SimpleTour>(newRootNodeForTour, SimpleTour(subTour)));
//                }
////                if (pTour->size() <= 1){ // move subTour into pTour, neighbor count should be 1
////                    std::cout << " SIZE of pTour => " << pTour->size() << " NODETOREMOVE => " << pNodeToRemove->getKey() << " NEIGHBORS " << pNodeToRemove->getTotalNeighbors() << std::endl;
////                    printList("pTour size 1 : ROOT " + std::to_string(pNodeToRemove->getRootNodeOfTour()), pTour);
////                }
//
//                // pTour should be rooted at NodeToRemove at this point
//                // should stay that way until pNodeToRemove is gone
//                pNodeToRemove->removeNeighbor(pNeighbor); // should be in balance after removing
//                // if no more pNodeToRemove in pTour, can't reroot to node to remove
//            }
//
////            std::cout << "     NODE TO REMOVE : " << pNodeToRemove->getKey() << " " << pNodeToRemove->getTotalNeighbors() << std::endl;
//
//            // remove pTour (should contain only a single node and be rooted to nodeToRemove)
////            if (pTour->size() > 1){
////                printList(" pTOUR size Greater than 1", pTour);
////                exit(0);
////            }
//
//            //resetRootNodesInSubTourOfTour(pTour); // pTour should always contain the node to remove as a single tour
//            simpleTours.erase(pNodeToRemove->getRootNodeOfTour());
////        }
//
//    } else { // no neighbors
//        //
//        // tours is a map => std::map<int, std::list< Node *> > tours; // key is the root of the tour
//        // list must be cleared first
////        SimpleTour * pTour = &(simpleTours.find(pNodeToRemove->getRootNodeOfTour())->second);
////        pTour = &tours[pNodeToRemove->getRootNodeOfTour()];
////        pTour->erase(pTour->begin(), pTour->end()); // erase invalidates existing iterators
////        pTour->clear();
//        simpleTours.erase(pNodeToRemove->getRootNodeOfTour());
//    }
//
//    // remove node
//    // remove the node from nodes
//    if (pNodeToRemove->getTotalNeighbors() == 0){
////        std::cout << " REMOVING NODE : " << pNodeToRemove->getKey() << " ROOTNODE :" << pNodeToRemove->getRootNodeOfTour() << std::endl;
////        if (nodes.find(pNodeToRemove->getKey()) == nodes.end()){
////            std::cout << "ERROR NOT FOUND in NODES LIST!  CANT REMOVE " << pNodeToRemove->getKey() << std::endl;
////            exit(0);
////        }
//        // nodes is a map => std::map<int, Node> nodes;
//        nodes.erase(pNodeToRemove->getKey());
//    } else {
//        std::cout << " NODE NEIGHBOR LIST NOT EMPTY PROBLEM WITH CODE?" << std::endl;
//        exit(0);
//    }
////    validateNodes("AFTER REMOVE NODES");
//    totalComponents = simpleTours.size();
//    return totalComponents;
//}





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

/**
 * Check that the node and its neighbors are present in node list
 * @param str
 * @return
 */
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

//        if ((*it).second.getPointerToTour()->size() == 0){
//            std::cout << (*it).first << " ZERO SIZE TOUR " << text << std::endl;
//            return true;
//        }

       ++count;
    } // use shared pointer? make instance on heap

    return false;
}


bool EulerTour::validateList(std::string text){

    std::cout << " VALIDATING LISTS IN TOURS :::::: " << text << std::endl;
    int count = 1;
    for(std::map<int, std::list< Node *> >::iterator it = tours.begin(); it != tours.end(); it++){

        int tourIndex = it->first; // tour index should be an existing node
        std::list<Node *> * pNodelist = &it->second;

        // check if the root of tour is in the nodes list
        if (nodes.find(tourIndex) == nodes.end()){
            std::cout << " @ " << text << " Node Does Not exist => " << tourIndex << " tour size : " << pNodelist->size()<< std::endl;
            return true;
        }

        // check if nodelist size violates limit, can only be 1 or greater than 2 such as a-b-a
        if (pNodelist->size() == 2){
            std::cout << " @ " << text << " LIST SIZE == 2 ************************************************* tourindex => " << tourIndex << std::endl;
            return true;
        }

        // check if nodelist contains a legitimate tour (every edge contains reciprocal pair)
        if (pNodelist->size() > 2 && !validateTour(pNodelist)){
            std::cout << " @ " << text << " INVALID LIST  " << tourIndex << std::endl;
            printList("INVALID LIST or TOUR ", pNodelist);
            return true;
        }

        // checks that first and last nodes of tour are the same
        if (*(pNodelist->begin()) != pNodelist->back()){
            printList(" DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR ", pNodelist);
            std::cout << " DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR " << std::endl;
            return false;
        }


        // go through each element of list and validate it is a node
        for(std::list<Node *>::iterator nodeIt=pNodelist->begin(); nodeIt != pNodelist->end(); nodeIt++){
            int key = (*nodeIt)->getKey();
            if (nodes.find(key) == nodes.end()){ // if not found we have an extra node
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>        tour MAP IND : " << tourIndex << " of " << tours.size() << std::endl;
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>       NODELIST SIZE : " << pNodelist->size() << std::endl;
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>                KEY  : " << key << std::endl;
                std::cout << "  =>>> NONSENSE KEY  : " << key << std::endl;
                printList("MESSED UP TOUR ", pNodelist);
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

    int oldNodeKey = subTour->front()->getRootNodeOfTour();
    int newRootNodeForTour = subTour->front()->getKey();

    if (oldNodeKey == newRootNodeForTour){
        for(std::list< Node * >::iterator it = subTour->begin(); it != subTour->end(); ++it) {
            auto pTempNode = &nodes.find((*it)->getKey())->second;
            pTempNode->setRootNodeOfTour(newRootNodeForTour);
        }
    } else {

        auto pTempTour = tours.find(newRootNodeForTour);
        if (pTempTour == tours.end()) {
            tours.insert(std::pair<int, std::list<Node *> >(newRootNodeForTour, std::list<Node *>()));  // key is the root of the tour
        } else {
            std::cout  << "resetRootNodesInSubTourOfTour SUBTOUR " << std::endl;
            exit(0);
        }

        // get the pointer to the new Node List
        std::list< Node *> * pNewRootNodeTour = &(tours.find(newRootNodeForTour)->second); // std::map<int, std::list< Node *> > tours;

        // empty tour at this point
        for(std::list< Node * >::iterator it = subTour->begin(); it != subTour->end(); ++it) {
            auto pTempNode = &nodes.find((*it)->getKey())->second;
            pTempNode->setRootNodeOfTour(newRootNodeForTour);
            pNewRootNodeTour->push_back(pTempNode);
        }

        // if oldNodeKey and newRootNodeForTour are not the same, then insert will not insert, must delete first and create new entry
        auto pOldKeyTour = tours.find(oldNodeKey);
        if (oldNodeKey != newRootNodeForTour && pOldKeyTour != tours.end()){
            pOldKeyTour->second.erase(pOldKeyTour->second.begin(), pOldKeyTour->second.end());
            pOldKeyTour->second.clear(); // clear the list
            tours.erase(oldNodeKey); // remove tour from MAP based on key.
        }
    }


//    auto pTempTour = tours.find(newRootNodeForTour);
//    if (pTempTour == tours.end()) {
//        tours.insert(std::pair<int, std::list<Node *> >(newRootNodeForTour, std::list<Node *>()));  // key is the root of the tour
//    }
//    else {
//std::cout  << "resetRootNodesInSubTourOfTour SUBTOUR " << std::endl;
//exit(0);
//       // pTempTour->second.clear();
//        //tours[newRootNodeForTour].clear();
//    }
//
//    // get the pointer to the new Node List
//    std::list< Node *> * pNewRootNodeTour = &(tours.find(newRootNodeForTour)->second); // std::map<int, std::list< Node *> > tours;
//
//    // empty tour at this point
//    for(std::list< Node * >::iterator it = subTour->begin(); it != subTour->end(); ++it) {
//        auto pTempNode = &nodes.find((*it)->getKey())->second;
//        pTempNode->setRootNodeOfTour(newRootNodeForTour);
//        pNewRootNodeTour->push_back(pTempNode);
////        pNewRootNodeTour->push_back(&nodes.find((*it)->getKey())->second);
////        (*it)->setPointerToTour(pNewRootNodeTour);
////        nodes.find((*it)->getKey())->second.setPointerToTour(pNewRootNodeTour);
////        nodes.find((*it)->getKey())->second.setRootNodeOfTour(newRootNodeForTour);
//    }
//
//    //tours[newRootNodeForTour] = *subTour;
//
//    // if oldNodeKey and newRootNodeForTour are not the same, then insert will not insert, must delete first and create new entry
//    auto pOldKeyTour = tours.find(oldNodeKey);
//    if (oldNodeKey != newRootNodeForTour && pOldKeyTour != tours.end()){
//        pOldKeyTour->second.erase(pOldKeyTour->second.begin(), pOldKeyTour->second.end());
//        pOldKeyTour->second.clear(); // clear the list
//        tours.erase(oldNodeKey); // remove tour from MAP based on key.
//    }
}


/**
 *
 */
void EulerTour::resetAccessed(std::vector<int> * checkedNodes){
    int totalChecked = checkedNodes->size();
    for(int i=0; i<totalChecked; i++){
        nodes.find((*checkedNodes)[i])->second.setAccessed(false);
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

            int neighborIndex = pNode->getPointerToNeighborByIndex(index);

            subTourToLoad->push_back(&nodes.find(neighborIndex)->second);

//            subTourToLoad->push_back(&nodes.find(pNode->getPointerToNeighborByIndex(index)->getKey())->second);
            //subTourToLoad->push_back(pNode->getPointerToNeighborByIndex(index));
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
        nodeToInsert = *(beginIt + i);

        nodes.insert ( std::pair<int, Node>(nodeToInsert, Node(nodeToInsert) ) );
        pNode = &(nodes.find(nodeToInsert)->second);

        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        it = pModel->getPointerToNeighborhood(pNode->getKey());
        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            neighbor = *(it+j);
            // if not found, itIndex will report last
            if (neighbor > -1 && nodes.find(neighbor) != nodes.end()){
                pNode->addNeighbor(&(nodes.find(neighbor)->second)); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == -1) {
                break;
            }
        }

        //std::cout << "  " << i << " ADDING NODE -> " << nodeToInsert << std::endl;
        addToTour(nodeToInsert);
//        validateList(std::to_string(i) + " FROM EulerTour::createInitialTour ");
        //validateNodes(std::to_string(i) + " FROM EulerTour::createInitialTour ");
        // create subtour
    } // end of adding beads

    //std::cout << " INITIAL TOURS SIZE " << tours.size() << std::endl;
    //validateList("FINAL FROM EulerTour::createInitialTour ");
    //validateNodes("FINAL FROM EulerTour::createInitialTour ");
    //totalComponents = simpleTours.size();
    //this->printTourSizes();
    totalComponents = tours.size();
}



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


// check against nodes list, should match
void EulerTour::printTourSizes(){

    for(auto &entry : tours) {
        auto const &index = entry.first;
        auto &base = entry.second;
        std::cout << " front bead : " << index << " tour size : " << base.size() << std::endl;
    }
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

