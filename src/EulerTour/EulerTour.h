//
// Created by Robert Rambo on 11/01/2017.
//

#ifndef IKETAMA_EULERTOUR_H
#define IKETAMA_EULERTOUR_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "Node.h"
#include <ctime>

class Model;
class Node;

/**
 *   Class is used to monitor the graph connectivity
 *   EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
 *   currentNumberOfComponents = eulerTour.getNumberOfComponents();
 */
class EulerTour {

    typedef std::map<int,Node>::iterator it_type;

    std::map<int, Node> nodes;
    std::vector<int>::iterator * pSelectedLattice;

    std::map<int, Node> backedUpNodes;
    std::map<int, std::list< Node *> > backedUpTours;

    std::map<int, std::list< Node *> > tours; // key is the root of the tour

    void createInitialTour(int workingLimit, Model *pModel);
    bool addToTour(int nodeToAdd);

    void createSubTour(Node * pNode, std::list< Node * > * subTourToLoad);
    void rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);
    void rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);

    std::list< Node * >::iterator getFirstOccurrenceInTour(Node * pNode, std::list<Node *> * pTour);
    std::list< Node * >::reverse_iterator getLastOccurrenceInTour(Node * pNode, std::list<Node *> * pTour);

    void testSetOne();

    void printList(std::string text, std::list< Node * > * list);
    void resetAccessed(std::vector<int> *checkedNodes);
    void resetRootNodesInSubTour(std::list<Node *> * subTour);
    void resetRootNodesInSubTourOfTour(std::list<Node *> * subTour);
    bool validateList(std::string note);
    bool validateTour(std::list<Node *> * tourtocheck);
    bool validateNodes();

    int totalComponents;

public:

    // pass in iterator to beginning of selected lattice points that are sorted upto WorkingLimit
    // want to use this Class to determine connectivity of a graph
    // Connectivity is the numberOfComponents
    //
    EulerTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel);
    ~EulerTour(){
        // delete pSelectedLattice;
    }

    int addNode(int latticePoint, Model *pModel);
    int removeNode(int indexOfNode);
    int getNumberOfComponents(){return totalComponents;}
    void createBackUp();

};

#endif //IKETAMA_EULERTOUR_H
