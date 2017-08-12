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
#include "SimpleTour.h"
#include <ctime>

class Model;
//class Node;

/**
 *   Class is used to monitor the graph connectivity
 *   EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
 *   currentNumberOfComponents = eulerTour.getNumberOfComponents();
 */


struct find_by_key : std::unary_function< Node *, bool>{
    find_by_key(int keyToFind) : key(keyToFind){}
    bool operator () (Node * p) { return p->getKey() == key; }
private:
    int key;
};


class EulerTour {

    typedef std::map<int,Node>::iterator it_type;

    int totalComponents;
    std::map<int, Node> nodes; // use shared pointer? make instance on heap
    //std::map<int, std::shared_ptr<Node> > nodes;
    //std::vector<int>::iterator * pSelectedLattice;

    std::map<int, SimpleTour> simpleTours;
    std::map<int, std::list< Node *> > tours; // key is the root of the tour
    // std::map<int, Tours> => tours is an object that contains List and Set

    void createInitialTour(int workingLimit, Model *pModel, std::vector<int>::iterator beginIt);
    bool addToTour(int nodeToAdd);
    bool addToSimpleTour(int nodeToAdd);

    void createSubTour(Node * pNode, std::list< Node * > * subTourToLoad);
    void rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);
    void rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);

    std::list< Node * >::iterator getFirstOccurrenceInTour(Node * pNode, std::list<Node *> * pTour);
    std::list< Node * >::reverse_iterator getLastOccurrenceInTour(Node * pNode, std::list<Node *> * pTour);

    void testSetOne();
    void printList(std::string text, std::list< Node * > * list);
    void resetAccessed(std::vector<int> *checkedNodes);
    //void resetRootNodesInSubTour(std::list<Node *> * subTour);
    void resetRootNodesInSubTourOfTour(std::list<Node *> * subTour);
    bool validateNodesAndTours(std::string text);
    bool validateTour(std::list<Node *> * tourtocheck);
    void removeListFromTour(int key);

public:

    // pass in iterator to beginning of selected lattice points that are sorted upto WorkingLimit
    // want to use this Class to determine connectivity of a graph
    // Connectivity is the numberOfComponents
    //
    EulerTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel);
    EulerTour();
    ~EulerTour(){

        // clear pointers to nodes
        for (std::map<int, std::list< Node *>>::iterator it=tours.begin(); it!=tours.end(); ++it){
            // delete points in tour
//            for(std::list<Node *>::iterator lit = it->second.begin(); lit != it->second.end(); ++lit){
//                it->second.erase(lit);
//            }
            it->second.clear();
        }
        // remove all nodes
        std::map<int,Node>::iterator it = nodes.begin();

        while(it != nodes.end()){
            it = nodes.erase(it);
        }

        tours.clear();
        nodes.clear();
    }

    int addNode(int latticePoint, Model *pModel);
    int removeNode(int indexOfNode);
    int removeNodeFromSimpleTour(int indexOfNode);
    void rerootActiveTour(Node * newRoot);

    int removeNodeOLD(int indexOfNode);
    int getNumberOfComponents(){return totalComponents;}
    //void createBackUp();

    int newTour(std::vector<int>::iterator beginIt, int workingLimit, Model *pModel);
    bool validateNodes(std::string st);
    void checkTourSize(std::string note);
    bool validateList(std::string note);
    bool checkNodesList(std::set<int> * beads_in_use);
    static bool deleteAll( Node * theElement ) { delete theElement; return true; }
};

#endif //IKETAMA_EULERTOUR_H
