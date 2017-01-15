//
// Created by Robert Rambo on 13/01/2017.
//

#ifndef IKETAMA_ANNEAL_H
#define IKETAMA_ANNEAL_H
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <Phase.h>
#include <Bead.h>

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/pending/disjoint_sets.hpp>


#ifdef __cplusplus
extern "C" {
#endif
#include <qhull_a.h>
#include <libqhull.h>
#ifdef __cplusplus
}
#endif

class Data;
class Model;
class Objective;
class Bead;

class Anneal {

    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertices_size_type VertexIndex;
    typedef VertexIndex* Rank;
    typedef Vertex* Parent;

    float highT;
    float highTempStartForCooling;
    int numberOfCoolingTempSteps;
    float percentAddRemove;
    float lowTempStop, totalTempStop;
    float expansionSlope, stepFactor;
    int lowerV, upperV, contactsPerBead;
    float highTempExchangeCutoff = 0.20;
    float highTempAcceptance;
    float contactCutOff;
    float interconnectivityCutOff; /*!< distance between lattice points to be counted as a neighbor */

    float expSlowCoolConstant = 0.87;
    float eta; // weight for compactness
    float lambda; // weight for connectivity
    float beta;
    float mu;
    int highTempRounds;
    int stepsPerTemp;
    std::string filenameprefix;

    int maxbin, totalBins;
    int totalNumberOfPhasesForSeededModeling=1;

    std::random_device randomObject;     // only used once to initialise (seed) engine

    void addLatticPositionToModel(std::vector<int>::iterator * pBeginIt,
                                  std::vector<int>::iterator * pEndIt,
                                  std::vector<int> * pBackUpState,
                                  int * pWorkingLimit,
                                  std::vector<int>::iterator * pItIndex);

    void restoreAddingFromBackUp(std::vector<int>::iterator * pBeginIt,
                                 std::vector<int> * pBackUpState,
                                 int * pWorkingLimit,
                                 std::vector<int> * pBinCountBackUp,
                                 std::vector<int>::iterator * pBinCountBegin);

    void removeLatticePositionToModel(std::vector<int>::iterator * pBeginIt,
                                      std::vector<int> & bead_indices,
                                      std::vector<int> & pBinCount,
                                      int * const pBin,
                                      int * pWorkingLimit,
                                      int totalBeadsInSphere,
                                      const int * pLatticePointToRemove);

    void removeLatticePositionToModelSym(std::vector<int>::iterator * pBeginIt,
                                         std::vector<int> & bead_indices,
                                         std::vector<int> & pBinCount,
                                         int * pWorkingLimit,
                                         const int * pLatticePointToRemove, Model * pModel, Data *pData);

    void restoreRemovingLatticePointFromBackUp(std::vector<int>::iterator * pBeginIt,
                                               int * pWorkingLimit,
                                               std::vector<int> * pBinCountBackUp,
                                               std::vector<int>::iterator * pBinCountBegin);

    float calculateLocalContactPotentialPerBead(std::set<int> *beads_in_use, Model *pModel, const int selectedIndex);

    float calculateLocalContactPotentialOfNeighborhood(std::set<int> *beads_in_use,
                                                       Model *pModel,
                                                       int const selectedIndex);
    bool checkIfWithinContact(int potentialPosition, std::set<int> * beads_in_use, Model *pModel);

    int numberOfContactsFromSet(std::set<int> *beads_in_use, Model *pModel, int const selectedIndex);
    int numberOfContacts(int &beadIndex, std::vector<int> *bead_indices, int &workingLimit, float &contactCutOff, Model *pModel, float * pDistance);

    int numberOfContactsFromSetSeeded(std::set<int> *beads_in_use,
                                      const std::set<int> *true_model,
                                      Model *pModel,
                                      int const selectedIndex);

    float connectivityPotentialSeeded(std::vector<int> *activeIndices, std::set<int> trueModel, std::set<int> * fixedPoints, int availableWorkingLimit, float *pDistances,
                                      int totalBeads, int &numberOfComponents);

    bool checkSet(int limit, std::vector<int> beadIndices, std::set<int> setvalues);

    float symViolationsPotential(int specifiedIndex, std::vector<int> * beads_indices, int workingLimit, Model * pModel);

    void fillPrBinsAndAssignTotalBin(int * const pBin, float * pDistance, unsigned long int totalDistancesInSphere, Data *pData);

public:

    Anneal(float highT, float percent, int lowerV, int upperV, int highTempRounds, int contacts, std::string prefix, int totalSASteps, float stepFactor, float eta, float lambda);

    ~Anneal(){
    }


    void createInitialModelSymmetry(Model *pModel, Data *pData);

    float calculateAverageDistance(float * pDistance, int *stopAt, std::vector<int> *bead_indices, Model * pModel);


    void setHighTempExchangeCutoff(float percent){ highTempExchangeCutoff = percent;}

    inline float calculateKLEnergy(std::vector<int> *bead_indices, std::vector<int> *bins, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData);

    inline float calculateKLEnergySymmetry(std::vector<int> *bead_indices, std::vector<int> *binCount, int beadIndiciesWorkingLimit, int totalBeadsInSphere, int &violation, Model *pModel, Data *pData);

    float calculateCVXHULLVolume(char * flags, std::vector<int> *bead_indices, int upTo, double *points, Model *pModel);

    void removeFromPr(int removeMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalBeadsInUniverse, std::vector<int> & prBins);

    void addToPr(int addMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalDistancesInSphere, std::vector<int> & prBins);

    float cvxHullCutoff(int temperature, int initial, float hill);


    void readPDB(Model *pModel, std::vector<int> * keptBeads, std::string filename);
    void printStatus(std::string text, int index, int wl, float energy, float kl, float dl);

    void medianFiltered(float &volumeAverage, float &volumeStdev, std::vector<int> &averageV, int &volumeCount);

    void beadToPoint(pointT *ptestPoint, Bead *pBead);

    int recalculateDeadLimit(int workingLimit, std::vector<int> &bead_indices, Model * pModel, int totalBeads );

    void refineCVXHull(std::vector<int> &bead_indices, std::vector<int> &active_indices,
                       int totalBeadsInSphere, int workingLimit, int *deadLimit,  Model *pModel);


    void enlargeCVXHull(
            std::vector<int> &bead_indices,
            std::vector<int> &active_indices,
            int totalBeadsInSphere,
            int workingLimit,
            int *deadLimit,
            Model *pModel,
            float &cbrtK);

    float calculateEnergy(float &klValue, float &alpha, float &volume, int &numberOfBeads, float &beadVolume, float &invBeadVolume, float &cbrtK);


    void bayesianUpdateMeanVariance(float *pPriorMean, float *pInvPriorVariance, float *pVolumeAverage,
                                    float *pVolumeStdev,
                                    int * pCount);


    std::string refineSymModel(Model *pModel, Data *pData, int index);

    bool checkForRepeats(std::vector<int> beads);

    void enlargeDeadLimit(std::vector<int> &vertexIndices,
                          int totalV,
                          std::vector<int> &bead_indices,
                          int workingLimit,
                          int *deadLimit,
                          int totalBeadsInSphere,
                          Model *pModel);


    void removeFromPrSym(int removeMeSubUnitIndex, std::vector<int> &beadsInUse, int workingLimit, std::vector<int> &prBins,
                         Model *pModel, Data *pData);

    void addToPrSym(int addMeSubUnitIndex, std::vector<int> &beadsInUse, int workingLimit, std::vector<int> &prBins,
                    Model *pModel, Data *pData);

    float calculateKLEnergySymmetryFast(std::vector<int> *subUnit_indices, std::vector<int> *binCount,
                                        int indicesWorkingLimit, int totalBeadsInSphere, int &violation, Model *pModel,
                                        Data *pData, float &klDivergence);

    void fillStepsArray(int * array, int stopLimit);

    int getContactsPerBead(){return contactsPerBead;}
    int getMaxBinUsedInAnnealing(){ return maxbin;}
    int gettotalBinsDerivedFromData(){ return totalBins;}
    int getHighTempRounds(){return highTempRounds;}
    float getHighTempStartForCooling(){return highTempStartForCooling;}
    float getLowTempStop(){return lowTempStop;}
    int getNumberOfCoolingSteps(){return numberOfCoolingTempSteps;}
    float getPercentAddRemove(){ return percentAddRemove;}
    float getLambda(){ return lambda;}
    float getMu(){ return mu;}
    float getEta(){ return eta;}
    float getBeta(){ return beta;}
    float getStepFactor(){ return stepFactor;}

    void setInterconnectivityCutoff(float value){ this->interconnectivityCutOff = value;}

    void refineMultiphase(Model *pModel, Objective *pObjective, std::vector<Phase *> * pPhases);

    float calculateKLEnergyMP(std::vector<int>::iterator *beginBeadIndices, std::vector<int> *pBinCount,
                              int totalBeadsInSphere, std::vector<int> *pDistancesConvertedToBinPerDataset, int totalDistancesInSphere, Data *pData, Model *pModel);

    int findPhases(std::vector<Phase *> *pPhases, int beadInUse);

    float calculateCVXHULLVolumeMP(char *flags, std::vector<int> *bead_indices, int startAt, int upTo, double *points,
                                   Model *pModel);

    float calculateEnergyMP(std::vector<float> *divKLs, float alpha, std::vector<float> *cvxHullVolumes,
                            std::vector<Phase *> *pPhases, float beadVolume, float invBeadVolume);

    int numberOfContactsMP(std::vector<int> *bead_indices, std::vector<Phase *> *pPhases, int phase1, int index1,
                           int phase2, int index2, float &contactCutOff, Model *pModel, float *pDistance);

    std::vector<int> createRandomConnectedGraph(std::vector<int> &active_indices, int startIndex, int upToIndex,
                                                int &totalAvailableIndices, Model *pModel);

    bool isConnectedComponent(std::vector<int> *activeIndices, int availableWorkingLimit, float *pDistances,
                              int totalDistances, int &numberOfComponents);

    void createRandomConnectedGraphSinglePhase(std::vector<int> &active_indices, int upToIndex, Model *pModel,
                                               const int dmax);


    std::string refineHomogenousBodyASACVX(Model *pModel, Data *pData, int iteration);
    std::string reAssignLatticeModel(std::string PDBFilename, Model *pModel, Data *pData);
    std::string refineHomogoenousBodyFromPDBSeed(Model *pModel, Data *pData, int iteration);

    void updateASATemp(int index, float evalMax, float acceptRate, float &temp, float &inv_temp);

    int createShortestConnectedPath(std::vector<int> &active_indices,  Model *pModel, const int dmax);


    void initialHighTempSearchSeeded(Model *pModel, Data *pData, std::string name);

    void populateLayeredDeadlimit(std::vector<int>::iterator iteratorBeadIndices, int workingLimit, int *pDeadLimit, Model *pModel, int totalBeads);

    void populateLayeredDeadlimitUsingSet(std::vector<int> * beadIndices, std::set<int> *beads_in_use, int workingLimit,
                                          int *pDeadLimit, Model *pModel);


    void adjustDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                Model *pModel, int totalWorkingBeads, const int selectedIndex);

    bool checkDKL(float currentKL, std::vector<int> *bead_indices, std::vector<int> *binCount, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData);

    void createSeedFromPDB(Model *pModel, Data *pData, std::string name, std::string PDBFilename, int totalPhases);
    std::string refineInputPDB(Model *pModel, Data *pData, std::string name, std::string PDBFilename, float temperature);

    float calculateKLDivergenceAgainstPDBPR(std::vector<int> &modelPR, std::vector<float> &targetPR);

    void populateLayeredDeadlimitSeed(std::vector<int>::iterator partialSetIndices, int workingLimit,
                                      std::vector<int> *trueModelIndices, int totalBeadsInTrue, int *deadLimit,
                                      Model *pModel, int totalWorkingSet);

    void adjustDeadLimitPerBeadSeed(std::vector<int>::iterator partialSetIndices, int workingLimit,
                                    std::vector<int> *trueModelIndices, int trueModelLimit, int *pDeadLimit, Model *pModel,
                                    int totalWorkingSet);

    void removeNeighborsFromDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit,
                                             int *pDeadLimit, Model *pModel, int totalWorkingBeads,
                                             const int selectedIndex);


    float calculateLocalContactPotentialPerBeadSeeded(std::set<int> *beads_in_use, const std::set<int> *trueModel,
                                                      Model *pModel,
                                                      int const selectedIndex);


    float calculateLocalContactPotentialOfNeighborhoodSeeded(std::set<int> *beads_in_use,
                                                             const std::set<int> *trueModel,
                                                             Model *pModel,
                                                             int const selectedIndex);

    float contactsPotential(float value);


    int changeToDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                 Model *pModel, int totalWorkingBeads, const int selectedInde);

    int additionToDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                   Model *pModel, int totalWorkingBeads, const int selectedIndex);


    float calculateTotalContactEnergy(std::vector<int> *bead_indices, const int workingLimit, Model *pModel, float *pDistance);

    float addToTotalContactEnergy(const int addMe, std::vector<int> *bead_indices, const int workingLimit, Model *pModel,
                                  int totalWorkingBeads, float *pDistance);

    void createInitialModelCVXHull(Model *pModel, Data *pData, std::string name);

    void calculateAverageNumberOfContacts(float *averageContacts, std::vector<int> *bead_indices, const int workingLimit,
                                          Model *pModel, float *pDistance);

    float connectivityPotential(int numberOfComponents);


    void populateNonSeedVector(std::vector<int> *nonSeed, int *totalNonSeed, std::vector<int> *activeIndices,
                               int availableWorkingLimit, std::set<int> *trueModel);

    void populateLayeredDeadlimitUsingFixedPoints(std::vector<int> *bead_indices, std::set<int> *beads_in_use, const std::set<int> *true_model,
                                                  std::set<int> *fixedPoints, const int workingLimit, int *pDeadLimit,
                                                  Model *pModel);

    bool isAnchorSafe(int swapPt, std::set<int> *beadsInUseTree, std::set<int> *reducedTrueModelTree, std::set<int> *fixedPoints,
                      int workingLimit, Model *pModel);
};

#include "Model.h"


/**
 *
 */
inline void Anneal::addLatticPositionToModel(std::vector<int>::iterator * pBeginIt,
                                             std::vector<int>::iterator * pEndIt,
                                             std::vector<int> * pBackUpState,
                                             int * pWorkingLimit,
                                             std::vector<int>::iterator * pItIndex){

    std::copy(*pBeginIt, *pEndIt, pBackUpState->begin());
    // make the swap at the border (workingLimit)
    std::iter_swap(*pBeginIt + *pWorkingLimit, *pItIndex);
    // increment workingLimit to include new position
    *pWorkingLimit += 1;
    std::sort(*pBeginIt, *pBeginIt + *pWorkingLimit);
}


/**
 * addMe must be present in beadsInUse
 * pBin is a pointer to the vector of binned distances in Model
 * prBin is the vector holding counts per bin
 * upperLimit is usually the working Limit of the set
 */
inline void Anneal::addToPr(int addMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalBeadsInUniverse, std::vector<int> & prBins ){

    // beadsInUse must be sorted
    int * row;
    unsigned long int row2;

    int i=0;
    // add down column
    while (beadsInUse[i] < addMe && i < upperLimit){
        row = &beadsInUse[i];
        row2 = (*row)*totalBeadsInUniverse - ((*row)*(*row+1)*0.5) - (*row) - 1;
        prBins[ *(pBin + row2 + addMe) ]++;
        i++;
    }

    // out of first loop, i = addMe
    i++;
    // Add across row
    row2 = addMe*totalBeadsInUniverse - addMe*(addMe+1)*0.5 - addMe - 1;
    while (i < upperLimit){
        prBins[ *(pBin + row2 + beadsInUse[i]) ]++;
        i++;
    }
}


inline void Anneal::beadToPoint(pointT *ptestPoint, Bead *pBead) {
    ptestPoint[0] = pBead->getX();
    ptestPoint[1] = pBead->getY();
    ptestPoint[2] = pBead->getZ();
}

/**
 * recalcualte P(r) distribution then compare against dataset for KL divergence (bead indices must be sorted)
 *
 */
inline float Anneal::calculateKLEnergy(std::vector<int> *bead_indices, std::vector<int> *binCount, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    fill(binCount->begin(), binCount->end(), 0.0);

    int row, row2;
    int * modelBin = pModel->getBins();

    // calculate P(r) for beads
    for(int m=0; m < upTo; m++){
        row = (*bead_indices)[m];
        row2 = row*totalBeadsInSphere - row*(row+1)*0.5 - row - 1;

        for(int n=(m+1); n < upTo; n++){
            //col = (*bead_indices)[n];
            (*binCount)[ *(modelBin + row2 + (*bead_indices)[n]) ]++;
        }
    }

    // calculate KLDivergence
    return pData->calculateKLDivergence(*binCount);
}

/**
 * Given the selectedIndex, how many of its neighbors are in use
 *
 * beads_in_use must not include the newly added lattice position given by selectedIndex
 */
inline float Anneal::calculateLocalContactPotentialOfNeighborhood(std::set<int> *beads_in_use,
                                                                  Model *pModel,
                                                                  int const selectedIndex) {

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    float neighborhoodEnergy = 0;

    for (int i=0; i< totalNeighbors; i++){

        int neighbor = *(it+i);

        if (beads_in_use->find(neighbor) != endOfSet){

            neighborhoodEnergy += calculateLocalContactPotentialPerBead(beads_in_use, pModel, neighbor);

        } else if (neighbor == -1) {
            break;
        }
    }

    return neighborhoodEnergy;
}

/**
 * for selected bead, determine how many neighbors it has in selected set of beads
 * beads_in_use -> pointer to selected set of beads
 * pModel -> pointer to model that contains, per lattice point, its respective neighborhood
 * selectedIndex -> bead to consider
 */
inline float Anneal::calculateLocalContactPotentialPerBead(std::set<int> *beads_in_use,
                                                           Model *pModel,
                                                           int const selectedIndex){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    float neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    for (int i=0; i< totalNeighbors; i++){

        int neighbor = *(it+i);

        if (beads_in_use->find(neighbor) != endOfSet){
            neighborContacts += 1;
        } else if (neighbor == -1) {
            break;
        }
    }

    return contactsPotential(neighborContacts);
}



inline float Anneal::contactsPotential(float neighborContacts){
    float value = 0;

    if (neighborContacts < (contactsPerBead - 4) ) { // 1
        value = 1000;
    } else if (  neighborContacts == (contactsPerBead - 4) ){ // 2
        value = 100;
    } else if ( (neighborContacts > (contactsPerBead - 4) ) && neighborContacts <= (contactsPerBead - 3) ){ // 3
        value = 0.1;
    } else if ( (neighborContacts > (contactsPerBead - 3) ) && neighborContacts <= (contactsPerBead - 2) ){ // 4
        value = 0.05;
    } else if ( (neighborContacts > (contactsPerBead - 2) ) && neighborContacts <= (contactsPerBead - 1) ){ // 5
        value = 0.001;
    } else if ( (neighborContacts > (contactsPerBead - 1) ) && neighborContacts < (contactsPerBead + 1) ){ // 0
        value = 0;
    } else if ((neighborContacts <= (contactsPerBead + 1) ) && neighborContacts < contactsPerBead + 2 ){ // 7
        value = 0.001;
    } else if (neighborContacts >= contactsPerBead + 2 ){ // 8
        value = 0.1;
    }

    return value;
}



inline void Anneal::fillPrBinsAndAssignTotalBin(int * const pBin, float * pDistance, unsigned long int totalDistancesInSphere, Data * pData){
    maxbin=0;
    for(int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    totalBins = pData->getShannonBins(); // set by experimental input Data file
    // maxbin and totalBins will be differenct
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
}


inline int Anneal::numberOfContactsFromSet(std::set<int> *beads_in_use,
                                           Model *pModel,
                                           int const selectedIndex){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    int neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    for (int i=0; i< totalNeighbors; i++){

        int neighbor = *(it+i);

        if (beads_in_use->find(neighbor) != endOfSet){
            neighborContacts += 1;
        } else if (neighbor == -1) {
            break;
        }
    }

    return neighborContacts;
}

/**
 *
 */
inline void Anneal::removeLatticePositionToModel(std::vector<int>::iterator * pBeginIt,
                                                 std::vector<int> & bead_indices,
                                                 std::vector<int> & pBinCount,
                                                 int * const pBin,
                                                 int * pWorkingLimit,
                                                 int totalBeadsInSphere,
                                                 const int * pLatticePointToRemove){

    std::vector<int>::iterator itIndex = find(*pBeginIt, *pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    removeFromPr(*pLatticePointToRemove, bead_indices, *pWorkingLimit, pBin, totalBeadsInSphere, pBinCount);
    // reduce the workingLimit
    // if wl = 10
    // 0 1 2 3 4 5 6 7 8 9 10
    // remove 4
    // 0 1 2 3 9 5 6 7 8 4 10
    //
    // 0 1 2 3 5 6 7 8 9 4 10
    //
    *pWorkingLimit -= 1;
    // swap selected point to the workingLimit
    std::iter_swap(itIndex, *pBeginIt + *pWorkingLimit);
    // still need to sort, swap changes the order
    // to save some time, sort should only be from swapped point to workingLimit
    std::sort(*pBeginIt, *pBeginIt + *pWorkingLimit);
}

/**
 * requires a sorted beadsInuse vector
 * removeMe is the index of the bead from origination
 */
inline void Anneal::removeFromPr(int removeMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalBeadsInUniverse, std::vector<int> & prBins ){

    // beads in Use is sorted
    int * row;
    unsigned long int row2;

    int i=0;
    // remove down column
    while (beadsInUse[i] < removeMe && i < upperLimit){
        row = &beadsInUse[i];
        row2 = (*row)*totalBeadsInUniverse - ((*row)*(*row+1)*0.5) - (*row) - 1;
        prBins[ *(pBin + row2 + removeMe) ]--;
        i++;
    }
    // when loop ends beadsInUse[i] => removeMe
    i++;
    // remove across row
    row2 = removeMe*totalBeadsInUniverse - removeMe*(removeMe+1)*0.5 - removeMe - 1;
    while (i < upperLimit){
        prBins[ *(pBin + row2 + beadsInUse[i]) ]--;
        i++;
    }

}

inline void Anneal::restoreAddingFromBackUp(std::vector<int>::iterator * pBeginIt,
                                            std::vector<int> * pBackUpState,
                                            int * pWorkingLimit,
                                            std::vector<int> * pBinCountBackUp,
                                            std::vector<int>::iterator * pBinCountBegin){

    std::copy(pBackUpState->begin(), pBackUpState->end(), *pBeginIt);
    *pWorkingLimit -= 1;
    std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
}

inline void Anneal::restoreRemovingLatticePointFromBackUp(std::vector<int>::iterator * pBeginIt,
                                                          int * pWorkingLimit,
                                                          std::vector<int> * pBinCountBackUp,
                                                          std::vector<int>::iterator * pBinCountBegin){

    *pWorkingLimit += 1;
    // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
    // sorting is n*log(n) with n = 200?  should be much smaller
    std::sort(*pBeginIt, *pBeginIt+*pWorkingLimit); // swapped index is at workingLimit
    std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
}

#endif //IKETAMA_ANNEAL_H
