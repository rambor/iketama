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
#include <Component.h>
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/pending/disjoint_sets.hpp>

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
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

    std::vector<double> connectivityPotentialTable;

    float highT;
    float highTempStartForCooling;

    float percentAddRemove;
    //float lowTempStop, totalTempStop;
    float totalTempStop;
    float expansionSlope;
    int lowerV, upperV;
    float highTempExchangeCutoff = 0.20;
    float alpha = 0.03;
    float highTempAcceptance;
    float contactCutOff, contactsPerBead;
    float interconnectivityCutOff; /*!< distance between lattice points to be counted as a neighbor */

    float expSlowCoolConstant = 0.87;
    float eta; // weight for compactness
    float lambda; // weight for connectivity
    float beta;
    float mu;
    int highTempRounds;
    int stepsPerTemp, ccmultiple;
    std::string filenameprefix;

    int maxbin, totalBins;
    int totalComponents=0;
    int totalNumberOfPhasesForSeededModeling=1;
    std::vector<Component> components;

    std::random_device randomObject;     // only used once to initialise (seed) engine

    void addLatticPositionToModel(std::vector<int> * pIndices,
                                  std::vector<int> * pBackUpState,
                                  int * pWorkingLimit,
                                  std::vector<int>::iterator * pItIndex);

    float calculateLocalContactPotentialPerBead(std::set<int> *beads_in_use, Model *pModel, const int selectedIndex);

    float calculateLocalContactPotentialOfNeighborhood(std::set<int> *beads_in_use,
                                                       Model *pModel,
                                                       int const selectedIndex);
    bool checkIfWithinContact(int potentialPosition, std::set<int> * beads_in_use, Model *pModel);


    bool checkSet(int limit, std::vector<int> beadIndices, std::set<int> setvalues);

    float connectivityPotentialSeeded(std::vector<int> *activeIndices,
                                      std::set<int> trueModel,
                                      std::set<int> * fixedPoints,
                                      int availableWorkingLimit,
                                      float *pDistances,
                                      int totalBeads,
                                      int &numberOfComponents);

    void createPlacesToCheck(int workingLimit,
                                      int average_number_of_contacts,
                                      int swap1,
                                      std::set<int> *beads_in_use,
                                      std::set<int> * returnMe,
                                      Model *pModel);

    void fillPrBinsAndAssignTotalBin(int * const pBin, float * pDistance, unsigned long int totalDistancesInSphere, Data *pData);


    int numberOfContactsFromSet(std::set<int> *beads_in_use, Model *pModel, int const selectedIndex);
    int numberOfContactsFromSetExclusive(std::set<int> *beads_in_use, Model *pModel, int const selectedIndex, int const excludeIndex);
    int numberOfContacts(int &beadIndex, std::vector<int> *bead_indices, int &workingLimit, Model *pModel, float * pDistance);
    int numberOfContactsExclusive(int &beadIndex, int excludeIndex, std::vector<int> *bead_indices, int &workingLimit, Model *pModel, float * pDistance);


    int numberOfContactsFromSetSeeded(std::set<int> *beads_in_use,
                                      const std::set<int> *true_model,
                                      Model *pModel,
                                      int const selectedIndex);

    void restoreAddingFromBackUp(std::vector<int> * pIndices,
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

    float symViolationsPotential(int specifiedIndex, std::vector<int> * beads_indices, int workingLimit, Model * pModel);

    void printContactList(std::vector<int> &bead_indices, std::set<int> * beads_in_use_tree, int workingLimit, Model * pModel);

    void printContactsFromSet(std::vector<int> &bead_indices, int workingLimit, std::set<int> *beads_in_use,
                              Model *pModel, int const selectedIndex);


    int getRandomNeighbor(int &locale, std::set<int> *beads_in_use, Model * pModel);

    double totalContactsPotential(int average);
    double averageContactsPotential(double average);

    int recalculateContactsSumAdd(std::set<int> *beads_in_use,
                                Model *pModel,
                                int const selectedIndex, int const currentSum);

    double recalculateContactsPotentialSumAdd(std::set<int> *beads_in_use,
                                                          Model *pModel,
                                                          int const selectedIndex,
                                                          double const currentSum);

    int recalculateContactsSumRemove(std::set<int> *beads_in_use,
                                     Model *pModel,
                                     int const selectedIndex,
                                     int const currentSum);

    double recalculateContactsPotentialSumRemove(std::set<int> *beads_in_use,
                                     Model *pModel,
                                     int const selectedIndex,
                                     double const currentSum);

    float connectivityPotentialPhases(int mainConnectivity);
    bool checkSetAndVector(int workingLimit, std::vector<int> * indices, std::set<int> * beads_in_use);
    int getComponentIndex(int index);

    void printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<int> * wl);
    void getExtremes(std::vector<int> & indices, int workingLimit, std::set<int> & extremes, Model * pModel);
public:

    Anneal(float highT,
           float percent,
           int lowerV,
           int upperV,
           int highTempRounds,
           float contacts,
           std::string prefix,
           int totalSASteps,
           float eta,
           float lambda,
           float alpha,
           float mu,
           int multiple);

    ~Anneal(){
    }

    bool createInitialModelSymmetry(Model *pModel, Data *pData);

    float calculateAverageDistance(float * pDistance, int *stopAt, std::vector<int> *bead_indices, Model * pModel);


    void setHighTempExchangeCutoff(float percent){ highTempExchangeCutoff = percent;}

    float calculateKLEnergy(std::vector<int> *bead_indices, std::vector<int> *bins, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData);

    float calculateKLEnergySymmetry(std::vector<int> *bead_indices, std::vector<int> *binCount, int beadIndiciesWorkingLimit, int totalBeadsInSphere, int &violation, Model *pModel, Data *pData);

    //float calculateCVXHULLVolume(char * flags, std::vector<int> *bead_indices, int upTo, double *points, Model *pModel);
    float calculateCVXHULLVolume(char * flags, std::vector<int> *bead_indices, int upTo, Model *pModel);

    void removeFromPr(int removeMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalBeadsInUniverse, std::vector<int> & prBins);

    void addToPr(int addMe, std::vector<int> & beadsInUse, int upperLimit, int * const pBin, int & totalDistancesInSphere, std::vector<int> & prBins);

    void sortMe(std::vector<int> * bead_indices, std::vector<int>::iterator &pSwap1, std::vector<int>::iterator &pSwap2, int * workinglimit);

    float cvxHullCutoff(int temperature, int initial, float hill);

    float calculateContactPotentialFromSum(float workingAverage);

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


    std::string refineSymModel(Model *pModel, Data *pData, std::string nameTo);

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
//    float getLowTempStop(){return lowTempStop;}
    //int getNumberOfCoolingSteps(){return numberOfCoolingTempSteps;}
    float getPercentAddRemove(){ return percentAddRemove;}
    float getLambda(){ return lambda;}
    float getMu(){ return mu;}
    float getEta(){ return eta;}
    float getBeta(){ return beta;}

    void setInterconnectivityCutoff(float value){ this->interconnectivityCutOff = value;}

    void refineMultiphase(Model *pModel, Objective *pObjective, std::vector<Phase *> * pPhases);

    float calculateKLEnergyMP(std::vector<int>::iterator *beginBeadIndices, std::vector<int> *pBinCount,
                              int totalBeadsInSphere, std::vector<int> *pDistancesConvertedToBinPerDataset, unsigned long int totalDistancesInSphere, Data *pData, Model *pModel);

    int findPhases(std::vector<Phase *> *pPhases, int beadInUse);

    float calculateCVXHULLVolumeMP(char *flags, std::vector<int> *bead_indices, int startAt, int upTo, double *points,
                                   Model *pModel);

    float calculateEnergyMP(std::vector<float> *divKLs, float alpha, std::vector<float> *cvxHullVolumes,
                            std::vector<Phase *> *pPhases, float beadVolume, float invBeadVolume);

    int numberOfContactsMP(std::vector<int> *bead_indices, std::vector<Phase *> *pPhases, int phase1, int index1,
                           int phase2, int index2, float &contactCutOff, Model *pModel, float *pDistance);

    bool isConnectedComponent(std::vector<int> *activeIndices, int availableWorkingLimit, float *pDistances,
                              int totalDistances, int &numberOfComponents);

    std::string refineHomogenousBodyASACVX(Model *pModel, Data *pData, std::string name);
    std::string refineHomogenousBodyASACVXSeeded(Model *pModel, Data *pData, std::string outputname);
    std::string reAssignLatticeModel(std::string PDBFilename, Model *pModel, Data *pData);
    std::string refineHomogoenousBodyFromPDBSeed(Model *pModel, Data *pData, int iteration);

    void updateASATemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp);

    int createShortestConnectedPath(std::vector<int> &active_indices,  Model *pModel, const int dmax);

    void initialHighTempSearchSeeded(Model *pModel, Data *pData, std::string name);

    void populateLayeredDeadlimit(std::vector<int>::iterator iteratorBeadIndices, int workingLimit, int *pDeadLimit, Model *pModel, int totalBeads);

    void populateLayeredDeadlimitUsingSet(std::vector<int> * beadIndices, std::set<int> *beads_in_use, int workingLimit,
                                          int *pDeadLimit, Model *pModel);

    void rePopulateLayeredDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use,
                                                    int * pDeadLimit, Model * pModel, int indexOfNewPosition);
    void removeFromdDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, const int workingLimit,
    int * pDeadLimit, Model * pModel, int indexOfRemovedPosition);

    void adjustDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                Model *pModel, int totalWorkingBeads, const int selectedIndex);

    bool checkDKL(float currentKL, std::vector<int> *bead_indices, std::vector<int> *binCount, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData);

    void createSeedFromPDB(Model *pModel, Data *pData, std::string name, std::string PDBFilename, int totalPhases);
    std::string refineInputPDB(Model *pModel, Data *pData, std::string name, std::string PDBFilename, float temperature);

    float calculateKLDivergenceAgainstPDBPR(std::vector<int> &modelPR, std::vector<double> &targetPR);

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

    int changeToDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                 Model *pModel, int totalWorkingBeads, const int selectedInde);

    int additionToDeadLimitPerBead(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit, int *pDeadLimit,
                                   Model *pModel, int totalWorkingBeads, const int selectedIndex);


    double calculateTotalContactSumPotential(std::set<int> *beads_in_use, Model *pModel);
    double calculateTotalContactSum(std::set<int> *beads_in_use, Model *pModel);

    float addToTotalContactEnergy(const int addMe, std::vector<int> *bead_indices, const int workingLimit, Model *pModel,
                                  int totalWorkingBeads, float *pDistance);

    bool createInitialModelCVXHull(Model *pModel, Data *pData, std::string name);
    bool createInitialModelCVXHullSeeded(Model *pModel, Data *pData, std::string name);

    void calculateAverageNumberOfContacts(float *averageContacts, std::vector<int> *bead_indices, const int workingLimit,
                                          Model *pModel, float *pDistance);

    float connectivityPotential(int numberOfComponents);


    void populateNonSeedVector(std::vector<int> *nonSeed, int *totalNonSeed, std::vector<int> *activeIndices,
                               int availableWorkingLimit, std::set<int> *trueModel);

    void populateLayeredDeadlimitUsingFixedPoints(std::vector<int> *bead_indices, std::set<int> *beads_in_use, const std::set<int> *true_model,
                                                  std::set<int> *fixedPoints, const int workingLimit, int *pDeadLimit,
                                                  Model *pModel);

    void populatedDeadLimitExcludingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, std::vector<int> * beads_in_component, int totalInComponent, const int workingLimit,
    int * pDeadLimit, Model * pModel);


    bool isAnchorSafe(int swapPt, std::set<int> *beadsInUseTree, std::set<int> *reducedTrueModelTree, std::set<int> *fixedPoints,
                      int workingLimit, Model *pModel);

    float getAlpha(){return alpha;}

    void populatePotential(int totalNeighbors);

    int getUseableNeighborFromSetDepth(std::set<int> *beads_in_use, Model *pModel, int selectedIndex);
    int getUseableNeighborFromSet(std::set<int> *beads_in_use, Model *pModel, int & selectedIndex);

    bool setAnchorPoints(std::string anchorFileName, std::string pdbFile, Model *pModel);

    bool inComponents(int index);

    bool canRemoveIfAnchor(int index);
};

#include "Model.h"


/**
 *
 */
inline void Anneal::addLatticPositionToModel(std::vector<int> * pIndices,
                                             std::vector<int> * pBackUpState,
                                             int * pWorkingLimit,
                                             std::vector<int>::iterator * pItIndex){


    std::copy(pIndices->begin(), pIndices->end(), pBackUpState->begin()); // make backup
    // make the swap at the border (workingLimit)
    std::iter_swap(pIndices->begin() + *pWorkingLimit, *pItIndex); // this swaps to working position and changes backup
    // increment workingLimit to include new position
    *pWorkingLimit += 1;
    std::sort(pIndices->begin(), pIndices->begin() + *pWorkingLimit);
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
        row2 = (*row)*(unsigned long int)totalBeadsInUniverse - ((*row)*(*row+1)*0.5) - (*row) - 1;
        prBins[ *(pBin + row2 + addMe) ]++;
        i++;
    }

    // out of first loop, i = addMe
    i++;
    // Add across row
    row2 = addMe*(unsigned long int)totalBeadsInUniverse - addMe*(addMe+1)*0.5 - addMe - 1;
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
 * recalcualte P(r) distribution then compare against dataset for KL divergence (bead indices must be sorted!)
 *
 */
inline float Anneal::calculateKLEnergy(std::vector<int> *bead_indices, std::vector<int> *binCount, int upTo, int totalBeadsInSphere, Model *pModel, Data *pData) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    int row;
    unsigned long int row2;
    int * modelBin = pModel->getPointerToBins();

    int totalbins = binCount->size();

    // calculate P(r) for beads
    for(int m=0; m < upTo; m++){ // set row
        row = (*bead_indices)[m];
        //row2 = row*(unsigned long int)totalBeadsInSphere - row*(row+1)*0.5 - row - 1;
        row2 = row*totalBeadsInSphere - row*(row+1)/2 - row - 1;

        for(int n=(m+1); n < upTo; n++){ // iterate over columns
            //col = (*bead_indices)[n];
            (*binCount)[ *(modelBin + row2 + (*bead_indices)[n]) ]++;
            if (*(modelBin + row2 + (*bead_indices)[n]) >= totalbins){
                std::cout << "EXCEEDED " << totalbins << " <= " << *(modelBin + row2 + (*bead_indices)[n]) << std::endl;
                exit(0);
            }

            if (row2 + (*bead_indices)[n] >= pModel->getTotalDistances()){
                std::cout << "EXCEEDED DISTANCES " << (row2 + (*bead_indices)[n]) << " <= " << pModel->getTotalDistances() << std::endl;
                exit(0);
            }
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

        if (neighbor > -1 && beads_in_use->find(neighbor) != endOfSet){

            neighborhoodEnergy += calculateLocalContactPotentialPerBead(beads_in_use, pModel, neighbor);

        } else if (neighbor == -1) {
            break;
        }
    }

    return neighborhoodEnergy;
}

/**
 * Given the selectedIndex, calculate the new sum of contacts for the selected Set defined in beads_in_use
 */
inline int Anneal::recalculateContactsSumRemove(std::set<int> *beads_in_use,
                                             Model *pModel,
                                             int const selectedIndex,
                                             int const currentSum) {

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    int sum = currentSum, totalNewContacts=0;

    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if (neighbor > -1 && beads_in_use->find(neighbor) != endOfSet){ // -1 will be endOfSet and also beads not in use
            sum -=1;
            totalNewContacts+=1;
        } else if (neighbor == -1) {
            break;
        }
    }

    return sum - totalNewContacts;
}


/**
 * Given the selectedIndex, calculate the new sum of contacts for the selected Set defined in beads_in_use
 */
inline double Anneal::recalculateContactsPotentialSumRemove(std::set<int> *beads_in_use,
                                                         Model *pModel,
                                                         int const selectedIndex,
                                                         double const currentSum) {

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    double sum = currentSum;
    int tempContactCount, totalNewContacts=0;

    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if ((neighbor > -1) && (beads_in_use->find(neighbor) != endOfSet)){ // -1 will be endOfSet and also beads not in use
            // get contacts for the neighbor
            // contact count is before removing the selectedIndex
            tempContactCount = numberOfContactsFromSet(beads_in_use, pModel, neighbor);
            sum += totalContactsPotential(tempContactCount-1) - totalContactsPotential(tempContactCount);
            //sum -= 1.0;
            totalNewContacts += 1.0; // counts contacts with respect to selectedIndex
        } else if (neighbor == -1) {
            break;
        }
    }

    // remove contributions of selectedIndex
    return (sum - totalContactsPotential(totalNewContacts));
    //return (sum - totalNewContacts);
}


/**
 * Given the selectedIndex, calculate the new sum of contacts for the selected Set defined in beads_in_use
 */
inline double Anneal::recalculateContactsPotentialSumAdd(std::set<int> *beads_in_use,
                                             Model *pModel,
                                             int const selectedIndex,
                                             double const currentSum) {

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    double sum = currentSum;
    int currentContactCountNeighbor, newContactsFromSelectedIndex=0;

    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if ((neighbor > -1) && beads_in_use->find(neighbor) != endOfSet){ // -1 will be endOfSet and also beads not in use
            // get contacts for the neighbor
            // since the selectedIndex isn't added yet, we add one to the count
            currentContactCountNeighbor = numberOfContactsFromSet(beads_in_use, pModel, neighbor);

            // adjust sum by removing old contribution and adding new one
            sum += totalContactsPotential(currentContactCountNeighbor+1) - totalContactsPotential(currentContactCountNeighbor);
            //sum += 1.0;
            newContactsFromSelectedIndex += 1.0;
        } else if (neighbor == -1) {
            break;
        }
    }

    return (sum + totalContactsPotential(newContactsFromSelectedIndex));
    //return (sum + newContactsFromSelectedIndex);
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

    return totalContactsPotential(neighborContacts);
}

/**
 *
 */
inline double Anneal::totalContactsPotential(int value){
    return connectivityPotentialTable[value];
}


inline double Anneal::averageContactsPotential(double average){
    double diff = 1.0-exp(-0.013*(average-contactsPerBead));
//        double diff = exp(abs(average-contactsPerBead));
//        double diff = (average-contactsPerBead); // harmonic potential
    return 1.0*(diff*diff);
//    return diff;
}



inline void Anneal::fillPrBinsAndAssignTotalBin(int * const pBin, float * pDistance, unsigned long int totalDistancesInSphere, Data * pData){
    maxbin=0;
    for(unsigned long int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    maxbin += 1;
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

        if ((neighbor > -1 ) && beads_in_use->find(neighbor) != endOfSet){
            neighborContacts += 1;
        } else if (neighbor == -1) {
            break;
        }
    }

    return neighborContacts;
}


/**
 * Find a possible empty neighbor directly incontact with a
 */
inline int Anneal::getUseableNeighborFromSetDepth(
        std::set<int> *beads_in_use,
        Model *pModel,
        int selectedIndex){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    int count=0;
    std::vector<int> possibleNeighbors(totalNeighbors);

    // create list of neighbors that are connected to selectedIndex
    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if ((neighbor > -1 ) && beads_in_use->find(neighbor) != endOfSet){
            possibleNeighbors[count] = neighbor;
            count++;
        } else if (neighbor == -1) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return -1;
    } else {
        // pick a neighbor and use
        int neighbor_to_use = possibleNeighbors[rand()%count];
        // find a neighbor
        std::vector<int> possibleNeighborsToSelect(totalNeighbors);
        std::vector<int>::iterator kit = pModel->getPointerToNeighborhood(neighbor_to_use);
        count=0;
        for (int i=0; i< totalNeighbors; i++){
            int neighbor = *(kit+i);
            if ((neighbor > -1 ) && (beads_in_use->find(neighbor) == endOfSet) && (neighbor != selectedIndex) ){
                possibleNeighborsToSelect[count]=neighbor;
                count++;
            } else if (neighbor == -1) {
                break;
            }
        }

        if (count==0){
            return -1;
        } else {
            return possibleNeighborsToSelect[rand()%count];
        }

    }
}


/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
inline int Anneal::getUseableNeighborFromSet(std::set<int> *beads_in_use,
                                           Model *pModel,
                                           int & selectedIndex){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<int> possibleNeighbors(totalNeighbors);

    int count=0;
    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if ((neighbor > -1 ) && beads_in_use->find(neighbor) == endOfSet){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            count++;
        } else if (neighbor == -1) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return -1;
    }

    return possibleNeighbors[rand()%count];
}





inline int Anneal::numberOfContactsFromSetExclusive(std::set<int> *beads_in_use,
                                           Model *pModel,
                                           int const selectedIndex, int const excludedIndex){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    int neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    for (int i=0; i< totalNeighbors; i++){

        int neighbor = *(it+i);

        if (neighbor != excludedIndex && beads_in_use->find(neighbor) != endOfSet){
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
        row2 = (*row)*(unsigned long int)totalBeadsInUniverse - ((*row)*(*row+1)*0.5) - (*row) - 1;
        prBins[ *(pBin + row2 + removeMe) ]--;
        i++;
    }
    // when loop ends beadsInUse[i] => removeMe
    i++;
    // remove across row
    row2 = removeMe*(unsigned long int)totalBeadsInUniverse - removeMe*(removeMe+1)*0.5 - removeMe - 1;
    while (i < upperLimit){
        prBins[ *(pBin + row2 + beadsInUse[i]) ]--;
        i++;
    }

}

inline void Anneal::restoreAddingFromBackUp(std::vector<int> * pIndices,
                                            std::vector<int> * pBackUpState,
                                            int * pWorkingLimit,
                                            std::vector<int> * pBinCountBackUp,
                                            std::vector<int>::iterator * pBinCountBegin){

    std::copy(pBackUpState->begin(), pBackUpState->end(), pIndices->begin());
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


inline float Anneal::connectivityPotentialPhases(int mainConnectivity){
    float currentConnectivityPotential = (mainConnectivity-1.0)*(mainConnectivity-1.0);
    float temp;
    for(int i=0; i < components.size(); i++) {
        //std::cout << i << " component tour => "<< components[i].getTotalNumberOfComponents() << std::endl;
        temp = components[i].getTotalNumberOfComponents()-1.0;
        currentConnectivityPotential += temp*temp;
    }
    return currentConnectivityPotential;
}

inline bool Anneal::inComponents(int index){

    for (int t=0; t<totalComponents; t++){
        if (components[t].inUse(index)){
            return true;
        }
    }
    return false;
}


inline void Anneal::sortMe(std::vector<int> * bead_indices, std::vector<int>::iterator &pOldValue, std::vector<int>::iterator &pNewValue, int * workinglimit){

    if (pOldValue > pNewValue){
      //  std::sort(bead_indices->begin(), bead_indices->begin()+)
        std::iter_swap(pOldValue, pNewValue);
        std::sort(bead_indices->begin(), pOldValue); // bead_indices needs to be sorted
    } else {
        std::iter_swap(pOldValue, pNewValue);
        std::sort(pOldValue, bead_indices->begin() + *workinglimit); // bead_indices needs to be sorted
    }
}


inline int Anneal::getComponentIndex(int index) {
    for(int i=0; i < totalComponents; i++) { // the selected set of beads that make up each component will be held by Component object
        if (components[i].inUse(index)){
            return i;
        }
    }
    return totalComponents;
}

inline float Anneal::calculateContactPotentialFromSum(float workingAverage){
    float diff = workingAverage - contactsPerBead;
    return diff*diff;
}

/**
 * define the set of points within quadrilateral and exclude from positional refinement, only consider points near the convex hull
 *
 */
inline void Anneal::getExtremes(std::vector<int> & indices, int workingLimit, std::set<int> & extremes, Model * pModel){

/**
 * FIND THE EXTREMES (which indices in indices)
 *
 * MAXIMIZE
 * (x-y+z) q1
 * (x-y-z) q2
 * (x+y+z) q3
 * (x+y-z) q4
 *
 * MINIMIZE
 * (x-y+z) q1
 * (x-y-z) q2
 * (x+y+z) q3
 * (x+y-z) q4
 */

    float tempMinQ1 = pModel->getBead(indices[0])->getQ1();
    float tempMinQ2 = pModel->getBead(indices[0])->getQ2();
    float tempMinQ3 = pModel->getBead(indices[0])->getQ3();
    float tempMinQ4 = pModel->getBead(indices[0])->getQ4();
    float tempMaxQ1 = pModel->getBead(indices[0])->getQ1();
    float tempMaxQ2 = pModel->getBead(indices[0])->getQ2();
    float tempMaxQ3 = pModel->getBead(indices[0])->getQ3();
    float tempMaxQ4 = pModel->getBead(indices[0])->getQ4();

    int indexMinQ1=indices[0];
    int indexMaxQ2=indexMinQ1;
    int indexMaxQ3=indexMinQ1, indexMinQ2=indexMinQ1, indexMaxQ1=indexMinQ1;
    int indexMaxQ4=indexMinQ1, indexMinQ3=indexMinQ1, indexMinQ4=indexMinQ1;


    for(int i=1; i<workingLimit; i++){

        int tempIndex = indices[i];
        Bead * pBead = pModel->getBead(tempIndex);
        bool kept = false;
        // min
        if(pBead->getQ1() < tempMinQ1){
            tempMinQ1 = pBead->getQ1();
            indexMinQ1 = tempIndex;
            kept = true;
        }

        if(pBead->getQ2() < tempMinQ2){
            tempMinQ2 = pBead->getQ2();
            indexMinQ2 = tempIndex;
            kept = true;
        }

        if(pBead->getQ3() < tempMinQ3){
            tempMinQ3 = pBead->getQ3();
            indexMinQ3 = tempIndex;
            kept = true;
        }

        if(pBead->getQ4() < tempMinQ4){
            tempMinQ4 = pBead->getQ4();
            indexMinQ4 = tempIndex;
            kept = true;
        }

        //max
        if(pBead->getQ1() > tempMaxQ1){
            tempMaxQ1 = pBead->getQ1();
            indexMaxQ1 = tempIndex;
            kept = true;
        }

        if(pBead->getQ2() > tempMaxQ2){
            tempMaxQ2 = pBead->getQ2();
            indexMaxQ2 = tempIndex;
            kept = true;
        }

        if(pBead->getQ3() > tempMaxQ3){
            tempMaxQ3 = pBead->getQ3();
            indexMaxQ3 = tempIndex;
            kept = true;
        }

        if(pBead->getQ4() > tempMaxQ4){
            tempMaxQ4 = pBead->getQ4();
            indexMaxQ4 = tempIndex;
            kept = true;
        }

    }

    // using min and max define outerset
    extremes.insert(indexMinQ1);
    extremes.insert(indexMinQ2);
    extremes.insert(indexMinQ3);
    extremes.insert(indexMinQ4);

    extremes.insert(indexMaxQ1);
    extremes.insert(indexMaxQ2);
    extremes.insert(indexMaxQ3);
    extremes.insert(indexMaxQ4);
}

#endif //IKETAMA_ANNEAL_H
