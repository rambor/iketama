//
// Created by Robert Rambo on 07/08/2017.
//

#ifndef IKETAMA_MULTIDATACOMPONENT_H
#define IKETAMA_MULTIDATACOMPONENT_H

#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <Data.h>
#include <regex>
#include <sys/stat.h>
#include <utility>
//#include <Component.h>
#include <Universe.h>

// create single universe using smallest possible lattice dimensions
// each model belongs to one or more distributions

struct find_by_string : std::unary_function< std::string, bool>{

    find_by_string(std::string keyToFind) : key(keyToFind){}

    bool operator () (std::string p) { //
        boost::to_upper(p); // string from read in line
        return (key == p);
    }
    private:
    std::string key;
};

class MultiDataComponent {

   // std::vector<Component> components;
   typedef std::map<std::string,Data>::iterator it_type;
    //Objective datasets;
    std::map<std::string, Data> datasets;
    std::map<std::string, float> weights;
    std::map<std::string, Universe> universes;
    std::set<std::string> phases;
    std::map<std::string, std::string> links;
    std::vector<int> link_potential;

    std::map<std::string, Universe * > baseUniversesForEachPhase; // maps phase to baseUniversie
   // std::map<std::string, Phase> phases;
   // std::map<float, std::string> orderOfUniverses;

    //boost::regex prDataFormat;

    //double smallestBinWidth =100.0f;
    double universalDmax, universalBinWidth;
    double inv_total_universes;

    double eta, lambda, alpha, mu, percentAddRemove;

    std::string setupFile;

    bool fileExists(const std::string & file);
    bool checkPhaseFormat();

    void addDataset(std::string filename, std::string name);

    double calculateVolumeEnergy();
    double calculateTotalKLEnergy();

    int highTempRounds;

    int totalSASteps;
    float contacts;
    double totalInvPhases;

    void makeCopyOfCurrentState();

    double calculateLinkPotential(double lambda);
    double calculateTotalComponentEnergy(double hLambda);
    double calculateContactPotential();
    void updateASATemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp);

    void printSummary();

public:

    MultiDataComponent(std::string setupfilename,
                       int highTempRounds,
                       float contacts,
                       int totalSASteps,
                       float eta,
                       float lambda,
                       float alpha,
                       float mu,
                       float percentAddRemove,
                       float contactsPerBead
    );

    bool createInitialModel();

    void printInfo(std::string text);
    void printError(std::string error);

    void initializeEulerTours();


    void anneal(int ccmulitple, float lowTempStop);
};


inline double MultiDataComponent::calculateVolumeEnergy(){

    double totalVolumeEnergy=0;

    // For each segment, minimize on volume per bead (add more beads or decrease spatial distribution)
    // if segment A is 100000 and segment B is 10000
    //
    // change has to be the same for each segment based on volume
    //
    //
    for(auto &entry : baseUniversesForEachPhase) {
        auto &uni = entry.second;
        totalVolumeEnergy += uni->getSegments(entry.first)->getNormalizedVolume();
    }

//    for(auto &entry : universes) {
//        if (entry.second.getTotalSegments() > 1){
//            totalVolumeEnergy += entry.second.getNormalizedVolume();
//        }
//    }
//
//    totalVolumeEnergy *= 0.2;


    return totalVolumeEnergy;
}


inline double MultiDataComponent::calculateContactPotential(){

    double temp=0;
//
    for(auto &entry : universes) {
        auto &uni = entry.second;
        if (uni.getTotalSegments() > 1){
            temp += uni.getContactsPotential();
        }
        //uni.testContactsUpdating();
    }

    temp *= 0.1; // weight

    // what about contacts potential of the segment
    for(auto &entry : baseUniversesForEachPhase) {
        temp += entry.second->getSegments(entry.first)->getTotalContactsEnergy();
    }

    return temp;
}



inline void MultiDataComponent::makeCopyOfCurrentState(){
    for(auto &entry : universes) {
        auto &uni = entry.second;
        uni.makeCopyOfCurrentState();
    }
}

/**
 *
 * Euler tour of the Universe is not initialized after building the initial models for each phase.
 * Must explicitly initialize tour using total_beads_in_use set that is created after making segments
 *
 * @return
 */
inline void MultiDataComponent::initializeEulerTours(){
    for(auto &entry : universes) {
        auto &uni = entry.second;
        std::cout << " " << uni.getNameOfUniverse() << std::endl;
        uni.initializeEulerTour();
    }
}


inline double MultiDataComponent::calculateTotalKLEnergy() {
    double currentTotalKLEnergy=0.0d;
    for(auto &entry : universes) {
        auto &uni = entry.second;
        currentTotalKLEnergy += uni.calculateKLEnergy();
    }
    return inv_total_universes*currentTotalKLEnergy;
}

inline double MultiDataComponent::calculateTotalComponentEnergy(double hLambda){

    double currentTotalEnergy=0.0d;
    for(auto &entry : universes) {
        auto &uni = entry.second;
        currentTotalEnergy += uni.calculateHarmonicComponentEnergy();
    }

    int diff;
    // for each segment -> get number of components
    for(auto &entry : baseUniversesForEachPhase) {
        auto const &phase = entry.first;
        auto &baseUniverse = entry.second;
        // make initial model of the phase
        diff = baseUniverse->getSegments(phase)->getEulerComponents() - 1;
        currentTotalEnergy += diff*diff;
    }

    return hLambda*currentTotalEnergy;
}






inline double MultiDataComponent::calculateLinkPotential(double lambda) {

    double value=0;

    if (links.size() > 0){

        for(auto &entry : links) {
            auto const &phase1 = entry.first;
            auto const &phase2 = entry.second;
            // make initial model of the phase
            Model * pModel = baseUniversesForEachPhase[phase1]->getModel();
            int totalNeighbors = pModel->getSizeOfNeighborhood();

            Segment * pPhase1 = baseUniversesForEachPhase[phase1]->getSegments(phase1);
            Segment * pPhase2 = baseUniversesForEachPhase[phase2]->getSegments(phase2);
            // need to check how many beads are linked between the two phases
            std::set<int>::iterator p1End = pPhase1->getPointerToBeadsInUseSet()->end();
            int count=0;
            for(std::set<int>::iterator p1It = pPhase1->getPointerToBeadsInUseSet()->begin(); p1It != p1End; ++p1It){
                // check if neighbor is present in pPhase2
                std::vector<int>::iterator pNeighbor = pModel->getPointerToNeighborhood(*p1It);
                for (int i=0; i< totalNeighbors; i++){
                    int neighbor = *(pNeighbor+i);
                    if ((neighbor > -1 ) && pPhase2->isSelected(neighbor)){ // if end of set, means not in use
                        count++;
                        //break; // count number of unique contacts, break moves to next bead position
                    } else if (neighbor == -1) {
                        break;
                    }
                }
            }
            value += (count < 4) ? link_potential[count] : 0;
        }

        return value;
    }
    return 0;
}


inline void MultiDataComponent::updateASATemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=0;

    if (stepEval < 0.15) {
        lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
    } else if (stepEval >= 0.15 && stepEval < 0.65){
        lamRate = 0.44;
    } else if (stepEval >= 0.65){
        lamRate = 0.44*pow(440, -(stepEval - 0.65)*2.857142857);
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
        changed=true;
    } else {
        temp = temp*1.001001001001;
        changed=true;
    }

    if (changed){
        inv_temp = 1.0/temp;
    }
}




#endif //IKETAMA_MULTIDATACOMPONENT_H
