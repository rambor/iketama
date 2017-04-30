//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_MODEL_H
#define IKETAMA_MODEL_H

#include <string>
#include <vector>
#include <bead.h>
#include <functions.h>
#include <Phase.h>
#include <Component.h>
#include <Data.h>
#include <math.h>
#include <regex>
#include <iostream>

// Forward declaration
class Data;
//class Component;
class Anneal;

class Model {

private:
    int number_of_beads, total_in_seed, total_in_reduced_seed, total_in_working_universe;
    int workingLimit, deadLimit, startingWorkingLimit, startingDeadLimit;
    float beadAverage, beadStDev;
    int sizeOfNeighborhood;
    float cutOffNeighbor;
    unsigned long int totalDistances;
    bool tooLarge;

    float limit;
    float bead_radius, inv_bead_radius;
    float bead_volume;
    float cvx_volume;
    float averageNumberOfContactsInModel;
    int numberOfSubUnits=1;
    int symmetryIndex=1;
    int lowestWorkingLimit;
    int lowerN, upperN; // number of beads
    std::string symmetry="C1";

    std::vector<Bead> beads;
    std::vector<vector3> beadsWithSymmetry;
    std::vector<float> distances; // linear array (n*(n-1)/2)
    std::vector<int> bins; // linear array (n*(n-1)/2)
    std::vector<int> bead_indices;
    std::vector<int> anchors;
    std::vector<int> neighbors; // 12 * bead_indices size
    std::vector<int> starting_set;
    std::vector<int> selected;
    std::vector<int> seed_indices;  // contains true model of PDB model converted to lattice model
    std::set<int> reduced_seed;
    std::vector<Component> components;

    std::vector<Data *> datasets;

    // setup for dynamically allocated arrays
    float * rThetaPhiAtomType = NULL;

    static constexpr float sqrt6 = 2.449489742783178;
    static constexpr float invsqrt6 = 0.4082482904638631;
    static constexpr float sqrt3 = 1.7320508075688772;
    static constexpr float inv3 = 0.3333333333333333;

    int getIndexInDistanceVector(int row, int col);
    void createDistancesAndConvertToSphericalCoordinates();
    void checkNeighborsList();


public:

    Model(float size, float bead_r, bool fastmode);
    Model(float size, float bead_r, bool fastmode,  std::string sym);

    ~Model(){
        delete[] rThetaPhiAtomType;
        rThetaPhiAtomType = NULL;

        for(std::vector<Data *>::iterator lit = datasets.begin(); lit != datasets.end(); ++lit){
            lit = datasets.erase(lit);
        }
        datasets.clear();

    }

    const float * getrThetaPhiAtomType() const {return rThetaPhiAtomType;}
    int getTotalNumberOfBeadsInUniverse() const {return number_of_beads;}

    //void createAmplitudes(int Lmax);
    void setNumberOfPoints(int lower, int upper){this->lowerN = lower; this->upperN = upper;}
    int getLowerNumberOfPoints() const {return lowerN;}
    int getUpperNumberOfPoints() const {return upperN;}
    // creates vectors for each phase => length is  number_of_beads
    // void initializePhases(std::vector<Phase *> &phases);

    void printSelectedBeads(int startAt, int stopAt, std::vector<int> &beadIDs);

    float getBeadVolume() const {return bead_volume;}
    unsigned long int getTotalDistances() const {return totalDistances;}

    /**
     * SizeOfNeighborhood depends on interconnectivity distance
     * 12 if defined by 2*bead_radius,
     * 18 if defined by 2*sqrt(2)*bead_radius,
     * 40 if defined by 3*sqrt(2)*beadadius
     */
    int getSizeOfNeighborhood(){ return sizeOfNeighborhood;}
    int * getPointerToBins(){ return &bins[0];}

    int getWorkingLimit(){ return workingLimit;}
    int getDeadLimit(){ return deadLimit;}


    void setModelBeginIterator(std::vector<int>::iterator &it);
    void setModelEndIterator(std::vector<int>::iterator &it);

    void setModelBeginStartingSetIterator(std::vector<int>::iterator &it);
    void setModelEndStartingSetIterator(std::vector<int>::iterator &it);


    void setBeadAverageAndStdev(float volume, float stdev);
    float getVolumeAverage(){return beadAverage;}
    float getVolumeStdev(){return beadStDev;}

    float * getDistanceAt(int i){ return &(distances[i]);}
    float distanceToPlane(float dx, float dy, float dz, vector3 plane, float rlength);
    float * getPointerToDistance(){ return &distances[0]; }

    int findIndex (const float &x, const float &y, const float &z);

    Bead * getBead(int i) { return &beads[i];}

    const int * getIndexOfBeadInGroup(int groupIndex, int subUnitIndex) const;

    float getBeadRadius(){ return bead_radius;}

    void makeCopyBeadIndices(std::vector<int> &copyOf) const { std::copy(bead_indices.begin(), bead_indices.end(), copyOf.begin());}
    void createCoordinatesOfSymBeadUniverse();
    void clearCoordinatesOfSymBeadUniverse(){ std::vector<vector3>().swap(beadsWithSymmetry);}

    void writeModelToFile(int workingLimit, std::vector<int> &beads, std::string nameOf, int steps);

    float getNumberOfSubUnits(){return numberOfSubUnits;}
    void getNeighbors(int indexOfBead, std::vector<int> &indices, int &totalNeighbors);

    std::string getSymmetryString(){ return symmetry;}

    void updateBeadIndices(int workingLimit, int deadlimit, std::vector<int> &indices);

    void transformCoordinatesBySymmetry(int subunitIndex, int workingLimit, int &startCount, std::vector<vector3> &coordinates);

    std::string writeSymModelToFile(float dkl, int workingLimit, std::vector<int> &selectedBeads, std::vector<int> &pofrModel, std::string nameOf, Anneal *annealedObject, Data *pData, int steps, float volume, float averageContacts);
    std::string writeModelToFile2(float dkl, int workingLimit, std::vector<int> &selectedBeads, std::vector<int> &pofrModel, std::string nameOf, Anneal *annealedObject, Data *pData, int steps, float volume, float averageContacts);

    void transformCoordinateBySym(int subunitIndex, vector3 &vec, vector3 &newVec);
    void setCVXHullVolume(float volume){ this->cvx_volume = volume; }
    float getCVXHullVolume () const { return this->cvx_volume; }

    void setAverageNumberOfContactsInModel(float number){ this->averageNumberOfContactsInModel = number; }
    float getAverageNumberOfContactsInModel() const { return this->averageNumberOfContactsInModel; }

    int mapSubUnitIndex(int index);

    vector3 *getVectorOfBeadCoordinateInSymBeadUniverse(int index);

    std::string createHeader(float dkl, Anneal *annealedObject, Data *pData, int steps, int workingNumber, float volume, float averageContacts);

    void writeSubModelToFile(int startIndex, int workingLimit, std::vector<int> &selectedBeads, std::string nameOf);

    void setStartingWorkingLimit(int limit){ startingWorkingLimit = limit;}
    int getStartingWorkingLimit(){return startingWorkingLimit;}

    void setStartingDeadLimit(int limit){ startingDeadLimit = limit;}
    int getStartingDeadLimit(){return startingDeadLimit;}


    /**
     *
     */
    void setStartingSet(std::vector<int> & beads){
        total_in_working_universe = beads.size();
        starting_set.resize(total_in_working_universe);
        std::copy(beads.begin(), beads.end(), starting_set.begin());
    }

    float getDistanceBetweenTwoBeads(int indexOne, int indexTwo);
    void setLowestWorkingLimit(int limit){this->lowestWorkingLimit = limit;}
    int getLowestWorkingLimit() const { return this->lowestWorkingLimit;}

    std::vector<int>::iterator getPointerToNeighborhood(int index);

    void createSeedFromPDB(std::string filename, Data * pData, int totalBins, std::vector<double> * pdbPr);

    void setReducedSeed(int limit, std::vector<int> &reduced);
    int getTotalInSeed(){ return this->total_in_seed; }
    int getTotalInReducedSeed(){ return this->total_in_reduced_seed; }
    int getTotalInWorkingUniverse(){ return this->total_in_working_universe; }
    const std::vector<int>::const_iterator getSeedBegin() const { return seed_indices.cbegin(); }
    const std::vector<int>::const_iterator getSeedEnd() const { return seed_indices.cend(); }


    const std::vector<int>::const_iterator getAnchorsIterator() const { return anchors.cbegin(); }
    const std::set<int>::const_iterator getReducedSeedBegin() const { return reduced_seed.cbegin(); }
    const std::set<int>::const_iterator getReducedSeedEnd() const { return reduced_seed.cend(); }


    void centerLatticeModel(int limit, std::vector<int> &reduced);
    bool inReducedSeed(int index);
    int getTotalAnchors(){ return anchors.size(); }


    //void estimatePointsPerComponent(float value);
    void printBeadsFromSet(std::set<int> &beadIDs);
    void printNeighborhood(int index);
    //Component * getComponentByIndex(int index);
    //int getEstimatedLatticePointsPerComponentByIndex(int index);
};


#endif //IKETAMA_MODEL_H
