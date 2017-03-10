//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_PDBMODEL_H
#define IKETAMA_PDBMODEL_H


#include <boost/regex.hpp>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <numeric>
#include <functional>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "functions.h"
#include <math.h> // debug
#include <fstream>
#include <sstream>
//#include "Waters.h"
#include "Coords.h"
#include <dirent.h>
#include <boost/lexical_cast.hpp>

const float invFourPI = 1/(4.0*M_PI);

class PDBModel {

private:
    std::string filename;
    bool ifRNA;
    bool ifWaters;
    int totalAtoms, totalResidues, waterCount;

    int watersPerResidue;  // maximum number of waters per residue, should match declaration in waters

    // Following vectors will all be same length
    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;
    // waters

    std::vector < int > resID;
    std::vector < float > x,y,z,atomVolume;

    Coords * keptWaters = NULL;

    std::set <int> atomList;

    float volume, dmax, cutoff;

    // pointer to a float array
    float * pconvertXYZ = NULL;

    // setup for dynamically allocated arrays
    float * rThetaPhiAtomType = NULL;
    float * atomicExpVolume = NULL;
    int * atomicNumbers = NULL;

    float * rWater = NULL;
    float * phiWater = NULL;
    float * waterCosTheta = NULL;

    float * centeredX = NULL;
    float * centeredY = NULL;
    float * centeredZ = NULL;

    void setWaterRPhi(int index, float rValue, float phiValue);

public:
    PDBModel(std::string filename, bool rna, bool discardWaters, float lower);
    ~PDBModel(){

        delete[] rThetaPhiAtomType;
        rThetaPhiAtomType = NULL;
        delete[] atomicExpVolume;
        atomicExpVolume = NULL;
        delete[] atomicNumbers;
        atomicNumbers= NULL;

        delete[] centeredX;
        centeredX = NULL;
        delete[] centeredY;
        centeredY = NULL;
        delete[] centeredZ;
        centeredZ = NULL;

        delete[] keptWaters;
        keptWaters = NULL;

        if (ifWaters){
            delete[] rWater;
            rWater=NULL;
            delete[] phiWater;
            phiWater=NULL;
            delete[] waterCosTheta;
            waterCosTheta=NULL;
        }
    } // destructor defined in header file


    void writeCenteredCoordinates() const;
    void writeKeptWaters() const;
    //void addWaters(Waters & waterModel);

    int getWaterCount() const { return waterCount;}

    int getTotalAtoms() const { return totalAtoms;}
    float getDmax() const { return dmax; }
    float getVolume() const { return volume; }
    const std::set <int> &getAtomList() const { return atomList;}

    const float * getrThetaPhiAtomType() const {return rThetaPhiAtomType;}
    const float * getAtomicExpVolume() const {return atomicExpVolume;}

    const float * getAtomVolume() const {return &atomVolume[0];}

    const int * getAtomicNumbers() const {return atomicNumbers;}
    const float * getCenteredX() const {return centeredX;}
    const float * getCenteredY() const {return centeredY;}
    const float * getCenteredZ() const {return centeredZ;}
    const float * getWaterR() const {return &rWater[0];}
    const float * getWaterPhi() const {return &phiWater[0];}
    const float * getWaterCosTheta() const {return &waterCosTheta[0];}

    // resID;
    const std::vector<int>::const_iterator getResIDIterator() const { return resID.cbegin(); }
    const std::vector<std::string>::const_iterator getAtomTypeIterator() const { return trimmedAtomType.cbegin(); }
    const std::vector<std::string>::const_iterator getChainIDIterator() const { return chainID.cbegin(); }

    //void setCutoff(float number){ cutoff = number*number;} // in angstroms
    const Coords * getKeptWaters() const {return &keptWaters[0];}

    void writeCenteredCoordinatesToFile(std::string name) const;

    std::string getAtomTypeByIndex(int index){ return atomType[index];}
};

#endif //IKETAMA_PDBMODEL_H
