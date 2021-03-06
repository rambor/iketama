//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_DATA_H
#define IKETAMA_DATA_H

#include <string>
#include <vector>
#include "Datum.h"
#include "math.h"
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include "Phase.h"

//class Model;

class Data {

private:
    std::string iofqfilename;
    std::string pofrfilename;

    std::vector<Datum> values;
    std::vector<int> pointsPerBin;
    std::vector<float> lowerBounds;

    std::vector<Datum> workingSet;
    std::vector <float> qvalues;

    // i indexes from 0...(probability_per_bin.size()-1)
    // lower => i*bin_width
    // upper => (i+1)*bin_width
    // r_value in bin => 0.5*bin_width + i*bin_width
    std::vector<float> probability_per_bin;
    std::vector<float> moore_coefficients;
    std::vector<float> bin_coefficients;
    std::vector<float> working_probability_per_bin;
    std::vector<Phase *> phases;

    void parseMooreCoefficients();
    void parseBins();

    void normalize(float norm);
    void normalizePrBins();
    float integrateMoore(float lower, float upper);
    //Partials amplitudes;

    struct RealSpace {
        float r, pr, sigma;
        RealSpace(float rvalue, float pvalue, float svalue) : r(rvalue), pr(pvalue), sigma(svalue){}
    };

    std::vector<RealSpace> pofr;

    // array of pointers to my Models being modeled
//    std::vector<Model *> models;

    int totalDataPoints, totalPhases, zeroBin;

    float qmin, qmax, ns_dmax, dmax;
    int lmax, id;
    float shannon_bins, weight=1.0;
    int data_bins;
    float volume;
    float rg;
    float rave;
    float bin_width;

    float integrateMooreSingle(float midpointrvalue, float upper);
    int findInArray(std::string value, std::vector<std::string> * strings);

public:

    Data(); // default constructor, needed for multidatacomponent
    Data(std::string filename);

    ~Data(){
//        for(std::vector<Phase *>::iterator lit = phases.begin(); lit != phases.end(); ++lit){
//            lit = phases.erase(lit);
//        }
//        phases.clear();
    }

    void addPofRData(std::string filename);

    // create CDF

    const float getQmax() const {return qmax;}
    const float getQmin() const {return qmin;}

    const float getDmax() const {return dmax;}

    const Datum getDatum(int index) const {return values[index];}

    float getQ(int index) {return values[index].getQ();}
    float getI(int index) {return values[index].getI();}

    /**
     * error associated with data
     */
    float getSigma(int index) {return values[index].getSigma();}
    float getVar(int index) {return values[index].getVar();}
    float getInvVar(int index) {return values[index].getInvVar();}

    float getVolume() {return volume;}
    float getBinRValue(int index);
    double getBinWidth(){return (double)bin_width;}

    int getShannonBins(){return (int)shannon_bins;}
    int getZeroBin(){return (int)zeroBin;}

    float convertToBinProbability(float distance);
    float convertBinToDistance(int bin);

    // probability_per_bin is the experimental distribution from dat file
    // maxbin must be greater than or equal to probability_per_bin.size
    // maxbin is based on the size of the bead universe
    void createWorkingDistribution(int maxBin){
        working_probability_per_bin.resize(maxBin);
        std::fill(working_probability_per_bin.begin(), working_probability_per_bin.end(), 0);
        //
        std::copy(probability_per_bin.begin(), probability_per_bin.end(), working_probability_per_bin.begin());

        //last nonzero bin
        zeroBin=0;
        for (int i=0; i<maxBin; i++){ // if maxbin is greater than probability_per_bin.size, we have an empty bin
            if (working_probability_per_bin[zeroBin] <= 0){
                break;
            }
            zeroBin++;
        }
    }

    /**
     * converts distances to bins, can be larger than dmax
     */
    inline int convertToBin(float distance){

        float ratio = distance/bin_width;
        int floored = floor(ratio);
        int binlocale;

        // less than or equal to
        // lower < distance <= upper
        // bin = 0 is the distance between two beads
        // due to floating point errors, sometimes distance will be slightly larger than 1
        // need to subtract 1 to make it zero
        //
        if ((ratio-floored) < 0.001 ){
            binlocale = (floored > 0) ? (floored - 1) : 0;
            //binlocale = 0;
//            if (binlocale > 0){
//                std::cout << bin_width <<  " ratio " << ratio << "  " << floored << " => " << binlocale << " distance " << distance << std::endl;
//            }
//            if (binlocale == 0){
//                std::cout << bin_width <<  " ratio " << ratio << "  " << floored << " => " << binlocale << " distance " << distance << std::endl;
//            }
        } else {
            binlocale = floored;
        }

        return binlocale;
    }


    float getProbabilityPerBin(int bin){return probability_per_bin[bin];}

    void setDataBinSize(int bins);
    void normalizePofR(int count);
    void normalizeMoorePofR();
    float calculatePofRUsingMoore(float rvalue);

    //void creatingWorkingSet(Model & model);
    void addPhase(Phase &phase);

    void printKLDivergence(std::vector<unsigned int> &modelPR);
    float calculateKLDivergence(std::vector<unsigned int> &modelPR);
    float calculateKLDivergenceContrast(std::vector<float> &modelPR);
    float calculateKLDivergenceMultiComponent(std::vector<float> &modelPR);
    void calculateRatioPr(std::vector<float> &modelPR);

    // HOLDS POINTERS TO MODELS
    //void addModel(Model & model) { models.push_back(&model); }

    const int getTotalData() const { return totalDataPoints; }

    const std::string getPofRFilename() const { return pofrfilename; }

    static bool checkPofRFile(std::string file);

    int getTotalPhases() const { return phases.size(); }
    Phase * getPhase(int locale){return phases[locale];}

    bool isPhasePresent(std::string id);

    void setWeight(float wt){ weight = wt;}
    float getWeight() const { return weight;}

    void setID(int value){ id=value;}
    int getID() const { return id;}

    float calculateVolumeFromDistribution(int totalBins, std::vector<int> * binCount);
};


#endif //IKETAMA_DATA_H
