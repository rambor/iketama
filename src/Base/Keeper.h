//
// Created by xos81802 on 01/12/2017.
//

#ifndef IKETAMA_KEEPER_H
#define IKETAMA_KEEPER_H


#include <vector>

class Keeper {

    std::vector<int> indices;
    int upperLimit;
    float d_kl;
    float contactsPotential;


public:
    float totalEnergy;

    Keeper (int upperLimit);

    bool operator< (const Keeper &other) const {
        return totalEnergy < other.totalEnergy;
    }

    void updateIndices(int upperLimit, std::vector<int> &vertexIndices, double d_kl, double contactsE, double totalEnergy);
    void recalculateEnergy(double constant);
    double getTotalEnergy(){ return totalEnergy;}

    int * getStartOfIndices(){ return &indices.front();}
    int getupperLimit(){ return upperLimit;}
    double getDkl() { return d_kl;}
};


#endif //IKETAMA_KEEPER_H
