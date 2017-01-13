//
// Created by Robert Rambo on 12/01/2017.
//

#ifndef IKETAMA_DATUM_H
#define IKETAMA_DATUM_H

class Datum {

    float q, iofq, sigma;
    float var, invvar;

public:
    Datum(float q, float iofq, float sigma);

    float getQ(){return q;}
    float getI(){return iofq;}
    float getSigma(){return sigma;}
    float getVar(){return var;}
    float getInvVar(){return invvar;}

};


#endif //IKETAMA_DATUM_H
