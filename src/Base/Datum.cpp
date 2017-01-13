//
// Created by Robert Rambo on 12/01/2017.
//


#include "Datum.h"

Datum::Datum(float q, float iofq, float sigma) {
    this->q = q;
    this->iofq = iofq;
    this->sigma = sigma;
    var = sigma*sigma;
    invvar = 1.0/var;
}