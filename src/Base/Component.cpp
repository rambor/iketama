//
// Created by Robert Rambo on 15/02/2017.
//
#include "Component.h"
//#include "../EulerTour/EulerTour.h"
#include "../Model.h"

using namespace std;

Component::Component(std::string id, float volume, Model *pfModel) {
    std::cout <<" CREATING NEW COMPONENT " << id << " VOL: " << volume << std::endl;
    this->id = id;
    this->volume = volume;
    this->pModel = pfModel;
}


void Component::addResid(int index, std::string chain) {
    std::cout << " ADDING RESID : " << index << " CHAIN => " << chain << endl;
    resids.push_back(index);
    chains.push_back(chain);
}

void Component::addCenteredAnchors(int index){
        centeredAnchors.insert(index);
}

/**
 * value should be in units of beads per volume
 */
void Component::setTargetNumberOfLatticePoints(float value){
    this->targetNumberOfBeads = (int)ceil(volume/value);
    if (anchors.size() > 0){
        this->targetNumberOfBeads += anchors.size();
    }
    // populate potential
    this->populatePotential(0.15);
    std::cout << " SET COMPONENT " << id << " TARGET NUMBER => " << targetNumberOfBeads << " " << (volume/value) << std::endl;
}

int Component::getTargetNumberOfLatticePoints(){
    std::cout << " GET COMPONENT " << id << " TARGET NUMBER => " << this->targetNumberOfBeads << " vol " << volume << std::endl;
    return this->targetNumberOfBeads;
}

// Euler tour should include the anchor points
int Component::addLatticePoint(int index){
    beads_in_use.insert(index);
    totalComponents = tour.addNode(index, pModel);
    if (isAnchor(index)){
        anchorCount++;
    }
    return totalComponents;
}

int Component::removeLatticePoint(int index){
    beads_in_use.erase(index);
    totalComponents = tour.removeNode(index);
    if (isAnchor(index)){
        anchorCount--;
    }
    return totalComponents;
}

// ADD ANCHOR OCCURS DURING INITIALIZATION FROM THE FILE
// IF NOT ANCHORS, THEN NOTHIHNG ADDED TO TOUR
void Component::addAnchor(int index){

    if (anchors.find(index) == anchors.end()){
        std::cout << "        ADDING ANCHOR " << index << " " << anchors.size() << endl;
        anchors.insert(index);
        anchorCount++;
        beads_in_use.insert(index);
        totalComponents = tour.addNode(index, pModel);
        empty=false;
    }
}

// if anchor point, returns true
bool Component::isAnchor(int index){
    return !(anchors.find(index) == anchors.end());
}


bool Component::isCenteredAnchor(int index){
    return !(centeredAnchors.find(index) == centeredAnchors.end());
}


void Component::printAnchors(){
    cout << anchors.size() << endl;
}

void Component::populatePotential(float percentage){

    percentageStep = 1.0/ceil(percentage*targetNumberOfBeads);

    potentialFunction.resize(5);
    potentialFunction[0] = 0;
    potentialFunction[1] = 0.0001;
    potentialFunction[2] = 0.0005;
    potentialFunction[3] = 0.001;
    potentialFunction[4] = 0.1;
}


float Component::potential() {

    int absdiff = (int)ceil(std::abs(beads_in_use.size() - targetNumberOfBeads)*percentageStep) - 1;

    if (absdiff > 4){
        absdiff = 4;
    }
    return potentialFunction[absdiff];

//    int diff = beads_in_use.size() - targetNumberOfBeads;
//    int tenPercent = (int)ceil(0.1*targetNumberOfBeads);
//
//    if (diff < -8){
//        return 0.02;
//    } else if (diff >= -8 && diff >= -7) {
//        return 0.007;
//    } else if (diff >= -6 && diff <= -5) {
//        return 0.005;
//    } else if (diff >= -4 && diff <= -3) {
//        return 0.0001;
//    } else if (diff >= -2 && diff <= 2) {
//        return 0;
//    } else if (diff >=3 && diff <= 4) {
//        return 0.0001;
//    } else if (diff >=5 && diff <= 6) {
//        return 0.0002;
//    } else if (diff >=7 && diff <= 8) {
//        return 0.0005;
//    } else if (diff > 8) {
//        return 0.0007;
//    }
}

float Component::calculateCVXVolume() {
    char flags[25];
    sprintf(flags, "qhull s FA");
    int numpoints = 3*beads_in_use.size();
    coordT points[numpoints];

    std::vector<int> active_indices(beads_in_use.size());

    int next, countIt=0;
    for (std::set<int>::iterator it = beads_in_use.begin(); it != beads_in_use.end(); ++it) {
        next = 3*countIt;
        Bead * pBead = pModel->getBead(*it);
        points[next] = pBead->getX();
        points[next+1] = pBead->getY();
        points[next+2] = pBead->getZ();
        active_indices[countIt] = *it;
        countIt++;
    }

    // needs to be optimized
    qh_new_qhull(3, countIt, points, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;
    cvxPoints.resize(totalV);

    for (int v = 0; v < totalV; v++) { //
        cvxPoints[v] = active_indices[qh_pointid(vertices->point)];
        vertices = vertices->next;
    }

    volume = qh totvol;
    qh_freeqhull(true);
    return volume;
}

void Component::printSet(){
    pModel->printBeadsFromSet(beads_in_use);
}

void Component::printBest(){
    pModel->printBeadsFromSet(best);
}


void Component::setBest(){
    best = beads_in_use;
}


void Component::copyBestToInUse(){

    beads_in_use = best;
    std::vector<int> temp;

    for (std::set<int>::iterator it = best.begin(); it != best.end(); ++it) {
        //beads_in_use.insert(*it);
        temp.push_back(*it);
    }

    std::sort(temp.begin(), temp.end());
    totalComponents = tour.newTour(temp.begin(), temp.size(), pModel);
    //tour = EulerTour(temp.begin(), temp.size(), pModel);
    //totalComponents = tour.newTour(temp.begin(), temp.size(), pModel);
}


/**
 * write subset of selected bead to file
 */
void Component::writeToFile(std::string nameOf){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    string residue_index;

    std::vector<int> indices;
    for(std::set<int>::iterator it = beads_in_use.begin(); it != beads_in_use.end(); ++it){
        indices.push_back(*it);
    }

    std::sort(indices.begin(), indices.end());
    fprintf(pFile, "REMARK 265         NUMBER LATTICE POINTS : %s\n", nameOf.c_str());
    fprintf(pFile, "REMARK 265                        TARGET : %i\n", targetNumberOfBeads );
    fprintf(pFile, "REMARK 265                        ACTUAL : %i\n", beads_in_use.size() );
    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    int residueCnt = 1;
    for (int i=0; i<indices.size(); i++){
        currentBead = pModel->getBead(indices[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}


/**
 * write subset of selected bead to file
 */
void Component::writeAnchorsToFile(std::string nameOf){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    string residue_index;

    std::vector<int> indices;
    for(std::set<int>::iterator it = anchors.begin(); it != anchors.end(); ++it){
        indices.push_back(*it);
    }

    std::sort(indices.begin(), indices.end());

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    int residueCnt = 1;
    for (int i=0; i<indices.size(); i++){
        currentBead = pModel->getBead(indices[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
}


void Component::printConstraints(){
    std::cout<< " TARGET NUMBER BEADS => " << targetNumberOfBeads << "( " << beads_in_use.size() << " )" << endl;
}


