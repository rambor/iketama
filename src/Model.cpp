//
// Created by Robert Rambo on 12/01/2017.
//

#include <PDBModel.h>
#include "Model.h"
#include "Anneal.h"

using namespace std;

Model::Model(float size, float bead_r, bool fastmode) {

    float radius = size*0.5;
    limit = radius + radius*0.1011;
    //fastslow = fastmode;

    bead_radius = bead_r;
    inv_bead_radius = 1.0/bead_r;
    bead_volume = 4.0/3.0*bead_radius*bead_radius*bead_radius*M_PI; // 4/3*PI*r^3

    sizeOfNeighborhood = (fastmode) ? 12 : 18;
    //sizeOfNeighborhood = 18;
    // change to 12 if neighborhood is defined by 2*bead_radius, 18 if 2*sqrt(2)*bead_radius, 40 if 3*sqrt(2)*bead radius
    if (sizeOfNeighborhood == 12){
        cutOffNeighbor = 2.001*bead_radius;
    } else if (sizeOfNeighborhood == 18){
        cutOffNeighbor = this->getBeadRadius()*2.8285; // 2*sqrt(2)*bead_radius
    } else {
        cutOffNeighbor = this->getBeadRadius()*3.464*1.001; // 3*sqrt(2)*bead radius
    }

    int klimit = (int) (limit*inv_bead_radius*3.0/2.0*invsqrt6);
    int count=0;

    float distance, dx, dy, dz;
    // float * pconvertXYZ = NULL;
    // positive first
    for (int k=-klimit; k<=klimit; k++){
        // for each k index over i and j
        dz = 2.0*inv3*sqrt6*k;

        if (dz*bead_radius > limit){
            break;
        }

        float inv3kmod2 = inv3*(k%2);

        for(int j=-klimit; j<=klimit; j++){
            dy = sqrt3*(j + inv3kmod2);

            if (dy*bead_radius <= limit){
                float jkmod2 = (j+k)%2;

                for(int i=-klimit; i<=klimit; i++){
                    // compute distance from center
                    dx = 2*i + jkmod2;

                    distance = bead_radius*sqrt(dx*dx + dy*dy + dz*dz);

                    if (distance <= limit){
                        // add bead to vector
                        beads.push_back(Bead(dx*bead_radius, dy*bead_radius, dz*bead_radius,1));
                        //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp( beads[count].getX(), beads[count].getY(), beads[count].getZ());
                        count++;
                    }
                } // end of i loop
            }
        } // end of j loop
    } // end of k loop


    number_of_beads = beads.size();
    cout << " NUMBER OF BEADS " <<  number_of_beads << endl;
    selected.resize(number_of_beads);

    rThetaPhiAtomType = new float[number_of_beads*5];

    // create spherical coordinates
    // string residue_index;
    this->createDistancesAndConvertToSphericalCoordinates();
    /*
    Bead * currentbead;
    vector3 diff;
    int locale, next;

    totalDistances = (int)(number_of_beads*(number_of_beads-1)/2.0);

    distances.resize(totalDistances);
    bins.resize(totalDistances);

    bead_indices.resize(number_of_beads);

    float root_dis;
    count=0;

    // for each bead, calculate SHE expansion
    for (int n=0; n < number_of_beads; n++) {

        currentbead = &(beads[n]);
        pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp( currentbead->getX(), currentbead->getY(), currentbead->getZ());

        locale = n*5;
        rThetaPhiAtomType[locale] = *pconvertXYZ;               // [0] r
        rThetaPhiAtomType[locale+1] = cosf(*(pconvertXYZ+1));   // [1] cos(theta)
        rThetaPhiAtomType[locale+2] = *(pconvertXYZ+2);         // [2] phi
        rThetaPhiAtomType[locale+3] = 1;             // [3] atomic number
        rThetaPhiAtomType[locale+4] = 1.0;                      // [4]
        //string residue_index = to_string(n+1);
        //printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), beads[n].getX(), beads[n].getY(), beads[n].getZ() );

        next = n+1;

        // populate distance matrix
        for(int m=next; m < number_of_beads; m++){
            diff = currentbead->getVec() - (&(beads[m]))->getVec();
            root_dis = diff.length();
            distances[count] = root_dis;

            count++;
        }

        bead_indices[n] = n;
    }
    */
}



Model::Model(float size, float bead_r, bool fastmode, string sym) : Model(size, bead_r, fastmode) {

    symmetry = sym;

    if (symmetry.substr(0,1) == "C" || symmetry.substr(0,1) == "c" ){ // rotation group
        symmetryIndex = atoi(sym.erase(0,1).c_str());
        numberOfSubUnits = symmetryIndex; // number of identical subunits
    } else if (symmetry.substr(0,1) == "D" || symmetry.substr(0,1) == "d" ) { // dihedral group
        symmetryIndex = atoi(sym.erase(0,1).c_str());
        numberOfSubUnits = 2*symmetryIndex; // number of identical subunits
    }

    if (number_of_beads*(number_of_beads-1)*0.5*numberOfSubUnits*2 < 250000000){ // number of bytes used to store P(r) as Shannon bins
        tooLarge = false;
    } else {
        tooLarge = true;
    }

    beadsWithSymmetry.resize(number_of_beads*numberOfSubUnits);
    this->createCoordinatesOfSymBeadUniverse();
}


/**
 * create coordinates of symmetry bead universe, destroy after use
 */
void Model::createCoordinatesOfSymBeadUniverse(){

    // for each bead
    int count=0;
    vector3 * tempVec;
    Bead * tempBead;

    for (int i=0; i< number_of_beads; i++){
        tempBead = &beads[i];
        tempVec = &beadsWithSymmetry[count];

        tempVec->x = tempBead->getX();
        tempVec->y = tempBead->getY();
        tempVec->z = tempBead->getZ();

        count++;

        for (int s=1; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector
            this->transformCoordinateBySym(s, *tempVec, beadsWithSymmetry[count]);
            count++;
        }
    }
}

/**
 *
 */
vector3 * Model::getVectorOfBeadCoordinateInSymBeadUniverse(int index){
    return &beadsWithSymmetry[index];
}


/**
 * map subunit index to beadsWithSymmetry vector
 */
int Model::mapSubUnitIndex(int index){
    return index*numberOfSubUnits;
}

/**
 *
 */
void Model::transformCoordinateBySym(int subunitIndex, vector3 &vec, vector3 &newVec){

    float invSymmetryIndex = 1.0/(float)symmetryIndex;
    float angle = M_PI*2.0*subunitIndex*invSymmetryIndex;
    float cosine = cos(angle);
    float sine = sin(angle);
    float x, y, z;

    if (subunitIndex < symmetryIndex){ // rotations for C and D symmetry groups

        x = vec.x;
        y = vec.y;
        z = vec.z;
        newVec.x = cosine*x - sine*y;
        newVec.y = sine*x + cosine*y;
        newVec.z = z;

    } else { // reflections for the D_group
        // rotate on x-axis first by 180 and then by
        // rotating on x-axis (x,y,z) -> (x,-y,-z)
        x = vec.x;
        y = -vec.y;
        z = -vec.z;

        //newVec.x = cosine*x + sine*y;
        //newVec.y = sine*x - cosine*y;
        newVec.x = cosine*x - sine*y;
        newVec.y = sine*x + cosine*y;
        newVec.z = z;
    }
}


/**
 * for each bead in the bead universe, calculate distances and convert to spherical coordinates
 * pre-populate a neighbors list for each bead based on a distance cutoff
 */
void Model::createDistancesAndConvertToSphericalCoordinates(){

    float * pconvertXYZ;
    Bead * currentbead;
    vector3 diff;
    int locale, next;

    totalDistances = ((unsigned long int)number_of_beads*(number_of_beads-1.0)*0.5);
    cout << " TOTAL DISTANCES " << totalDistances << " ( MAX MEMORY => " << bead_indices.max_size() << " )" << endl;
    cout << "       MAX INDEX " << std::numeric_limits<int>::max() << endl;

    try{
        if (totalDistances > bead_indices.max_size()){
            throw std::invalid_argument( "PHYSICAL SYSTEM IS TOO SMALL < NOT ENOUGH MEMEORY => REDUCE RESOLUTION: \n");
        } else {
            cout << " RESIZING DISTANCES VECTOR " << endl;
            distances.resize(totalDistances);
            cout << "      RESIZING BINS VECTOR " << endl;
            bins.resize(totalDistances);

            bead_indices.resize(number_of_beads);

            neighbors.resize(sizeOfNeighborhood*number_of_beads);
            starting_set.resize(number_of_beads);

            std::fill(neighbors.begin(), neighbors.end(), -1);
            float root_dis;
            unsigned long int discount=0;

            for (int n=0; n < number_of_beads; n++) {

                currentbead = &(beads[n]);
                pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp( currentbead->getX(), currentbead->getY(), currentbead->getZ());

                locale = n*5;
                rThetaPhiAtomType[locale] = *pconvertXYZ;               // [0] r
                rThetaPhiAtomType[locale+1] = cosf(*(pconvertXYZ+1));   // [1] cos(theta)
                rThetaPhiAtomType[locale+2] = *(pconvertXYZ+2);         // [2] phi
                rThetaPhiAtomType[locale+3] = 1;                        // [3] atomic number
                rThetaPhiAtomType[locale+4] = 1.0;                      // [4]
                //string residue_index = to_string(n+1);
                //printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), beads[n].getX(), beads[n].getY(), beads[n].getZ() );

                next = n+1;
                // populate distance matrix as 1-D
                for(int m=next; m < number_of_beads; m++){
                    diff = currentbead->getVec() - (&(beads[m]))->getVec();
                    root_dis = diff.length();
                    distances[discount] = root_dis;
                    discount++;
                }

                bead_indices[n] = n;
            }

            cout << " POPULATING NEIGHBORS " << endl;
            // populate neighbors list

            for (int n=0; n < number_of_beads; n++){

                int count=0;
                // going down rows (n is fixed column position) m < n
                for (int m=0; m < n ; m++){
                    //if (distances[(int)(n*number_of_beads - 0.5*n*(n+1)) - n - 1 + m] < cutOffNeighbor){
                    if (distances[((unsigned long int)m*number_of_beads - 0.5*m*(m+1)) - m - 1 + n] < cutOffNeighbor){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }

                // fixed row, going across columns n < m
                for (int m=(n+1); m < number_of_beads ; m++){

                    if (distances[((unsigned long int)n*number_of_beads - 0.5*n*(n+1)) - n - 1 + m] < cutOffNeighbor){
                        neighbors[sizeOfNeighborhood*n + count] = m;
                        count++;
                    }
                }

            }

            //this->checkNeighborsList();
            cout << " FINISHED NEIGHBORS " << endl;
        }
    } catch (exception &err) {
        cerr<<"Caught "<<err.what()<<endl;
        cerr<<"Type "<<typeid(err).name()<<endl;
        exit(0);
    }

}


// error is in distance or index
void Model::checkNeighborsList(){

    std::vector<int>::iterator it;

    for (int n=0; n < number_of_beads; n++){

        int anchor = n;
        Bead * currentbead = &(beads[n]);
        it = this->getPointerToNeighborhood(anchor);
        for(int neigh =0; neigh<sizeOfNeighborhood; neigh++){
            int neighbor = *(it + neigh);
            if (neighbor > -1){
                // calculate distance
                vector3 diff = currentbead->getVec() - (&(beads[neighbor]))->getVec();
                double root_dis = diff.length();
                if (root_dis > cutOffNeighbor){
                    cout << n << " TOO FAR " << neighbor << endl;
                }
            }
            if (neighbor == 0){
                cout << n << " neighbors of zero " << neighbor << endl;
            }
        }
    }
}



/**
 * indices of lattice points that comprise neighborhood of lattice point at index
 */
std::vector<int>::iterator Model::getPointerToNeighborhood(int index){
//    if (neighbors[sizeOfNeighborhood*index] != *(neighbors.begin() + sizeOfNeighborhood*index)){
//        cout << "NOT EQUAL Q!!!!!" << endl;
//    }
//    cout << index << " " << neighbors[sizeOfNeighborhood*index] << " " << this->getDistanceBetweenTwoBeads(index, neighbors[sizeOfNeighborhood*index]) <<  endl;
    return (neighbors.begin() + sizeOfNeighborhood*index);
}


void Model::printNeighborhood(int location){

    string residue_index;

    for(int i=0; i<sizeOfNeighborhood; i++){
        int index = *(neighbors.begin() + sizeOfNeighborhood*location + i);
        cout<<"INDEX " << index <<endl;
        residue_index = to_string(location);
        printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), beads[index].getX(), beads[index].getY(), beads[index].getZ() );
        //cout << index << " => "  << i << " NEIGHBOR " << *(neighbors.begin() + sizeOfNeighborhood*index + i) << endl;
    }
}


/**
 * calculate distance to plane where point is given by dx, dy, dz
 * plane is from the cross product of two vectors
 * distance can be negative or positive depending on how cross product was taken.
 */
inline float Model::distanceToPlane(float dx, float dy, float dz, vector3 plane, float rlength){
    return (plane.x*dx + plane.y*dy + plane.z*dz)*rlength;
}


/**
 * find bead index corresponding to coordinates
 */
int Model::findIndex(const float &x, const float &y, const float &z)  {

    Bead * temp;
    int found;

    for(int i=0; i < number_of_beads; i++){
        temp = &beads[i];

        if ((temp->getX() == x) && (temp->getY() == y) && (temp->getZ() == z)) {
            found = i;
            break;
        } else if (abs(temp->getX()-x) + abs(temp->getY()-y) + abs(temp->getZ()-z) < 1){
            cout << "USING DIFF " << endl;
            found = i;
            break;
        }
    }

    return found;
}


/**
 * calculate distance between two lattice points (unordered)
 */
float Model::getDistanceBetweenTwoBeads(int indexOne, int indexTwo){
    int first = indexOne;
    int second = indexTwo;

    if (indexOne > indexTwo){
        first = indexTwo;
        second = indexOne;
    }

    return distances[(unsigned long int)first*number_of_beads - (first*(first+1)*0.5) - first - 1 + second];
}


void Model::getNeighbors(int indexOfBead, std::vector<int> &indices, int &totalNeighbors){

    // beadsInUse must be sorted
    // float contactCutOff = bead_radius*3 + 0.02*bead_radius;
    float contactCutOff = 3.5*bead_radius;
    float lowerCutOff = 2.01*bead_radius;
    float dist;
    int sizeOfInputArray = 60;

    unsigned long int row2;

    int i=0, toKeep=0;
    // add down column (column is beadIndex) over all beads
    while(i < indexOfBead){
        // row => i;
        // row2 = row*number_of_beads - (row*(row+1)*0.5) - row - 1;
        // column => indexOfBead
        dist = distances[(unsigned long int)i*number_of_beads - (i*(i+1)*0.5) - i - 1 + indexOfBead ];
        if (lowerCutOff < dist && dist < contactCutOff) {
            //if (dist < lowerCutOff) {
            indices[toKeep] = i;
            toKeep++;
        }
        i++;
    }

    // out of first loop, i = addMe
    i++;
    // Add across row
    if (toKeep < sizeOfInputArray){
        row2 = (unsigned long int)indexOfBead*number_of_beads - indexOfBead*(indexOfBead+1)*0.5 - indexOfBead - 1;
        while (i < number_of_beads && toKeep < sizeOfInputArray){
            dist = distances[row2 + i ];
            if (lowerCutOff < dist && dist < contactCutOff) {
                //if (dist < lowerCutOff) {
                indices[toKeep] = i;
                toKeep++;
            }
            i++;
        }
    }

    totalNeighbors = toKeep;
}

/*
void Model::initializePhases(std::vector<Phase *> &phases) {

    int totalPhases = phases.size();

    for(int i=0; i<totalPhases; i++){
        phases[i]->setStatusVector(number_of_beads);
    }
}
 */


void Model::printBeadsFromSet(std::set<int> &beadIDs){

    std::string residue_index;
    for(std::set<int>::iterator it=beadIDs.begin(); it!=beadIDs.end(); ++it){
        int n = *it;
        residue_index = to_string(n+1);
        printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), beads[n].getX(), beads[n].getY(), beads[n].getZ() );
    }
}


void Model::printSelectedBeads(int startAt, int stopAt, vector<int> &beadIDs){

    string residue_index;
    for(int i=startAt; i<stopAt; i++){
        int n = beadIDs[i];
        residue_index = to_string(n+1);
        printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), beads[n].getX(), beads[n].getY(), beads[n].getZ() );
    }
}


int Model::getIndexInDistanceVector(int row, int col) {

    int locale;

    if (col > row){
        locale = row*number_of_beads - row*(row+1)*0.5 + (col - row - 1);
    } else {
        locale = col*number_of_beads - col*(col+1)*0.5 + (row - col - 1);
    }

    return locale;
}

/**
 * holds the number of beads, volume is detemined by multiplying by bead volume
 */
void Model::setBeadAverageAndStdev(float number_of_beads, float stdev) {
    this->beadAverage = number_of_beads;
    this->beadStDev = stdev;
}


/**
 * write subset of selected bead to file
 */
void Model::writeSubModelToFile(int startIndex, int workingLimit, vector<int> &selectedBeads, string nameOf){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    string residue_index;

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    int residueCnt = 1;
    for (int i=startIndex; i<workingLimit; i++){
        currentBead = this->getBead(selectedBeads[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(residueCnt);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, " CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1," CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        residueCnt++;
    }

    fclose(pFile);

}


void Model::writeModelToFile(int workingLimit, vector<int> &selectedBeads, string nameOf, int steps){
    FILE * pFile;

    //const char *outputFileName;
    nameOf = nameOf + ".pdb";
    //const char * outputFileName = nameOf.c_str() ;
    pFile = fopen(nameOf.c_str(), "w");

    Bead * currentBead;
    std::string residue_index;

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );

    for (int i=0; i < workingLimit; i++){
        currentBead = this->getBead(selectedBeads[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(i + 1);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1," CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fclose(pFile);
}


std::string Model::writeModelToFile2(float dkl, int workingNumber, std::vector<int> &selectedBeads, std::vector<int> &pofrModel, std::string nameOf, Anneal *annealedObject, Data *pData, int steps, float volume, float averageContacts){
    FILE * pFile;

    const char *outputFileName;
    nameOf = nameOf + ".pdb";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    Bead * currentBead;
    std::string residue_index;

    std::string temp = createHeader(dkl, annealedObject, pData, steps, workingNumber, volume, averageContacts);
    fprintf(pFile, temp.c_str());

    // Add P(r) distributions
    // final D_kl, volume and energy
    float totalCounts = 0.0;
    float prob;
    int totalm = pofrModel.size();

    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    for (int i=0; i<totalm; i++){
        totalCounts += pofrModel[i];
    }

    int shannon_bins = annealedObject->gettotalBinsDerivedFromData();
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0

    float binwidth = pData->getBinWidth(), r_value;
    for (int i=0; i < shannon_bins; i++){
        prob = pData->getProbabilityPerBin(i);  // bounded by experimental Shannon Number
        r_value = i*binwidth + 0.5*binwidth;
        fprintf(pFile, "REMARK 265  %8.3f  %.4E  %.4E  %-.4E\n", r_value, prob, (pofrModel[i]/totalCounts), (prob * log(prob/pofrModel[i]*totalCounts)));
    }

    fprintf(pFile, "REMARK 265\n");
    // write coordinates
    for (int i=0; i<workingNumber; i++){
        currentBead = this->getBead(selectedBeads[i]);
        //residue_index = std::to_string(selectedBeads[i]);
        residue_index = std::to_string(i + 1);
        //residue_index = std::to_string(selectedBeads[i]);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", i+1," CA ", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    }

    fprintf(pFile, "END\n");
    fclose(pFile);
    return nameOf;
}


void Model::updateBeadIndices(int workingLimit, int deadlimit, std::vector<int> &indices) {

    this->workingLimit = workingLimit;
    this->deadLimit = deadlimit;

    std::copy(indices.begin(), indices.end(), bead_indices.begin());
}


void Model::setModelBeginStartingSetIterator(std::vector<int>::iterator &it){
    it = starting_set.begin();
}

void Model::setModelEndStartingSetIterator(std::vector<int>::iterator &it){
    it = starting_set.end();
}

void Model::setModelBeginIterator(std::vector<int>::iterator &it){
    it = bead_indices.begin();
}

void Model::setModelEndIterator(std::vector<int>::iterator &it){
    it = bead_indices.end();
}


/**
 * transform subunit contained within coordinates limited to workingLimit
 */
void Model::transformCoordinatesBySymmetry(int subunitIndex, int workingLimit, int &startCount, std::vector<vector3> &coordinates){
    vector3 * transformed, * nextVec;
    // if subunit index > sym_index, implies D_n, otherwise, its all rotations

    float angle = M_PI*2.0*subunitIndex/(float)symmetryIndex;
    float cosine = cos(angle);
    float sine = sin(angle);
    float x, y, z;

    if (subunitIndex < symmetryIndex){ // rotations for C and D symmetry groups

        for(int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
            x = (*transformed).x;
            y = (*transformed).y;
            z = (*transformed).z;

            nextVec = &coordinates[startCount];
            (*nextVec).x = cosine*x - sine*y;
            (*nextVec).y = sine*x + cosine*y;
            (*nextVec).z = z;

            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }

    } else { // reflections for the D_group
        // rotate on x-axis first by 180 and then by
        // rotating on x-axis (x,y,z) -> (x,-y,-z)
        for(int i=0; i<workingLimit; i++){
            transformed = &coordinates[i];
            x = (*transformed).x;
            y = -(*transformed).y;
            z = -(*transformed).z;

            nextVec = &coordinates[startCount];
//            (*nextVec).x = cosine*x + sine*y;
//            (*nextVec).y = sine*x - cosine*y;
            (*nextVec).x = cosine*x - sine*y;
            (*nextVec).y = sine*x + cosine*y;
            (*nextVec).z = z;
            //cout << startCount << " " << angle << " SUBUNIT " << subunitIndex << " " << (*nextVec).x << " " << (*nextVec).y << " " << (*nextVec).z << endl;
            startCount++;
        }
    }
}


string Model::writeSymModelToFile(float dkl, int workingLimitS, vector<int> &beads, vector<int> &pofrModel, string name, Anneal *annealedObject, Data *pData, int totalSteps, float volume, float averageContacts) {

    char alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ;
    std::vector<char> alphabet( alpha, alpha+sizeof(alpha)-1 ) ;

    FILE * pFile;

    const char *outputFileName;
    std::string newName = name + ".pdb";
    outputFileName = newName.c_str() ;
    pFile = fopen(outputFileName, "w");

    std::string temp = createHeader(dkl, annealedObject, pData, totalSteps, workingLimitS, volume, averageContacts);
    fprintf(pFile, temp.c_str());

    std::string residue_index;

    int totalCoordinates = numberOfSubUnits*workingLimitS;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create sym partners and add to Pr
    // calculate P(r) for subunit
    for (int i=0; i<workingLimitS; i++){
        tempBead = this->getBead(beads[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    int count = workingLimitS;
    for (int s=1; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector
        this->transformCoordinatesBySymmetry(s, workingLimitS, count, coordinates);
    }

    string chain;
    int index=0;
    for (int s=0; s<numberOfSubUnits; s++){  // create sym related subunits and add to coordinates vector

        chain = alphabet[s];

        for (int i=0; i<workingLimitS; i++){
            // convert coordinate to
            residue_index = std::to_string(i+1);
            fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", i+1, "CA", "ALA", chain.c_str(), residue_index.c_str(), coordinates[index].x, coordinates[index].y, coordinates[index].z );
            index++;
        }
        fprintf(pFile, "TER\n");
    }
    fprintf(pFile, "END\n");
    fclose(pFile);
    return newName;
}

std::string Model::createHeader(float dkl, Anneal * annealedObject, Data *pData, int totalSteps, int workingNumber, float volume, float averageContacts){

    string pofrFileName = pData->getPofRFilename();
    char buffer[80];
    int cstring;

    string tempHeader = "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";

    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";
    tempHeader += "REMARK 265              RADIATION SOURCE : BENDING MAGNET\n";
    tempHeader += "REMARK 265             SYNCHROTRON (Y/N) : Y\n";
    tempHeader += "REMARK 265                      BEAMLINE : B21 DIAMOND LIGHT SOURCE (example)\n";
    tempHeader += "REMARK 265          TEMPERATURE (KELVIN) : 298 (example)\n";
    tempHeader += "REMARK 265                            PH : 6.9 (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 20 mM HEPES (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 150 mM KCl (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 2 mM TCEP (example)\n";
    tempHeader += "REMARK 265              BUFFER COMPONENT : 5 mM KNitrate (example)\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265       DATA REDUCTION SOFTWARE : SCATTER (v3.x)\n";
    tempHeader += "REMARK 265               SOFTWARE AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265        DATA ANALYSIS SOFTWARE : SCATTER (v3.x)\n";
    tempHeader += "REMARK 265               SOFTWARE AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 MODELING METHOD\n";
    tempHeader += "REMARK 265                      SOFTWARE : IKETAMA (eKay-TaMa) v 1.9\n";
    tempHeader += "REMARK 265                        METHOD : Cross-Entropy Minimization Density Modeling\n";
    tempHeader += "REMARK 265                        TARGET : REAL-SPACE P(r)-distribution\n";
    tempHeader += "REMARK 265                        AUTHOR : RP RAMBO\n";
    tempHeader += "REMARK 265 \n";

    cstring = snprintf(buffer, 80, "REMARK 265  EXPERIMENTAL REAL SPACE FILE : %s\n", pofrFileName.c_str());
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 SEARCH AND REFINEMENT CONSTRAINTS\n";


    cstring = snprintf(buffer, 80, "REMARK 265                   EDGE RADIUS : %.3f\n", this->bead_radius);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265    CONTACTS PER LATTICE POINT : %i\n", annealedObject->getContactsPerBead());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265      AVG PTS HIGH TEMP SEARCH : %.0f\n", this->beadAverage);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265                         STDEV : %.0f\n", this->beadStDev);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265             SEARCH SPACE DMAX : %.1f Angstrom\n", 2*limit);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265        EXPERIMENTAL P(r) DMAX : %i Angstrom\n", pData->getDmax());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265        EMPIRICAL POROD VOLUME : %i Angstrom^3\n", (int)pData->getVolume());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265   SEARCH SPACE SHANNON NUMBER : %i \n", annealedObject->getMaxBinUsedInAnnealing());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265   EXPERIMENTAL SHANNON NUMBER : %i\n", annealedObject->gettotalBinsDerivedFromData());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265                      SYMMETRY : %s\n", symmetry.c_str());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265             TEMP RANGE (HIGH) : %.4E \n", annealedObject->getHighTempStartForCooling());
    tempHeader.append(buffer);
//    cstring = snprintf(buffer, 80, "REMARK 265              TEMP RANGE (LOW) : %.4E\n", annealedObject->getLowTempStop());
//    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265       TOTAL TEMPERATURE STEPS : %i\n", totalSteps);
    tempHeader.append(buffer);
//    cstring = snprintf(buffer, 80, "REMARK 265        TEMP STEP SCALE FACTOR : %.4f\n", annealedObject->getStepFactor());
//    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265         LAMBDA (CONNECTIVITY) : %.4E\n", annealedObject->getLambda());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265              MU (COMPACTNESS) : %.4E\n", annealedObject->getMu());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265     BETA (LATTICE VIOLATIONS) : %.4E\n", annealedObject->getBeta());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265                 ETA (CONTACT) : %.4E\n", annealedObject->getEta());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265 ALPHA (COMPACTNESS POTENTIAL) : %.4E\n", annealedObject->getAlpha());
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265            PERCENT ADD REMOVE : %.2f\n", annealedObject->getPercentAddRemove());
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\n";

    tempHeader += "REMARK 265 REFINED VALUES\n";
    cstring = snprintf(buffer, 80, "REMARK 265 TOTAL LATTICE POINTS IN MODEL : %i \n", workingNumber);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265          TOTAL LATTICE VOLUME : %.0f Angstrom^3\n", workingNumber*bead_volume);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265               CVX HULL VOLUME : %i Angstrom^3\n", (int)volume);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265     AVERAGE CONTACTS PER BEAD : %.2f \n", averageContacts);
    tempHeader.append(buffer);
    cstring = snprintf(buffer, 80, "REMARK 265   KULLBACK LEIBLER DIVERGENCE : %.5E\n", dkl);
    tempHeader.append(buffer);

    tempHeader += "REMARK 265\nREMARK 265\tFollowing table describes model-data agreement\nREMARK 265\t<r> is the average r-value per Shannon bin\n";
    tempHeader += "REMARK 265\nREMARK 265      <r>\t   P(r)_OBS   P(r)_MODEL   CROSS-ENTROPY\n";

//cout << tempHeader << endl;
    return tempHeader;
}


//void Model::estimatePointsPerComponent(float value){
//
//    for(std::vector<Component>::iterator it = components.begin(); it != components.end(); ++it) {
//        Component temp = *it;
//        cout << temp.getID() << " SETTING COMPONENT POINTS " << value << " SIZE: " << components.size() << endl;
//        temp.setTargetNumberOfLatticePoints(value);
//    }
//}


/**
 * returns sorted indices of lattice positions that match input PDB model
 */
void Model::createSeedFromPDB(string filename, Data * pData, int totalBins, std::vector<double> * pdbPr){

    PDBModel pdbModel(filename, true, true, this->bead_radius); // coordinates are centered

    const int totalBeads = this->getTotalNumberOfBeadsInUniverse();

    Bead * currentBead;
    float diffx, diffy, diffz, bx, by, bz;
    const int totalAtoms = pdbModel.getTotalAtoms();

    float beadradius = this->getBeadRadius();
    float b2 = beadradius*beadradius*1.000037;


    std::vector<int> indicesToCheck(totalAtoms);

    int * ptr = (totalAtoms != 0) ? &indicesToCheck.front() : NULL;

    for(int i = 0; i < totalAtoms; i++) {
        ptr[i] = i;
    }
    int limit = totalAtoms;

    cout << "** CONVERT TO LATTICE MODEL => please wait" << endl;
    int tempCnt=0;
    seed_indices.resize(totalBeads);
    for(int i=0; i < totalBeads && limit>0; i++){ // iterate over each bead in Universe

        currentBead = this->getBead(i);
        bx = currentBead->getX();
        by = currentBead->getY();
        bz = currentBead->getZ();

        for (int t=0; t < limit; t++){

            diffx = bx - *(pdbModel.getCenteredX() + ptr[t]);
            diffy = by - *(pdbModel.getCenteredY() + ptr[t]);
            diffz = bz - *(pdbModel.getCenteredZ() + ptr[t]);

            if ((diffx*diffx + diffy*diffy + diffz*diffz) <= b2){
                seed_indices[tempCnt] = i; // if lattice point is within distance, push index into Array and break
                //limit--;
                //iter_swap(indicesToCheck.begin()+t, indicesToCheck.begin()+limit);
                tempCnt++;
                break;
            }
        }
    }

    seed_indices.resize(tempCnt);

    // make PDB P(r) for refining beadmodel
    double sum=0.0;
    for (int i=0; i<totalAtoms; i++){

        bx = *(pdbModel.getCenteredX()+i);
        by = *(pdbModel.getCenteredY()+i);
        bz = *(pdbModel.getCenteredZ()+i);

        int next = i+1;
        for (int j=next; j<totalAtoms; j++){
            diffx = bx - *(pdbModel.getCenteredX()+j);
            diffy = by - *(pdbModel.getCenteredY()+j);
            diffz = bz - *(pdbModel.getCenteredZ()+j);
            float dis = sqrt((diffx*diffx + diffy*diffy + diffz*diffz));
            // convert to bin
            (*pdbPr)[pData->convertToBin(dis)]++;
            sum += 1.0;
        }
    }

    //normalize
    double invSum = 1.0/sum;
    for (int i=0; i<totalBins; i++){
        cout << i << " " << (*pdbPr)[i] << endl;
        (*pdbPr)[i] *= invSum;
    }
    // plot should integrate to 1

    std::sort(seed_indices.begin(), seed_indices.end());
    total_in_seed = seed_indices.size();
    //pModel->printSelectedBeads(0, totalKept, *keptBeads);
    pdbModel.writeCenteredCoordinatesToFile("centeredPDBSeed");
    this->writeSubModelToFile(0, total_in_seed, seed_indices, "centeredLatticeModel");
}

void Model::setReducedSeed(int limit, std::vector<int> &reduced){
    //reduced_seed.resize(limit);
    total_in_reduced_seed = limit;

    //for(int i=0; i<limit; i++){
    //    reduced_seed[i] = reduced[i];
    //}
    reduced_seed.clear();
    for(int i=0; i<limit; i++){
        reduced_seed.insert(reduced[i]);
    }
}

/**
 * Search reduced seed model to see if bead index is within selected set of indices
 */
bool Model::inReducedSeed(int index){

    if (reduced_seed.find(index) != reduced_seed.end()){
        return true;
    }
    return false;
}

/**
 * center the input set of indices, reassigning the bead indices
 *
 */
void Model::centerLatticeModel(int workingLimit, std::vector<int> & indices){

    vector3 totalVec;// = this->getBead(indices[0])->getVec();

    for(int i=0; i<workingLimit; i++){
        totalVec += this->getBead(indices[i])->getVec();
    }

    float invTotal = 1.0/(float)workingLimit;
    totalVec.x = totalVec.x*invTotal;
    totalVec.y = totalVec.y*invTotal;
    totalVec.z = totalVec.z*invTotal;

    std::vector<int> temp_bead_indices(bead_indices.begin(), bead_indices.end());
    std::vector<int>::iterator tempStartIt = temp_bead_indices.begin();

    int searchLimit = 0, selectedIndex;
    float distance, mindistance;

    for(int i=0; i<workingLimit; i++){

        this->getBead(indices[i])->getVec();
        vector3 newVec = totalVec - this->getBead(indices[i])->getVec(); // translate the bead
        // determine which bead is closest to this new vector

        selectedIndex=0;
        mindistance = 100;

        for (int j=searchLimit; j<number_of_beads; j++){

            distance = (newVec - this->getBead(temp_bead_indices[j])->getVec()).length();

            if (distance < mindistance){ // find the bead that is closest to translated vector
                selectedIndex = j;
                mindistance = distance;
            }
        }

        std::iter_swap(tempStartIt+ searchLimit, tempStartIt+selectedIndex);
        searchLimit++;
        // how does it handle two beads that map to same location?
    }

    std::copy(tempStartIt, temp_bead_indices.end(), indices.begin());
    this->writeModelToFile(workingLimit, temp_bead_indices, "centered_model", 0);
}

// return pointer to component
//Component * Model::getComponentByIndex(int index){
//        return &(components[index]);
//}
//
//int Model::getEstimatedLatticePointsPerComponentByIndex(int index){
//    Component * temp = &(components[index]);
//    return temp->getTargetNumberOfLatticePoints();
//}
