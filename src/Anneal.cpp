//
// Created by Robert Rambo on 13/01/2017.
//

#include "Anneal.h"
#include "Model.h"
#include "PDBModel.h"

using namespace std;

Anneal::Anneal(float highT,
               float percent,
               int lowerV,
               int upperV,
               int highTempRounds,
               float contactsPerBead,
               string fileprefix,
               int totalSteps,
               float eta,
               float lambda,
               float alpha,
               int multiple) {

    this->highT = highT;
    this->percentAddRemove = percent;
    this->lowerV = lowerV;
    this->upperV = upperV;
    this->highTempRounds = highTempRounds;
    this->stepsPerTemp = 5;
    this->contactsPerBead = contactsPerBead;
    this->totalTempStop = totalSteps;

    this->highTempStartForCooling = 0.00001; //0.00001
    this->filenameprefix = fileprefix;

//    this->numberOfCoolingTempSteps = (int) ceil(log(lowTempStop/highTempStartForCooling)/(log(expSlowCoolConstant)));
    this->eta = eta;
    this->lambda = lambda;
    this->alpha = alpha;
    this->ccmultiple = multiple;

    // expansionSlope = (stepsPerTemp - 5)/numberOfCoolingTempSteps;
    // printf("  SETTING EXPANSION SLOPE => %.2f (%.0f -> %.0f)\n", expansionSlope, highT, lowTempStop);
}


/**
 * Random add/remove to generate initial model for seeded modeling
 * Try to find the minimial set of lattice points that agree with atomistic P(r)
 * contactCutoff determines
 */
//void Anneal::createSeedFromPDB(Model *pModel, Data *pData, string name, string PDBFilename, int numberOfUniqueConnectedPhases){
//
//    totalNumberOfPhasesForSeededModeling = numberOfUniqueConnectedPhases;
//    this->lowTempStop = highTempStartForCooling;
//    contactCutOff = interconnectivityCutOff;
//
//    populatePotential(pModel->getSizeOfNeighborhood());
//
//    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
//    float * pDistance = pModel->getPointerToDistance();
//    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class
//
//    // convert distances within the large search space to bins based on input P(R)-DATA file
//    this->fillPrBinsAndAssignTotalBin( pBin,  pDistance,  totalDistancesInSphere,  pData);
//
//    std::vector<float> prPDB(maxbin);
//    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
//
//    // create working observed probability distribution that encompasses search sphere
//    pData->createWorkingDistribution(maxbin);
//
//    const std::vector<int>::const_iterator trueModelBeginIt = pModel->getSeedBegin();
//    //const std::vector<int>::const_iterator trueModelEndIt = pModel->getSeedEnd();
//    int workingLimit = pModel->getTotalInSeed();
//
//    // copy trueModel into BeadIndices
//    // number of shannon bins for the model is calculated over the Universe (not the data)
//    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
//    std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
//    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50
//
//    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
//    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
//    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
//    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;
//
//    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
//    int minWorkingLimit = 0.17*workingLimit;
//
//    int deadLimit = workingLimit;; // as bead indices are discarded, set upper limit of vector
//    std::vector<int> bead_indices(workingLimit); // large vector ~1000's
//    std::vector<int> lowest_bead_indices(workingLimit); // large vector ~1000's
//    std::vector<int> active_indices(workingLimit); // large vector ~1000's
//    std::vector<int> backUpState(workingLimit);
//
//    std::clock_t start;
//    // c-style is slightly faster for large vector sizes
//    const int num = workingLimit;
//    int * ptr = (num != 0) ? &bead_indices.front() : NULL;
//    for(int i = 0; i < num; i++) {
//        ptr[i] = *(trueModelBeginIt+i);
//    }
//    // prepare bead_indices by copying in truemodel
//    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
//    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
//
//    std::random_device rd;
//    std::mt19937 gen(rd());
//
//    float inv_kb_temp = 1.0/0.000001;
//
//    int lowerN = 0.1*workingLimit, upperN = workingLimit;
//    cout << " TOTAL LATTICE IN SEED : " << workingLimit << endl;
//    cout << "        LATTICE LIMITS : " << lowerN << " <= N <= " << upperN << endl;
//    // randomize and take the workingLength as first set, shuffling takes about 10x longer than copy and sort
//
//    int tempNumberOfComponents, currentComponents;
//    bool isConnected;
//    // bead_indices contains only the indices that relate to the input PDB
//    isConnected = isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, currentComponents);
//
//    cout << " CHECKING CONNECTIVITY  " << endl;
//    cout << "        IS CONNECTED ? : " << isConnected << endl;
//    cout << "  NUMBER OF COMPONENTS : " << currentComponents << endl;
//
//    const int components = tempNumberOfComponents;
//    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
//    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
//    std::vector<int>::iterator endBinCount = binCount.end();
//
//    // calculate KL against known
//    // create custom D_KL function supplying Pr_of_PDB as target
//    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
//    currentKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
//
//    //float tempTotalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
//
//    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
//    cout << "          INITIAL D_KL : " << currentKL << endl;
//
//    float energy_prior, testKL;
//    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
//
//    int original, addMe;
//    std::uniform_real_distribution<float> distribution(0.0,1.0);
//
//    std::vector<int>::iterator beginIt = bead_indices.begin();
//    std::vector<int>::iterator endIt = bead_indices.end();
//
//    float lowestE = currentKL;
//    int lowestWorkingLimit, lowestDeadLimit = deadLimit;
//    lowestWorkingLimit = workingLimit;
//    std::copy(beginIt, endIt, lowest_bead_indices.begin());
//
//    int counter=1;
//    int priorWorkingLimit;
//
//    float inv_divideBy, average_x=0.5*workingLimit, stdev=0.2*workingLimit;
//    float sum_x_squared=0, sum_x=0, divideBy=0, acceptRate = 0.5, inv500 = 1.0/500.0;
//    //output for plotting
//    bool isUpdated = false;
//
//    int seedHighTempRounds = 2*highTempRounds;
//    //int seedHighTempRounds = 5000;
//
//    float startContactAVG=0.0;
//    for (int i=0; i<workingLimit; i++){
//        startContactAVG += (float)numberOfContacts(bead_indices[i], &bead_indices, workingLimit, pModel, pDistance);
//    }
//
//    double startContactsPotential = calculateTotalContactSum(&beads_in_use_tree, workingLimit, pModel)/(float)workingLimit;
//
//    startContactAVG *= 1.0/(float)workingLimit;
//
//    float startKL = currentKL;
//
//    int high;
//    for (high=0; high < seedHighTempRounds; high++){ // iterations during the high temp search
//
//        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
//        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
//        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
//
//        if (distribution(gen) >= 0.5){ // ADD BEAD?
//            cout << "*******************                  ADD                   *******************" << endl;
//
//            printf("     ADDME => %i \n", 1);
//            if (workingLimit < deadLimit){
//                // find a bead within workinglimit -> deadLimit
//                double afterAdding;
//                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
//                addMe = bead_indices[randomSpot];
//
//                itIndex = bead_indices.begin() + randomSpot;
//                //itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.begin() + deadLimit, addMe);
//                //make the swap at the border (workingLimit)
//                std::iter_swap(bead_indices.begin() + workingLimit, itIndex);
//                workingLimit++;
//                //std::copy(bead_indices.begin(), bead_indices.begin() + workingLimit, backUpState.begin());
//                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
//
//                addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
//
//                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
//                // since neighbor list is already determined, don't need to check for contacts
//                if (testKL < currentKL && numberOfContacts(addMe, &bead_indices, workingLimit, pModel, pDistance) >= 1) {
//                    currentKL = testKL;
//                    currentComponents = tempNumberOfComponents;
//                    isUpdated = true;
//                } else if (exp(-(testKL - currentKL) * inv_kb_temp) > distribution(gen) && numberOfContacts(addMe, &bead_indices, workingLimit, pModel, pDistance) >= 1) {
//                    currentKL = testKL;
//                    currentComponents = tempNumberOfComponents;
//                    isUpdated = true;
//                } else { // undo changes (rejecting)
//                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
//                    workingLimit--;
//                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//                }
//            }
//
//        } else { // REMOVE BEADS?
//            cout << "*******************                 REMOVE                 *******************" << endl;
//            // test for deletion
//            priorWorkingLimit = 1;
//            printf("     REMOVE => %i\n", priorWorkingLimit);
//
//            // collect indices that are not part of PDBModel
//            int randomSpot = rand() % workingLimit;
//            original = bead_indices[randomSpot];
//
//            itIndex = bead_indices.begin() + randomSpot;
//            //itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, original);
//            // remove original from P(r)
//            std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); //backup binCount
//            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
//            removeFromPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
//            // make the swap
//            workingLimit--;
//            std::iter_swap(itIndex, bead_indices.begin() + workingLimit);
//            // still need to sort, swap changes the order
//            std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
//            isConnected = isConnectedComponent(&bead_indices, workingLimit, pDistance, pModel->getTotalNumberOfBeadsInUniverse(), tempNumberOfComponents);
//            // don't remove lattice point if model is not interconnected
//            //isConnected = true;
//            if (isConnected){
//                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
//                if (testKL < currentKL) {
//                    currentKL = testKL;
//                    currentComponents = tempNumberOfComponents;
//                    isUpdated = true;
//                } else if (exp(-(testKL - currentKL)*inv_kb_temp) > distribution(gen)){
//                    currentKL = testKL;
//                    currentComponents = tempNumberOfComponents;
//                    isUpdated = true;
//                } else { // undo changes and move to next bead (rejecting)
//                    workingLimit++;
//                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
//                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//                }
//            } else { // undo changes and move to next bead (rejecting)
//                workingLimit++;
//                //sort(beginIt, beginIt+workingLimit); // swapped index is at workingLimit
//                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
//                std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//            }
//
//        }
//
//        if (currentKL < lowestE){
//            lowestE = currentKL;
//            lowestWorkingLimit = workingLimit;
//            lowestDeadLimit = deadLimit;
//            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
//            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_");
//        }
//
//        sum_x_squared += workingLimit*workingLimit;
//        sum_x += workingLimit;
//        divideBy += 1.0;
//
//        counter++;
//
//        cout << "*******************                                        *******************" << endl;
//        printf("       TEMP => %-.8f \n     INVKBT => %.4f\n", lowTempStop, inv_kb_temp);
//        printf("   MAXSTEPS => %i (%4i) \n", seedHighTempRounds, high);
//        printf("      GRAPH => %i \n", currentComponents);
//        printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", average_x, stdev);
//        printf("LIMIT: %5i (>= MIN: %i) DEADLIMIT: %5i D_KL: %.4E \n", workingLimit, minWorkingLimit, deadLimit, currentKL);
//
//
//        if (isUpdated){
//            acceptRate = inv500*(499*acceptRate+1);
//            isUpdated = false;
//        } else {
//            acceptRate = inv500*(499*acceptRate);
//        }
//
//        updateASATemp(high, seedHighTempRounds, acceptRate, lowTempStop, inv_kb_temp);
//    } // end of HIGH TEMP EQUILIBRATION
//
//    // average number of contacts per bead
//
//    float contactSum=0.0;
//    for (int i=0; i<lowestWorkingLimit; i++){
//        contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, pModel, pDistance);
//    }
//    std::set<int> lowest_beads_in_use_tree(lowest_bead_indices.begin(), lowest_bead_indices.begin() + lowestWorkingLimit);
//    double runningContactsSum = calculateTotalContactSum(&lowest_beads_in_use_tree, lowestWorkingLimit, pModel);
//
//
//    cout << "KL DIVERGENCE : " << endl;
//    cout << "  INITIAL D_KL => " << startKL << endl;
//    cout << "   LOWEST D_KL => " << lowestE << endl;
//    cout << "AVERAGE CONTACTS : (per lattice point)" << endl;
//    cout << "       INITIAL => " << startContactAVG << " ENERGY => " << startContactsPotential << endl;
//    cout << "         FINAL => " << contactSum/(float)lowestWorkingLimit << " ENERGY => " << runningContactsSum/(float)lowestWorkingLimit  << endl;
//
//    cout << " Contacts Per Bead " << contactsPerBead << endl;
//
//    // this is fixed model for initial high temp search?
//    // set this as the seed
//    //    pModel->setBeadAverageAndStdev(average_x, stdev);
//    //    pModel->setReducedSeed(workingLimit, bead_indices);
//    //    pModel->writeModelToFile(workingLimit, bead_indices, name);
//    //    pModel->setBeadAverageAndStdev(average_x, stdev);
//
//    // using lowestEnergy model, organize beads of the entire Model search space.
//    bead_indices.resize(totalBeadsInSphere);
//    beginIt = bead_indices.begin();
//    endIt = bead_indices.end();
//
//    ptr = &bead_indices.front();
//    for(int i = 0; i < totalBeadsInSphere; i++) {
//        ptr[i] = i;
//    }
//
//    std::vector<int>::iterator it;
//    for (int i=0; i<lowestWorkingLimit; i++){
//        it = std::find(beginIt, endIt, lowest_bead_indices[i]); // if itTrueIndex == endTrue, it means point is not within the set
//        std::iter_swap(beginIt+i, it);
//    }
//
//    //pModel->setStartingSet(lowest_bead_indices);
//    pModel->setStartingSet(bead_indices);
//
//    pModel->setStartingWorkingLimit(lowestWorkingLimit);
//    pModel->setStartingDeadLimit(lowestDeadLimit);
//
//    cout << "LOWEST WORKING LIMIT : " << lowestWorkingLimit << endl;
//    // set seed model to be invariant during reconstruction
//    pModel->setReducedSeed(lowestWorkingLimit, lowest_bead_indices);
//    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);
//}

float Anneal::calculateAverageDistance(float * pDistance, int *stopAt, vector<int> *bead_indices, Model * pModel){

    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    int row, row2, count = 0;
    float distanceSum=0.0;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    for(int m=0; m < *stopAt; m++){
        row = (*bead_indices)[m];
        row2 = row * (totalBeads) - row*(row+1)*0.5 - row - 1;

        // parallel
        for(int n=m+1; n < *stopAt; n++){
            distanceSum +=  *(pDistance + row2 + (*bead_indices)[n]);
            count++;
        }
    }
    return distanceSum/count;
}


float Anneal::calculateCVXHULLVolume(char *flags, std::vector<int> *bead_indices, int upTo, coordT * points, Model *pModel) {

    // points should already be set to 3*upTo
//    int dim = 3;
//    char flagCVX[25];
//    sprintf(flagCVX, "qhull s FA");
    //coordT pointSet[3*upTo];

    for (int i=0; i<upTo; i++){
        beadToPoint(&points[i*3], pModel->getBead((*bead_indices)[i]));
        //  beadToPoint(&(pointSet[i*3]), pModel->getBead((*bead_indices)[i]));
    }

    qh_new_qhull(3, upTo, points, 0, flags, NULL, NULL);

    int volume_test = qh totvol;
    qh_freeqhull(true);
    return volume_test;
}


/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSum(std::set<int> *beads_in_use, Model *pModel){

    //calculate contacts per bead for selectedIndex
    //int limit = workingLimit;

    double sum=0;
//    for(int i=0; i<workingLimit; i++){
//        //sum += totalContactsPotential( numberOfContacts((*bead_indices)[i], bead_indices, limit, pModel, pDistance) );
//        sum += totalContactsPotential(numberOfContactsFromSet(beads_in_use, pModel,(*bead_indices)[i]));
//        //sum += numberOfContacts((*bead_indices)[i], bead_indices, limit, pModel, pDistance);
//    }

    std::set<int>::iterator it;
    for (it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {
        //int test = numberOfContactsFromSet(beads_in_use, pModel, *it);
        sum += totalContactsPotential(numberOfContactsFromSet(beads_in_use, pModel, *it));
    }

    return sum;
}


/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSumComponent(std::set<int> *beads_in_use, std::set<int> * beads_in_component, int const workingLimit,
                                        Model *pModel){

    //calculate contacts per bead for selectedIndex
    double sum=0;

    std::set<int>::iterator it;
    for (it = beads_in_component->begin(); it != beads_in_component->end(); ++it) {
        //int test = numberOfContactsFromSet(beads_in_use, pModel, *it);
        sum += totalContactsPotential(numberOfContactsFromSet(beads_in_use, pModel, *it));
    }

    return sum;
}



// need to return the selected index
void Anneal::createPlacesToCheck(int workingLimit,
                                          int average_number_of_contacts,
                                          int swap1,
                                          std::set<int> * beads_in_use,
                                          std::set<int> * returnMe,
                                          Model * pModel){

    std::vector<int>::iterator primaryNeighborhood, secondaryNeighborhood;
    int secondaryN;
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    // select first element from randomized active_set
    // check if available neighbor can be added
    primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);

    for (int n=0; n < totalNeighbors; n++){ // run into problem of check
        int neighbor = *(primaryNeighborhood + n);

        if (neighbor > -1){
            set<int>::iterator inSet = beads_in_use->find(neighbor);

            if (inSet == beads_in_use->end()){ // not in use, so its in deadlimit
                // if number of contacts at new position (less current) is 0, skip
                returnMe->insert(neighbor);

            } else { // its already a neighbor within WorkSet so check its neighborhood
                secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);

                for (int s=0; s < totalNeighbors; s++){
                    // adding to any of these already insures at least 1 contact
                    // check if secondary neighbor is in beads_in_use
                    secondaryN = *(secondaryNeighborhood + s);
                    inSet = beads_in_use->find(secondaryN);
                    if (secondaryN > -1 && inSet == beads_in_use->end()){ // if .end(), means not in use
                        returnMe->insert(secondaryN);
                    } else if (secondaryN == -1){
                        break;
                    }
                }
            }
        } else if (neighbor == -1) {
            break;
        }
    }


//    for (int i=0; i < workingLimit; i++) {
//
//        swap1 = *(activeIt + i);
//        //swap1 = *active_indices[i];
//
//        numberContacts = numberOfContactsFromSet(beads_in_use, pModel, swap1);
//        if (numberContacts < average_number_of_contacts){
//            // select first element from randomized active_set
//            // check if available neighbor can be added
//            primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);
//
//            for (int n=0; n < totalNeighbors; n++){ // run into problem of check
//                int neighbor = *(primaryNeighborhood + n);
//
//                if (neighbor > -1){
//                    set<int>::iterator inSet = beads_in_use->find(neighbor);
//
//                    if (inSet == beads_in_use->end()){ // not in use, so its in deadlimit
//                        // if number of contacts at new position (less current) is 0, skip
//                        placesToCheck.insert(neighbor);
//
//                    } else { // its already a neighbor within WorkSet so check its neighborhood
//                        secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);
//
//                        for (int s=0; s < totalNeighbors; s++){
//                            // adding to any of these already insures at least 1 contact
//                            // check if secondary neighbor is in beads_in_use
//                            secondaryN = *(secondaryNeighborhood + s);
//                            inSet = beads_in_use->find(secondaryN);
//                            if (secondaryN > -1 && inSet == beads_in_use->end()){ // if .end(), means not in use
//                                placesToCheck.insert(secondaryN);
//                            } else if (secondaryN == -1){
//                                break;
//                            }
//                        }
//                    }
//                } else if (neighbor == -1) {
//                    break;
//                }
//            }
//            break;
//        }
//    }
}

bool Anneal::checkForRepeats(std::vector<int> beads) {
    bool state = false;
    int beadSize = beads.size();

    std::set<int> testSet(beads.begin(), beads.end());
    cout << "______________________________________________________________________________" << endl;
    cout << "*******************                 TEST                   *******************" << endl;
    cout << "*******************              -----------               *******************" << endl;
    cout << " TEST SET " << testSet.size() << " vector set " << beadSize << endl;
    if (testSet.size() != beadSize){
        for(int i=0; i<10; i++){
            cout << "                 !!!!!!!!!!!!!!!!DUPLICATE ENTRIES FOUND! " << testSet.size() << " " << beadSize <<  endl;
        }
        state = true;
    }
    return state;
}

/*!
 *
 */
float Anneal::connectivityPotential(int numberOfComponents){

    switch(numberOfComponents) {
        case 1:
            return 0.f;
        case 2:
            return 100.f;
        case 3:
            return 1000.f;
        case 4:
            return 10000.f;
        case 5:
            return 100000.f;
        case 6:
            return 1000000.f;
        default:
            return 1000000.f*numberOfComponents;
    }
}

/*
 * go through each vertex of hull and move one unit further away
 */
void Anneal::enlargeDeadLimit(std::vector<int> &vertexIndices,
                              int totalV,
                              std::vector<int> &bead_indices,
                              int workingLimit,
                              int *deadLimit,
                              int totalBeadsInSphere,
                              Model *pModel) {

    int * pIndexOfHullpt;
    vector<int>::iterator beginIt = bead_indices.begin(), it;
    vector<int>::iterator endIt = bead_indices.end();

    float sum_x=0.0, sum_y=0.0, sum_z=0.0, countf=0.0;
    Bead * currentBead;
    // calculate center of the hull
    for (int v = 0; v < totalV; v++) {
        currentBead = pModel->getBead(vertexIndices[v]);
        sum_x += currentBead->getX();
        sum_y += currentBead->getY();
        sum_z += currentBead->getZ();
        countf += 1.0;
    }

    float ave_x, ave_y, ave_z, distance_from_center, temp_distance, diff_x, diff_y, diff_z;
    // calculate center of HULL based on vertices
    float inv_count = 1.0/countf;
    ave_x = sum_x*inv_count;
    ave_y = sum_y*inv_count;
    ave_z = sum_z*inv_count;
    // end calculate center of hull

    std::vector<int> neighbors(60);
    int totalNeighbors, replaceWith, distanceTo;
    int modifiedWorkingLImit = workingLimit;

    // Search region between workingLimit and DeadLimit for a better point
    for (int v = 0; v < totalV; v++) { // SWAP VERTICES OF HULL TO BETTER SPOT

        pIndexOfHullpt = &vertexIndices[v];
        // find all contact neighboring beads
        // for through each one and accept first one that increase hull volume
        pModel->getNeighbors(*pIndexOfHullpt, neighbors, totalNeighbors);

        currentBead = pModel->getBead(*pIndexOfHullpt);

        diff_x = currentBead->getX()-ave_x;
        diff_y = currentBead->getY()-ave_y;
        diff_z = currentBead->getZ()-ave_z;
        distance_from_center = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z; // current distance from center of HULL

        replaceWith = *pIndexOfHullpt; // if nothing greater use same vertex

        // PARALLELIZE
        for (int i=0; i < totalNeighbors; i++){ // not all beads are
            currentBead = pModel->getBead(neighbors[i]);
            diff_x = currentBead->getX()-ave_x;
            diff_y = currentBead->getY()-ave_y;
            diff_z = currentBead->getZ()-ave_z;
            temp_distance = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z;

            if (temp_distance > distance_from_center){
                // add bead to workinglist
                distance_from_center = temp_distance;
                replaceWith = neighbors[i];
            }
        }

        it = find(beginIt, endIt, replaceWith);
        distanceTo = distance(beginIt, it);

        if ( distanceTo >= *deadLimit){ // outside feasible region
            std::iter_swap(beginIt + *deadLimit, it);
            std::iter_swap(beginIt + *deadLimit, beginIt + modifiedWorkingLImit);
            (*deadLimit)++;
            modifiedWorkingLImit++;
        } else if (distanceTo >= workingLimit) { // point is within Feasible region
            std::iter_swap(it, beginIt + modifiedWorkingLImit);
            modifiedWorkingLImit++;
        }
    }

    // add these points to region within workinglimit and deadlimit (working zone)
    // add points, increase working limit
    *deadLimit = recalculateDeadLimit(modifiedWorkingLImit, bead_indices, pModel, totalBeadsInSphere);
    // calculate convex HULL with this set of points
}

/**
 * has to be sorted
 * &nonSeed, nonSeedCount, pDistances, totalBeads, numberOfComponentsToBeModeled
 */
bool Anneal::isConnectedComponent(std::vector<int> *activeIndices, int availableWorkingLimit, float *pDistances,
                                  int totalBeads, int &numberOfComponents) {

    if (!std::is_sorted(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit)){
        cout << " NOT SORTED FROM isConnectedComponent " << endl;
        return false;
    }

    bool test=false;
    std::vector<int> tempIndices(availableWorkingLimit);

    std::copy(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit, tempIndices.begin());
    std::sort(tempIndices.begin(), tempIndices.end());

    // int lastBeadIndex =tempIndices[availableWorkingLimit-1];
    //Graph graph(lastBeadIndex+1);  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    Graph graph;  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    boost::graph_traits<Graph>::edge_descriptor edge;
    bool flag;
    float distance;
    // boost::tie(edge, flag) = add_edge(0, 1, graph);  // tie is a reference based tuple
    // ds.union_set(0,1);
    // cout << " SIZE OF : " << availableWorkingLimit << endl;

    int count=0;

    // create edges within cutoff
    std::vector< Vertex> vertices(availableWorkingLimit); ////typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    for(int i=0; i < availableWorkingLimit; i++){
        vertices[i] = boost::add_vertex(graph);
    }

    Vertex * pOne, * pTwo;

    // create edges between points that are within cutoff
    for(int i=0; i < availableWorkingLimit; i++){

        int firstBead = tempIndices[i];
        pOne = &vertices[i];
        int next = i+1;

        for(int j=next; j<availableWorkingLimit; j++){
            int secondBead = tempIndices[j];
            // is secondBead a neighbor of firstBead
            // distance = *(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead);
            if ((*(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead)) <= interconnectivityCutOff){
                pTwo = &vertices[j];
                //cout << " j " << j << " " << secondBead << " " << distance << endl;
                //boost::tie(edge, flag) = boost::add_edge(firstBead, secondBead, graph);
                boost::tie(edge, flag) = boost::add_edge(*pOne, *pTwo, graph);
                //cout << " j " << j << " EDGE ADDED " << endl;
                //ds.union_set(firstBead, secondBead);
                count++;
            }

        }
    }

    std::vector<VertexIndex> rank(num_vertices(graph));
    std::vector<Vertex> parent(num_vertices(graph));

    boost::disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);

    initialize_incremental_components(graph, ds);
    incremental_components(graph, ds);

    boost::component_index<VertexIndex> components(parent.begin(), parent.end());

    // one single component should be
    // parent.size() - availableWorkingLimit + 1
    //cout << " Components Size " << components.size() << " Parent size " << parent.size() << " Working limit " << availableWorkingLimit << endl;
    numberOfComponents = components.size();
/*
    BOOST_FOREACH(VertexIndex current_index, components) {
                    std::cout << "component " << current_index << " contains: ";
                    // std::cout << boost::num_vertices(components[current_index].first) << endl;
                    // Iterate through the child vertex indices for [current_index]
                    BOOST_FOREACH(VertexIndex child_index,
                                  components[current_index]) {
                                    std::cout << child_index << " ";
                                }

                    std::cout << std::endl;
                }
*/
    // map index of vertex back to active indices via sorted temp indices
    // go through each compnent and redistribute to largest to make connected
    if (numberOfComponents == 1){
        test = true;
    }

    return test;
}

/**
 * contacts calculation must be performed on a sorted list and includes the beadIndex in the list
 */
int Anneal::numberOfContacts(int &beadIndex, vector<int> *bead_indices, int &workingLimit, Model *pModel, float * pDistance){

    int count=0;
    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();
    // beadsInUse must be sorted
    int row;
    unsigned long int row2;

    int i=0;
    // add down column (column is beadIndex
    while ((*bead_indices)[i] < beadIndex && i < workingLimit){
        row = (*bead_indices)[i];
        row2 = row*totalBeads - (row*(row+1)*0.5) - row - 1;
        if (*(pDistance + row2 + beadIndex) < contactCutOff) {
            count++;
        }
        i++;
    }
    // out of first loop, assume bead_indices[i] == beadIndex at this point
    i++;
    // Add across row (constant row)
    row2 = beadIndex*totalBeads - beadIndex*(beadIndex+1)*0.5 - beadIndex - 1;
    while (i < workingLimit){
        if (*(pDistance + row2 + (*bead_indices)[i]) < contactCutOff) {
            count++;
        }
        i++;
    }
    return count;
}



int Anneal::getRandomNeighbor(int &locale, std::set<int> *beads_in_use, Model * pModel){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(locale);
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<int> temp(totalNeighbors);
    std::set<int>::iterator endOfSet = beads_in_use->end();

    int count=0; // possible beads to use
    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if (beads_in_use->find(neighbor) == endOfSet && (neighbor > -1)){ // positions not in use should not be found
            temp[count] = neighbor;
            count++;
        } else if (neighbor == -1) {
            break;
        }
    }
    // if count == 0, no neighbors to use
    if (count == 0){
        return 0;
    }else {
        return (temp[rand()%count]);
    }
}


/**
 * contacts calculation must be performed on a sorted list and includes the beadIndex in the list
 */
int Anneal::numberOfContactsExclusive(int &beadIndex, int excludeIndex, vector<int> *bead_indices, int &workingLimit, Model *pModel, float * pDistance){

    int count=0;
    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();
    // beadsInUse must be sorted
    int row;
    unsigned long int row2;

    int i=0;
    // add down column (column is beadIndex
    while ((*bead_indices)[i] < beadIndex && i < workingLimit){
        row = (*bead_indices)[i];
        row2 = row*totalBeads - (row*(row+1)*0.5) - row - 1;
        if (*(pDistance + row2 + beadIndex) < contactCutOff) {
            count++;
        }
        i++;
    }

    // Add across row (constant row)
    row2 = beadIndex*totalBeads - beadIndex*(beadIndex+1)*0.5 - beadIndex - 1;
    while (i < workingLimit){
        if (*(pDistance + row2 + (*bead_indices)[i]) < contactCutOff) {
            count++;
        }
        i++;
    }
    return count;
}

/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within deadLimit
 */
void Anneal::populateLayeredDeadlimit(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit,
                                      int *pDeadLimit, Model *pModel, const int totalBeads) {
    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex;
    int distance;
    int neighbor;

    for (int i = 0; i<workingLimit; i++){

        it = pModel->getPointerToNeighborhood(*(iteratorBeadIndices + i));

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);

            itIndex = std::find(iteratorBeadIndices+(*pDeadLimit), iteratorBeadIndices + totalBeads, neighbor);
            distance = (itIndex - iteratorBeadIndices); // distance from beginning of vector
            // if not found, itIndex will report last

            if ((neighbor > -1) && (distance >= *pDeadLimit) && (distance < totalBeads)){
                std::iter_swap(iteratorBeadIndices + (*pDeadLimit), itIndex);
                (*pDeadLimit)++;
            } else if (neighbor == -1){
                break;
            }

//            if ( distance >= *pDeadLimit && (neighbor != -1) && (distance < totalBeads)) {
//                //cout << j << " SELECTED INDEX: " << *(iteratorBeadIndices + i) << " NEIGHHBOR: "  << neighbor << " " << distance << " DL " << *pDeadLimit << endl;
//                std::iter_swap(iteratorBeadIndices + (*pDeadLimit), itIndex);
//                (*pDeadLimit)++;
//            } else if (neighbor == -1) {
//                break;
//            }
        }
    }
}

/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within deadLimit
 */
void Anneal::populateLayeredDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, const int workingLimit,
                                              int * pDeadLimit, Model * pModel) {
    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex, iteratorBeadIndices = bead_indices->begin(), endIt = bead_indices->end();
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int neighbor;

    for (int i = 0; i < workingLimit; i++){ // excludes true model positions

        it = pModel->getPointerToNeighborhood((*bead_indices)[i]);

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);

            if ((neighbor > -1) && (beads_in_use->find(neighbor) == endOfSet)){ // if true, add neighbor to search space

                itIndex = std::find( (bead_indices->begin() + (*pDeadLimit)), endIt, neighbor);

                // if neighbor was already added from previous point, it will not be found
                if (itIndex != endIt){
                    std::iter_swap( (bead_indices->begin() + (*pDeadLimit)), itIndex);
                    (*pDeadLimit)++;
                }
            } else if (neighbor == -1) {
                break;
            }
        }
    }
}


int Anneal::recalculateDeadLimit(int workingLimit, vector<int> &bead_indices, Model * pModel, int totalBeadsInSphere ){

    int deadLimit = totalBeadsInSphere;

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[25];
    sprintf(flags, "qhull s FA");

    coordT hullPoints2[3*workingLimit];

    // can be threaded
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints2[i*3], pModel->getBead(bead_indices[i]));
    }

    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;

    // UPDATE DEADZONE : move points not selected that are outside of hull to deadZone
    std::vector<int>::iterator beginIt = bead_indices.begin();

    std::vector<int> inside(totalBeadsInSphere-workingLimit);
    std::vector<int> outside(totalBeadsInSphere-workingLimit);
    int insideCount=0, outsideCount=0;
    int * beadPosition;

    for(int i=workingLimit; i < deadLimit; i++){
        beadPosition = &bead_indices[i];
        beadToPoint(testPoint, pModel->getBead(*beadPosition));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (isoutside){
            outside[outsideCount] = *beadPosition;
            outsideCount++;
        } else {
            inside[insideCount] = *beadPosition;
            insideCount++;
        }
    }

    qh_freeqhull(true);

    std::copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
    deadLimit = workingLimit + insideCount;
    std::copy(outside.begin(), outside.begin()+outsideCount, beginIt + deadLimit);

    return deadLimit;
}


/*
 *
 */
void Anneal::populatedDeadLimitExcludingSet(std::vector<int> * bead_indices,
                                            std::set<int> * beads_in_use,
                                            std::vector<int> * beads_in_component,
                                            int totalInComponent,
                                            const int workingLimit,
                                            int * pDeadLimit,
                                            Model * pModel){

    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex, endIt = bead_indices->end();
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int neighbor;

    for (int i = 0; i < totalInComponent; i++){ // excludes true model positions

        it = pModel->getPointerToNeighborhood((*beads_in_component)[i]);

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);
            // need to check if neighbor is found in use, if not add it to possible
            if ((neighbor > -1) && (beads_in_use->find(neighbor) == endOfSet)){ // if true, add neighbor to search space

                itIndex = std::find( (bead_indices->begin() + (*pDeadLimit)), endIt, neighbor);
                // if neighbor was already added from previous point in beads_in_component, it will not be found
                if (itIndex != endIt){
                    std::iter_swap( (bead_indices->begin() + (*pDeadLimit)), itIndex);
                    (*pDeadLimit)++;
                }
            } else if (neighbor == -1) {
                break;
            }
        }
    }
}

void Anneal::refineCVXHull(std::vector<int> &bead_indices,
                           std::vector<int> &active_indices,
                           int totalBeadsInSphere,
                           int workingLimit,
                           int *pDeadLimit,
                           Model *pModel){

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[25];
    sprintf(flags, "qhull s FA");

    coordT hullPoints[3*workingLimit];

    // create set to determine CVX Hull
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints[i*3], pModel->getBead(bead_indices[i]));
        active_indices[i] = bead_indices[i]; // copy of the indices before any sorting up to workingLimit
    }

    // needs to be optimized
    qh_new_qhull(3, workingLimit, hullPoints, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;

    pModel->setCVXHullVolume(qh totvol);

    // UPDATE DEADZONE : move points not selected that are outside of hull to deadZone
    *pDeadLimit = totalBeadsInSphere;

    std::vector<int>::iterator beginIt;
    std::vector<int> inside(totalBeadsInSphere-workingLimit);
    std::vector<int> outside(totalBeadsInSphere-workingLimit);
    int insideCount=0, outsideCount=0;
    int * beadPosition;

    // sort unselected beads into those inside and outside CVX Hull
    for(int i=workingLimit; i < *pDeadLimit; i++){
        beadPosition = &bead_indices[i];
        beadToPoint(testPoint, pModel->getBead(*beadPosition));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (isoutside){
            outside[outsideCount] = *beadPosition;
            outsideCount++;
        } else {
            inside[insideCount] = *beadPosition;
            insideCount++;
        }
    }

    beginIt = bead_indices.begin();
    // copy sorted indices back to bead_indices
    std::copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
    *pDeadLimit = workingLimit + insideCount;
    std::copy(outside.begin(), outside.begin()+outsideCount, beginIt + *pDeadLimit);

    // make copy of points from vertices and then free qh_hull
    std::vector<int> vertexIndices(totalV);

    for (int v = 0; v < totalV; v++) { // SWAP VERTICES OF HULL TO BETTER SPOT
        vertexIndices[v] = active_indices[qh_pointid(vertices->point)];
        vertices = vertices->next;
    }

    qh_freeqhull(true);
    enlargeDeadLimit(vertexIndices, totalV, bead_indices, workingLimit, pDeadLimit, totalBeadsInSphere, pModel);
    // return vertexIndices to see if they are at the edges
    //cout <<"TOTAL POINTS IN HULL : " << totalV  << endl;
}

void Anneal::updateASATemp(int index, float evalMax, float acceptRate, float &temp, float &inv_temp){

    bool changed = false;
    float stepEval = index/evalMax;
    float lamRate;

    if (stepEval < 0.15) {
        lamRate = 0.44+0.56*pow(560, -stepEval*6.6666666666666667);
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


/**
 * Diagnostic for calculating potential around each lattice in the model
 */
void Anneal::printContactList(std::vector<int> &bead_indices, std::set<int> * beads_in_use_tree, int workingLimit, Model * pModel){
    cout << "CONTACT LIST" << endl;
    for(int i=0; i<workingLimit; i++){
        int index = bead_indices[i];
        float value = calculateLocalContactPotentialPerBead(beads_in_use_tree, pModel, index);
        int cc = numberOfContactsFromSet(beads_in_use_tree, pModel, index);
        cout<< i << " " << index << " => " << cc <<  "  POTENTIAL => " << value << " <=> " << totalContactsPotential(cc) <<  endl;
    }
}


/**
 * Treat the connectivity potential as a look up table.
 *
 */
void Anneal::populatePotential(int totalNeighbors){

    connectivityPotentialTable.resize(totalNeighbors+1);

    for(int i=0; i<(totalNeighbors+1); i++){
        double diff = 1.0-exp(-alpha*(i-contactsPerBead));
        connectivityPotentialTable[i]=1000.0*(diff*diff);
    }

    connectivityPotentialTable[0] *= 100.0;
    connectivityPotentialTable[1] *= 1.2;
}


void Anneal::modPotential(float factor){

//    float inv = 1.0/factor;
//    int upper = (int)ceil(contactsPerBead);
//    int next = upper-1;
//    connectivityPotentialTable[0] *= factor;
//    connectivityPotentialTable[1] *= factor;
//    connectivityPotentialTable[2] *= 1.2*inv;

//    for(int i=0; i< next-1; i++){
//        connectivityPotentialTable[i] *= factor;
//    }
//
//    connectivityPotentialTable[next-1] *= inv;

//    for(int i=next; i< upper+1; i++){
//        connectivityPotentialTable[i] *= inv;
//    }
//
//    next = upper+1;
//    int total = connectivityPotentialTable.size();
//    for(int i=next; i<total; i++){
//        connectivityPotentialTable[i] *= 1.7*factor;
//    }
}



void Anneal::printContactsFromSet(std::vector<int> &bead_indices, int workingLimit, std::set<int> *beads_in_use,
                                           Model *pModel,
                                           int const selectedIndex){

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> tempVec(workingLimit);
    std::copy(bead_indices.begin(), bead_indices.begin()+workingLimit, tempVec.begin());
    std::shuffle(tempVec.begin(), tempVec.end(), gen);


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

}



/**
 * modelPR and targetPR are the same size
 * targetPR is derived from PDB
 */
float Anneal::calculateKLDivergenceAgainstPDBPR(vector<int> &modelPR, vector<float> &targetPR){

    float totalCounts = 0.0;
    float kl=0.0, prob, * value;
    int totalm = modelPR.size();
    std::vector<float> modelPR_float(modelPR.begin(), modelPR.end());
    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    //last nonzero bin
    int last=0;
    for (int i=0; i<totalm; i++){
        if (targetPR[last] <= 0){
            break;
        }
        last++;
    }

    for (int i=0; i<totalm; i++){
        totalCounts += modelPR_float[i];
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    for (int i=0; i < last; i++){
        prob = targetPR[i];  // bounded by experimental Shannon Number
        //tempPR = modelPR[i];
        value = &modelPR_float[i];
        if (prob > 0 && *value > 0){
            kl += prob * log(prob/(*value) * totalCounts);
        } else if (prob > 0 && *value <= 0){ // severely penalize any model bin that is zero
            kl += 10000000000000000;
        }
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    // severely penalize choices that make a bin zero for values < dmax
    /*
    int last = totalm-1;
    for (int i=0; i < last; i++){

        if (modelPR_float[i] <= 0){
            kl += 10000000000000000000;
        } else {
            prob = targetPR[i];  //
            kl += prob * log(prob/(modelPR_float[i])*totalCounts);
        }
    }

    // assume if last bin in dataset is 0, then 0*log0 = 0
    if (probability_per_bin[last] > 0){
        prob = targetPR[last];
        kl += prob * log(prob/(modelPR_float[last])*totalCounts);
        //kl += prob * log(prob/(modelPR_float[last])*totalCounts);
    }
     */

    return kl;  // returns value per bin
}

// convert PDB model into lattice model
void Anneal::readPDB(Model *pModel, vector<int> * keptBeads, string filename){

    PDBModel pdbModel(filename, true, true, pModel->getBeadRadius());

    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    Bead * currentBead;
    float diffx, diffy, diffz, bx, by, bz;
    int totalAtoms = pdbModel.getTotalAtoms();
    float beadradius = pModel->getBeadRadius();
    float b2 = beadradius*beadradius;

    for(int i=0; i<totalBeads; i++){

        currentBead = pModel->getBead(i);
        bx = currentBead->getX();
        by = currentBead->getY();
        bz = currentBead->getZ();

        for (int t=0; t<totalAtoms; t++){
            diffx = bx - *(pdbModel.getCenteredX()+t);
            diffy = by - *(pdbModel.getCenteredY()+t);
            diffz = bz - *(pdbModel.getCenteredZ()+t);

            if ((diffx*diffx + diffy*diffy + diffz*diffz) < b2){
                keptBeads->push_back(i);
                break;
            }
        }
    }

    int totalKept = keptBeads->size();
    //pModel->printSelectedBeads(0, totalKept, *keptBeads);
    pModel->writeSubModelToFile(0, totalKept, *keptBeads, "ideal");
}


bool Anneal::setAnchorPoints(std::string anchorFileName, std::string pdbFile, Model *pModel){

    PDBModel pdbModel(pdbFile, true, true, pModel->getBeadRadius()); // centered Coordinates

    // if anchor points are in pdbFile, return true, else return false
    // CHAIN, RESIDUE NUMBER, ATOM?
    // ATOM     54  O   GLY A   8
    const int totalAtoms = pdbModel.getTotalAtoms();
    string line;
    int acceptedLines = 0;

    ifstream anchorFile (anchorFileName.c_str());
    boost::regex pdbStart("ATOM");
    boost::regex residue("RESID");
    boost::regex lineFormat("\\w+\\s+[0-9]+\\s+\\w+[A-Z0-9]+", boost::regex::icase);
    boost::regex component_id("COMPONENT_ID");
    boost::regex volume("VOLUME");
    boost::regex chain("CHAIN");
    boost::regex wat("HOH");
    boost::regex hash("#");

    std::vector<int>::const_iterator pdbResIDs = pdbModel.getResIDIterator();
    std::vector<string>::const_iterator pdbAtomTypes = pdbModel.getAtomTypeIterator();
    std::vector<string>::const_iterator pdbChainIds = pdbModel.getChainIDIterator();
    const int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    Bead * currentBead;

    // find closest non-seed bead position!
    // format of Anchor file
    std::vector<std::string> tempLine;
    std::vector<std::string> splitLine;
    std::vector<int> resids;
    std::vector<float> volumes;
    std::vector<std::string> ids;

    int currentResidID;
    std::string currentComponentID;

    // get lines in the file
    if (anchorFile.is_open()) {
        while(!anchorFile.eof()) {
            getline(anchorFile, line);
            boost::algorithm::trim(line);
            tempLine.push_back(line);
        }
    }
    anchorFile.close();

    // get componentIDs and volumes
    try {
        for(std::vector<std::string>::iterator it = tempLine.begin(); it != tempLine.end(); ++it) {

            boost::split(splitLine, *it, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((*it).size() > 0 && boost::regex_search(splitLine[0], component_id) && splitLine.size() == 4 && boost::regex_search(*it, volume)){
                components.push_back( Component(splitLine[1], stof(splitLine[3]), pModel) );
                volumes.push_back(stof(splitLine[3]));
            } else if ( boost::regex_search(splitLine[0], component_id) && splitLine[0].size() > 0) {
                throw std::invalid_argument( "COMPONENT ID or VOLUME NOT SPECIFIED : \n\t" + *it  + " \n");
            }
        }

    } catch (exception &err) {
        cerr<<"Caught "<<err.what()<<endl;
        cerr<<"Type "<<typeid(err).name()<<endl;
        exit(0);
    }

    // for each component add the resids
    try {
        for(std::vector<std::string>::iterator it = tempLine.begin(); it != tempLine.end(); ++it) {

            boost::split(splitLine, *it, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((*it).size() > 0 && boost::regex_search(splitLine[0], residue) && splitLine.size() == 6 && boost::regex_search(*it, component_id) && boost::regex_search(*it, chain)){

                // component_ID must be in the list, if not throw exception
                std::string tempId = splitLine[5];
                auto fit = std::find_if(components.begin(), components.end(), [&tempId](const Component& obj) {return obj.getID() == tempId;});

                if (fit != components.end()) {
                    // found element. it is an iterator to the first matching element.
                    // if you really need the index, you can also get it:
                    auto index = std::distance(components.begin(), fit);
                    int tempResid = stoi(splitLine[1]);
                    if (tempResid > 1){ // check if RESID is in PDB model
                        (*fit).addResid(tempResid, splitLine[3]);
                    } else {
                        throw std::invalid_argument( "IMPROPER RESID: \n\t" + *it  + " RESID \n" + std::to_string(tempResid) + " \n");
                    }
                } else {
                    throw std::invalid_argument( "COMPONENT ID MISSING OR INCORRECT: \n\t" + *it  + " \n");
                }

            } else if ( (*it).size() > 0 && boost::regex_search(splitLine[0], residue) && splitLine.size() < 6 && boost::regex_search(*it, component_id) ) {
                throw std::invalid_argument( "COMPONENT ID or RESID NOT SPECIFIED : \n\t" + *it  + " \n");
            }
        }
    } catch (exception &err) {
        cerr<<"Caught "<<err.what()<<endl;
        cerr<<"Type "<<typeid(err).name()<<endl;
        exit(0);
    }

    // for each Component, get resids
    // map resid to structure, for each resid, grab all the atoms and calculate average
    float xpos, ypos, zpos;
    float diffx, diffy, diffz;
    for(std::vector<Component>::iterator it = components.begin(); it != components.end(); ++it) {
        //Component * temp = it;
        float dis2;
        for(int r=0; r< (it->getTotalResids()); r++){
            xpos=0;
            ypos=0;
            zpos=0;
            int atomCounter=0;
            float min = 10000;

            for (int i=0; i < totalAtoms; i++){ // calculate average position of residue
                if ( it->getResidByIndex(r) == *(pdbResIDs + i) && (it->getChainsByIndex(r).compare(*(pdbChainIds + i)) == 0) ) {
                    //if ( temp.getResidByIndex(r) == *(pdbResIDs + i)  ) {
                    xpos += *(pdbModel.getCenteredX() + i);
                    ypos += *(pdbModel.getCenteredY() + i);
                    zpos += *(pdbModel.getCenteredZ() + i);
                    atomCounter++;
                } else if (it->getResidByIndex(r) < *(pdbResIDs + i)) {
                    break;
                }
            }
            // calculate average position
            float inv = 1.0/(float)atomCounter;
            int keeper;
            xpos *= inv;
            ypos *= inv;
            zpos *= inv;
            // find in bead universe the bead that is closest
            for(int b=0; b < totalBeads; b++){ // iterate over each bead in Universe

                currentBead = pModel->getBead(b);
                diffx = currentBead->getX() - xpos;
                diffy = currentBead->getY() - ypos;
                diffz = currentBead->getZ() - zpos;
                dis2 =(diffx*diffx + diffy*diffy + diffz*diffz);

                if (dis2 <= min){
                    min = dis2;
                    keeper = b;
                }
            }
            it->addAnchor(keeper);
        }
    }

    anchorFile.close();
    anchorFile.clear();

    totalComponents = components.size();
    if (components.size() == 0){
        return false;
    }
    return true;
}