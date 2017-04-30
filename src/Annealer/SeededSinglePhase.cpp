//
// Created by Robert Rambo on 15/01/2017.
//
#include "../Anneal.h"
#include "Data.h"
#include "../Model.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"

using namespace std;

/**
 * Random add/remove to generate initial model for seeded modeling
 * Try to find the minimial set of lattice points that agree with atomistic P(r)
 * contactCutoff determines
 */
void Anneal::createSeedFromPDB(Model *pModel, Data *pData, string name, string PDBFilename, int numberOfUniqueConnectedPhases){

    totalNumberOfPhasesForSeededModeling = numberOfUniqueConnectedPhases;
    contactCutOff = interconnectivityCutOff;
    double lowTempStop = (double)highTempStartForCooling;

    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    // convert distances within the large search space to bins based on input P(R)-DATA file
    this->fillPrBinsAndAssignTotalBin( pBin,  pDistance,  totalDistancesInSphere,  pData);

    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size

    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    const std::vector<int>::const_iterator trueModelBeginIt = pModel->getSeedBegin();
    int workingLimit = pModel->getTotalInSeed();

    // copy trueModel into BeadIndices
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    // create set of anchors if exists
    std::set<int> excludeAnchorsList;
    for(int i=0; i < totalComponents; i++) {
        // the selected set of beads that make up each component will be held by Component object
        // add randomly selected indices to Component
        Component * component = &(components[i]);
        if (!component->anchorsEmpty()){
            for(std::set<int>::iterator it = component->getAnchors()->begin(); it != component->getAnchors()->end(); ++it) {  // find the anchor
                excludeAnchorsList.insert(*it);
            }
            component->writeAnchorsToFile("anchors_"+std::to_string(i+1));
        }
    }

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    int minWorkingLimit = 0.17*workingLimit;

    std::vector<int> bead_indices(workingLimit); // large vector ~1000's
    std::vector<int> lowest_bead_indices(workingLimit); // large vector ~1000's
    std::vector<int> backUpState(workingLimit);

    std::clock_t start;
    // c-style is slightly faster for large vector sizes
    const int num = workingLimit;
    int * ptr = (num != 0) ? &bead_indices.front() : NULL;
    for(int i = 0; i < num; i++) {
        ptr[i] = *(trueModelBeginIt+i);
    }
    // prepare bead_indices by copying in truemodel
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    // check for anchors, if anchors are present, move them out of the core model (anchors are only considered in the component)
//    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << endl;
//    std::cout << "              INITIAL WORKINGLIMIT   " << workingLimit << endl;
//    if(excludeAnchorsList.size() > 0){
//        for(int i=0; i<workingLimit; i++){
//            if (excludeAnchorsList.find(bead_indices[i]) != excludeAnchorsList.end()){
//                // move to workingLimit and decrement
//                int decrement=1;
//                while( !( excludeAnchorsList.find(*(bead_indices.begin()+(workingLimit-decrement))) == excludeAnchorsList.end() ) ){
//                    decrement--;
//                }
//                cout << " FOUND ANCHOR AND SWAPPING TO END " << bead_indices[i] << endl;
//                std::iter_swap(bead_indices.begin() + i, bead_indices.begin() + (workingLimit-decrement));
//                workingLimit -= decrement;
//            }
//        }
//    }
//    std::cout << " (AFTER ANCHOR CHECK) WORKINGLIMIT  " << workingLimit << endl;
    int deadLimit = workingLimit;; // as bead indices are discarded, set upper limit of vector
    //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.end());

    std::random_device rd;
    std::mt19937 gen(rd());

    double inv_kb_temp = 1.0/0.000001;

    int lowerN = 0.1*workingLimit, upperN = workingLimit;
    cout << " TOTAL LATTICE IN SEED : " << workingLimit << endl;
    cout << "        LATTICE LIMITS : " << lowerN << " <= N <= " << upperN << endl;
    // randomize and take the workingLength as first set, shuffling takes about 10x longer than copy and sort

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isConnected = true;
    // bead_indices contains only the indices that relate to the input PDB
    //isConnected = isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

    cout << " CHECKING CONNECTIVITY  " << endl;
    cout << "        IS CONNECTED ? : " << isConnected << endl;
    cout << "  NUMBER OF COMPONENTS : " << currentNumberOfComponents << endl;

    //const int components = tempNumberOfComponents;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();

    // calculate KL against known
    // calculateKLEnergy populates binCount based on model in bead_indices
    // create custom D_KL function supplying Pr_of_PDB as target
    float testKL, currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    std::cout << "   DIRECT INITIAL D_KL : " << currentKL << " (AGAINST DATA) " << std::endl;
    currentKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
    //float tempTotalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    std::cout << "          INITIAL D_KL : " << currentKL << std::endl;
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::vector<int>::iterator beginIt = bead_indices.begin();
    std::vector<int>::iterator endIt = bead_indices.end();

    int lowestWorkingLimit = workingLimit, lowestDeadLimit = deadLimit;
    std::copy(beginIt, endIt, lowest_bead_indices.begin());

    int counter=1;
    int priorWorkingLimit;

    float average_x=0.5*workingLimit, stdev=0.2*workingLimit;
    float sum_x_squared=0, sum_x=0, divideBy=0, acceptRate = 0.5, inv500 = 1.0/500.0;
    //output for plotting
    bool isUpdated = false;

    int seedHighTempRounds = 2*highTempRounds;
    int deadUpdate = std::ceil(seedHighTempRounds*0.091);

    populatePotential(pModel->getSizeOfNeighborhood());
    float startContactSum=0.0;
    for (int i=0; i<lowestWorkingLimit; i++){
        startContactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }

    double lowestRunningContactSum, runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, pModel);

    double etaConstant = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) + 4 );
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);
    double etaFactor = 1.0/std::pow(10000, 1.0/(seedHighTempRounds/(float)deadUpdate));

    startContactSum = startContactSum/(float)workingLimit;
    double newSum;
    float this_energy, startKL = currentKL;
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    float lowestE = current_energy;
    char addRemoveText[50];

    int high;
    for (high=0; high < seedHighTempRounds; high++){ // iterations during the high temp search
        beginIt = bead_indices.begin();
        endIt = bead_indices.end();
        std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
        cout << "******************************************************************************" << endl;
        if (distribution(gen) >= 0.5){ // ADD BEAD?

           // std::copy(beginIt+workingLimit, beginIt + deadLimit, active_indices.begin()+workingLimit);
           // std::shuffle(active_indices.begin()+workingLimit, active_indices.begin() + deadLimit, gen);
           // std::shuffle(beginIt+workingLimit, endIt, gen);

            printf("     ADDME => %i \n", 1);
            if (workingLimit < deadLimit){
                cout << "*******************                  ADD                   *******************" << endl;
                sprintf(addRemoveText, "");
                double afterAdding;
                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                original = bead_indices[randomSpot];

//                while( !( excludeAnchorsList.find(original) == excludeAnchorsList.end() ) ){ // do not add anchor point
//                    randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
//                    original = bead_indices[randomSpot];
//                }

                if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) > 0){
                    itIndex = bead_indices.begin() + randomSpot;
                    // select first element from randomized active_set
                    // check if available neighbor can be added
                    newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);
                    afterAdding = etaConstant*newSum/(double)(workingLimit+1);
                    // make the swap at the border (workingLimit)
                    addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);
                    addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                    //testKL = pData->calculateKLDivergence(binCount);
                    testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
                    // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                    tempNumberOfComponents = 1;
                    this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);
                    beads_in_use_tree.insert(original);
                    tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                    if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterAdding;
                        runningContactsSum = newSum;
                        isUpdated = true;
                        eulerTour.addNode(original, pModel);
                        sprintf(addRemoveText, "     ADD => %i", 1);
                    } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterAdding;
                        runningContactsSum = newSum;
                        isUpdated = true;
                        eulerTour.addNode(original, pModel);
                        sprintf(addRemoveText, "     ADD => %i", 1);
                    } else { // undo changes (rejecting)
                        beads_in_use_tree.erase(original);
                        restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

            //sort(beginIt, beginIt + workingLimit);
        } else { // REMOVE BEADS?
            cout << "*******************                 REMOVE                 *******************" << endl;
            // test for deletion
            priorWorkingLimit = 1;
            printf("     REMOVE => %i\n", priorWorkingLimit);
            sprintf(addRemoveText, "     REMOVE => %i", priorWorkingLimit);

            double afterRemoving;
            int randomSpot = rand() % workingLimit;
            original = bead_indices[randomSpot];
            tempNumberOfComponents = eulerTour.removeNode(original);

            if (tempNumberOfComponents == 1){
                // grab from randomized active_indices list
                newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                afterRemoving = etaConstant*newSum/(double)(workingLimit-1);
                removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                // still need to sort, swap changes the order
                //testKL = pData->calculateKLDivergence(binCount);
                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
                this_energy = testKL + lambda*connectivityPotential(1);

                beads_in_use_tree.erase(original);
                tempTotalContactEnergy = (afterRemoving - totalContactEnergy);
                if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    isUpdated = true;
                } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    isUpdated = true;
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.addNode(original, pModel);
                }
            } else {
                eulerTour.addNode(original, pModel);
            }
        }


        if (current_energy < lowestE){
            lowestE = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            lowestRunningContactSum = runningContactsSum;
            std::copy(beginIt, endIt, lowest_bead_indices.begin());
        }

        sum_x_squared += workingLimit*workingLimit;
        sum_x += workingLimit;
        divideBy += 1.0;

        cout << "*******************                                        *******************" << endl;
        printf("       TEMP => %-.8f \n     INVKBT => %.4f\n", lowTempStop, inv_kb_temp);
        printf("   MAXSTEPS => %i (%4i) \n", seedHighTempRounds, high);
        printf("      GRAPH => %i \n", currentNumberOfComponents);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f \n", totalContactEnergy, etaConstant, runningContactsSum);
        printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", average_x, stdev);
        printf("LIMIT: %5i (>= MIN: %i) DEADLIMIT: %5i D_KL: %.4E \n", workingLimit, minWorkingLimit, deadLimit, currentKL);

        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;
        } else {
            acceptRate = inv500*(499*acceptRate);
        }

        updateASATemp(high, seedHighTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        if (counter % deadUpdate == 0){ // if counter too small, add/remove may not sample sufficiently
            etaConstant *= etaFactor;
            totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);
        }
        counter++;

    } // end of HIGH TEMP EQUILIBRATION

    // average number of contacts per bead

    float contactSum=0.0;
    //for (int i=0; i<lowestWorkingLimit; i++){
    for (int i=0; i<workingLimit; i++){
        //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
        contactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }

    float average_number_of_contacts = contactSum/(float)workingLimit;
    double unscaled = runningContactsSum / (double)workingLimit;

    cout << "KL DIVERGENCE : " << endl;
    cout << "                   INITIAL D_KL => " << startKL << endl;
    cout << "                    LOWEST D_KL => " << lowestE << endl;
    cout << "   UNSCALED RUNNING CONTACT SUM => " << unscaled << endl;
    cout << "AVERAGE CONTACTS : (per lattice point)" << endl;
    cout << "                        INITIAL => " << startContactSum << endl;
    cout << "                          FINAL => " << average_number_of_contacts << " " << totalContactEnergy << endl;
    cout << " Contacts Per Bead " << contactsPerBead << endl;

    // this is fixed model for initial high temp search?
    // set this as the seed
    //    pModel->setBeadAverageAndStdev(average_x, stdev);
    //    pModel->setReducedSeed(workingLimit, bead_indices);
    //    pModel->writeModelToFile(workingLimit, bead_indices, name);
    //    pModel->setBeadAverageAndStdev(average_x, stdev);

    // remap to bead universe for storage and printing
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    bead_indices.resize(totalBeadsInSphere);
    beginIt = bead_indices.begin();
    endIt = bead_indices.end();
    ptr = &bead_indices.front();
    for(int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::vector<int>::iterator it;
    for (int i=0; i<workingLimit; i++){
        it = std::find(beginIt, endIt, lowest_bead_indices[i]); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(beginIt+i, it);
    }

    //pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setStartingDeadLimit(deadLimit);

    cout << "Starting lowest wl " << lowestWorkingLimit << endl;
    // set seed model to be invariant during reconstruction

    char flags[25];
    sprintf(flags, "qhull s FA");
    //int numpoints = 3*totalBeadsInSphere;
    //coordT points[numpoints];
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // REMOVE ALL ANCHORS FROM REDUCED MODEL
    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << endl;
    std::cout << "              INITIAL WORKINGLIMIT   " << workingLimit << endl;
    if(excludeAnchorsList.size() > 0){
        for(std::set<int>::iterator it = excludeAnchorsList.begin(); it != excludeAnchorsList.end(); ++it){
            std::vector<int>::iterator found = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, *it);
            int dis = std::distance(bead_indices.begin(), found);
            if (dis < workingLimit){ // move it
                cout << *it <<  " FOUND ANCHOR AND SWAPPING TO END " << " " << workingLimit << endl;
                std::iter_swap(found, bead_indices.begin() + workingLimit-1);
                workingLimit--;
            }
        }
    }

    std::cout << "                      WORKINGLIMIT   " << workingLimit << endl;
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    pModel->setReducedSeed(workingLimit, bead_indices);

    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            name,
            this,
            pData,
            high,
            cvx,
            average_number_of_contacts);
}

// create set of true lattice points that will remain during the search
// add anchor points
// look for a connected space that only includes anchors from true set
//
// initial high temp search is find connected tour that includes anchor(s)


bool Anneal::createInitialModelCVXHullSeeded(Model *pModel, Data *pData, std::string name) {
    
    // do short high temp to create Euler tours for each component
    srand(time(0));
    // convert distances to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class
    // faster if pointers to above is reused - conversion happens only once
    // pBin will be reused many times
    maxbin=0;
    for(int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);
    std::vector<int>::iterator beginIt = bead_indices.begin(), endIt = bead_indices.end();
    //std::clock_t start;
    // c-style is slightly faster for large vector sizes
    const int num = totalBeadsInSphere;
    int * ptr = (num != 0) ? &bead_indices.front() : NULL;
    for(int i = 0; i < num; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    float invBeadVolume = 1.0/pModel->getBeadVolume();
    // connectivity potential
    // volume potential
    // pick random number between lowerN and upperN
    pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "universe", 0);
    // create initial model using lattice points from seed
    // add additional points randomly based on the volume of each component
    const std::set<int> true_model_tree(pModel->getSeedBegin(), pModel->getSeedEnd());            // complete lattice model of PDB
    std::set<int> reduced_model_tree(pModel->getReducedSeedBegin(), pModel->getReducedSeedEnd()); // reducedModel
    int totalInReducedSeed = pModel->getTotalInReducedSeed();

    // fill bead indices
    int locale=0;
    for(std::set<int>::iterator it = pModel->getReducedSeedBegin(); it != pModel->getReducedSeedEnd(); ++it) {
        std::vector<int>::iterator fit = std::find(beginIt, endIt, *it); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(beginIt+locale, fit);
        locale++;
    }

    // estimate volume per reduced seed, use this to estimate beads per Component
    //const int totalComponents = totalComponents;
    cout << " TOTAL NUMBER OF COMPONENTS " << totalComponents << endl;
    //pModel->estimatePointsPerComponent(pModel->getBeadVolume()*true_model_tree.size()/(float)totalInReducedSeed);
    // for each component in pModel, get target number of beads based on estimate
    int workingLimit = totalInReducedSeed;
    std::shuffle(bead_indices.begin()+workingLimit, bead_indices.end(), gen);
    std::set<int> anchorsInUse;

    for(int i=0; i < totalComponents; i++) { // the selected set of beads that make up each component will be held by Component object
        // add randomly selected indices to Component
        Component * component = &(components[i]);
        component->setTargetNumberOfLatticePoints(pModel->getBeadVolume()*true_model_tree.size()/(float)totalInReducedSeed);
        int toAdd = component->getTargetNumberOfLatticePoints();
        cout << i <<" ESTIMATED LATTICE POINTS " << toAdd << endl;

        for(int j=0; j<toAdd; j++){
            component->addLatticePoint(bead_indices[workingLimit]); // grab random set of points starting from anchor?
            workingLimit++;
        }

        if (!component->anchorsEmpty()){
            // add all anchor points to initial model
            for(std::set<int>::iterator it = component->getAnchors()->begin(); it != component->getAnchors()->end(); ++it) {  // find the anchor
                std::vector<int>::iterator fit = std::find(beginIt, endIt, *it);  // if itTrueIndex == endTrue, it means point is not within the set
                // remove anchor from core if it is there and put into component
                int dis = std::distance(beginIt, fit);
                if (dis >= workingLimit){
                    cout << "ANCHOR NOT FOUND ADDING TO INITIAL MODEL " << *it << endl;
                    std::iter_swap(beginIt+workingLimit, fit);
                    workingLimit++;
                }
                cout << "ANCHOR ADDED TO INUSE " << *it << endl;
                //anchorsInUse.insert(*it);
            }
            // set centered Anchors
            for(std::set<int>::iterator it = component->getCenteredAnchors()->begin(); it != component->getCenteredAnchors()->end(); ++it) {  // find the anchor
                std::vector<int>::iterator fit = std::find(beginIt, endIt, *it);  // if itTrueIndex == endTrue, it means point is not within the set
                // remove anchor from core if it is there and put into component
                anchorsInUse.insert(*it);
            }
        }
    }

    int lowerN = workingLimit - 0.315*workingLimit;
    int upperN = (int)(workingLimit + 0.315*workingLimit);
    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN );

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    // 1. create euler tour of selected set
    // 2. populate search space
    // setup parameters for hull
    char flags[25];
    sprintf(flags, "qhull s FA");
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];

    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int alterMe, tempNumberOfComponents;
    // isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();
    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    lambda = 0.01;
    beta = 0.00001;
    mu = 0.000001; // scale this to volume of search space
    // the number of components in each Component object must be looked at independently
    float componentPotential=0;
    cout << " TOTAL COMPONENT(S) CHECK => " << totalComponents << endl;
    for(int i=0; i < totalComponents; i++) {
        // the selected set of beads that make up each component will be held by Component object
        componentPotential+=components[i].potential();
    }
    cout << "   CONNECTIVITY POTENTIAL => " << connectivityPotentialPhases(eulerTour.getNumberOfComponents()) << endl;
    float testKL, test_energy, current_energy = currentKL + lambda*connectivityPotentialPhases(eulerTour.getNumberOfComponents()) + mu*current_volume/(float)workingLimit  + beta*componentPotential;
    float lowest_energy = current_energy;

    int testIndex;
    float inv_kb_temp = 1.0/(float)highT, tempConnectivityPotential, currentNumberOfComponents;
    inv_kb_temp = 1.0/0.0001;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    cout << " POPULATING DEAD LIMIT" << endl;
    int deadLimit; // as bead indices are discarded, set upper limit of vector
    populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere);
    cout << " DEAD LIMIT SET => " << deadLimit << endl;

    int lowestWorkingLimit, lowestDeadLimit;
    lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    lowestDeadLimit = deadLimit;

    bool updateCVX =false;

    cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< endl;
    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            "initial_random_model",
            this,
            pData,
            0,
            0,
            0);

    int high=0, componentIndex, terminateCount=0;
    bool startCount=false;

while (terminateCount < 10000){
    //for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy

        if (distribution(gen) < 0.57){

            alterMe = number_of_beads_to_use(gen);  // uniform distribution

            beginIt = bead_indices.begin();
            endIt = bead_indices.end();
            componentIndex = totalComponents;

            if (alterMe > workingLimit){
                cout << "*******************                  ADD                   ******************* "  << endl;
                // pick component
                // int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                int randomSpot = bead_indices[rand() % workingLimit];
                // find neighbor of randomSpot that is available
                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomSpot);
                // add only to non-reduced model indices
                while( !( (reduced_model_tree.find(randomSpot) == reduced_model_tree.end()) && (testIndex > -1) ) ){
                    randomSpot =  bead_indices[rand() % workingLimit];
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomSpot);
                }

                // which component does it belong to?
                for(int i=0; i < totalComponents; i++) {
                    // the selected set of beads that make up each component will be held by Component object
                    if (components[i].inUse(randomSpot)){
                        componentIndex=i;
                        break;
                    }
                }

                itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);
                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = pData->calculateKLDivergence(binCount);
                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                // figure out which component it belongs too
                // need to calculate the number of components
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                if (componentIndex < totalComponents)
                    components[componentIndex].addLatticePoint(testIndex);

                componentPotential=0;
                for(int i=0; i < totalComponents; i++)
                    componentPotential+=components[i].potential();

                tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);
                test_energy = testKL + lambda*tempConnectivityPotential + mu*test_volume/(float)workingLimit + beta*componentPotential;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempConnectivityPotential;
                    beads_in_use_tree.insert(testIndex);
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempConnectivityPotential;
                    beads_in_use_tree.insert(testIndex);
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                    if (componentIndex < totalComponents){
                        components[componentIndex].removeLatticePoint(testIndex);
                    }
                }

            } else if (alterMe < workingLimit) {

                cout << "*******************                 REMOVE                 ******************* "  << endl;

                if (workingLimit > lowerN){
                    // go through entire list
                    testIndex = bead_indices[rand() % workingLimit];

                    // do not remove special positions (reduced model and anchors)
                    // if testIndex is not in reduced_model, find will return end
                    // maintain at least one anchorInUse

                    while( !((reduced_model_tree.find(testIndex) == reduced_model_tree.end()) && (canRemoveIfAnchor(testIndex)) ) ){
                        testIndex = bead_indices[rand() % workingLimit];
                    }

                    for(int i=0; i < totalComponents; i++) {
                        // the selected set of beads that make up each component will be held by Component object
                        if (components[i].inUse(testIndex)){
                            componentIndex=i;
                            break;
                        }
                    }
                    // grab from randomized active_indices list

                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = pData->calculateKLDivergence(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    if (componentIndex < totalComponents)
                        components[componentIndex].removeLatticePoint(testIndex);

                    componentPotential=0;
                    for(int i=0; i < totalComponents; i++) { // the selected set of beads that make up each component will be held by Component object
                        componentPotential+=components[i].potential();
                    }

                    tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);
                    test_energy = testKL + lambda*tempConnectivityPotential + mu*test_volume/(float)workingLimit + beta*componentPotential;

                    if (test_energy < current_energy ) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempConnectivityPotential;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else if (exp(-(test_energy - current_energy)*inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempConnectivityPotential;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        if (componentIndex < totalComponents){
                            components[componentIndex].addLatticePoint(testIndex);
                        }
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

        } else {

            cout << "*******************               POSITIONAL               *******************" << endl;

            int swap1, bb;
            std::vector<int>::iterator pSwap2;
            // pick component index
            // calculate CVX hull
            // only move non-anchor point
            componentIndex = rand() % totalComponents;
            int countIt=0;
            for (int i = 0; i < workingLimit; i++) {
                int toCheck = bead_indices[i];
                // excluded reduced model but include anchors
                if ( components[componentIndex].inUse(toCheck) ){ // includes anchor points
                    beadToPoint(&points[countIt*3], pModel->getBead(toCheck));
                    active_indices[countIt] = toCheck;
                    countIt++;
                }
            }

            // needs to be optimized
            qh_new_qhull( 3, countIt, points, 0, flags, NULL, NULL);
            vertexT * vertices = qh vertex_list;
            int totalV = qh num_vertices;

            // make copy of points from vertices and then free qh_hull
            bool isSwapped;

            // only move CVX hull points
            std::vector<int> indices_to_check(totalV); // large vector ~1000's
            // I want to move a bead that is not part of reduced seed
            // want to only move beads that are part of the newly defined component
            countIt=0;
            for (int v = 0; v < totalV; v++) {
                indices_to_check[countIt] = active_indices[qh_pointid( vertices->point)];
                countIt++;
                vertices = vertices->next;
            }
//
            qh_freeqhull(true);

            populatedDeadLimitExcludingSet(&bead_indices,
                                           &beads_in_use_tree,
                                           &indices_to_check,
                                           countIt,
                                           workingLimit,
                                           &deadLimit,
                                           pModel);

            std::shuffle(bead_indices.begin()+workingLimit, bead_indices.begin()+deadLimit, gen);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.begin()+countIt, gen);

            componentPotential=0;
            for(int i=0; i < totalComponents; i++) {
                // the selected set of beads that make up each component will be held by Component object
                componentPotential+=components[i].potential();
            }

            if (countIt > 3){
                for (int v = 0; v < 1; v++) { // try 3 times
                    swap1 = indices_to_check[rand()%countIt];
                    while (!(anchorsInUse.find(swap1) == anchorsInUse.end())){
                        swap1 = indices_to_check[rand()%countIt];
                    }

                    isSwapped = false;
                    // find bead to swap in active set
                    itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                    // remove selected index from P(r)
                    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
                    removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                    std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // make copy of altered P(r)
                    // find better position
                    eulerTour.removeNode(swap1);
                    beads_in_use_tree.erase(swap1);
                    components[componentIndex].removeLatticePoint(swap1);

                    //for (bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
                    for (bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
                        // if neighbor is not in use, then find it in bead_indices and assign to pSwap2
                        int neighbor = bead_indices[bb];
                        //cout << componentIndex << " Neighbor " << neighbor << endl;
                        // make the swap, sort and update P(r)
                        //pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), neighbor);
                        pSwap2 = bead_indices.begin() + bb;

                        std::iter_swap(itIndex, pSwap2);
                        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                        addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        testKL = pData->calculateKLDivergence(workingBinCount);

                        tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);
                        //cout << "MAIN CONNECTIVITY " << tempNumberOfComponents << endl;
                        components[componentIndex].addLatticePoint(neighbor);

                        //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                        tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);

                        test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                        test_energy = testKL + lambda*tempConnectivityPotential + mu*test_volume/(float)workingLimit + beta*componentPotential;

                        if (test_energy < current_energy) {
                            beads_in_use_tree.insert(neighbor);
                            //move swapped position to the deadLimit
                            std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
                            currentNumberOfComponents = tempConnectivityPotential;
                            currentKL = testKL;
                            current_volume = test_volume;
                            current_energy = test_energy;
                            isSwapped = true;
                            updateCVX = true;
                            break;
                        } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                            beads_in_use_tree.insert(neighbor);
                            std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin());
                            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
                            currentNumberOfComponents = tempConnectivityPotential;
                            currentKL = testKL;
                            current_volume = test_volume;
                            current_energy = test_energy;
                            isSwapped = true;
                            updateCVX = true;
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, swap1); // find swapped value in the sorted array
                        eulerTour.removeNode(neighbor);
                        components[componentIndex].removeLatticePoint(neighbor);

                    } // end of neighborhood for-loop
                    // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                    if (!isSwapped){ // itIndex already reverses on exit from loop
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                        eulerTour.addNode(swap1, pModel);
                        components[componentIndex].addLatticePoint(swap1);
                        beads_in_use_tree.insert(swap1);
                    }
                }
            }
        }

        if (updateCVX){
            //deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel,  totalBeadsInSphere);
            //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
            populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
            //populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere);
            updateCVX=false;
            if (startCount){
                terminateCount++;
            }
        }


        cout << "*******************                                        *******************" << endl;
        printf("   MAXSTEPS => %i (%5i) \n", high, terminateCount);
        printf("  GRAPH POT => %.2f\n", currentNumberOfComponents);
        printf("      UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("     VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, mu, mu*current_volume);
        printf("LIMIT: %5i DEADLIMIT: %5i D_KL: %.4E ENRGY: %.4E \n", workingLimit, deadLimit, currentKL, current_energy);
        cout << "*******************                                        *******************" << endl;
        cout << "*******************              CONDENSATION              *******************" << endl;
        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return false;
//        }

        if (currentNumberOfComponents == 0 && current_energy < lowest_energy){
            lowest_energy = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            for(int i=0; i<totalComponents; i++){
                components[i].setBest();
            }

            if (!startCount) {
              startCount = true;
            }
        }

        high++;

        // do reset if not converged after 10000 steps
        if (high > 20000 && !startCount){
          return false;
        }
}

    pModel->printBeadsFromSet(anchorsInUse);

//    pModel->writeModelToFile(workingLimit, bead_indices, name, high);
    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);
//    pModel->writeModelToFile(lowestDeadLimit, lowest_bead_indices, "initial_search_layer", high);
    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    //pModel->setStartingWorkingLimit(lowestDeadLimit);
    pModel->setStartingDeadLimit(lowestDeadLimit);
    pModel->setBeadAverageAndStdev(1.1*pModel->getStartingWorkingLimit(), 0.17*1.1*pModel->getStartingWorkingLimit());

    cout << "*******************                                        *******************" << endl;
    cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << endl;
    printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    cout << "*******************                                        *******************" << endl;

    if (currentNumberOfComponents==0 ) {
        return true;
    } else {
        cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << endl;
        cout << "INCREASE highTempRounds, g" << endl;
        pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "failed", high);
        //reset tour;
        return false;
    }
}

/**
 *
 */
std::string Anneal::refineHomogenousBodyASACVXSeeded(Model *pModel, Data *pData, std::string outputname){
    cout << "########################<<<<<>>>>>>############################## " << endl;
    cout << "#                                                               # " << endl;
    cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << endl;
    cout << "#                                                               # " << endl;
    cout << "########################<<<<<>>>>>>############################## " << endl;

    //float oldN = pModel->getVolumeAverage();     // need to reset this for modeling via rounds
    //float oldStdev = pModel->getVolumeStdev();
    float runningAverage = pModel->getVolumeAverage();
    float runningVariance = pModel->getVolumeStdev();

    //populatePotential(pModel->getSizeOfNeighborhood());
    double lowTempStop = (double)highTempStartForCooling;
    int priorWorkingLimit;
    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));

    // make copy of bead_indices
    std::vector<int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);
    std::vector<int>::iterator beginIt, endIt, itIndex, pSwap2;

    // copy Starting_Set from initial model
    pModel->setModelBeginStartingSetIterator(beginIt);
    pModel->setModelEndStartingSetIterator(endIt);
    std::copy(beginIt, endIt, bead_indices.begin()); // complete model, includes the components
    std::set<int> reduced_model_tree(pModel->getReducedSeedBegin(), pModel->getReducedSeedEnd()); // reducedModel
    std::vector<int> reduced_model_vec(reduced_model_tree.size());
    std::copy(reduced_model_tree.begin(), reduced_model_tree.end(), reduced_model_vec.begin());

    int workingLimit = pModel->getStartingWorkingLimit();

    // reset components
    std::set<int> anchorsInUse;
    for(std::vector<Component>::iterator it = components.begin(); it != components.end(); ++it){
        it->copyBestToInUse();
        //it->printBest();
        if (!it->anchorsEmpty()){
            // add any anchor points, anchor points can never be deleted and belong to the component via Euler tour
            // from random selection above, anchor points
            for(std::set<int>::iterator sit = it->getCenteredAnchors()->begin(); sit != it->getCenteredAnchors()->end(); ++sit) {  // find the anchor
                anchorsInUse.insert(*sit);
            }
        }
    }
    // reset iterators to internal bead_indices
    // pModel->centerLatticeModel(workingLimit, bead_indices);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // set deadLimit of the selected set
    int deadLimit;
    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);

    // convert distances in Search Sphere to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    this->fillPrBinsAndAssignTotalBin(pBin,  pDistance,  totalDistancesInSphere,  pData);
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS: " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS: " << maxbin << endl;
    cout << "              BINWIDTH: " << pData->getBinWidth() << endl; // lattice should be limited by binwidth

    std::vector<int>::iterator beginBinCount = binCount.begin();
    std::vector<int>::iterator endBinCount = binCount.end();

    float testKL, currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

    //CVX HULL STUFF
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    char flags[25];
    sprintf(flags, "qhull s FA");

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float tempAverageContacts;

    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;

    float step_limit = (updateCount < 10000) ? 27000 : updateCount;
    int deadUpdate = std::ceil(updateCount*0.091);

    std::vector<float> tempDuringRun(step_limit);
    std::vector<float> divergenceDuringRun(step_limit);
    std::vector<int> workingLimitDuringRun(step_limit);

    float * pTempDuringRun = &tempDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float sum_x_squared=0, sum_x=0, divideBy=0, average_x, stdev, acceptRate = 0.5, inv500 = 1.0/500.0;

    int counter=1;
    double inv_kb_temp = 1.0/lowTempStop;
    float this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    // calculate initial/current energy of the system
    double runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, pModel);
    float componentPotential=0;
    double componentContactSum=0;
    for(int i=0; i < totalComponents; i++) {
        componentPotential+=components[i].potential();
        components[i].setTotalContactSum(calculateTotalContactSum(components[i].getBeadsInUse(), pModel));
        componentContactSum += components[i].getTotalContactSum()/(double)components[i].getBeadsInUse()->size();
    }

    double etaConstant = pow(10, std::floor(log10(currentKL) - log10(runningContactsSum / (double)workingLimit + componentContactSum)) + 4 );
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit + componentContactSum);
    double etaFactor = 1.0/std::pow(10000, 1.0/(step_limit/(float)deadUpdate));

    float low_temp_limit = step_limit*0.91;
    // want D_KL to settle < 10^-5 to 10^-6
    beta=0.001;
    percentAddRemove = 0.1175;

    float currentComponentPotential = componentPotential;
    float current_energy = currentKL + lambda*connectivityPotentialPhases(currentNumberOfComponents) + beta*componentPotential;
    std::normal_distribution<float> volumeGen(pModel->getVolumeAverage(), pModel->getVolumeStdev());

    std::clock_t startTime;
    int attempts=0, failures=0;
    double newSum;
    double runtime;
    int numberOfCoolingTempSteps, componentIndex;
    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
        beginIt = bead_indices.begin();
        endIt = bead_indices.end();
        componentIndex = totalComponents;
        componentContactSum = 0;
        if ( distribution(gen) < percentAddRemove ) { //add or remove bead within working Set (exclude deadzone)
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            int alterMe = (int) volumeGen(gen);
            // build a list of indices within defined region of convex hull
//            if (alterMe > workingLimit){ // ADD BEAD?
            if ( distribution(gen) < 0.4985761 ) {
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                std::sprintf(addRemoveText, "");
                startTime = std::clock();
                double afterAdding;

                //pick random bead => should favor the component
                std::set<int> * tree = &beads_in_use_tree;
                int randomSpot = bead_indices[rand() % workingLimit];
                int toAdd = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomSpot);
                // NEED TO INCLUDE THE ANCHORS
                if (distribution(gen) < 0.071173) { // add to a non component - can't be in anchor list
                    while (  !((anchorsInUse.find(randomSpot) == anchorsInUse.end()) && (toAdd > -1)) || inComponents(randomSpot)  ) {
                        randomSpot = bead_indices[rand() % workingLimit];
                        toAdd = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomSpot);
                    }
                } else { // add to existing component(s)
                    while (  !( inComponents(randomSpot) && (toAdd > -1) ) ) { // find a new position
                        randomSpot = bead_indices[rand() % workingLimit]; //find random index in components
                        toAdd = getUseableNeighborFromSet(&beads_in_use_tree, pModel, randomSpot);
                    }
                    componentIndex= getComponentIndex(randomSpot);
                    tree = components[componentIndex].getBeadsInUse();
                }


                if (numberOfContactsFromSet(tree, pModel, toAdd) > 1){ // find neighbor of randomSpot that is available
                    itIndex = std::find(bead_indices.begin(), bead_indices.end(), toAdd);
                    // select first element from randomized active_set
                    // check if available neighbor can be added
                    newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, toAdd, runningContactsSum);
                    // if added to Component, must recalculate contactSum
                    if (componentIndex < totalComponents) {
                        Component * comp = &components[componentIndex];
                        comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, toAdd, comp->getTotalContactSum()));
                        comp->addLatticePoint(toAdd);
                    }
                    // if added to Component, must recalculate componentPotential
                    componentPotential=0;
                    for(int i=0; i < totalComponents; i++) {
                        Component * comp = &components[i];
                        componentPotential+=comp->potential(); // constraints the number of beads
                        componentContactSum += comp->getTotalContactSum()/(double)comp->getBeadsInUse()->size();
                    }

                    afterAdding = etaConstant*(newSum/(double)(workingLimit+1) + componentContactSum);
                    // make the swap at the border (workingLimit)
                    addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);
                    addToPr(toAdd, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                    testKL = pData->calculateKLDivergence(binCount);

                    // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                    tempNumberOfComponents = 1;
                    this_energy = testKL + lambda*connectivityPotentialPhases(tempNumberOfComponents) + beta*componentPotential;;

                    beads_in_use_tree.insert(toAdd);
                    tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                    if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterAdding;
                        runningContactsSum = newSum;
                        currentComponentPotential=componentPotential;
                        isUpdated = true;
                        eulerTour.addNode(toAdd, pModel);
                        std::sprintf(addRemoveText, "     ADD => %i", 1);
                    } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterAdding;
                        runningContactsSum = newSum;
                        currentComponentPotential=componentPotential;
                        isUpdated = true;
                        eulerTour.addNode(toAdd, pModel);
                        std::sprintf(addRemoveText, "     ADD => %i", 1);
                    } else { // undo changes (rejecting)
                        beads_in_use_tree.erase(toAdd);
                        restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);

                        if (componentIndex < totalComponents) {
                            Component * comp = &components[componentIndex];
                            comp->setTotalContactSum(recalculateContactsPotentialSumRemove(comp->getBeadsInUse(), pModel, toAdd, comp->getTotalContactSum()));
                            comp->removeLatticePoint(toAdd);
                        }
                    }
                }
                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

            } else { // REMOVE BEADS?
                cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
                priorWorkingLimit = 1;
                std::sprintf(addRemoveText, "     REMOVE => %i", priorWorkingLimit);
                double afterRemoving;
                startTime = std::clock();

                int original = bead_indices[rand() % workingLimit];
                // do not remove special positions (reduced model and anchors)
                // if testIndex is not in reduced_model, find will return end
                while( !((reduced_model_tree.find(original) == reduced_model_tree.end()) && (canRemoveIfAnchor(original) )) ){
                    original = bead_indices[rand() % workingLimit];
                }

                componentIndex = getComponentIndex(original);
                // check euler tour
                tempNumberOfComponents = eulerTour.removeNode(original);
                if (componentIndex < totalComponents) {
                    Component * comp = &components[componentIndex];
                    comp->setTotalContactSum(recalculateContactsPotentialSumRemove(comp->getBeadsInUse(), pModel, original, comp->getTotalContactSum()));
                    comp->removeLatticePoint(original); // recalculate the totalContactSum
                }

                if (connectivityPotentialPhases(tempNumberOfComponents) == 0){
                    // grab from randomized active_indices list
                    newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                    componentPotential=0;
                    for(int i=0; i < totalComponents; i++) {
                        Component * comp = &components[i];
                        componentPotential+=comp->potential(); // constraints the number of beads
                        componentContactSum += comp->getTotalContactSum()/(double)comp->getBeadsInUse()->size();
                    }

                    afterRemoving = etaConstant*(newSum/(double)(workingLimit-1) + componentContactSum);
                    beginIt = bead_indices.begin();
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                    // still need to sort, swap changes the order
                    testKL = pData->calculateKLDivergence(binCount);
                    this_energy = testKL + lambda*connectivityPotentialPhases(tempNumberOfComponents) + beta*componentPotential;

                    beads_in_use_tree.erase(original);
                    tempTotalContactEnergy = (afterRemoving - totalContactEnergy);

                    if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
                        currentComponentPotential=componentPotential;
                        isUpdated = true;
                    } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
                        currentComponentPotential=componentPotential;
                        isUpdated = true;
                    } else { // undo changes and move to next bead (rejecting)
                        beads_in_use_tree.insert(original);
                        beginIt = bead_indices.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.addNode(original, pModel);
                        if (componentIndex < totalComponents) {
                            Component * comp = &components[componentIndex];
                            comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, original, comp->getTotalContactSum()));
                            comp->addLatticePoint(original);
                        }
                    }
                } else {
                    eulerTour.addNode(original, pModel);
                    if (componentIndex < totalComponents) {
                        Component * comp = &components[componentIndex];
                        comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, original, comp->getTotalContactSum()));
                        comp->addLatticePoint(original);
                    }
                }

                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
            }
        } else { // positional refinement
            // flipping/swap mechanism
            // pick a non-component lattice point near anchor
            // only search within deadLimit, no need to recalculate at end
            bool isSwapped;
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            startTime = std::clock(); // randomly select an index to move
            int originalSwap2Value, position = rand() % workingLimit;
            int swap1 = bead_indices[ position ];
            // only use positions that are not in reduced model and anchor
            while( !( inComponents(swap1) && canRemoveIfAnchor(swap1)  ) ){ // only pick component
                position = rand() % workingLimit;
                swap1 = bead_indices[ position ];
            }

            std::vector<int> indices_to_check; //
            int totalToCheck;
            std::set<int> * tree = &beads_in_use_tree;
            for(int i=0; i < totalComponents; i++) { // the selected set of beads that make up each component will be held by Component object
                Component * comp = &components[i];
                if (comp->inUse(swap1)){
                    componentIndex=i;
                    totalToCheck = comp->getBeadsInUse()->size();
                    indices_to_check.resize(totalToCheck);
                    // I want to move a bead that is not part of reduced seed
                    // want to only move beads that are part of the newly defined component
                    std::copy(comp->getBeadsInUse()->begin(), comp->getBeadsInUse()->end(), indices_to_check.begin());
                    tree = comp->getBeadsInUse();
                    break;
                }
            }

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin()); // backup current state
            // check euler tour of main and components
            componentPotential=0;
            for(int i=0; i < totalComponents; i++) {
                componentPotential+=components[i].potential();
            }

            //remove nodes
            tempNumberOfComponents = eulerTour.removeNode(swap1);
            if (componentIndex < totalComponents) {
                Component * comp = &components[componentIndex];
                comp->setTotalContactSum(recalculateContactsPotentialSumRemove(comp->getBeadsInUse(), pModel, swap1, comp->getTotalContactSum()));
                comp->removeLatticePoint(swap1);
            }

            if (connectivityPotentialPhases(tempNumberOfComponents) == 0){
                // prefer to move index if it has too few contacts
                // remove contribution of swap1
                //int wldl = (numberOfCoolingTempSteps > low_temp_limit) ? deadLimit : (workingLimit + (int)((deadLimit-workingLimit)*0.31));
                double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);
                isSwapped = false;
                itIndex = bead_indices.begin() + position;
                // remove selected index from P(r)
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
                beads_in_use_tree.erase(swap1); // remove from in_use tree

                double newPotential, invWorkingLimit = 1.0/(double)(workingLimit); // break on first success

                int retry;
                for (int bb = 0; bb < 3; bb++) { // go through and find a better position within unused beads
                    // make the swap, sort and update P(r)
                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, swap1);
                    retry = 0;
                    while (  !(originalSwap2Value > -1) ) { // find a new position
                        originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, swap1);
                        retry++;
                        if (retry > 11){
                            break;
                        }
                    }

                    if (retry > 11){ // in case we have picked a point that has no available neighbors
                        break;
                    }

                    pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), originalSwap2Value);
                    //pSwap2 = (bead_indices.begin() + bb); // swapped point belongs to a component specified by componentIndex
                    //originalSwap2Value = *pSwap2;

                    if (numberOfContactsFromSet(tree, pModel, originalSwap2Value) > 0){ // should make this component specific?
                        newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
                        newPotential = etaConstant*newSum*invWorkingLimit;

                        beads_in_use_tree.insert(*pSwap2); // add new lattice

                        std::iter_swap(itIndex, pSwap2);
                        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                        addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        testKL = pData->calculateKLDivergence(workingBinCount);

                        tempNumberOfComponents=1;
                        if (componentIndex < totalComponents){
                            Component * comp = &components[componentIndex];
                            comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, originalSwap2Value, comp->getTotalContactSum()));
                            comp->addLatticePoint(originalSwap2Value);
                        }

                        componentContactSum=0;
                        for(int i=0; i < totalComponents; i++) {
                            Component * comp = &components[i];
                            componentContactSum += comp->getTotalContactSum()/(double)comp->getBeadsInUse()->size();
                        }
                        newPotential += etaConstant*componentContactSum;

                        this_energy = testKL + lambda * connectivityPotentialPhases(tempNumberOfComponents) + beta*componentPotential;
                        tempTotalContactEnergy = newPotential-totalContactEnergy;

                        if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            currentComponentPotential=componentPotential;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel); // adding node will always be in contact (so connectivity doesn't change)
                            std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                            break;
                        } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {
                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            currentComponentPotential=componentPotential;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
                            std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = bead_indices.begin() + position;
                        beads_in_use_tree.erase(originalSwap2Value);
                        if (componentIndex < totalComponents){
                            Component * comp = &components[componentIndex];
                            comp->setTotalContactSum(recalculateContactsPotentialSumRemove(comp->getBeadsInUse(), pModel, originalSwap2Value, comp->getTotalContactSum()));
                            comp->removeLatticePoint(originalSwap2Value);
                        }
                    }
                } // end of for loop

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                    std::sprintf(addRemoveText, "");
                    eulerTour.addNode(swap1, pModel);

                    if (componentIndex < totalComponents) {
                        Component * comp = &components[componentIndex];
                        comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, swap1, comp->getTotalContactSum()));
                        comp->addLatticePoint(swap1);
                        //cout << "IN COMP CONTACT SUM : " << comp->getTotalContactSum() << " " << calculateTotalContactSum(comp->getBeadsInUse(), pModel) << endl;
                    }
                }
            } else {

                eulerTour.addNode(swap1, pModel);
                if (componentIndex < totalComponents) {
                    Component * comp = &components[componentIndex];
                    comp->setTotalContactSum(recalculateContactsPotentialSumAdd(comp->getBeadsInUse(), pModel, swap1, comp->getTotalContactSum()));
                    comp->addLatticePoint(swap1);
                    //cout << "OUT COMP CONTACT SUM : " << comp->getTotalContactSum() << " " << calculateTotalContactSum(comp->getBeadsInUse(), pModel) << endl;
                }
            }
            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)  ){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

//        if (connectivityPotentialPhases(eulerTour.getNumberOfComponents()) > 0){
//            cout << " STOPPED TOO MANY COMPONENTS " << currentNumberOfComponents << " " << eulerTour.getNumberOfComponents() << endl;
//            string name = "failed_" + std::to_string(numberOfCoolingTempSteps);
//            //pModel->writeModelToFile(workingLimit, bead_indices, name);
//            return "Stopped";
//        }

        printf("       TEMP => %-.8f \n     ACCEPT => %.4f  FAILURES => %i\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop, acceptRate, failures, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
        printf("  UPDATECNT => %7i TIME => %.5f (SECONDS) \n", deadUpdate, runtime);
        printf("      GRAPH => %3i COMP POTENTIAL %.5f \n", currentNumberOfComponents, currentComponentPotential);
        printf(" LATTCE AVG => %.0f        STDEV => %0.3f  %s\n", runningAverage, runningVariance, addRemoveText);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f \n", totalContactEnergy, etaConstant, runningContactsSum);
        printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, (current_energy+totalContactEnergy));

        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;
        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;
            if (currentKL < lowestKL )
                lowestKL = currentKL;

            failures=0;
        } else {
            acceptRate = inv500*(499*acceptRate);
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        //update for running average
        sum_x_squared += workingLimit*workingLimit;
        sum_x += workingLimit;
        divideBy += 1.0;

        if (counter % deadUpdate == 0){ // if counter too small, add/remove may not sample sufficiently
            float inv_divideBy = 1.0/divideBy;
            average_x = sum_x*inv_divideBy;

            float sum_x_square_inv_divideBy = sum_x_squared*inv_divideBy;
            float average_x2 = average_x*average_x;
            if (sum_x_square_inv_divideBy > (1+average_x2)){
                stdev = sqrt(sum_x_square_inv_divideBy - average_x*average_x);
            } else {
                stdev = 2;
            }

            sum_x_squared = 0.0;
            sum_x = 0.0;
            divideBy= 0.0;

            runningAverage = 0.5*(average_x + runningAverage);
            runningVariance = stdev;
            // recalculate eta, increase last term to decrease weight of eta
            etaConstant *= etaFactor;
            componentContactSum=0;
            for(int i=0; i < totalComponents; i++) {
                Component * comp = &components[i];
                componentContactSum += comp->getTotalContactSum()/(double)comp->getBeadsInUse()->size();
            }
            totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit + componentContactSum);
            //populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
            std::normal_distribution<float> volumeGen(runningAverage, runningVariance);
        }
        counter++;
    } // end of steps

    // At end of each temp, update a probability model for volume?  Use this to select
    cout << "------------------------------------------------------------------------------" << endl;
    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    cout << "------------------------------------------------------------------------------" << endl;

    tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << " POSITIONAL SUCCESS/FAILURE " << ((double)failures/(double)attempts)*100 << endl;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    // perform positional refinement until delta E stabilizes?
    int totalLessThanAverage=0, numberContacts;
    for (int i=0; i<workingLimit; i++){
        numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        if (numberContacts < average_number_of_contacts ){
            totalLessThanAverage +=1;
        }
    }

    cout << " TOTAL LATTICE POINTS LESS THAN AVERAGE " << totalLessThanAverage << endl;
    for(int i=0; i < totalComponents; i++) {
        components[i].printConstraints(); // constraints the number of beads
    }
    //this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);

    //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);

    // if multithread, must put a lock on these two step
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    //pModel->setCVXHullVolume(calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel));

    pData->printKLDivergence(binCount);
    // write components
    for (int i=0; i<totalComponents; i++){
        components[i].writeToFile("component_" + std::to_string(i+1));
        cout << i << " TARGET NUMBER OF POINTS => " << components[i].getTargetNumberOfLatticePoints() << " ( " << components[i].getBeadsInUse()->size() << " )" << endl;
    }

    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            outputname,
            this,
            pData,
            numberOfCoolingTempSteps,
            cvx,
            average_number_of_contacts);


    return nameOfModel;
}
