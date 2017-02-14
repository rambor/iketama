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
    this->lowTempStop = highTempStartForCooling;

    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    // convert distances within the large search space to bins based on input P(R)-DATA file
    this->fillPrBinsAndAssignTotalBin( pBin,  pDistance,  totalDistancesInSphere,  pData);

    std::vector<float> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size

    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    const std::vector<int>::const_iterator trueModelBeginIt = pModel->getSeedBegin();
    //const std::vector<int>::const_iterator trueModelEndIt = pModel->getSeedEnd();
    int workingLimit = pModel->getTotalInSeed();

    // copy trueModel into BeadIndices
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    int minWorkingLimit = 0.17*workingLimit;

    int deadLimit = workingLimit;; // as bead indices are discarded, set upper limit of vector
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
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.end());

    std::random_device rd;
    std::mt19937 gen(rd());

    float inv_kb_temp = 1.0/0.000001;

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

    const int components = tempNumberOfComponents;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();

    // calculate KL against known
    // create custom D_KL function supplying Pr_of_PDB as target
    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    cout << "   DIRECT INITIAL D_KL : " << currentKL << endl;
    currentKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
    //float tempTotalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    cout << "          INITIAL D_KL : " << currentKL << endl;

    float energy_prior, testKL;
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::vector<int>::iterator beginIt = bead_indices.begin();
    std::vector<int>::iterator endIt = bead_indices.end();


    int lowestWorkingLimit = workingLimit, lowestDeadLimit = deadLimit;
    std::copy(beginIt, endIt, lowest_bead_indices.begin());

    int counter=1;
    int priorWorkingLimit;

    float inv_divideBy, average_x=0.5*workingLimit, stdev=0.2*workingLimit;
    float sum_x_squared=0, sum_x=0, divideBy=0, acceptRate = 0.5, inv500 = 1.0/500.0;
    //output for plotting
    bool isUpdated = false;

    int seedHighTempRounds = 2*highTempRounds;
    int deadUpdate = std::ceil(seedHighTempRounds*0.091);
    populatePotential(pModel->getSizeOfNeighborhood());
    float startContactSum=0.0;
    for (int i=0; i<lowestWorkingLimit; i++){
        startContactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        //startContactSum += numberOfContacts(bead_indices[i], &bead_indices, workingLimit, contactCutOff, pModel, pDistance);
    }

    double runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, workingLimit, pModel);

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

        if (distribution(gen) >= 0.5){ // ADD BEAD?
            cout << "*******************                  ADD                   *******************" << endl;

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
                if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) > 0){
                    itIndex = bead_indices.begin() + randomSpot;
                    // select first element from randomized active_set
                    // check if available neighbor can be added
                    newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);
                    afterAdding = etaConstant*newSum/(double)(workingLimit+1);
                    // make the swap at the border (workingLimit)
                    addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);
                    addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                    testKL = pData->calculateKLDivergence(binCount);
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
                        restoreAddingFromBackUp(&beginIt, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }

                // find a bead within workinglimit -> deadLimit
//                for (int d = workingLimit; d < deadLimit; d++) {
//                    //
//                    addMe = active_indices[d]; // remove from active_indices if used
//                    itIndex = find(beginIt + workingLimit, beginIt + deadLimit, addMe);
//                    //make the swap at the border (workingLimit)
//                    std::iter_swap(beginIt + workingLimit, itIndex);
//                    workingLimit++;
//                    std::copy(beginIt, beginIt + workingLimit, backUpState.begin());
//                    sort(beginIt, beginIt + workingLimit);
//
//                    addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
//
//                    testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
//                    // since neighbor list is already determined, don't need to check for contacts
//                    if (testKL < currentKL && numberOfContacts(addMe, &bead_indices, workingLimit, contactCutOff, pModel, pDistance) >= 1) {
//                        currentKL = testKL;
//                        isUpdated = true;
//                        break;
//                    } else if (exp(-(testKL - currentKL) * inv_kb_temp) > distribution(gen) && numberOfContacts(addMe, &bead_indices, workingLimit, contactCutOff, pModel, pDistance) >= 1) {
//                        currentKL = testKL;
//                        isUpdated = true;
//                        break;
//                    } else { // undo changes (rejecting)
//                        std::copy(backUpState.begin(), backUpState.begin() + workingLimit, beginIt);
//                        workingLimit--;
//                        std::iter_swap(beginIt + workingLimit, itIndex); //reverse swap
//                        copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//                    }
//                }
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
                testKL = pData->calculateKLDivergence(binCount);

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

//            std::copy(beginIt, beginIt + workingLimit, active_indices.begin());
//            shuffle(active_indices.begin(), active_indices.begin() + workingLimit, gen);
//            // collect indices that are not part of PDBModel
//            while (priorWorkingLimit > 0){
//                // go through entire list
//                for (int i=0; i < workingLimit; i++) { // removing bead always decreases total volume?
//                    // grab from randomized active_indices list
//                    original = active_indices[i];
//                    itIndex = find(beginIt, beginIt + workingLimit, original);
//                    // remove original from P(r)
//                    copy(beginBinCount, endBinCount, binCountBackUp.begin()); //backup binCount
//                    copy(beginIt, beginIt + deadLimit, backUpState.begin());
//                    removeFromPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
//                    // make the swap
//                    workingLimit--;
//                    std::iter_swap(itIndex, beginIt + workingLimit);
//                    // still need to sort, swap changes the order
//                    sort(beginIt, beginIt+workingLimit);
//                    isConnected = isConnectedComponent(&bead_indices, workingLimit, pDistance, pModel->getTotalNumberOfBeadsInUniverse(), tempNumberOfComponents);
//                    // don't remove lattice point if model is not interconnected
//                    //cout << "IsConnected? " << isConnected << " " << tempNumberOfComponents << endl;
//                    //isConnected = true;
//                    if (isConnected){
//                        testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB);
//                        if (testKL < currentKL) {
//                            currentKL = testKL;
//                            //cout << "REMOVED " << endl;
//                            //std::iter_swap(active_indices.begin() + i, active_indices.begin() + workingLimit); // move
//                            isUpdated = true;
//                            break;
//                        } else if (exp(-(testKL - currentKL)*inv_kb_temp) > distribution(gen)){
//                            //cout << "REMOVED SA" << endl;
//                            currentKL = testKL;
//                            //std::iter_swap(active_indices.begin() + i, active_indices.begin() + workingLimit);
//                            isUpdated = true;
//                            break;
//                        } else { // undo changes and move to next bead (rejecting)
//                            workingLimit++;
//                            //sort(beginIt, beginIt+workingLimit); // swapped index is at workingLimit
//                            copy(backUpState.begin(), backUpState.begin() + deadLimit, beginIt);
//                            copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//                        }
//
//                    } else { // undo changes and move to next bead (rejecting)
//                        workingLimit++;
//                        //sort(beginIt, beginIt+workingLimit); // swapped index is at workingLimit
//                        copy(backUpState.begin(), backUpState.begin() + deadLimit, beginIt);
//                        copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); //copy to bin count
//                    }
//
//                } // decrement if I make it through whole list or found a success
//                priorWorkingLimit--;
//            }
        }

//          isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,tempNumberOfComponents);  // 10x longer than contactEnergy calculation
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices) || tempNumberOfComponents != 1 ){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            cout << " COMP " << tempNumberOfComponents << endl;
//            return;
//        }
//        if (tempNumberOfComponents > 1){
//            cout << " STOPPED TOO MANY COMPONENTS " << currentNumberOfComponents << endl;
//            string name = "failed_" + std::to_string(high);
//           // pModel->writeModelToFile(workingLimit, bead_indices, name);
//            return;
//        }

        if (current_energy < lowestE){
        //if (currentKL < lowestE){
            //lowestE = currentKL;
            lowestE = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            std::copy(beginIt, endIt, lowest_bead_indices.begin());
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_");
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


    cout << "KL DIVERGENCE : " << endl;
    cout << "  INITIAL D_KL => " << startKL << endl;
    cout << "   LOWEST D_KL => " << lowestE << endl;
    cout << "AVERAGE CONTACTS : (per lattice point)" << endl;
    cout << "       INITIAL => " << startContactSum << endl;
    cout << "         FINAL => " << average_number_of_contacts << " " << totalContactEnergy << endl;
    cout << " Contacts Per Bead " << contactsPerBead << endl;

    // this is fixed model for initial high temp search?
    // set this as the seed
    //    pModel->setBeadAverageAndStdev(average_x, stdev);
    //    pModel->setReducedSeed(workingLimit, bead_indices);
    //    pModel->writeModelToFile(workingLimit, bead_indices, name);
    //    pModel->setBeadAverageAndStdev(average_x, stdev);

    // using lowestEnergy model, organize beads of the entire Model search space.
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    bead_indices.resize(totalBeadsInSphere);
    beginIt = bead_indices.begin();
    endIt = bead_indices.end();

    // remap to bead universe for storage and printing
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
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);

    pModel->setReducedSeed(workingLimit, bead_indices);
    //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name);

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