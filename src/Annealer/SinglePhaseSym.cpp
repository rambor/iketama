//
// Created by Robert Rambo on 13/01/2017.
//
#include "../Anneal.h"
#include "Data.h"
#include "../Model.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"

using namespace std;

/**
 * create initial model of complete object
 */
bool Anneal::createInitialModelSymmetry(Model *pModel, Data *pData) {

    srand(time(0));
    contactCutOff = interconnectivityCutOff;

    //float inv_kb_temp = 1.0/highT;
    float inv_kb_temp = 1.0/0.0001;

    // convert distances to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins();

    maxbin=0;
    int violations;
    for(int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    // adjust this part
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();
    int deadLimit; // as bead indices are discarded, set upper limit of vector
    std::vector<int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);

    const int num = totalBeadsInSphere;
    int * ptr = (num != 0) ? &subUnit_indices.front() : NULL;
    for(int i = 0; i < num; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    float invBeadVolume = 1.0/pModel->getBeadVolume();

    //float inv_kb = totalBeadsInSphere*beadVolume/(pData->getShannonBins());
    float inv_sym = 1.0/pModel->getNumberOfSubUnits();

    int lowerN = round(lowerV*invBeadVolume*inv_sym);
    int upperN = round(upperV*invBeadVolume*inv_sym);

    int alterMe;

    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN);   // number of beads in ASU

    int subUnitWorkingLimit = number_of_beads_to_use(gen);
    pModel->writeModelToFile(totalBeadsInSphere, subUnit_indices, "universe", 0);

    cout << "Bead Search Limited to: " << lowerN << " <= N <= " << upperN << endl;
    cout << "      INITIAL MODEL WL: " << subUnitWorkingLimit << endl;
    cout << "              SYMMETRY: " << pModel->getSymmetryString() << endl;
    cout << "        TOTAL SUBUNITS: " << pModel->getNumberOfSubUnits() << endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    std::vector<int>::iterator beginIt = subUnit_indices.begin();
    std::vector<int>::iterator endIt = subUnit_indices.end();

    cout << "        CREATING INITIAL RANDOM MODEL " << endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    std::set<int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);

    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    //vertexT * vertices;  // not sure if this needs to be explicitly deleted
    //int indexOfHullpt, count, potentialContacts, index;
    char flags[25];
    sprintf(flags, "qhull s FA");
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, points, pModel);
    //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
    populateLayeredDeadlimit(subUnit_indices.begin(), subUnitWorkingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    // calculate starting energy

    std::clock_t start;
    start = std::clock();

    // fill binCount for first time
    float currentKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
    lambda = 0.001;
    mu = 0.000001; // scale this to volume of search space

    cout << " THREAD TIME "  << (std::clock() - start)/(double) CLOCKS_PER_SEC << " seconds " << endl;
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");
    // BIN calculation is  0.000865
    // HULL calculation is 0.0015 seconds

    float testKL, test_energy; // sets alpha as a constant during high temp eq

    int lowestWorkingLimit, lowestDeadLimit;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    int tempNumberOfComponents, currentNumberOfComponents;
    EulerTour eulerTour(subUnit_indices.begin(), subUnitWorkingLimit, pModel);
    currentNumberOfComponents = eulerTour.getNumberOfComponents();
//    isConnectedComponent(&subUnit_indices, subUnitWorkingLimit, pDistance, totalBeadsInSphere, currentNumberOfComponents, pModel);

    float current_energy = currentKL + lambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + mu*current_volume/(float)subUnitWorkingLimit;

    cout << " INITIAL => ENERGY :  " << current_energy << endl;
    cout << " INITIAL =>   D_KL : " <<  currentKL << endl;
    cout << " INITIAL => VOLUME : " << current_volume << endl;
    cout << "        VIOLATIONS : " << violations << endl;
    cout << " INITIAL        WL : " << subUnitWorkingLimit << " DL : " << deadLimit << endl;

    float lowest_energy = current_energy;

    std::vector<int>::iterator itIndex;
    std::vector<int>::iterator beginBinCount = binCount.begin();
    std::vector<int>::iterator endBinCount = binCount.end();

    //TEST
    //std::uniform_int_distribution<> randomTest (0, subUnitWorkingLimit-1);
    //for (int i=0; i<10; i++){
    //    int testI = randomTest(gen);
    //    removeFromPrSym(subUnit_indices[testI], subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
    //    addToPrSym(subUnit_indices[testI], subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
    //    testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
    //    if (pData->calculateKLDivergence(binCount) != testKL){
    //        cout << i  << " FAiled " << pData->calculateKLDivergence(binCount) << " != " << testKL << endl;
    //        cout << subUnitWorkingLimit << " => " << "index " << testI << " => " << subUnit_indices[testI] << endl;
    //        return;
    //    }
    //}

    char addRemoveText[50];
    std::vector<int> averageV(highTempRounds);
    std::vector<int>::iterator pSwap2;
    int volumeCount=0, testIndex, originalSwap2Value;
    float volumeSum=0, workingLimitSum=0;

    float tempViolations, totalViolations=0;
    beta = 0.0001;
    beta = 0;
    for(int i=0; i<subUnitWorkingLimit; i++){
        totalViolations += symViolationsPotential(subUnit_indices[i], &subUnit_indices, subUnitWorkingLimit, pModel);
    }
    current_energy += beta * totalViolations;

    cout << "Total Violations " << totalViolations << endl;
    bool updateCVX = false;

    int high;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
        std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy

        if (distribution(gen) < 0.31){
            // create points from workSet to determine HULL
            float invSubUnitworkingLimit = mu/(float)subUnitWorkingLimit;
            bool isSwapped;

            sprintf(addRemoveText, "POSITIONAL");

            std::vector<int>::iterator pSwap2;
            // recalculate
            for (int i = 0; i < subUnitWorkingLimit; i++) {
                beadToPoint(&points[i*3], pModel->getBead(subUnit_indices[i]));
                active_indices[i] = subUnit_indices[i];
            }

            // needs to be optimized
            qh_new_qhull(3, subUnitWorkingLimit, points, 0, flags, NULL, NULL);
            vertexT * vertices = qh vertex_list;
            int totalV = qh num_vertices;
            // only move CVX hull points
            std::vector<int> indices_to_check(totalV); // large vector ~1000's
            for (int v = 0; v < totalV; v++) { //
                indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                vertices = vertices->next;
            }
            qh_freeqhull(true);

            std::shuffle(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.begin()+deadLimit, gen); // randomize possible swap space
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            int swapCount = 0, swap1, bb;

            for (int v = 0; v < 1; v++) { // try 3 times

                swap1 = indices_to_check[v];
                isSwapped = false;
                // find bead to swap in active set
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);
                // remove selected index from P(r)
                std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
                //removeFromPr(swap1, subUnit_indices, subUnitWorkingLimit, pBin, totalBeadsInSphere, binCount);
                removeFromPrSym(swap1, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
                std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // make copy of altered P(r)
                beads_in_use_tree.erase(swap1);

                for (bb = 0; bb < sizeOfNeighborhood; bb++) {
                    int neighbor = subUnit_indices[subUnitWorkingLimit+bb];
                    // make the swap, sort and update P(r)
                    pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.begin()+deadLimit, neighbor);
                    std::iter_swap(itIndex, pSwap2);
                    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit); // bead_indices needs to be sorted

                    addToPrSym(neighbor, subUnit_indices, subUnitWorkingLimit, workingBinCount, pModel, pData);
                    testKL = pData->calculateKLDivergence(workingBinCount);
                    //testKL = calculateKLEnergySymmetry(&subUnit_indices, &workingBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
                    //tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
                    isConnectedComponent(&subUnit_indices, subUnitWorkingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                    test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, points, pModel);
                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + test_volume*invSubUnitworkingLimit;

                    if (test_energy < current_energy) {
                        beads_in_use_tree.insert(neighbor);
                        //move swapped position to the deadLimit
                        std::iter_swap(pSwap2, subUnit_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                        std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        current_energy = test_energy;
                        isSwapped = true;
                        updateCVX = true;
                        swapCount++;
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        beads_in_use_tree.insert(neighbor);
                        std::iter_swap(pSwap2, subUnit_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin());
                        std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        current_energy = test_energy;
                        isSwapped = true;
                        updateCVX = true;
                        swapCount++;
                        break;
                    }
                    // reverse changes; find index in sorted list and replace
                    std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                    std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                    itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit, swap1); // find swapped value in the sorted array
                }

                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    //eulerTour.addNode(swap1, pModel);
                    beads_in_use_tree.insert(swap1);
                }
            }


        } else { // add or remove

            // dead limit is set by convex hull
            alterMe = number_of_beads_to_use(gen); // pick number between lower and upper bounds
            beginIt = subUnit_indices.begin();
            endIt = subUnit_indices.end();

            if (subUnitWorkingLimit > alterMe && alterMe > lowerN) { // REMOVE beads from sorted list into useable range < deadLimit
                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                sprintf(addRemoveText, "  REMOVE  ");

                int randomSpot = rand() % subUnitWorkingLimit;
                testIndex = subUnit_indices[randomSpot];

                tempViolations = (totalViolations - symViolationsPotential(testIndex, &subUnit_indices, subUnitWorkingLimit, pModel)); // pointless since calculation is a removal
                removeLatticePositionToModelSym(&beginIt, subUnit_indices, binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);

                testKL = pData->calculateKLDivergence(binCount);
                //testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);
                test_volume = calculateCVXHULLVolume(flags,
                                                     &subUnit_indices,
                                                     subUnitWorkingLimit,
                                                     points,
                                                     pModel);

                //tempNumberOfComponents = eulerTour.removeNode(testIndex);
                isConnectedComponent(&subUnit_indices, subUnitWorkingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                test_energy = testKL + lambda * (tempNumberOfComponents - 1) * (tempNumberOfComponents - 1) +
                              mu * test_volume / (float) subUnitWorkingLimit + beta*tempViolations;

                if (test_energy < current_energy) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = tempViolations;
                    beads_in_use_tree.erase(testIndex);
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = tempViolations;
                    beads_in_use_tree.erase(testIndex);
                } else { // undo changes and move to next bead (rejecting)
                    //currentNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                    restoreRemovingLatticePointFromBackUp(&beginIt, &subUnitWorkingLimit, &binCountBackUp,
                                                          &beginBinCount);
                }

            } else { // ADD beads
                // randomly add from deadLimit
                sprintf(addRemoveText, "    ADD   ");
                int randomSpot = rand() % (deadLimit - subUnitWorkingLimit) + subUnitWorkingLimit;

                testIndex = subUnit_indices[randomSpot];
                itIndex = beginIt + randomSpot;

                addLatticPositionToModel(&beginIt, &endIt, &backUpState, &subUnitWorkingLimit, &itIndex);
                addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);

                tempViolations = (totalViolations + symViolationsPotential(testIndex, &subUnit_indices, subUnitWorkingLimit, pModel));
                testKL = pData->calculateKLDivergence(binCount);

                test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, points, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                //tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                isConnectedComponent(&subUnit_indices, subUnitWorkingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1)
                              + mu*test_volume/(float)subUnitWorkingLimit
                              + beta*tempViolations;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = tempViolations;
                    beads_in_use_tree.insert(testIndex);
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = tempViolations;
                    beads_in_use_tree.insert(testIndex);
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&beginIt, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                    //currentNumberOfComponents = eulerTour.removeNode(testIndex);
                }
            }

        } // END OF ADD?REMOVE

        if (updateCVX){
            populateLayeredDeadlimitUsingSet(&subUnit_indices, &beads_in_use_tree, subUnitWorkingLimit, &deadLimit, pModel);
            updateCVX=false;
        }
        //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
        //populateLayeredDeadlimitUsingSet(&subUnit_indices, &beads_in_use_tree, subUnitWorkingLimit, &deadLimit, pModel);
        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts

        printf("*******************             %s                 ******************* \n", addRemoveText);
        printf("   MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        printf("      GRAPH => %3i\n", currentNumberOfComponents);
        printf("      UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("     VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, mu, mu*current_volume);
        printf(" VIOLATIONS => %f  BETA => %.4E  \n", totalViolations, beta);
        printf("LIMIT: %5i DEADLIMIT: %5i D_KL: %.4E ENRGY: %.4E \n", subUnitWorkingLimit, deadLimit, currentKL, current_energy);
        cout << "*******************                                        *******************" << endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            //std::string name = "initial_sym_search_" + std::to_string(high);
            //pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, name, high);
            //pModel->writeModelToFile(deadLimit, subUnit_indices, name+"_skin", high);
            //pModel->writeSymModelToFile(subUnitWorkingLimit, subUnit_indices, name+"sym_");
            //counter++;
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;
            lowest_energy = current_energy;

            workingLimitSum += subUnitWorkingLimit;
            volumeSum += current_volume;
            volumeCount++;
        }

//        if (checkForRepeats(subUnit_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << subUnitWorkingLimit << " D_KL " << currentKL << " <=> " << endl;
//            return;
//        }
//        testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (currentKL != testKL){
//            // if this fails, update of P(r) is wrong or subUnit_indices is corrupted
//            cout << high  << "           Failed " << pData->calculateKLDivergence(binCount) << " != " << testKL << endl;
//            cout << high  << " Failed currentKL " << currentKL << " != " << testKL << endl;
//            return;
//        }

    } // end of HIGH TEMP EQUILIBRATION

    pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage, volumeStdev;

    volumeAverage = workingLimitSum/(float)volumeCount;
    volumeStdev = 0.17*volumeAverage;
    // remove points close to hull
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, "initial_CVX_sym_" + filenameprefix, high);
    pModel->setStartingSet(lowest_subUnit_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);

    pModel->setBeadAverageAndStdev(1.37*pModel->getStartingWorkingLimit(), 0.17*1.37*pModel->getStartingWorkingLimit());

    //std::normal_distribution<float> volumeGen(beadAverage, beadStDev);
    cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << endl;
    cout << "*******************                                        *******************" << endl;
    cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << endl;
    printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    cout << "*******************                                        *******************" << endl;

    if (currentNumberOfComponents ==1 ) {
        return true;
    } else {
        cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << endl;
        cout << "INCREASE highTempRounds, g" << endl;
        return false;
    }
}


/**
 *
 */
string Anneal::refineSymModel(Model *pModel, Data *pData, std::string nameTo){

    srand(time(0));
    cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY" << endl;
    //pModel->setBeadAverageAndStdev(pModel->getStartingWorkingLimit(), 0.17*pModel->getStartingWorkingLimit());
    float oldN = pModel->getVolumeAverage();  // need to reset this for additional modeling
    float oldStdev = pModel->getVolumeStdev();
    float runningAverage = pModel->getVolumeAverage();
    float runningVariance = pModel->getVolumeStdev();
    std::normal_distribution<float> volumeGen(runningAverage, runningVariance);

    this->lowTempStop = highTempStartForCooling;
    int swap1, alterMe, totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());

    // make copy of bead_indices
    std::vector<int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);

    std::vector<int>::iterator beginIt, endIt, itIndex, pSwap2;

    // copy Starting_Set from initial model
    pModel->setModelBeginStartingSetIterator(beginIt);
    pModel->setModelEndStartingSetIterator(endIt);
    std::copy(beginIt, endIt, bead_indices.begin());

    int workingLimit = pModel->getStartingWorkingLimit();

    // reset iterators to internal bead_indices
    beginIt = bead_indices.begin();
    endIt = bead_indices.end();

    //pModel->centerLatticeModel(workingLimit, bead_indices);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);


    int deadLimit;
    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    //int deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
    //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);

    // convert distances in Search Sphere to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    maxbin=0;
    for(int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

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

    int violations;
    float testKL, currentKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );

    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

    //CVX HULL STUFF
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];

    //int dim = 3;
    char flags[25];
    sprintf(flags, "qhull s FA");

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;
    int original, addMe;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float inv_kb_temp, relative_diff, tempAverageContacts;
    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;
    float step_limit = (updateCount < 10000) ? 27000 : updateCount;
    int deadUpdate = std::ceil(updateCount*0.091);
    float low_temp_limit = step_limit*0.91;

    std::vector<float> tempDuringRun(step_limit);
    std::vector<float> divergenceDuringRun(step_limit);
    std::vector<int> workingLimitDuringRun(step_limit);

    float * pTempDuringRun = &tempDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();
    //bool isConnected = isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, interconnectivityCutOff, tempNumberOfComponents, pModel);
    bool isUpdated = false;
    float sum_x_squared=0, sum_x=0, divideBy=0, average_x, stdev, acceptRate = 0.5, inv500 = 1.0/500.0;

    int counter=1, randomSpot;

    inv_kb_temp = 1.0/this->lowTempStop ;
    int deltaECount=0, failures=0;
    float this_energy, lowestKL = currentKL;
    char addRemoveText[50];
    char titleText[50];

    double newSum, runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, workingLimit, pModel);

    double etaConstant = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) + 4 );
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);
    double etaFactor = 1.0/std::pow(10000, 1.0/(step_limit/(float)deadUpdate));

    double percentAddRemoveSwitch = 0.51;

    // want D_KL to settle < 10^-5 to 10^-6
    // int span = (int)(std::floor(log10f(currentKL)) - std::floor(log10f(totalContactEnergy/(float)workingLimit))) + 3;

    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    //
    //

    float tempViolations, totalViolations=0, afterAdding, afterRemoving;
    beta = 0.000001; // violations
    beta=0;
    for(int i=0; i<workingLimit; i++){
        totalViolations += symViolationsPotential(bead_indices[i], &bead_indices, workingLimit, pModel);
    }
    current_energy += beta*totalViolations;

    int numberOfCoolingTempSteps;
    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        beginIt = bead_indices.begin();
        endIt = bead_indices.end();
        std::copy(beginIt, endIt, backUpState.begin());
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

        if (distribution(gen) < percentAddRemoveSwitch ){ //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            sprintf(titleText, "ADD?REMOVE");
            alterMe = (int) volumeGen(gen);
            // build a list of indices within defined region of convex hull
            if (alterMe > workingLimit){ // ADD BEAD?

                sprintf(addRemoveText, "     ADD => %i", 1);
                randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                //
                // rather than try to find a bead to add, see if the one random one is acceptable
                //
                addMe = bead_indices[randomSpot]; // remove from active_indices if used
                // calculate local potential for beads in that neighbor
                //itIndex = std::find(beginIt + workingLimit, beginIt + deadLimit, addMe);
                itIndex = bead_indices.begin() + randomSpot;
                // select first element from randomized active_set
                // check if available neighbor can be added
                newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, addMe, runningContactsSum);
                afterAdding = etaConstant*newSum/(double)(workingLimit+1);
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);
                addToPrSym(addMe, bead_indices, workingLimit, binCount, pModel, pData);

                tempViolations = totalViolations + symViolationsPotential(addMe, &bead_indices, workingLimit, pModel);
                testKL = pData->calculateKLDivergence(binCount);

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                //tempNumberOfComponents = eulerTour.addNode(original, pModel);
                tempNumberOfComponents = 1;
                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*tempViolations;

                beads_in_use_tree.insert(addMe);
                tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    totalViolations = tempViolations;
                    isUpdated = true;
                    eulerTour.addNode(addMe, pModel);
                } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    totalViolations = tempViolations;
                    isUpdated = true;
                    eulerTour.addNode(addMe, pModel);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(addMe);
                    restoreAddingFromBackUp(&beginIt, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                // test for deletion
                sprintf(addRemoveText, "     REMOVE     ");

                randomSpot = std::rand() % workingLimit;
                original = bead_indices[randomSpot];
                tempNumberOfComponents = eulerTour.removeNode(original);

                if (tempNumberOfComponents == 1){

                    newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                    afterRemoving = etaConstant*newSum/(double)(workingLimit-1);
                    removeLatticePositionToModelSym(&beginIt, bead_indices, binCount, &workingLimit, &original, pModel, pData);

                    tempViolations = totalViolations - symViolationsPotential(original, &bead_indices, workingLimit, pModel);
                    testKL = pData->calculateKLDivergence(binCount);
                    this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*tempViolations;

                    beads_in_use_tree.erase(original);
                    tempTotalContactEnergy = (afterRemoving - totalContactEnergy);

                    if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
                        totalViolations = tempViolations;
                        isUpdated = true;
                    } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
                        totalViolations = tempViolations;
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

        } else { // positional refinement

            // only search within deadLimit, no need to recalculate at end
            sprintf(titleText, "POSITIONAL");
            sprintf(addRemoveText, "");
            bool isSwapped;
            std::shuffle(beginIt+workingLimit, beginIt + deadLimit, gen); // randomizes order of beads to select from
            std::copy(beginIt, endIt, backUpState.begin());
            //float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
            // find bead to swap in active set
            // only refine a single position with 2 or less contacts
            int originalSwap2Value;
            int wldl = (numberOfCoolingTempSteps > low_temp_limit) ? deadLimit : (workingLimit + (int)((deadLimit-workingLimit)*0.31));

            // randomly select an index to move
            int position = std::rand() % workingLimit;
            swap1 = bead_indices[ position ];

            float currentViolations;

            if (eulerTour.removeNode(swap1) == 1){
                double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);
                isSwapped = false;
                // swap
                itIndex = bead_indices.begin() + position;
                //itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1); // locate lattice within workingLimit
                // remove selected index from P(r)
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
                // violations in current position
                currentViolations = symViolationsPotential(swap1, &bead_indices, workingLimit, pModel);

                beads_in_use_tree.erase(swap1);
                double newPotential, invWorkingLimit = 1.0/(double)(workingLimit);

                // find better position
                for (int bb = workingLimit; bb < wldl; bb++) {
                    pSwap2 = (bead_indices.begin() + bb);
                    originalSwap2Value = *pSwap2;

                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, originalSwap2Value) > 0){
                        newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
                        newPotential = etaConstant*newSum*invWorkingLimit;

                        beads_in_use_tree.insert(*pSwap2); // add new lattice

                        std::iter_swap(itIndex, pSwap2);
                        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                        addToPrSym(originalSwap2Value, bead_indices, workingLimit, workingBinCount, pModel, pData);

                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        tempViolations = totalViolations + symViolationsPotential(originalSwap2Value, &bead_indices, workingLimit, pModel) - currentViolations;

                        testKL = pData->calculateKLDivergence(workingBinCount);
                        tempNumberOfComponents=1;
                        this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents) + beta* tempViolations;
                        tempTotalContactEnergy = newPotential-totalContactEnergy;

                        if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
                            totalViolations = tempViolations;
                            break;
                        } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
                            totalViolations = tempViolations;
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = bead_indices.begin() + position;
                        beads_in_use_tree.erase(originalSwap2Value);
                    }
                } // end for loop

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                    sprintf(addRemoveText, "");
                    eulerTour.addNode(swap1, pModel);
                }
            } else {
                sprintf(addRemoveText, "EULER TOUR CRITICAL PT");
                eulerTour.addNode(swap1, pModel);
            }


//            for (int i=0; i<0*workingLimit; i++){
//
//                swap1 = active_indices[i];
//                numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);
//
//                if (numberContacts < 3 ){
//                    isSwapped = false;
//                    // swap
//                    itIndex = find(beginIt, beginIt + workingLimit, swap1); // locate lattice within workingLimit
//
//                    // remove selected index from P(r)
//                    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
//                    removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
//                    std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
//
//                    // violations in current position
//                    currentViolations = symViolationsPotential(swap1, &bead_indices, workingLimit, pModel);
//                    // current local potential
//                    localPotentialAtOldPosition = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1) +
//                                                       calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, swap1));
//                    beads_in_use_tree.erase(swap1);
//                    localPotentialAtOldPositionWithOutBead = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1));
//
//                    // find better position
//                    for (int bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
//                        // make the swap, sort and update P(r)
//                        pSwap2 = &bead_indices[bb]; // located in workingZone
//                        originalSwap2Value = *pSwap2;
//
//                        localPotentialAtNewPosition = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, *pSwap2));
//                        beads_in_use_tree.insert(*pSwap2);
//
//                        if ( numberOfContactsFromSet(&beads_in_use_tree, pModel, *pSwap2) > 0){
//
//                            localPotentialAtNewPositionWithBead = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, *pSwap2) +
//                                                                       calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, *pSwap2));
//
//                            //*itIndex = *pSwap2;         // swapped here (at this point bead_indices contains duplicate values
//                            std::iter_swap(itIndex, pSwap2);
//                            std::sort(beginIt, beginIt + workingLimit); // bead_indices needs to be sorted
//
//                            addToPrSym(originalSwap2Value, bead_indices, workingLimit, workingBinCount, pModel, pData);
//
//                            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
//                            tempViolations = totalViolations + symViolationsPotential(originalSwap2Value, &bead_indices, workingLimit, pModel) - currentViolations;
//
//                            testKL = pData->calculateKLDivergence(workingBinCount);
//
//                            isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,
//                                                 tempNumberOfComponents
//                            );
//                            // 10x longer than contactEnergy calculation
//                            // new state => localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead
//                            // old state => localPotentialAtOldPosition + localPotentialAtNewPosition
//                            this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents) + beta* tempViolations;
//
//                            tempTotalContactEnergy = localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead - localPotentialAtOldPosition - localPotentialAtNewPosition;
//
//                            if ((this_energy + localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead) < (current_energy + localPotentialAtOldPosition + localPotentialAtNewPosition)) {
//
//                                copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
//                                currentKL = testKL;
//                                current_energy = this_energy;
//                                currentNumberOfComponents = tempNumberOfComponents;
//                                totalContactEnergy = tempTotalContactEnergy;
//                                isUpdated = true;
//                                isSwapped = true;
//                                totalViolations = tempViolations;
//                                break;
//
//                            } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {
//
//                                copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
//                                currentKL = testKL;
//                                current_energy = this_energy;
//                                currentNumberOfComponents = tempNumberOfComponents;
//                                totalContactEnergy = tempTotalContactEnergy;
//                                isUpdated = true;
//                                isSwapped = true;
//                                totalViolations = tempViolations;
//                                break;
//                            }
//
//                            // reverse changes; find index in sorted list and replace
//                            // pSwap2 points to &bead_indices[bb]
//                            // itIndex
//                            // binCount is P(r) distribution without swap1 value
//                            std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
//                            std::copy(backUpState.begin(), backUpState.begin() + (bb+1), beginIt);
//                            itIndex = find(beginIt, beginIt + workingLimit, swap1);
//                        }
//                        beads_in_use_tree.erase(originalSwap2Value);
//                    } // inner for loop exit on break
//
//                    // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
//                    if (!isSwapped){ // itIndex already reverses on exit from loop
//                        std::sort(beginIt, beginIt + workingLimit);
//                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
//                        beads_in_use_tree.insert(swap1); // add back swap one to tree
//                    }
//
//                    break;
//                }
//
//            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;
//        float testKL1 = calculateKLEnergySymmetry(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

        cout << "______________________________________________________________________________" << endl;
        printf("*******************               %s               *******************\n", titleText);
        cout << "*******************                                        *******************" << endl;
        printf("       TEMP => %-.8f \n     ACCEPT => %.4f  FAILURES => %i\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop , acceptRate, failures, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
        printf("      GRAPH => %3i\n", currentNumberOfComponents);
        printf("  UPDATECNT => %7i \n", deadUpdate);
        printf(" LATTCE AVG => %.0f        STDEV => %0.3f  %s\n", runningAverage, runningVariance, addRemoveText);
        printf("   CONTACTE => %.4E ETA => %.4E \n", totalContactEnergy, etaConstant);
        printf(" VIOLATIONS => %f  BETA => %.4E  \n", totalViolations, beta);
        printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, current_energy);

        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;
            if (currentKL < lowestKL ){
                lowestKL = currentKL;
            }
//            populateLayeredDeadlimit(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
            populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
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
            if (sum_x_square_inv_divideBy > average_x2){
                stdev = sqrt(sum_x_squared*inv_divideBy - average_x*average_x);
            } else {
                stdev = 2;
            }

            sum_x_squared = 0.0;
            sum_x = 0.0;
            divideBy= 0.0;

            runningAverage = 0.5*(average_x + runningAverage);
            runningVariance = stdev;

            etaConstant *= etaFactor;
            totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);

            percentAddRemoveSwitch = percentAddRemove;
            std::normal_distribution<float> volumeGen(runningAverage, runningVariance);
        }
        counter++;
        // check for early termination?
    } // end of steps

    // At end of each temp, update a probability model for volume?  Use this to select

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (int i=0; i<workingLimit; i++){
        tempAverageContacts += numberOfContacts(bead_indices[i], &bead_indices, workingLimit, pModel, pDistance);
    }
    tempAverageContacts =tempAverageContacts/(float)workingLimit;

    // perform positional refinement until delta E stabilizes?
    // final round to move any points close to body
    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
    //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);

    pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");

    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);

    //pModel->writeModelToFile(workingLimit, bead_indices, "refined_");
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");

    //  pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_" + filenameprefix);
    //  pModel->writeSymModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_sym_");

    return nameOfModel;
}



/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 *
 */
inline float Anneal::calculateKLEnergySymmetryFast(vector<int> *subUnit_indices, std::vector<int> *binCount, int indicesWorkingLimit, int totalBeadsInSphere, int &violation, Model *pModel, Data *pData, float &klDivergence) {


    // calculate KLDivergence
    //float squared = subUnits*subUnits*indicesWorkingLimit*indicesWorkingLimit;
    //float ratio = (float)violation/divisor;
    //cout << " RATIO " << ratio << " " << violation << endl;

    klDivergence = pData->calculateKLDivergence(*binCount);
    return klDivergence;
}



/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
inline float Anneal::calculateKLEnergySymmetry(vector<int> *subUnit_indices, std::vector<int> *binCount, int indicesWorkingLimit, int totalBeadsInSphere, int &violation, Model *pModel, Data *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    float limit = (pModel->getBeadRadius())*2.0;
    violation=0;
    //float divisor = totalBeadsInSphere*(totalBeadsInSphere-1)*0.5;

    int subUnits = (int)pModel->getNumberOfSubUnits();
    int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create first subunit from selected indices
    for (int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    int count = indicesWorkingLimit;
    for (int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetry(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    int next_i;
    for (int i = 0; i < totalCoordinates; i++) {
        tempVec1 = &coordinates[i];
        for (next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();

            (*binCount)[pData->convertToBin(distance_to)]++; // some distances will exceed dmax

            if (distance_to < limit){
                violation++;
            }
        }
    }
    // calculate KLDivergence
    //float squared = subUnits*subUnits*indicesWorkingLimit*indicesWorkingLimit;
    //float ratio = (float)violation/divisor;
    //cout << " RATIO " << ratio << " " << violation << endl;
    return pData->calculateKLDivergence(*binCount);
}


/**
 * beadsInUse must be sorted
 */
inline void Anneal::addToPrSym(int addMeSubUnitIndex, std::vector<int> & beadsInUse, int workingLimit, std::vector<int> & prBins, Model *pModel, Data *pData){
    // each bead has a sym related partner in pModel
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    float totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * refVec, * targetVec;
    int indexInMasterBeadArray = pModel->mapSubUnitIndex(addMeSubUnitIndex);
    int next;
    int stopAt=0;

    int * ptr = &beadsInUse.front();
    for(int i=0; i<workingLimit; i++){
        if (ptr[i] == addMeSubUnitIndex){
            stopAt = i;
            break;
        }
    }

    //
    // add interdomain distances
    //
    refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);
    // for selected bead (given by removeMe), calculate distance to all other selected beads
    next = 0;
    int bin;

    while(next < stopAt && (next < workingLimit)){ // stops when addMeSubUnitIndex is reached
        // for each index in subunit
        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[next]) );
        bin = pData->convertToBin((*refVec - *targetVec).length());

        prBins[ bin ] += totalSubUnits;
        //totalAdded += totalSubUnits;
        //for(int ss=0; ss<totalSubUnits; ss++){
        //    prBins[ bin ]++;
        //}
        next++;
    }

    next++; // move past the addMeSubUnitIndex by incrementing
    for(int b=next; b<workingLimit; b++){

        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]));
        bin = pData->convertToBin((*refVec - *targetVec).length());

        prBins[ bin ] += totalSubUnits;
        //totalAdded += totalSubUnits;
        //for(int ss=0; ss<totalSubUnits; ss++){
        //    prBins[ bin ]++;
        //}
    }

    //
    // add intra domain distances (between subunits)
    //
    int subUnitLimit = totalSubUnits-1;
    int outer;
    for(int s=0; s < subUnitLimit; s++) {

        // iterate over each position in outer subunit
        outer=0;
        for(; outer < stopAt; outer++){ // go over all the atoms in the subunit

            refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
            for(int ss=(s+1); ss < totalSubUnits; ss++) {
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
                //  totalAdded++;
            }
        }

        // outer is pointing to new position in outer subunit
        refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
        for(int ss=(s+1); ss < totalSubUnits; ss++) {

            for(int subIndex=0; subIndex < workingLimit; subIndex++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex(beadsInUse[subIndex]));
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
                //  totalAdded++;
            }
        }

        outer++;
        for(; outer < workingLimit; outer++){ // go over all the atoms in the subunit

            refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[outer]) + s);
            for(int ss=(s+1); ss < totalSubUnits; ss++) {
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
                // totalAdded++;
            }
        }
    }


// old method for updating
// total number of beads
    /*
    for(int s=0; s < totalSubUnits; s++){  // go through each subUnit


        currentIndex = s + indexInMasterBeadArray;
        refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(currentIndex);
        // for selected bead (given by removeMe), calculate distance to all other selected beads


        next = 0;
        while(next < stopAt && (next < workingLimit)){ // stops when addMeSubUnitIndex is reached
            // for each index in subunit
            indexOfTarget = pModel->mapSubUnitIndex(beadsInUse[next]);

            for(int ss=0; ss<totalSubUnits; ss++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexOfTarget);
                // calculate distance and convert to PrBin
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
                totalAdded++;
            }
            next++;
        }

        next++; // move past the addMeSubUnitIndex by incrementing
        for(int b=next; b<workingLimit; b++){
            indexOfTarget = pModel->mapSubUnitIndex(beadsInUse[b]);

            for(int ss=0; ss<totalSubUnits; ss++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexOfTarget);
                // calculate distance and convert to PrBin
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
                totalAdded++;
            }
        }


        // need to include distances between sym partners
        for(int ss=(s+1); ss<totalSubUnits; ss++){
            targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
            prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]++;
            totalAdded++;
        }

    }
*/
    //cout << workingLimit << " TOTAL ADDED => " << totalAdded << " === " << (totalSubUnits*(workingLimit-1) + (totalSubUnits-1)*workingLimit)<< endl;
}






/**
 * beadsInUse does not have to be sorted
 *
 * prBins is the model P(r) distribution
 */
inline void Anneal::removeFromPrSym(int removeMeSubUnitIndex, vector<int> & beadsInUse, int workingLimit, vector<int> & prBins, Model *pModel, Data *pData){

    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    float totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * refVec, * targetVec;
    int indexInMasterBeadArray = pModel->mapSubUnitIndex(removeMeSubUnitIndex);
    int indexOfTarget;
    int currentIndex, next;
    //int totalRemoved=0;
    int stopAt=0;

    // find index of bead to remove
    int * ptr = &beadsInUse.front();
    for(int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            stopAt = i;
            break;
        }
    }


/*
    for(int s=0; s<totalSubUnits; s++){ // remove bead and sym mates from P(r)
        currentIndex = s + indexInMasterBeadArray;
        refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(currentIndex);
        // for selected bead (given by removeMe), calculate distance to all other selected beads
        next = 0;
        //while((beadsInUse[next] != removeMeSubUnitIndex) && (next < workingLimit)){
        while(next < stopAt){
            // for each index in subunit
            indexOfTarget = pModel->mapSubUnitIndex(beadsInUse[next]);

            for(int ss=0; ss<totalSubUnits; ss++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexOfTarget);
                // calculate distance and convert to PrBin
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
                //totalRemoved++;
            }
            next++;
        }

        next++;
        for(int b=next; b<workingLimit; b++){
            indexOfTarget = pModel->mapSubUnitIndex(beadsInUse[b]);

            for(int ss=0; ss<totalSubUnits; ss++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexOfTarget);
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
                //totalRemoved++;
            }
        }

        // need to include distances between sym partners
        for(int ss=(s+1); ss<totalSubUnits; ss++){
            targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
            prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
            //totalRemoved++;
        }

    }
*/

    //
    // remove interdomain distances
    //
    refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);
    // for selected bead (given by removeMe), calculate distance to all other selected beads
    next = 0;
    int bin;

    while(next < stopAt && (next < workingLimit)){ // stops when addMeSubUnitIndex is reached
        // for each index in subunit
        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[next]) );
        bin = pData->convertToBin((*refVec - *targetVec).length());

        prBins[ bin ] -= totalSubUnits;
        //for(int ss=0; ss<totalSubUnits; ss++){
        //    prBins[ bin ]++;
        //}
        next++;
    }

    next++; // move past the addMeSubUnitIndex by incrementing
    for(int b=next; b<workingLimit; b++){

        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]));
        bin = pData->convertToBin((*refVec - *targetVec).length());

        prBins[ bin ] -= totalSubUnits;
        //for(int ss=0; ss<totalSubUnits; ss++){
        //    prBins[ bin ]++;
        //}
    }

    //
    // add intra domain distances (between subunits)
    //
    int subUnitLimit = totalSubUnits-1;
    int outer;
    for(int s=0; s < subUnitLimit; s++) {

        // iterate over each position in outer subunit
        outer=0;
        for(; outer < stopAt; outer++){ // go over all the atoms in the subunit

            refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
            for(int ss=(s+1); ss < totalSubUnits; ss++) {
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
            }
        }

        // outer is pointing to new position in outer subunit
        refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
        for(int ss=(s+1); ss < totalSubUnits; ss++) {

            for(int subIndex=0; subIndex < workingLimit; subIndex++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex(beadsInUse[subIndex]));
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
            }
        }

        outer++;
        for(; outer < workingLimit; outer++){ // go over all the atoms in the subunit

            refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[outer]) + s);
            for(int ss=(s+1); ss < totalSubUnits; ss++) {
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
                prBins[ pData->convertToBin((*refVec - *targetVec).length()) ]--;
            }
        }
    }
}

/**
 *
 *
 */
inline void Anneal::removeLatticePositionToModelSym(std::vector<int>::iterator * pBeginIt,
                                                    std::vector<int> & bead_indices,
                                                    std::vector<int> & modelPrBins,
                                                    int * pWorkingLimit,
                                                    const int * pLatticePointToRemove, Model * pModel, Data *pData){

    std::vector<int>::iterator itIndex = find(*pBeginIt, *pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    removeFromPrSym(*pLatticePointToRemove, bead_indices, *pWorkingLimit, modelPrBins, pModel, pData);
    // reduce the workingLimit
    // if wl = 10
    // 0 1 2 3 4 5 6 7 8 9 10
    // remove 4
    // 0 1 2 3 9 5 6 7 8 4 10
    //
    // 0 1 2 3 5 6 7 8 9 4 10
    //
    *pWorkingLimit -= 1;
    // swap selected point to the workingLimit
    std::iter_swap(itIndex, *pBeginIt + *pWorkingLimit);
    // still need to sort, swap changes the order
    // to save some time, sort should only be from swapped point to workingLimit
    std::sort(*pBeginIt, *pBeginIt + *pWorkingLimit);
}


/**
 * for the specified lattice index (position), count number of violations when symmetry mates are created
 * basically trying to miniminize overlaps of points within symmetry with symmetry related
 */
inline float Anneal::symViolationsPotential(int specifiedIndex, std::vector<int> * beads_indices, int workingLimit, Model * pModel){

    float violations=0;
    float totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * refVec, * targetVec;

    int indexInMasterBeadArray = pModel->mapSubUnitIndex(specifiedIndex);
    refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);

    int * ptr = &(beads_indices->front());

    float cutOff = 2*pModel->getBeadRadius();

    for(int ss=1; ss < totalSubUnits; ss++) {

        for(int subIndex=0; subIndex < workingLimit; subIndex++){
            targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex( ptr[subIndex] ));
            if ((*refVec - *targetVec).length() < cutOff){ /*!< distance between lattice points to be counted as a neighbor */
                violations++;
            }
        }
    }

    return violations/(float)workingLimit;
}