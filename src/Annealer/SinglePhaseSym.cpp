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
    for(unsigned long int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    // adjust this part
    maxbin += 2;
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
//    int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();
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
    std::vector<int>::iterator beginIt;

    std::cout << "        CREATING INITIAL RANDOM MODEL " << endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    std::set<int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::uniform_int_distribution<int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    //vertexT * vertices;  // not sure if this needs to be explicitly deleted
    //int indexOfHullpt, count, potentialContacts, index;
    char flags[25];
    sprintf(flags, "qhull s FA");
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
    //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
    populateLayeredDeadlimit(subUnit_indices.begin(), subUnitWorkingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    // calculate starting energy

    std::clock_t start = std::clock();

    // fill binCount for first time
    float currentKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
    float hlambda = 0.001;
    float muConstant = mu*currentKL/((1.0-mu)*current_volume);


    std::cout << " THREAD TIME "  << (std::clock() - start)/(double) CLOCKS_PER_SEC << " seconds " << std::endl;
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

    float current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + muConstant*current_volume;

    std::cout << " INITIAL => ENERGY :  " << current_energy << endl;
    std::cout << " INITIAL =>   D_KL : " <<  currentKL << endl;
    std::cout << " INITIAL => VOLUME : " << current_volume << endl;
    std::cout << "        VIOLATIONS : " << violations << endl;
    std::cout << " INITIAL        WL : " << subUnitWorkingLimit << " DL : " << deadLimit << endl;

    float lowest_energy = current_energy;

    std::vector<int>::iterator itIndex;
    std::vector<int>::iterator beginBinCount = binCount.begin();
    std::vector<int>::iterator endBinCount = binCount.end();

    //TEST
//    std::uniform_int_distribution<> randomTest (0, subUnitWorkingLimit-1);
//    for (int i=0; i<300; i++){
//        int testI = randomTest(gen);
//        removeFromPrSym(subUnit_indices[testI], subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
//        addToPrSym(subUnit_indices[testI], subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
//        testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (pData->calculateKLDivergence(binCount) != testKL){
//            std::cout << i  << " Failed " << pData->calculateKLDivergence(binCount) << " != " << testKL << std::endl;
//            std::cout << subUnitWorkingLimit << " => " << "index " << testI << " => " << subUnit_indices[testI] << std::endl;
//            return "failed";
//        }
//    }

    char addRemoveText[50];
    std::vector<int> averageV(highTempRounds);
    float sum_x_squared=0;
    int volumeCount=0, testIndex;
    float volumeSum=0, workingLimitSum=0;

    float totalViolations=violations/(double)subUnitWorkingLimit;
    beta = 0;
    current_energy += beta * totalViolations;

    bool updateCVX = false;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)

    int high;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) > 0.91){
            // create points from workSet to determine HULL
            float invSubUnitworkingLimit = muConstant;
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

            //std::shuffle(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.begin()+deadLimit, gen); // randomize possible swap space
            //std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            int swapCount = 0, swap1;

            for (int v = 0; v < 1; v++) { // try n times

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
                eulerTour.removeNode(swap1);
                int randind;

                for (int bb = 0; bb < 5; bb++) {
                    //int neighbor = subUnit_indices[subUnitWorkingLimit+bb];

                    randind = randomIndex(gen);
                    int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randind]);
                    while (neighbor == -1 || neighbor == swap1){
                        randind = randomIndex(gen);
                        neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randind]);
                    }

                    // make the swap, sort and update P(r)
                    pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);
                    std::iter_swap(itIndex, pSwap2);
                    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit); // bead_indices needs to be sorted

                    addToPrSym(neighbor, subUnit_indices, subUnitWorkingLimit, workingBinCount, pModel, pData);
                    testKL = pData->calculateKLDivergence(workingBinCount); // 100x faster than calculateKLEnergySymmetry
                    //testKL = calculateKLEnergySymmetry(&subUnit_indices, &workingBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );

                    tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

                    test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + test_volume*invSubUnitworkingLimit;

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
                    eulerTour.removeNode(neighbor);
                }

                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    eulerTour.addNode(swap1, pModel);
                    beads_in_use_tree.insert(swap1);
                }
            }


        } else { // add or remove

            // dead limit is set by convex hull
            alterMe = number_of_beads_to_use(gen); // pick number between lower and upper bounds
            beginIt = subUnit_indices.begin();

            if (subUnitWorkingLimit > alterMe && alterMe > lowerN) { // REMOVE beads from sorted list into useable range < deadLimit
                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                sprintf(addRemoveText, "  REMOVE  ");

                testIndex = subUnit_indices[randomIndex(gen)];

                removeLatticePositionToModelSym(&beginIt, subUnit_indices, binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);

//                testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);

                test_volume = calculateCVXHULLVolume(flags,
                                                     &subUnit_indices,
                                                     subUnitWorkingLimit,
                                                     pModel);

                tempNumberOfComponents = eulerTour.removeNode(testIndex);

                test_energy = testKL + hlambda * (tempNumberOfComponents - 1) * (tempNumberOfComponents - 1) +
                              muConstant * test_volume;

                if (test_energy < current_energy) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = violations;
                    beads_in_use_tree.erase(testIndex);
                    updateCVX=true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = violations;
                    beads_in_use_tree.erase(testIndex);
                    updateCVX=true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    beginIt = subUnit_indices.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &subUnitWorkingLimit, &binCountBackUp,
                                                          &beginBinCount);
                }

            } else { // ADD beads
                // randomly add from deadLimit
                std::sprintf(addRemoveText, "    ADD   ");
//                int randomSpot = rand() % (deadLimit - subUnitWorkingLimit) + subUnitWorkingLimit;
//                testIndex = subUnit_indices[randomSpot];
//                itIndex = subUnit_indices.begin() + randomSpot;

                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                while ( testIndex == -1) { // find a new position
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                }
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.end(), testIndex);

                addLatticPositionToModel(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &itIndex);
                addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);
                //testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);

                test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1)
                              + muConstant*test_volume;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = violations;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = violations;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                    currentNumberOfComponents = eulerTour.removeNode(testIndex);
                }
            }
        } // END OF ADD?REMOVE

        if (updateCVX){
            //populateLayeredDeadlimitUsingSet(&subUnit_indices, &beads_in_use_tree, subUnitWorkingLimit, &deadLimit, pModel);
            std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<int>(0,subUnitWorkingLimit-1); // guaranteed unbiased
            updateCVX=false;
        }

        //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
        //populateLayeredDeadlimitUsingSet(&subUnit_indices, &beads_in_use_tree, subUnitWorkingLimit, &deadLimit, pModel);
        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts

        printf("*******************             %s                 ******************* \n", addRemoveText);
        printf("         MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        printf("            GRAPH => %3i\n", currentNumberOfComponents);
        printf("            UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("           VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, muConstant, muConstant*current_volume);
        printf("       VIOLATIONS => %f  BETA => %.4E  \n", totalViolations, beta);
        printf("        LIMIT : %5i DEADLIMIT : %5i D_KL : %.4E ENRGY: %.4E \n", subUnitWorkingLimit, deadLimit, currentKL, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            //std::string name = "initial_sym_search_" + std::to_string(high);
            //pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, name, high);
            //pModel->writeModelToFile(deadLimit, subUnit_indices, name+"_skin", high);
            //pModel->writeSymModelToFile(subUnitWorkingLimit, subUnit_indices, name+"sym_");
            //counter++;
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;

            muConstant = mu*currentKL/((1.0-mu)*current_volume);

            current_energy = currentKL + muConstant*current_volume;

            lowest_energy = current_energy;
            workingLimitSum += subUnitWorkingLimit;
            sum_x_squared += subUnitWorkingLimit*subUnitWorkingLimit;
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
//            return "failed";
//        }

    } // end of HIGH TEMP EQUILIBRATION

    pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage, volumeStdev;

    volumeAverage = workingLimitSum/(float)volumeCount;
    volumeStdev = volumeAverage*volumeAverage - sum_x_squared/(float)volumeCount;
    // remove points close to hull
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, "initial_CVX_sym_" + filenameprefix, high);
    pModel->setStartingSet(lowest_subUnit_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);

    pModel->setBeadAverageAndStdev(workingLimitSum/(float)volumeCount, volumeStdev);

    std::cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    printf("    LOWEST => %d \n", lowestWorkingLimit);
    std::cout << "*******************                                        *******************" << std::endl;

    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }
}


/**
 *
 */
std::string Anneal::refineSymModel(Model *pModel, Data *pData, std::string nameTo){

    srand(time(0));
    cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY" << endl;
    //pModel->setBeadAverageAndStdev(pModel->getStartingWorkingLimit(), 0.17*pModel->getStartingWorkingLimit());
//    float oldN = pModel->getVolumeAverage();  // need to reset this for additional modeling
//    float oldStdev = pModel->getVolumeStdev();
//    float runningAverage = pModel->getVolumeAverage();

    double lowTempStop = (double)highTempStartForCooling;
    int swap1, totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());

    // make copy of bead_indices
    std::vector<int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    //std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
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

    int deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
    enlargeDeadLimit(bead_indices, &deadLimit, pModel);
    // set deadLimit of the selected set
//    int deadLimit;
//    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    //int deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
    //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
//    std::shuffle(bead_indices.begin() + workingLimit, bead_indices.begin() + deadLimit, gen); // randomize search space

    // convert distances in Search Sphere to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    maxbin=0;
    for(unsigned long int  i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    maxbin += 2;
    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    pData->createWorkingDistribution(maxbin); // create working observed probability distribution that encompasses search sphere

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    std::cout << "    TOTAL EXP N_S BINS: " << totalBins << std::endl;
    std::cout << "    MAX MODEL N_S BINS: " << maxbin << std::endl;
    std::cout << "              BINWIDTH: " << pData->getBinWidth() << std::endl; // lattice should be limited by binwidth

    std::vector<int>::iterator beginBinCount = binCount.begin();
    std::vector<int>::iterator endBinCount = binCount.end();

    int violations;
    float testKL, currentKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
    char flags[] = "qhull FA"; //CVX HULL STUFF
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    std::shuffle(bead_indices.begin()+workingLimit, bead_indices.end(), gen);
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;
    int original;//
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float tempAverageContacts;//, updateSteps = 7.0;;
    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;
    float step_limit = (updateCount < 30000) ? 30000 : (float)updateCount;

    std::vector<float> tempDuringRun(step_limit);
    std::vector<float> divergenceDuringRun(step_limit);
    std::vector<int> workingLimitDuringRun(step_limit);

    float * pTempDuringRun = &tempDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float acceptRate = 0.5, inv500 = 1.0/500.0;
    float inv500slash499 = 499.0/500.0;
    int counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    float this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    //double newSum, runningContactsSum = calculateTotalContactSumPotential(&beads_in_use_tree, pModel);
    double newSum, runningContactsSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);

    double muAt = mu;
    double baseFactor = eta;
    //double etaFactor = std::pow(eta/baseFactor, 1.0/(updateSteps));
    double invWL = 1.0/(double)workingLimit;

    double etaConstant = baseFactor*currentKL/((1.0-baseFactor)*runningContactsSum*invWL);
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum*invWL);

    // want D_KL to settle < 10^-5 to 10^-6
    // int span = (int)(std::floor(log10f(currentKL)) - std::floor(log10f(totalContactEnergy/(float)workingLimit))) + 3;
    float muConstant = mu*currentKL/((1.0-mu)*current_volume*invWL);
    float currentVolE = muConstant*current_volume*invWL;
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    //float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);

    float tempViolations, afterAdding, afterRemoving;
    beta = 0; // violations
    float totalViolations = violations/(double)workingLimit;
    current_energy += beta*totalViolations;

    float lowest_energy = current_energy + totalContactEnergy;

    std::clock_t startTime;
    std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<int> randomDeadIndex (workingLimit, deadLimit-1);

    std::vector<int> selections(workingLimit);
    for(unsigned long int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    int numberOfCoolingTempSteps=0;
    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

        if (distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)

            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            startTime = std::clock();
            if (distribution(gen) < 0.5){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                //
                if (distribution(gen) > acceptRate ){
                    // get lattice point not in use but touching a position in use, pick random bead to add to
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while(original == -1){
                        original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }
                } else {
                    original = bead_indices[randomDeadIndex(gen)];
                }

                // rather than try to find a bead to add, see if the one random one is acceptable
                itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);
                // check if available neighbor can be added

               // newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);
               // afterAdding = etaConstant*newSum/(double)(workingLimit+1);
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);

                beads_in_use_tree.insert(original);
                newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit+1, pModel);
                afterAdding = etaConstant*newSum/(double)(workingLimit+1);

                addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);
//                testKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//                tempViolations = violations/(double)(workingLimit+1);

                tempNumberOfComponents = eulerTour.addNode(original, pModel);

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);// + beta*tempViolations;
                tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                   // beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
//                    totalViolations = tempViolations;
                    isUpdated = true;
                    //eulerTour.addNode(original, pModel);
                    currentNumberOfComponents = tempNumberOfComponents;
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                  //  beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
//                    totalViolations = tempViolations;
                    isUpdated = true;
                    //eulerTour.addNode(original, pModel);
                    currentNumberOfComponents = tempNumberOfComponents;
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(original);

                    eulerTour.removeNode(original);
                    sprintf(addRemoveText, "     ADD => %i", 0);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                // test for deletion
                original = bead_indices[randomIndex(gen)];

                beginIt = bead_indices.begin();
                removeLatticePositionToModelSym(&beginIt, bead_indices, binCount, &workingLimit, &original, pModel, pData);

                beads_in_use_tree.erase(original);
                newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit-1, pModel);

                afterRemoving = etaConstant*newSum/(double)(workingLimit-1);
                testKL = pData->calculateKLDivergence(binCount);

//                testKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//                tempViolations = violations/(double)(workingLimit-1);

                tempNumberOfComponents = eulerTour.removeNode(original);
                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);// + beta*tempViolations;

                    tempTotalContactEnergy = (afterRemoving - totalContactEnergy);

                    if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
//                        beads_in_use_tree.erase(original);
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
//                        totalViolations = tempViolations;
                        isUpdated = true;
                    } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
//                        beads_in_use_tree.erase(original);
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
//                        totalViolations = tempViolations;
                        isUpdated = true;
                    } else { // undo changes and move to next bead (rejecting)
                        beads_in_use_tree.insert(original);
                        beginIt = bead_indices.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.addNode(original, pModel);
                    }
            }

        } else { // positional refinement

// only search within deadLimit, no need to recalculate at end
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            //float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
            // find bead to swap in active set
            // only refine a single position with 2 or less contacts
            int originalSwap2Value;
            startTime = std::clock();

            // randomly select an index to move
            // select only node I can move?
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            int position = selections[0];

            swap1 = bead_indices[ position ];
            for(int i=1;i<workingLimit; i++){
                if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                    break;
                } else { // if last node, will break out
                    eulerTour.addNode(swap1, pModel);
                    position = selections[i];
                    swap1 = bead_indices[position]; // what happens if I fail all the way?
                }
            }

//            double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);
            // isSwapped = false;
            // swap
            itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
            removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
            std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

            beads_in_use_tree.erase(swap1);
            double newPotential, invMuWorkingLimit = muConstant*invWL;

            // find better position
            for(int i=0;i<workingLimit; i++){
                int newPosition = bead_indices[selections[i]];
                if (newPosition != swap1){
                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, newPosition);
                    if (originalSwap2Value != -1 && originalSwap2Value != swap1){
                        break;
                    }
                }
            }

            // get available neighbor position
            pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);

//                newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
//                newPotential = etaConstant*newSum*invWorkingLimit;

            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPrSym(originalSwap2Value, bead_indices, workingLimit, workingBinCount, pModel, pData);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->calculateKLDivergence(workingBinCount);

            beads_in_use_tree.insert(originalSwap2Value);
            newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);
            newPotential = etaConstant*newSum*invWL;

//            double testKL1 = calculateKLEnergySymmetry(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//            if (testKL != testKL1){
//                std::cout << " Not equal " << testKL << " " << testKL1 << std::endl;
//                exit(0);
//            }
//            tempViolations = violations/(double)(workingLimit);

            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);// + beta*tempViolations;
            tempTotalContactEnergy = newPotential-totalContactEnergy;

            if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
//                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
//                    eulerTour.addNode(originalSwap2Value, pModel);
                std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                currentKL = testKL;
                current_energy = this_energy;
                currentNumberOfComponents = tempNumberOfComponents;
                totalContactEnergy = newPotential;
                runningContactsSum = newSum;
                isUpdated = true;
//                    totalViolations = tempViolations;
                sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {
//                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
//                    eulerTour.addNode(originalSwap2Value, pModel);
                std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                currentKL = testKL;
                current_energy = this_energy;
                currentNumberOfComponents = tempNumberOfComponents;
                totalContactEnergy = newPotential;
                runningContactsSum = newSum;
                isUpdated = true;
//                    totalViolations = tempViolations;
                sprintf(addRemoveText, "  SA SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                beads_in_use_tree.erase(originalSwap2Value);

                beads_in_use_tree.insert(swap1); // add back swap one to tree
                sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                eulerTour.addNode(swap1, pModel);
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;
//        float testKL1 = calculateKLEnergySymmetry(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

        printf("       TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("     ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        printf("       TIME : %.5f (SECONDS)  %s\n", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC), addRemoveText);
        printf("      LIMIT : %5i  VOL => %-4.3E AVG => %.3f \n", workingLimit, current_volume, binCount[0]/(float)workingLimit);
//        printf(" VIOLATIONS => %f  BETA => %.4E  \n", totalViolations, beta);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f BASE => %.2f\n", totalContactEnergy, etaConstant, runningContactsSum, muAt);
        printf("       D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, lowestKL, (current_energy+totalContactEnergy));

        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            invWL = 1.0/(double)workingLimit;
            updated++;

            if (updated % coupons == 0 && numberOfCoolingTempSteps > 1){
                current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
                enlargeDeadLimit(bead_indices, &deadLimit, pModel);
                //pModel->writeModelToFile(deadLimit, bead_indices, "dead", updated);
                //pModel->writeModelToFile(workingLimit, bead_indices, "cooling", updated);
            }

            randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
            randomDeadIndex = std::uniform_int_distribution<int>(workingLimit, deadLimit-1);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            selections.resize(workingLimit);
            int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned long int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);

        //update for running average
        if ((current_energy + totalContactEnergy) < lowest_energy){ // if counter too small, add/remove may not sample sufficiently

            if ( totalContactEnergy/current_energy < eta) {
                etaConstant = current_energy*eta/((1.0 - eta)*runningContactsSum*invWL);
                totalContactEnergy = etaConstant*runningContactsSum*invWL;
            }

            if (currentKL < lowestKL){
                lowestKL = currentKL;
            }

            lowest_energy = current_energy + totalContactEnergy;
        }
        counter++;
    } // end of steps


    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (int i=0; i<workingLimit; i++){
        tempAverageContacts += numberOfContacts(bead_indices[i], &bead_indices, workingLimit, pModel, pDistance);
    }

    tempAverageContacts =tempAverageContacts*invWL;
    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }
    // perform constant temp
    this->asaAcceptanceRate = 0.001;
    complementASAAcceptanceRate = 1.0 - this->asaAcceptanceRate;
    intASAAcceptanceRate = (int)(1000*this->asaAcceptanceRate);
    intComplementASAAcceptanceRate = (int)(1000*complementASAAcceptanceRate);

    lowTempStop *= 0.01;
    inv_kb_temp = 1.0/lowTempStop;
    int originalSwap2Value, finalCoupons = 3*(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    updated = 0;
    double newPotential, invWorkingLimit = etaConstant*invWL;

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < 0*finalCoupons; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

        std::sprintf(addRemoveText, " ");
        std::cout << "______________________________________________________________________________" << std::endl;
        std::cout << "*******************        CONSTANT TEMP POSITIONAL        *******************" << std::endl;
        std::cout << "*******************                                        *******************" << std::endl;

        // find bead to swap in active set
        // only refine a single position with 2 or less contacts

        startTime = std::clock();

        // randomly select an index to move
        // select only node I can move?
        bool tourtest = true;
        int position = randomIndex(gen);
        swap1 = bead_indices[ position ];

        while(tourtest){
            if (eulerTour.removeNode(swap1) == 1){
                tourtest = false;
            } else {
                eulerTour.addNode(swap1, pModel);
                position = randomIndex(gen);
                swap1 = bead_indices[ position ];
            }
        }

        // old above
        double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);
        // swap
        itIndex = bead_indices.begin() + position;
        // remove selected index from P(r)
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
        removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
        // violations in current position
        float currentViolations = symViolationsPotential(swap1, &bead_indices, workingLimit, pModel);

        beads_in_use_tree.erase(swap1);
        // find better position
        for(int i=0;i<workingLimit; i++){
            int newPosition = bead_indices[selections[i]];
            if (newPosition != swap1){
                originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, newPosition);
                if (originalSwap2Value != -1 && originalSwap2Value != swap1){
                    break;
                }
            }
        }

        // get available neighbor position
        pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);

        newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
        newPotential = newSum*invWorkingLimit;

        std::iter_swap(itIndex, pSwap2);
        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

        addToPrSym(originalSwap2Value, bead_indices, workingLimit, workingBinCount, pModel, pData);

        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
        tempViolations = totalViolations + symViolationsPotential(originalSwap2Value, &bead_indices, workingLimit, pModel) - currentViolations;
        testKL = pData->calculateKLDivergence(workingBinCount);
//                tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                tempVolE = invMuWorkingLimit * tempVolume;
//                this_energy = testKL + tempVolE + beta*tempViolations;
        this_energy = testKL + beta*tempViolations;
        tempTotalContactEnergy = newPotential-totalContactEnergy;

        if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            eulerTour.addNode(originalSwap2Value, pModel);

            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
            currentKL = testKL;
            current_energy = this_energy;
//                    current_volume = tempVolume;
//                    currentVolE = tempVolE;
            totalContactEnergy = newPotential;
            runningContactsSum = newSum;
            isUpdated = true;
            totalViolations = tempViolations;
            std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
        } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {
            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            eulerTour.addNode(originalSwap2Value, pModel);

            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
            currentKL = testKL;
            current_energy = this_energy;
//                    current_volume = tempVolume;
//                    currentVolE = tempVolE;
            totalContactEnergy = newPotential;
            runningContactsSum = newSum;
            isUpdated = true;
            totalViolations = tempViolations;
            std::sprintf(addRemoveText, "  SA SWAPPED => %i to %i", swap1, originalSwap2Value);
        } else {
            std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
            std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
            beads_in_use_tree.insert(swap1); // add back swap one to tree
            std::sprintf(addRemoveText, "      FAILED => %i", swap1);
//            eulerTour.removeNode(originalSwap2Value);
            eulerTour.addNode(swap1, pModel);
        }

        printf("       TEMP : %-.4E MAXSTEPS => %i (%4i) \n", lowTempStop, finalCoupons, numberOfCoolingTempSteps);
        printf("      LIMIT : %i UPDATES : %6i  FAILURES => %i  \n", workingLimit, updated, failures);
        printf("       TIME : %.5f (SECONDS)  %s\n", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC), addRemoveText);
        printf("       VOLE => %-5.4E ( %.0f ) AVG => %.3f \n", currentVolE, current_volume, binCount[0]*invWL);
        printf(" VIOLATIONS => %f  BETA => %.4E  \n", totalViolations, beta);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f BASE => %.2f\n", totalContactEnergy, etaConstant, runningContactsSum, muAt);
        printf("       D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, lowestKL, (current_energy+totalContactEnergy));

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;
            if (currentKL < lowestKL ){
                lowestKL = currentKL;
            }

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            failures=0;
            updated++;
        } else {
            acceptRate = inv500*(499*acceptRate);
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
    }

    // At end of each temp, update a probability model for volume?  Use this to select
    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (int i=0; i<workingLimit; i++){
        tempAverageContacts += numberOfContacts(bead_indices[i], &bead_indices, workingLimit, pModel, pDistance);
    }
    tempAverageContacts =tempAverageContacts/(float)workingLimit;

    // perform positional refinement until delta E stabilizes?
    // final round to move any points close to body
    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    // pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    // pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");

//    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_final_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);
//    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);

    //pModel->writeModelToFile(workingLimit, bead_indices, "refined_");
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");

    //  pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_" + filenameprefix);
    //  pModel->writeSymModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_sym_");

    return nameOfModel;
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

    float limit = (pModel->getBeadRadius())*0.51;
    violation=0;
    //float divisor = totalBeadsInSphere*(totalBeadsInSphere-1)*0.5;

    int subUnits = (int)pModel->getNumberOfSubUnits();
    int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;
    // create first subunit from selected indices and populate coordinates
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


    violation /= 2;

    // calculate KLDivergence
    //float squared = subUnits*subUnits*indicesWorkingLimit*indicesWorkingLimit;
    //float ratio = (float)violation/divisor;
    //cout << " RATIO " << ratio << " " << violation << endl;
//    double value = pData->calculateKLDivergence(*binCount);
//    if (value < 0){
//        std::cout << " neg" << std::endl;
//        int ct = 0;
//        for(std::vector<int>::iterator it = binCount->begin(); it != binCount->end(); ++it){
//            std::cout << ct << " => " << *it << std::endl;
//            ct++;
//        }
//        exit(0);
//    }
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
    int next;
    int stopAt=0;

    // find index of bead to remove
    int * ptr = &beadsInUse.front();
    for(int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            stopAt = i;
            break;
        }
    }
    //
    // remove interdomain distances
    //
    refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);
    // for selected bead (given by removeMe), calculate distance to all other selected beads
    next = 0;

    while(next < stopAt && (next < workingLimit)){ // stops when addMeSubUnitIndex is reached
        // for each index in subunit
        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[next]) );
        prBins[ pData->convertToBin((*refVec - *targetVec).length()) ] -= totalSubUnits;

        next++;
    }

    next++; // move past the addMeSubUnitIndex by incrementing
    for(int b=next; b<workingLimit; b++){

        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]));
        prBins[ pData->convertToBin((*refVec - *targetVec).length()) ] -= totalSubUnits;
    }

    //
    // remove intra domain distances (between subunits)
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

    std::vector<int>::iterator itIndex = std::find(*pBeginIt, *pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
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

    float cutOff = 1.001*pModel->getBeadRadius();

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

//need to determine, for each bead, number of contacts it makes

/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSumPotentialSym(std::set<int> *beads_in_use, std::vector<int> * beads_indices, int workingLimit, Model *pModel){

    //calculate contacts per bead for selectedIndex
    double sum=0;
    float totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * refVec, * targetVec;

    int * ptr = &(beads_indices->front());
    float cutOff = 2.001*pModel->getBeadRadius();


    std::set<int>::iterator it;
    for (it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {

        int subTotal = numberOfContactsFromSet(beads_in_use, pModel, *it);

        // now include contacts from symmetry
        int indexInMasterBeadArray = pModel->mapSubUnitIndex(*it);
        refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray); // calculate distances from ref to other

        for(int ss=1; ss < totalSubUnits; ss++) {

            for(int subIndex=0; subIndex < workingLimit; subIndex++){
                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex( ptr[subIndex] ));
                if ((*refVec - *targetVec).length() <= cutOff){ /*!< distance between lattice points to be counted as a neighbor */
                    subTotal++;
                }
            }
        }

        sum += totalContactsPotential(subTotal);
    }

    return sum;
}


/**
 * Given the selectedIndex, calculate the new sum of contacts for the selected Set defined in beads_in_use
 */
inline double Anneal::recalculateContactsPotentialSumAddSym(std::set<int> *beads_in_use,
                                                            Model *pModel,
                                                            int const selectedIndex,
                                                            std::vector<int> * beads_indices,
                                                            int workingLimit,
                                                            double const currentSum) {

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    double sum = currentSum;
    int currentContactCountNeighbor, newContactsFromSelectedIndex=0;

    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if ((neighbor > -1) && beads_in_use->find(neighbor) != endOfSet){ // -1 will be endOfSet and also beads not in use
            // get contacts for the neighbor
            // since the selectedIndex isn't added yet, we add one to the count
            currentContactCountNeighbor = numberOfContactsFromSet(beads_in_use, pModel, neighbor);

            // adjust sum by removing old contribution and adding new one
            sum += totalContactsPotential(currentContactCountNeighbor+1) - totalContactsPotential(currentContactCountNeighbor);
            //sum += 1.0;
            newContactsFromSelectedIndex += 1.0;
        } else if (neighbor == -1) {
            break;
        }
    }


    // check symmetry partners
    float totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * targetVec;

    vector3 * refVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(selectedIndex));

    int * ptr = &(beads_indices->front());

    float cutOff = 2.001*pModel->getBeadRadius();

    for(int ss=1; ss < totalSubUnits; ss++) {

        for(int subIndex=0; subIndex < workingLimit; subIndex++){
            targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex( ptr[subIndex] ));
            if ((*refVec - *targetVec).length() < cutOff){ /*!< distance between lattice points to be counted as a neighbor */
                newContactsFromSelectedIndex += 1.0;
            }
        }
    }

    return (sum + totalContactsPotential(newContactsFromSelectedIndex));
    //return (sum + newContactsFromSelectedIndex);
}


