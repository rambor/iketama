//
// Created by Robert Rambo on 13/01/2017.
//
#include "../Anneal.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"

using namespace std;

void Anneal::createInitialModelCVXHull(Model *pModel, Data *pData, std::string name) {

    srand(time(0));

    contactCutOff = interconnectivityCutOff;
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

    int deadLimit; // as bead indices are discarded, set upper limit of vector
    std::vector<int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);
    std::vector<int>::iterator beginIt = bead_indices.begin(), endIt = bead_indices.end();
    //std::clock_t start;
    // c-style is slightly faster for large vector sizes
    // start = std::clock();
    const int num = totalBeadsInSphere;
    int * ptr = (num != 0) ? &bead_indices.front() : NULL;
    for(int i = 0; i < num; i++) {
        ptr[i] = i;
    }
    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;

    std::random_device rd;
    std::mt19937 gen(rd());

    float invBeadVolume = 1.0/pModel->getBeadVolume();

    int lowerN = round(lowerV*invBeadVolume);
    int upperN = round(upperV*invBeadVolume);

    // connectivity potential
    // volume potential
    // pick random number between lowerN and upperN
    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN);
    int workingLimit = number_of_beads_to_use(gen);
    pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "universe");
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // setup parameters for hull
    char flags[25];
    sprintf(flags, "qhull s FA");
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];

    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int alterMe, tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();
    // isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();

    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    //lambda = 0.001;
    mu = 0.000001; // scale this to volume of search space

    float testKL, test_energy, current_energy = currentKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*current_volume/(float)workingLimit;

    float lowest_energy = current_energy;

    int testIndex;
    float inv_kb_temp = 1.0/(float)highT;
    inv_kb_temp = 1.0/0.0001;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    cout << " POPULATING DEAD LIMIT" << endl;
    //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
    populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere);
    cout << " DEAD LIMIT SET => " << deadLimit << endl;

    int lowestWorkingLimit, lowestDeadLimit;
    lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    lowestDeadLimit = deadLimit;
    double sumIt = 0.0, counter=0.0;
    bool updateCVX =false;

    std::clock_t start;
    int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();

    cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< endl;
    for (int high=0; high < highTempRounds; high++) { // iterations during the high temp search

        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy

        if (distribution(gen) < 0.67){

            alterMe = number_of_beads_to_use(gen);  // uniform distribution

            beginIt = bead_indices.begin();
            endIt = bead_indices.end();

            if (alterMe > workingLimit){
                cout << "*******************                  ADD                   ******************* "  << endl;
                // randomly add from deadLimit
                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;

                testIndex = bead_indices[randomSpot];
                itIndex = bead_indices.begin() + randomSpot;

                //make the swap at the border (workingLimit)
                addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);

                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->calculateKLDivergence(binCount);
                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                //tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&beginIt, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    //eulerTour.removeNode(testIndex);
                }

            } else if (alterMe < workingLimit) {

                cout << "*******************                 REMOVE                 ******************* "  << endl;

                if (workingLimit > lowerN){
                    // go through entire list
                    int randomSpot = rand() % workingLimit;

                    // grab from randomized active_indices list
                    testIndex = bead_indices[randomSpot];
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = pData->calculateKLDivergence(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);

                    //tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy ) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                    } else if (exp(-(test_energy - current_energy)*inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                    } else { // undo changes and move to next bead (rejecting)
                        // std::copy(backUpState.begin(), backUpState.begin()+workingLimit+1, beginIt);
                        //eulerTour.addNode(testIndex, pModel);
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

        } else {

            cout << "*******************               POSITIONAL               *******************" << endl;

            int swap1, bb;
            std::vector<int>::iterator pSwap2;
            // recalculate
            // calculate convex hull and get hull points
            // can be threaded
            for (int i = 0; i < workingLimit; i++) {
                beadToPoint(&points[i*3], pModel->getBead(bead_indices[i]));
                active_indices[i] = bead_indices[i];
            }

            // needs to be optimized
            qh_new_qhull(3, workingLimit, points, 0, flags, NULL, NULL);
            vertexT * vertices = qh vertex_list;
            int totalV = qh num_vertices;

            // make copy of points from vertices and then free qh_hull
            bool isSwapped;

            // only move CVX hull points
            std::vector<int> indices_to_check(totalV); // large vector ~1000's
            for (int v = 0; v < totalV; v++) { //
                indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                vertices = vertices->next;
            }
            qh_freeqhull(true);

            std::shuffle(bead_indices.begin()+workingLimit, bead_indices.begin()+deadLimit, gen);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            int swapCount = 0;

            for (int v = 0; v < 3; v++) { // try 3 times

                swap1 = indices_to_check[v];
                isSwapped = false;
                // find bead to swap in active set
                itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                // remove selected index from P(r)
                std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // make copy of altered P(r)
                // find better position
                //eulerTour.removeNode(swap1);
                beads_in_use_tree.erase(swap1);
                //for (bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
                for (bb = 0; bb < sizeOfNeighborhood; bb++) { // go through and find a better position within unused beads
                    // if neighbor is not in use, then find it in bead_indices and assign to pSwap2
                    //int neighbor = bead_indices[rand()%(deadLimit-workingLimit) + workingLimit ];
                    int neighbor = bead_indices[workingLimit+bb];

                    // make the swap, sort and update P(r)
                    pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.begin()+deadLimit, neighbor);

                    std::iter_swap(itIndex, pSwap2);
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                    addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                    // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                    testKL = pData->calculateKLDivergence(workingBinCount);
                    //tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);
                    isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy) {
                        beads_in_use_tree.insert(neighbor);
                        //move swapped position to the deadLimit
                        std::iter_swap(pSwap2, bead_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
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
                        std::iter_swap(pSwap2, bead_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin());
                        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
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
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    itIndex = std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, swap1); // find swapped value in the sorted array
                    //eulerTour.removeNode(neighbor);

                } // end of neighborhood for-loop

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                // if (bb >= deadLimit && !isSwapped){ // itIndex already reverses on exit from loop
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    //std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    //eulerTour.addNode(swap1, pModel);
                    beads_in_use_tree.insert(swap1);
                }
            }
            //cout << " SWAPPED POINTS " << swapCount << " <=> " << totalV <<  endl;
        }


        if (updateCVX){
            //deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel,  totalBeadsInSphere);
            //start = std::clock();
            //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
            populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
            //populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere);
            //cout << "DEADLIMIT TIME : " << ((std::clock() - start)/(double) CLOCKS_PER_SEC) << endl;
            updateCVX=false;
//            string name = "search_" + std::to_string(high);
//            pModel->writeModelToFile(workingLimit, bead_indices, name);
//            name = "search_skin" + std::to_string(high);
//            pModel->writeModelToFile(deadLimit, bead_indices, name);
        }

        cout << "*******************                                        *******************" << endl;
        printf("   MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        printf("      GRAPH => %3i\n", currentNumberOfComponents);
        printf("      UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("     VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, mu, mu*current_volume);
        printf("LIMIT: %5i DEADLIMIT: %5i D_KL: %.4E ENRGY: %.4E \n", workingLimit, deadLimit, currentKL, current_energy);
        cout << "*******************                                        *******************" << endl;

        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
        //float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
        //if (currentKL != testKL1 || checkForRepeats(bead_indices)){
        //    cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
        //    return;
        //}

//        string name = "search_" + std::to_string(high);
//        pModel->writeModelToFile(workingLimit, bead_indices, name);
//        name = "search_skin" + std::to_string(high);
//        pModel->writeModelToFile(deadLimit, bead_indices, name);

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowest_energy = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }

    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name);
    pModel->writeModelToFile(lowestDeadLimit, lowest_bead_indices, "initial_search_layer");
    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    //pModel->setStartingWorkingLimit(lowestDeadLimit);
    pModel->setStartingDeadLimit(lowestDeadLimit);

    double av = 1.0/counter*sumIt;
    pModel->setBeadAverageAndStdev(av, 0.2*av);
    cout << "*******************                                        *******************" << endl;
    cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << endl;
    printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    cout << "*******************                                        *******************" << endl;
}


string Anneal::refineHomogenousBodyASACVX(Model *pModel, Data *pData, int iteration){
    cout << "########################<<<<<>>>>>>############################## " << endl;
    cout << "#                                                               # " << endl;
    cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << endl;
    cout << "#                                                               # " << endl;
    cout << "########################<<<<<>>>>>>############################## " << endl;

    pModel->setBeadAverageAndStdev(1.37*pModel->getStartingWorkingLimit(), 0.17*1.37*pModel->getStartingWorkingLimit());
    float oldN = pModel->getVolumeAverage();     // need to reset this for modeling via rounds
    float oldStdev = pModel->getVolumeStdev();

    populatePotential(pModel->getSizeOfNeighborhood());

    this->lowTempStop = highTempStartForCooling;
    int priorWorkingLimit;
    int swap1;

    int alterMe;
    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));

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
    //pModel->centerLatticeModel(workingLimit, bead_indices);
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
    //int dim = 3;
    char flags[25];
    sprintf(flags, "qhull s FA");

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;
    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float inv_kb_temp, tempAverageContacts;

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
    inv_kb_temp = 1.0/this->lowTempStop;
    int deltaECount =0;
    float energy_i, this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    double runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, workingLimit, pModel);

    //float tempTotalContactEnergy, totalContactEnergy = eta*(totalContactsPotential(runningContactsSum / (float)workingLimit));
    //eta = pow(10, ceil(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) + 1);
    //eta = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) - 2.5 );
    eta = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) + 4 );
    double tempTotalContactEnergy, totalContactEnergy = eta*(runningContactsSum / (double)workingLimit);
    //double etaFactor = 1.0/std::pow(100000, 1.0/(step_limit/(float)deadUpdate));
    double etaFactor = 1.0/std::pow(100000, 1.0/(step_limit/(float)deadUpdate));

    float low_temp_limit = step_limit*0.91;
    // want D_KL to settle < 10^-5 to 10^-6
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    std::normal_distribution<float> volumeGen(pModel->getVolumeAverage(), pModel->getVolumeStdev());

    std::clock_t startTime;
    int attempts=0, failures=0;
    double newSum;
    double runtime;
    int contactsLimit = (int)ceil(contactsPerBead);

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
        beginIt = bead_indices.begin();
        endIt = bead_indices.end();

        if ( distribution(gen) < percentAddRemove ) { //add or remove bead within working Set (exclude deadzone)
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               ADD?REMOVE               *******************" << endl;
            alterMe = (int) volumeGen(gen);
            // build a list of indices within defined region of convex hull
            if (alterMe > workingLimit){ // ADD BEAD?
                cout << "*******************                  ADD                   *******************" << endl;
                sprintf(addRemoveText, "");
                startTime = std::clock();
                double afterAdding;
                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                original = bead_indices[randomSpot];
                itIndex = bead_indices.begin() + randomSpot;
                // select first element from randomized active_set
                // check if available neighbor can be added
                newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);
                afterAdding = eta*newSum/(double)(workingLimit+1);
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);
                addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = pData->calculateKLDivergence(binCount);

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                //tempNumberOfComponents = eulerTour.addNode(original, pModel);
                tempNumberOfComponents = 1;
                //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

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
                    //currentNumberOfComponents = eulerTour.removeNode(original);
                }

                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

            } else { // REMOVE BEADS?
                cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
                priorWorkingLimit = 1;
                sprintf(addRemoveText, "     REMOVE => %i", priorWorkingLimit);

                double afterRemoving;
                startTime = std::clock();
                int randomSpot = rand() % workingLimit;
                original = bead_indices[randomSpot];
                tempNumberOfComponents = eulerTour.removeNode(original);

                if (tempNumberOfComponents == 1){
                    // grab from randomized active_indices list

                    newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                    afterRemoving = eta*newSum/(double)(workingLimit-1);
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                    // still need to sort, swap changes the order
                    testKL = pData->calculateKLDivergence(binCount);

                    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                    this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);

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

                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
            }
        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end

            bool isSwapped;
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::shuffle(bead_indices.begin() + workingLimit, bead_indices.begin() + deadLimit, gen); // randomize search space
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin()); // backup current state
            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               POSITIONAL               *******************" << endl;
            cout << "*******************                                        *******************" << endl;
            int originalSwap2Value;
            int wldl = (workingLimit + (int)((deadLimit-workingLimit)*0.31));

            if (numberOfCoolingTempSteps > low_temp_limit){
                wldl = deadLimit;
            }

            startTime = std::clock();
            // randomly select an index to move
            swap1 = bead_indices[ rand() % workingLimit];

            if (eulerTour.removeNode(swap1) == 1){
                // prefer to move index if it has too few contacts
                // remove contribution of swap1
                double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);

                isSwapped = false;
                // swap
                itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1); // locate lattice within workingLimit
                // remove selected index from P(r)
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
                beads_in_use_tree.erase(swap1); // remove from in_use tree

                double newPotential, invWorkingLimit = 1.0/(double)(workingLimit);
                // possible multi-thread the loop
                // break on first success

                for (int bb = workingLimit; bb < wldl; bb++) { // go through and find a better position within unused beads
                    // make the swap, sort and update P(r)
                    pSwap2 = (bead_indices.begin() + bb);
                    originalSwap2Value = *pSwap2;

                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, originalSwap2Value) > 0){
                        newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
                        newPotential = eta*newSum*invWorkingLimit;

                        beads_in_use_tree.insert(*pSwap2); // add new lattice

                        std::iter_swap(itIndex, pSwap2);
                        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                        addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere,
                                workingBinCount);

                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        testKL = pData->calculateKLDivergence(workingBinCount);

//                isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,tempNumberOfComponents);  // 10x longer than contactEnergy calculation
                        tempNumberOfComponents=1;
                        this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents);
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
                            sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                            break;
                        } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) >
                                   distribution(gen)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
                            sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                        beads_in_use_tree.erase(originalSwap2Value);
                    }
                }

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                    sprintf(addRemoveText, "");
                    eulerTour.addNode(swap1, pModel);
                }
            } else {
                eulerTour.addNode(swap1, pModel);
            }
            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

//        isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,tempNumberOfComponents);  // 10x longer than contactEnergy calculation
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices) || tempNumberOfComponents != 1 ){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            cout << " COMP " << tempNumberOfComponents << endl;
//            return "stopped";
//        }
//        if (tempNumberOfComponents > 1){
//            cout << " STOPPED TOO MANY COMPONENTS " << currentNumberOfComponents << endl;
//            string name = "failed_" + std::to_string(numberOfCoolingTempSteps);
//            pModel->writeModelToFile(workingLimit, bead_indices, name);
//            return "Stopped";
//        }

        printf("       TEMP => %-.8f \n     ACCEPT => %.4f  FAILURES => %i\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop, acceptRate, failures, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
        printf("  UPDATECNT => %7i TIME => %.5f (SECONDS) \n", deadUpdate, runtime);
        printf("      GRAPH => %3i \n", currentNumberOfComponents);
        printf(" LATTCE AVG => %.0f        STDEV => %0.3f  %s\n", pModel->getVolumeAverage(), pModel->getVolumeStdev(), addRemoveText);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f \n", totalContactEnergy, eta, runningContactsSum);
        printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, (current_energy+totalContactEnergy));

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
            // I am only moving/removing one bead at a time, should this be updated every 10 steps?
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
            if (sum_x_square_inv_divideBy > (1+average_x2)){
                stdev = sqrt(sum_x_square_inv_divideBy - average_x*average_x);
            } else {
                stdev = 2;
            }

            sum_x_squared = 0.0;
            sum_x = 0.0;
            divideBy= 0.0;
            pModel->setBeadAverageAndStdev(0.5*(average_x + pModel->getVolumeAverage()), stdev);
            //recalculate eta, increase last term to decrease weight of eta
            // if too much eta, shape is determined by contact potential
            // need right balance between KL and contactPotential
            //eta = pow(10, ceil(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) - 0.911);
            //this->modPotential(1.2);
            //runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, workingLimit, pModel);
            eta *= etaFactor;
            totalContactEnergy = eta*(runningContactsSum / (double)workingLimit);

            string name = "model_" + std::to_string(numberOfCoolingTempSteps);
            pModel->writeModelToFile(workingLimit, bead_indices, name);
            //name = "newSkin_" + std::to_string(numberOfCoolingTempSteps);
            //pModel->writeModelToFile(deadLimit, bead_indices, name);
            std::normal_distribution<float> volumeGen(pModel->getVolumeAverage(), pModel->getVolumeStdev());
        }
        counter++;
    } // end of steps

    // At end of each temp, update a probability model for volume?  Use this to select
    // string tempName = "rough";
    // pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, tempName, this, pData);

    pModel->writeModelToFile(workingLimit, bead_indices, "before_final");
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
    deltaECount =0;
    energy_i = current_energy;

    int totalLessThanAverage=0, numberContacts;
    for (int i=0; i<workingLimit; i++){
        numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        if (numberContacts < average_number_of_contacts ){
            totalLessThanAverage +=1;
        }
    }

    int finalRounds = (totalLessThanAverage*std::log((double)totalLessThanAverage) + 0.5772156649*totalLessThanAverage + 0.5);
    cout << " TOTAL LATTICE POINTS LESS THAN AVERAGE " << totalLessThanAverage << endl;

    int swappedCounter =0;
    int round=0;
    double beforeFinalKL = currentKL, relative_diff;
    double beforeFinalEnergy = current_energy;
    for (; round < 0*finalRounds; round++){

        std::shuffle(beginIt + workingLimit, beginIt + deadLimit, gen);
        std::copy(beginIt, beginIt + workingLimit, active_indices.begin());
        std::copy(beginIt, endIt, backUpState.begin());
        std::shuffle(active_indices.begin(), active_indices.begin() + workingLimit, gen); // randomizes order of beads to select
        cout << "______________________________________________________________________________" << endl;
        cout << "*******************                 FINAL                  *******************" << endl;
        cout << "*******************               POSITIONAL               *******************" << endl;

        float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
        int originalSwap2Value;
        bool isSwapped;

        int numberContacts, contactsToConsider;
        // only perform positional refinement on beads that do not satisfy the contact limit
        beginIt = bead_indices.begin();

//        for (int i=0; i<workingLimit; i++){
//            // for each lattice point, check if better local position
//            swap1 = active_indices[i];
//            numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);
//
//            if (numberContacts < contactsLimit){
//                cout << "------------------------------------------------------------------------------" << endl;
//                printf("      CHECKING INDEX => %i               CONTACTS: %i (%.1f)\n", swap1, numberContacts, average_number_of_contacts);
//                //swap1 defined during setting of placesToCheck
//                isSwapped = false;
//                itIndex = std::find(beginIt, beginIt + workingLimit, swap1); // locate lattice within workingLimit
//                // remove selected index from P(r)
//                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
//                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
//                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
//                // current local potential
//                localPotentialAtOldPosition = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1) +
//                                                   calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, swap1));
//                beads_in_use_tree.erase(swap1);
//                localPotentialAtOldPositionWithOutBead = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1));
//                // find better position within neighborhood
//                std::vector<int>::iterator it = pModel->getPointerToNeighborhood(swap1);
//
//                for (int n=0; n< totalNeighbors; n++){ // run into problem of check
//
//                    int neighbor = *(it+n);
//                    if ((neighbor > -1) && beads_in_use_tree.find(neighbor) == beads_in_use_tree.end()){
//                        // if neighbor is not in use, it won't be found in beads_in_use_tree
//                        pSwap2 = std::find(beginIt, endIt, neighbor);
//                        localPotentialAtNewPosition = eta * (calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, *pSwap2));
//
//                        beads_in_use_tree.insert(*pSwap2);
//
//                        localPotentialAtNewPositionWithBead = eta * (calculateLocalContactPotentialOfNeighborhood(
//                                                                      &beads_in_use_tree, pModel, *pSwap2) +
//                                                               calculateLocalContactPotentialPerBead(
//                                                                       &beads_in_use_tree, pModel, *pSwap2));
//
//                        std::iter_swap(itIndex, pSwap2);
//                        std::sort(beginIt, beginIt + workingLimit); // bead_indices needs to be sorted
//
//                        addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere,
//                                workingBinCount);
//
//                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
//                        testKL = pData->calculateKLDivergence(workingBinCount);
//
//                        isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,
//                                             tempNumberOfComponents);  // 10x longer than contactEnergy calculation
//
//                        this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents);
//
//                        if ((this_energy + localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead) <
//                            (current_energy + localPotentialAtOldPosition + localPotentialAtNewPosition)) {
//
//                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
//                            currentKL = testKL;
//                            current_energy = this_energy;
//                            currentNumberOfComponents = tempNumberOfComponents;
//                            isSwapped = true;
//                            swappedCounter++;
//                            printf("*==>   SWAPPED INDEX => %i               CONTACTS: %i (%.1f)\n", swap1, contactsToConsider, average_number_of_contacts);
//                            break;
//                        }
//                        // reverse changes; find index in sorted list and replace
//                        // pSwap2 points to &bead_indices[bb]
//                        // itIndex
//                        // binCount is P(r) distribution without swap1 value
//                        std::copy(beginBinCount, endBinCount,
//                                  workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
//                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
//                        itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
//                        beads_in_use_tree.erase(neighbor);
//
//                    } else if (neighbor == -1) {
//                        break;
//                    }
//                }
//
//                if (!isSwapped){ // itIndex already reverses on exit from loop
//                    std::sort(beginIt, beginIt + workingLimit);
//                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
//                    beads_in_use_tree.insert(swap1); // add back swap one to tree
//                }
//            } // end Contacts check of selected bead
        //}

        for (int i=0; i < workingLimit; i++){ // find a lattice point that has too few contacts

            swap1 = active_indices[i];
            numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);

            if (numberContacts < contactsLimit){

                cout << "------------------------------------------------------------------------------" << endl;
                printf("      CHECKING INDEX => %i               CONTACTS: %i (%.1f)\n", swap1, numberContacts, average_number_of_contacts);


                //swap1 defined during setting of placesToCheck
                isSwapped = false;
                // swap
                itIndex = std::find(beginIt, beginIt + workingLimit, swap1); // locate lattice within workingLimit
                // remove selected index from P(r)

                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

                // current local potential
                localPotentialAtOldPosition = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1) +
                                                   calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, swap1));

                beads_in_use_tree.erase(swap1);
                localPotentialAtOldPositionWithOutBead = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, swap1));
                // find better position

                for (int bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads

                    // make the swap, sort and update P(r)
                    pSwap2 = (bead_indices.begin() + bb);
                    //contactsToConsider = numberOfContactsFromSet(&beads_in_use_tree, pModel, *pSwap2);
                    contactsToConsider = numberOfContactsFromSetExclusive(&beads_in_use_tree, pModel, *pSwap2, *itIndex);
                    // if ( (contactsToConsider > 1) && (contactsToConsider < contactsLimit) ) { // only move to location with few contacts
                        originalSwap2Value = *pSwap2;

                        localPotentialAtNewPosition = eta * (calculateLocalContactPotentialOfNeighborhood(
                                &beads_in_use_tree, pModel, *pSwap2));

                        beads_in_use_tree.insert(*pSwap2);

                        localPotentialAtNewPositionWithBead = eta *
                                                              (calculateLocalContactPotentialOfNeighborhood(
                                                                      &beads_in_use_tree, pModel, *pSwap2) +
                                                               calculateLocalContactPotentialPerBead(
                                                                       &beads_in_use_tree, pModel, *pSwap2));

                        std::iter_swap(itIndex, pSwap2);
                        std::sort(beginIt, beginIt + workingLimit); // bead_indices needs to be sorted

                        addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere,
                                workingBinCount);

                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        testKL = pData->calculateKLDivergence(workingBinCount);

                        isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere,
                                             tempNumberOfComponents);  // 10x longer than contactEnergy calculation

                        this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents);

                        if ((this_energy + localPotentialAtOldPositionWithOutBead +
                             localPotentialAtNewPositionWithBead) <
                            (current_energy + localPotentialAtOldPosition + localPotentialAtNewPosition)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            isSwapped = true;
                            swappedCounter++;
                            printf("*==>   SWAPPED INDEX => %i               CONTACTS: %i (%.1f)\n", swap1, contactsToConsider, average_number_of_contacts);
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // pSwap2 points to &bead_indices[bb]
                        // itIndex
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount,
                                  workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = std::find(beginIt, beginIt + workingLimit, swap1);
                        //}
                        beads_in_use_tree.erase(originalSwap2Value);
                    //}
                } // inner for loop exit on break

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(beginIt, beginIt + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                }

                // search whole space and breaking on first point is slow
                //break;
            }
        }

//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)) {
//                cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " <<
//                testKL1 << endl;
//                return "stopped";
//        }

        relative_diff = abs(current_energy - energy_i)/(0.5*(current_energy + energy_i));
        energy_i = current_energy;

        if (relative_diff <= 0.0001){
            deltaECount++;

            if (deltaECount > 31){
                break;
            }
        } else {
            deltaECount = 0;
        }
        cout << "------------------------------------------------------------------------------" << endl;
        printf("      COUNT => %i (%i)              delta_D_KL: %.7f\n", round, finalRounds, relative_diff);
        cout << "------------------------------------------------------------------------------" << endl;
    }

    float average_number_of_contacts_before = tempAverageContacts/(float)workingLimit;
    tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }
    average_number_of_contacts = tempAverageContacts/(float)workingLimit;

    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " FINAL LOW TEMP REFINEMENT" << endl;
    printf("    SWAPPED => %i             \n", swappedCounter);
    printf("     ROUNDS => %i (%i)            \n", round, finalRounds);
    cout << "              KL DIVERGENCE   ENERGY " << endl;
    printf("     BEFORE =>    %.4E           %.4E  \n", beforeFinalKL, beforeFinalEnergy);
    printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    cout << "    AVERAGE" << endl;
    printf("     BEFORE =>    %.2f  \n", average_number_of_contacts_before);
    printf("      AFTER =>    %.2f  \n", average_number_of_contacts);
    cout << "------------------------------------------------------------------------------" << endl;

    this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);


    // final round to move any points close to body
    pModel->setCVXHullVolume(calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel));
    pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    string nameTo = filenameprefix + "_annealed_" + std::to_string(iteration);
    pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");
    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, nameTo, this, pData);

    return nameOfModel;
}

