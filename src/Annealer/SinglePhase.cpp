//
// Created by Robert Rambo on 13/01/2017.
//
#include "../Anneal.h"
#include "Data.h"
#include "../Model.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"

using namespace std;

void Anneal::createInitialModelCVXHull(Model *pModel, Data *pData, std::string name) {

    srand(time(0));

    //interconnectivityCutOff = pModel->getBeadRadius()*2.8285;
    //interconnectivityCutOff = pModel->getBeadRadius()*2.0001;
    //interconnectivityCutOff = pModel->getBeadRadius()*3.464*1.001;
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

    float beadVolume = pModel->getBeadVolume();
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
    std::set<int> beads_in_use_tree(beginIt, beginIt + workingLimit);

    // setup parameters for hull
    char flags[25];
    sprintf(flags, "qhull s FA");
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];

    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    populateLayeredDeadlimit(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?

    int tempNumberOfComponents, alterMe;
    isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

    EulerTour eulerTour(beginIt, workingLimit, pModel);
    tempNumberOfComponents = eulerTour.getNumberOfComponents();
    int currentNumberOfComponents = tempNumberOfComponents;

    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();

    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    lambda = 0.001;
    mu = 0.00001; // scale this to volume of search space

    float testKL, test_energy, current_energy = currentKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*current_volume/(float)workingLimit;

    float lowest_energy = current_energy;

    int testIndex;
    float inv_kb_temp = 1.0/(float)highT;
    inv_kb_temp = 1.0/0.0001;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    cout << " POPULATING DEAD LIMIT" << endl;
    //populateLayeredDeadlimit(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere);
    refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
    cout << " DEAD LIMIT SET => " << deadLimit << endl;

    int lowestWorkingLimit, lowestDeadLimit;
    lowestWorkingLimit = workingLimit;
    std::copy(beginIt, endIt, lowest_bead_indices.begin());
    lowestDeadLimit = deadLimit;
    bool updateCVX =false;

    std::clock_t start;

    cout << " STARTING CONSTANT TEMP SEARCH " << endl;
    for (int high=0; high < highTempRounds; high++) { // iterations during the high temp search

        std::sort(beginIt, beginIt + workingLimit);
        std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy

        if (distribution(gen) < 0.95){

            alterMe = number_of_beads_to_use(gen);  // uniform distribution

            if (alterMe > workingLimit){
                cout << "*******************             ADD                 ******************* "  << endl;
                // randomly add from deadLimit

                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;

                testIndex = bead_indices[randomSpot];
                itIndex = beginIt + randomSpot;

                //make the swap at the border (workingLimit)
                addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);

                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->calculateKLDivergence(binCount);
                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1

                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    //copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
                    //std::iter_swap(active_indices.begin() + workingLimit - 1, active_indices.begin() + d);
                    //adjustDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                    //copy(beginIt, endIt, backUpState.begin());
                    beads_in_use_tree.insert(testIndex);
                    //copy(beginIt + workingLimit, beginIt + deadLimit, active_indices.begin() + workingLimit);
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {

                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    //copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
                    //std::iter_swap(active_indices.begin() + workingLimit - 1, active_indices.begin() + d);
                    //adjustDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                    //copy(beginIt, endIt, backUpState.begin());
                    beads_in_use_tree.insert(testIndex);
                    //copy(beginIt + workingLimit, beginIt + deadLimit, active_indices.begin() + workingLimit);
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&beginIt, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    currentNumberOfComponents = eulerTour.removeNode(testIndex);
                }

            } else if (alterMe < workingLimit) {

                cout << "*******************             REMOVE                 ******************* "  << endl;

                if (workingLimit > lowerN){
                    // go through entire list
                    int randomSpot = rand() % workingLimit;

                    // grab from randomized active_indices list
                    testIndex = bead_indices[randomSpot];
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = pData->calculateKLDivergence(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy ) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;

                        //std::iter_swap(active_indices.begin() + i, active_indices.begin() + workingLimit); // move swap position
                        //removeNeighborsFromDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                        beads_in_use_tree.erase(testIndex);
                    } else if (exp(-(test_energy - current_energy)*inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;

                        //std::iter_swap(active_indices.begin() + i, active_indices.begin() + workingLimit);
                        //removeNeighborsFromDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                        beads_in_use_tree.erase(testIndex);
                    } else { // undo changes and move to next bead (rejecting)
                        // std::copy(backUpState.begin(), backUpState.begin()+workingLimit+1, beginIt);
                        currentNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

            cout << "*******************             ADD/REMOVE                 *******************" << endl;
            printf("   MAXSTEPS => %i (%4i) \n", highTempRounds, high);
            printf("      GRAPH => %3i\n", currentNumberOfComponents);
            printf("      UPPER => %i LOWER => %i \n", upperN, lowerN);
            printf("     VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, mu, mu*current_volume);
            printf("LIMIT: %5i DEADLIMIT: %5i D_KL: %.4E ENRGY: %.4E \n", workingLimit, deadLimit, currentKL, current_energy);
            cout << "*******************                                        *******************" << endl;

        } else {

            cout << "*******************             POSITIONAL                 *******************" << endl;
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
            int swapCount=1;
            int swap1, bb;
            int * pSwap2;
            bool isSwapped;

            // only move CVX hull points
            std::vector<int> indices_to_check(totalV); // large vector ~1000's
            for (int v = 0; v < totalV; v++) { //
                indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                vertices = vertices->next;
            }
            qh_freeqhull(true);

            //pModel->writeModelToFile(deadLimit, bead_indices, "positional"+ std::to_string(high));
            //pModel->writeModelToFile(totalV, indices_to_check, "hull"+ std::to_string(high));

            //int randomSpot = rand() % totalV;
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);

            for (int v = 0; v < 3; v++) { // try 3 times

                swap1 = indices_to_check[v];
                isSwapped = false;
                // find bead to swap in active set
                itIndex = find(beginIt, beginIt + workingLimit, swap1);
                // remove selected index from P(r)
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)
                // find better position
                eulerTour.removeNode(swap1);

                for (bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
                    // make the swap, sort and update P(r)
                    pSwap2 = &bead_indices[bb]; // located in workingZone
                    *itIndex = *pSwap2;         // swapped here (at this point bead_indices contains duplicate values)

                    std::sort(beginIt, beginIt + workingLimit); // bead_indices needs to be sorted

                    addToPr(*pSwap2, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                    // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                    testKL = pData->calculateKLDivergence(workingBinCount);

                    tempNumberOfComponents = eulerTour.addNode(*pSwap2, pModel);
                    //   isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, interconnectivityCutOff, tempNumberOfComponents, pModel);

                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel);
                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy) {
                        *pSwap2 = swap1; // commit the swap in deadzone
                        std::iter_swap(beginIt + bb, beginIt + deadLimit - swapCount); // move the swapped lattice point to the end
                        std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        current_energy = test_energy;
                        isSwapped = true;
                        updateCVX = true;
                        swapCount++;
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        *pSwap2 = swap1; // commit swap in deadzone
                        std::iter_swap(beginIt + bb, beginIt + deadLimit - swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
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
                    // std::copy(backUpState.begin(), backUpState.begin()+workingLimit, beginIt);
                    std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                    itIndex = find(beginIt, beginIt+workingLimit, *pSwap2); // find swapped value in the sorted array
                    *itIndex = swap1; // replace with original value, beginIt is unsorted
                    currentNumberOfComponents=eulerTour.removeNode(*pSwap2);
                }
                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (bb >= deadLimit && !isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(beginIt, beginIt + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount);
                    currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                }
            }

            cout << " SWAPPED POINTS " << swapCount << " <=> " << totalV <<  endl;
        }

        //start = std::clock();

        if (updateCVX){
            //refineCVXHull(bead_indices, active_indices, totalBeadsInSphere, workingLimit, &deadLimit, pModel);
            populateLayeredDeadlimit(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere);
            updateCVX=false;
        }

        //std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        float currentKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        cout << "########################<<<<<>>>>>>############################## " << endl;
//        cout << currentKL << " <=> " << currentKL1 << endl;
//        string name = "search_" + std::to_string(high);
//        pModel->writeModelToFile(workingLimit, bead_indices, name);
//        name = "search_skin" + std::to_string(high);
//        pModel->writeModelToFile(deadLimit, bead_indices, name);

        if (current_energy < lowest_energy){
            lowest_energy = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            std::copy(beginIt, endIt, lowest_bead_indices.begin());
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }

    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name);
    pModel->writeModelToFile(lowestDeadLimit, lowest_bead_indices, "initial_search_layer");
    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setStartingDeadLimit(lowestDeadLimit);

}


string Anneal::refineHomogenousBodyASACVX(Model *pModel, Data *pData, int iteration){
    cout << "########################<<<<<>>>>>>############################## " << endl;
    cout << "#                                                               # " << endl;
    cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << endl;
    cout << "#                                                               # " << endl;
    cout << "########################<<<<<>>>>>>############################## " << endl;

    pModel->setBeadAverageAndStdev(pModel->getStartingWorkingLimit(), 0.17*pModel->getStartingWorkingLimit());
    float oldN = pModel->getVolumeAverage();     // need to reset this for modeling via rounds
    float oldStdev = pModel->getVolumeStdev();

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
    beginIt = bead_indices.begin();
    endIt = bead_indices.end();

    pModel->centerLatticeModel(workingLimit, bead_indices);
    std::set<int> beads_in_use_tree(beginIt, beginIt + workingLimit);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);

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
    int original, addMe;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float inv_kb_temp, relative_diff;
    float tempAverageContacts;

    // coupon collector's problem
    int updateCount = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);

    float step_limit = (updateCount < 259) ? 5900 : 7*updateCount;
    cout << " WORKING LIMIT => " << workingLimit << " " << step_limit << endl;
    int deadUpdate = updateCount*0.1;

    std::vector<float> tempDuringRun(step_limit);
    std::vector<float> divergenceDuringRun(step_limit);
    std::vector<int> workingLimitDuringRun(step_limit);

    float * pTempDuringRun = &tempDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(beginIt, workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float sum_x_squared=0, sum_x=0, divideBy=0, average_x, stdev, acceptRate = 0.5, inv500 = 1.0/500.0;

    int counter=1;
    inv_kb_temp = 1.0/this->lowTempStop;
    int deltaECount =0;
    float energy_i, this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    lambda = 0.0001;
    beta = 0.0;

    float tempTotalContactEnergy, totalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
    // want D_KL to settle < 10^-5 to 10^-6
    //int span = (int)(std::floor(log10f(currentKL)) - std::floor(log10f(totalContactEnergy/(float)workingLimit))) + 3;

    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    //
    std::normal_distribution<float> volumeGen(pModel->getVolumeAverage(), pModel->getVolumeStdev());

    std::clock_t startTime;
    int attempts=0, failures=0;
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    double runtime;
    bool stopMe = false;

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

        if (distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            beginIt = bead_indices.begin();
            endIt = bead_indices.end();

            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               ADD?REMOVE               *******************" << endl;
            alterMe = (int) volumeGen(gen);
            // build a list of indices within defined region of convex hull
            if (alterMe > workingLimit){ // ADD BEAD?

                cout << "*******************                  ADD                   *******************" << endl;
                //int start = workingLimit; // find first bead I can add, then break
                sprintf(addRemoveText, "     ADD => %i", 1);

                float beforeAdding, afterAdding;
                int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                // go through each selected bead, determine number of contacts, if less than average, consider adding
                startTime = std::clock();

                // select first element from randomized active_set
                // check if available neighbor can be added
                std::vector<int>::iterator it = pModel->getPointerToNeighborhood(bead_indices[randomSpot]);
                for (int i=0; i< totalNeighbors; i++){ // run into problem of check
                    int neighbor = *(it+i);
                    if ((neighbor > -1) && beads_in_use_tree.find(neighbor) == beads_in_use_tree.end()){
                        // if neighbor is not in use, it won't be found in beads_in_use_tree
                        itIndex = std::find(beginIt + workingLimit, beginIt + deadLimit, neighbor);

                        if (itIndex != (beginIt + deadLimit)){  // if itIndex is not within workignLimit skip

                            beforeAdding = eta*calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, neighbor);

                            // make the swap at the border (workingLimit)
                            addLatticPositionToModel(&beginIt, &endIt, &backUpState, &workingLimit, &itIndex);

                            addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                            testKL = pData->calculateKLDivergence(binCount);

                            isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);

                            beads_in_use_tree.insert(neighbor);
                            afterAdding = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, neighbor) +
                                               calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, neighbor));

                            tempTotalContactEnergy = afterAdding - beforeAdding;

                            if ( (this_energy + afterAdding) < (current_energy + beforeAdding) ) {
                                currentKL = testKL;
                                current_energy = this_energy;
                                currentNumberOfComponents = tempNumberOfComponents;
                                totalContactEnergy = tempTotalContactEnergy;
                                isUpdated = true;
                                break;
                            } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                                currentKL = testKL;
                                current_energy = this_energy;
                                currentNumberOfComponents = tempNumberOfComponents;
                                totalContactEnergy = tempTotalContactEnergy;
                                isUpdated = true;
                                break;
                            } else { // undo changes (rejecting)
                                beads_in_use_tree.erase(neighbor);
                                restoreAddingFromBackUp(&beginIt, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                            }
                        }

                    } else if (neighbor == -1) {
                        break;
                    }
                }

                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

            } else { // REMOVE BEADS?
                cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
                priorWorkingLimit = 1;
                sprintf(addRemoveText, "     REMOVE => %i", priorWorkingLimit);

                beginIt = bead_indices.begin();

                float beforeRemoving, afterRemoving;
                startTime = std::clock();

                // grab from randomized active_indices list
                int randomSpot = rand() % workingLimit;
                original = bead_indices[randomSpot];

                beforeRemoving = eta*(calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, original) +
                                      calculateLocalContactPotentialPerBead(&beads_in_use_tree, pModel, original));

                removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                // still need to sort, swap changes the order
                // don't remove lattice point if model is not interconnected
                testKL = pData->calculateKLDivergence(binCount);

                //tempNumberOfComponents = eulerTour.removeNode(original);
                isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);

                beads_in_use_tree.erase(original);
                afterRemoving = eta*calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, original);

                tempTotalContactEnergy = afterRemoving - beforeRemoving;

                if ((this_energy + afterRemoving) < (current_energy + beforeRemoving) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    //tempNumberOfComponents = eulerTour.addNode(original, pModel);
                }

                runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
            }

            printf("       TEMP => %-.8f \n     ACCEPT => %.4f\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop , acceptRate, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
            printf("  UPDATECNT => %7i TIME => %f SECONDS\n", deadUpdate, runtime);
            printf("      GRAPH => %3i\n", currentNumberOfComponents);
            printf(" LATTCE AVG => %.0f        STDEV => %0.3f  %s\n", pModel->getVolumeAverage(), pModel->getVolumeStdev(), addRemoveText);
            printf(" DELTA_D_KL => %.7E\n", relative_diff);
            printf("   CONTACTE => %.4E ETA => %.4E \n", totalContactEnergy, eta);
            printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, current_energy);

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            bool isSwapped;
            attempts +=1;

            std::shuffle(bead_indices.begin()+workingLimit, bead_indices.begin() + deadLimit, gen);           // randomizes order of positions outside workingLimit
            std::copy(bead_indices.begin(), bead_indices.begin() + workingLimit, active_indices.begin());     // copy working set into active indices
            std::shuffle(active_indices.begin(), active_indices.begin() + workingLimit, gen); // randomizes order of beads to move

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               POSITIONAL               *******************" << endl;
            cout << "*******************                                        *******************" << endl;

            float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
            //
            // find bead to swap in active set
            // only refine a single position with 2 or less contacts
            //
            int originalSwap2Value, numberContacts, secondaryN;
            int contactsLimit = contactsPerBead + 1;
            //std::vector<int>::iterator primaryNeighborhood, secondaryNeighborhood;
            //std::set<int> placesToCheck;

            startTime = std::clock();

            // find first lattice point with too few contacts
            for (int i=0; i < workingLimit; i++) {
                swap1 = active_indices[i];
                numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);
                if (numberContacts < contactsLimit){
                    // select first element from randomized active_set
                    // check if available neighbor can be added
//                    primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);
//
//                    for (int n=0; n < totalNeighbors; n++){ // run into problem of check
//                        int neighbor = *(primaryNeighborhood + n);
//
//                        if (neighbor > -1){
//                            set<int>::iterator inSet = beads_in_use_tree.find(neighbor);
//
//                            if (inSet == beads_in_use_tree.end()){ // not in use, so its in deadlimit
//                                // if number of contacts at new position (less current) is 0, skip
//                                placesToCheck.insert(neighbor);
//
//                            } else { // its already a neighbor within WorkSet so check its neighborhood
//                                secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);
//
//                                for (int s=0; s < totalNeighbors; s++){
//                                    // adding to any of these already insures at least 1 contact
//                                    // check if secondary neighbor is in beads_in_use
//                                    secondaryN = *(secondaryNeighborhood + s);
//                                    inSet = beads_in_use_tree.find(secondaryN);
//                                    if (secondaryN > -1 && inSet == beads_in_use_tree.end()){ // if .end(), means not in use
//                                        placesToCheck.insert(secondaryN);
//                                    } else if (secondaryN == -1){
//                                        break;
//                                    }
//                                }
//                            }
//                        } else if (neighbor == -1) {
//                            break;
//                        }
//                    }
                    break;
                }
            }

            beginIt = bead_indices.begin();

            // too local of a search ???

            // if (numberContacts < contactsLimit ){
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

            int contactsToConsider;

            for (int bb = workingLimit; bb < deadLimit; bb++) { // go through and find a better position within unused beads
                // make the swap, sort and update P(r)
                pSwap2 = (bead_indices.begin() + bb);
                contactsToConsider = numberOfContactsFromSet(&beads_in_use_tree, pModel, *pSwap2);

                if ( (contactsToConsider > 1) && (contactsToConsider < contactsLimit) ) { // only move to location with few contacts

                    originalSwap2Value = *pSwap2;

                    localPotentialAtNewPosition = eta * (calculateLocalContactPotentialOfNeighborhood(&beads_in_use_tree, pModel, *pSwap2));
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

                    tempTotalContactEnergy =
                            localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead -
                            localPotentialAtOldPosition - localPotentialAtNewPosition;

                    if ((this_energy + localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead) <
                        (current_energy + localPotentialAtOldPosition + localPotentialAtNewPosition)) {

                        copy(workingBinCount.begin(), workingBinCount.end(),
                             beginBinCount); // final copy so next round make backup
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = tempTotalContactEnergy;
                        isUpdated = true;
                        isSwapped = true;
                        break;
                    } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) >
                               distribution(gen)) {

                        copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = tempTotalContactEnergy;
                        isUpdated = true;
                        isSwapped = true;
                        break;
                    }
                    // reverse changes; find index in sorted list and replace
                    // pSwap2 points to &bead_indices[bb]
                    // itIndex
                    // binCount is P(r) distribution without swap1 value
                    std::copy(beginBinCount, endBinCount,
                              workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                    std::copy(backUpState.begin(), backUpState.begin() + (bb + 1), beginIt);
                    itIndex = std::find(beginIt, beginIt + workingLimit, swap1);
                    //}
                    beads_in_use_tree.erase(originalSwap2Value);
                }
            } // inner for loop exit on break

            // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
            if (!isSwapped){ // itIndex already reverses on exit from loop
                std::sort(beginIt, beginIt + workingLimit);
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                failures +=1;
            }


            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
            printf("       TEMP => %-.8f \n     ACCEPT => %.4f\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop, acceptRate, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
            printf("  UPDATECNT => %7i TIME => %f SECONDS\n", deadUpdate, runtime);
            printf("      GRAPH => %3i\n", currentNumberOfComponents);
            printf(" LATTCE AVG => %.0f        STDEV => %0.3f  \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
            printf(" DELTA_D_KL => %.7E\n", relative_diff);
            printf("   CONTACTE => %.4E ETA => %.4E \n", totalContactEnergy, eta);
            printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, current_energy);
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

//        cout << "______________________________________________________________________________" << endl;
//        cout << "*******************                 TEST                   *******************" << endl;
//        cout << "*******************              -----------               *******************" << endl;
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

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
        } else {
            acceptRate = inv500*(499*acceptRate);
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

            pModel->setBeadAverageAndStdev(average_x, stdev);
            string name = "newAverage_" + std::to_string(counter);
            pModel->writeModelToFile(workingLimit, bead_indices, name);
            //name = "newSkin_" + std::to_string(counter);
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
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
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

    for (int round = 0; round < 0*finalRounds; round++){

        std::shuffle(beginIt+workingLimit, beginIt + deadLimit, gen); // randomizes order of beads to select
        std::copy(beginIt, beginIt + workingLimit, active_indices.begin());
        std::copy(beginIt, endIt, backUpState.begin());
        std::shuffle(active_indices.begin(), active_indices.begin() + workingLimit, gen); // randomizes order of beads to select
        cout << "______________________________________________________________________________" << endl;
        cout << "*******************                 FINAL                  *******************" << endl;
        cout << "*******************               POSITIONAL               *******************" << endl;

        float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
        int originalSwap2Value;
        bool isSwapped;

        int numberContacts, secondaryN, contactsToConsider;
        int contactsLimit = contactsPerBead + 1;
        std::vector<int>::iterator primaryNeighborhood, secondaryNeighborhood;
//        std::set<int> placesToCheck;
//        for (int i=0; i < workingLimit; i++) {
//            swap1 = active_indices[i];
//            numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);
//            if (numberContacts < average_number_of_contacts){
//                // select first element from randomized active_set
//                // check if available neighbor can be added
//                primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);
//
//                for (int n=0; n < totalNeighbors; n++){ // run into problem of check
//                    int neighbor = *(primaryNeighborhood + n);
//
//                    if (neighbor > -1){
//                        set<int>::iterator inSet = beads_in_use_tree.find(neighbor);
//
//                        if (inSet == beads_in_use_tree.end()){ // not in use, so its in deadlimit
//                            // if number of contacts at new position (less current) is 0, skip
//                            placesToCheck.insert(neighbor);
//
//                        } else { // its already a neighbor within WorkSet so check its neighborhood
//                            secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);
//
//                            for (int s=0; s < totalNeighbors; s++){
//                                // adding to any of these already insures at least 1 contact
//                                // check if secondary neighbor is in beads_in_use
//                                secondaryN = *(secondaryNeighborhood + s);
//                                inSet = beads_in_use_tree.find(secondaryN);
//                                if (secondaryN > -1 && inSet == beads_in_use_tree.end()){ // if .end(), means not in use
//                                    placesToCheck.insert(secondaryN);
//                                } else if (secondaryN == -1){
//                                    break;
//                                }
//                            }
//                        }
//                    } else if (neighbor == -1) {
//                        break;
//                    }
//                }
//                break;
//            }
//        }
//
//        int sizeOfPlaces = placesToCheck.size();
        beginIt = bead_indices.begin();

        for (int i=0; i < workingLimit; i++){

            swap1 = active_indices[i];
            numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1);

            if (numberContacts < average_number_of_contacts){
                cout << "------------------------------------------------------------------------------" << endl;
                printf("      CHECKING INDEX => %i               CONTACTS: %i (%.1f)\n", swap1, numberContacts, average_number_of_contacts);
                cout << "------------------------------------------------------------------------------" << endl;

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
                    contactsToConsider = numberOfContactsFromSet(&beads_in_use_tree, pModel, *pSwap2);

                    if ( (contactsToConsider > 1) && (contactsToConsider < contactsLimit) ) { // only move to location with few contacts
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

                        tempTotalContactEnergy =
                                localPotentialAtOldPositionWithOutBead + localPotentialAtNewPositionWithBead -
                                localPotentialAtOldPosition - localPotentialAtNewPosition;

                        if ((this_energy + localPotentialAtOldPositionWithOutBead +
                             localPotentialAtNewPositionWithBead) <
                            (current_energy + localPotentialAtOldPosition + localPotentialAtNewPosition)) {

                            copy(workingBinCount.begin(), workingBinCount.end(),
                                 beginBinCount); // final copy so next round make backup
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = tempNumberOfComponents;
                            totalContactEnergy = tempTotalContactEnergy;
                            isSwapped = true;
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // pSwap2 points to &bead_indices[bb]
                        // itIndex
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount,
                                  workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.begin() + (bb + 1), beginIt);
                        itIndex = std::find(beginIt, beginIt + workingLimit, swap1);
                        //}
                        beads_in_use_tree.erase(originalSwap2Value);
                    }
                } // inner for loop exit on break

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(beginIt, beginIt + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                }

                // search whole space and breaking on first point is slow
                break;
            }
        }

        //float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
        //if (currentKL != testKL1 || checkForRepeats(bead_indices)){
        //if (currentKL != testKL1){
        //    cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
        //    return "stopped";
        //}

        relative_diff = abs(current_energy - energy_i)/(0.5*(current_energy + energy_i));
        energy_i = current_energy;

        if (relative_diff <= 0.0001){
            deltaECount++;

            if (deltaECount > 9){
                break;
            }
        } else {
            deltaECount = 0;
        }
        cout << "------------------------------------------------------------------------------" << endl;
        printf("      COUNT => %i (%i)              delta_D_KL: %.7f\n", round, finalRounds, relative_diff);
        cout << "------------------------------------------------------------------------------" << endl;
    }

    // final round to move any points close to body
    pModel->setCVXHullVolume(calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel));
    pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    string nameTo = filenameprefix + "_annealed_" + std::to_string(iteration);
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");
    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, nameTo, this, pData);

    //for (int i=0; i<workingLimit; i++){
    //    cout << i << "  " << (numberOfContacts(bead_indices[i], &bead_indices, workingLimit, contactCutOff, pModel, pDistance)) << endl;
    //}

    return nameOfModel;
}