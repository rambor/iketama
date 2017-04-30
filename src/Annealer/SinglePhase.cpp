//
// Created by Robert Rambo on 13/01/2017.
//
#include "../Anneal.h"
#include "PDBModel.h"
#include "../EulerTour/EulerTour.h"

using namespace std;

bool Anneal::createInitialModelCVXHull(Model *pModel, Data *pData, std::string name) {

    srand(time(0));
    contactCutOff = interconnectivityCutOff;
    // convert distances to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * pBin = pModel->getPointerToBins(); // initialized as empty in Model class
    // faster if pointers to above is reused - conversion happens only once
    // pBin will be reused many times
    maxbin=0;
    for(unsigned long int i=0; i < totalDistancesInSphere; i++){
        *(pBin+i) = pData->convertToBin(*(pDistance + i)); // some distances will exceed dmax
        if (*(pBin+i) > maxbin){
            maxbin = *(pBin+i);
        }
    }

    maxbin += 1;

    std::cout << "MAX BIN : " << maxbin << std::endl;

    int capacity=11;
    int averaging[capacity];
    int averagingIndex=0;

    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin); // smallish vector, typically < 50
    //std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);  // smallish vector, typically < 50


    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    //std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);
    std::vector<int>::iterator beginIt;
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
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "universe", 0);
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    cout << " POPULATING DEAD LIMIT" << endl;
    int deadLimit; // as bead indices are discarded, set upper limit of vector
    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    //populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    cout << "    DEADLIMIT SET => " << deadLimit << endl;
    cout << " WORKINGLIMIT SET => " << workingLimit << endl;

    // setup parameters for hull
    char flags[] = "qhull FA";
    //int numpoints = 3*totalBeadsInSphere;
    coordT points[3*(upperN+(int)(0.30*upperN))];

//    float lowest_volume, test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    float initialVolume = current_volume;
//
//    // area of the model P(r) distribution
//    cout << "       CVX VOLUME => " << current_volume  << endl;

    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int alterMe, tempNumberOfComponents, currentNumberOfComponents  = eulerTour.getNumberOfComponents();
    //int currentNumberOfComponents;
    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, currentNumberOfComponents);
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<int>::iterator beginBinCount = binCount.begin(), itIndex;
    std::vector<int>::iterator endBinCount = binCount.end();

    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    float lowest_volume, test_volume, current_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);//calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float initialVolume = current_volume;

    // area of the model P(r) distribution
    cout << "       CVX VOLUME => " << current_volume  << " WL => " << workingLimit << endl;
    //lambda = 0.001;
    mu = 0.000001; // scale this to volume of search space

    float testKL, test_energy, current_energy = currentKL + lambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + mu*current_volume/(float)workingLimit;

    cout << " STARTING KL => " << currentKL << endl;
    float lowest_energy = current_energy;

    int testIndex;
    //float inv_kb_temp = 1.0/(float)highT;
    float inv_kb_temp = 1.0/0.0001;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    int lowestWorkingLimit, lowestDeadLimit;
    lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    lowestDeadLimit = deadLimit;
    double sumIt = 0.0, counter=0.0;

    bool updateCVX =false;

    //int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();  seems on average there are less than 12 available to look at

    cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< endl;

    int high;

    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
        //std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
        std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased

        if (distribution(gen) < 0.73){

            alterMe = number_of_beads_to_use(gen);  // uniform distribution

            beginIt = bead_indices.begin();

            if (alterMe > workingLimit){
                cout << "*******************                  ADD                   ******************* "  << endl;
                // randomly add from deadLimit
                //int randomSpot = rand() % (deadLimit - workingLimit) + workingLimit;
                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while ( testIndex == -1) { // find a new position
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

                itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
//                testIndex = bead_indices[randomSpot];
//                itIndex = bead_indices.begin() + randomSpot;
                //make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);

                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->calculateKLDivergence(binCount);

                //test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                test_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX = true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX = true;
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                }

            } else if (alterMe < workingLimit) {

                cout << "*******************                 REMOVE                 ******************* "  << endl;

                if (workingLimit > lowerN){
                    // go through entire list
                    int randomSpot = randomIndex(gen);

                    // grab from randomized active_indices list
                    testIndex = bead_indices[randomSpot];
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = pData->calculateKLDivergence(binCount);
                    //test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                    test_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy ) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else if (exp(-(test_energy - current_energy)*inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else { // undo changes and move to next bead (rejecting)
                        // std::copy(backUpState.begin(), backUpState.begin()+workingLimit+1, beginIt);
                        eulerTour.addNode(testIndex, pModel);
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

        } else {

            cout << "*******************               POSITIONAL               *******************" << endl;

            int swap1, bb;
            std::vector<int>::iterator pSwap2;

//            std::set<int> extremes;
//            this->getExtremes(bead_indices, workingLimit, extremes, pModel);
//            std::shuffle(extremes.begin(), extremes.end(), gen);
//            std::cout << " NUMBER EXTREMES " << extremes.size() << std::endl;
//            std::vector<int> indices_to_check(extremes.size());
//            std::copy(extremes.begin(), extremes.end(), indices_to_check.begin());
            // recalculate
            // calculate convex hull and get hull points
            // can be threaded
            std::vector<int> active_indices(workingLimit);
            for (int i = 0; i < workingLimit; i++) {
                beadToPoint(&points[i*3], pModel->getBead(bead_indices[i]));
                active_indices[i] = bead_indices[i];
            }

            // needs to be optimized
            qh_new_qhull(3, workingLimit, points, 0, flags, NULL, NULL);
            vertexT * vertices = qh vertex_list;
            int totalV = qh num_vertices;
//
//            // make copy of points from vertices and then free qh_hull
            bool isSwapped;
//
//            // only move CVX hull points
            std::vector<int> indices_to_check(totalV); // large vector ~1000's
            for (int v = 0; v < totalV; v++) { //
                indices_to_check[v] = active_indices[qh_pointid( vertices->point)];
                vertices = vertices->next;
            }

            qh_freeqhull(true);


            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            //int swapCount = 0;
            int neighbor;

//            for(std::set<int>::iterator it=extremes.begin(); it != extremes.end(); ++it){
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
                beads_in_use_tree.erase(swap1);
                eulerTour.removeNode(swap1);
                int randind;

                for (bb = 0; bb < 5; bb++) { // go through and find a better position within unused beads
                    // find a bead in use, and see if it has an empty neighbor to use
                    randind = randomIndex(gen);
                    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randind]);
                    while (neighbor == -1 || neighbor == swap1){
                        randind = randomIndex(gen);
                        neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randind]);
                    }
                    //std::cout << v << " bb: "  << bb << " SWAP1 " << swap1 << " NEIGH => " << neighbor << " RND " << randind <<  " WL " << workingLimit << std::endl;
                    // if neighbor is not in use, then find it in bead_indices and assign to pSwap2
                    // make the swap, sort and update P(r)
                    //pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), neighbor);
                    pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), neighbor); // potential problem

                    std::iter_swap(itIndex, pSwap2);
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                    addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                    // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                    testKL = pData->calculateKLDivergence(workingBinCount);
                    tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);
                    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);

                    //test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                    test_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);
                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + mu*test_volume/(float)workingLimit;

                    if (test_energy < current_energy) { //move swapped position to the deadLimit
                        beads_in_use_tree.insert(neighbor);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
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
                        currentNumberOfComponents = tempNumberOfComponents;
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
                } // end of neighborhood for-loop

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    eulerTour.addNode(swap1, pModel);
                    beads_in_use_tree.insert(swap1);
                }
            }
        }

        if (updateCVX){
            updateCVX=false;
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
//            string name = "search_" + std::to_string(high);
//            pModel->writeModelToFile(workingLimit, bead_indices, name);
//            name = "search_skin" + std::to_string(high);
//            pModel->writeModelToFile(deadLimit, bead_indices, name);
        }

//        if (currentNumberOfComponents == 1){
//            break;
//        }

        cout << "*******************                                        *******************" << endl;
        printf("   MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        printf("      GRAPH => %3i\n", currentNumberOfComponents);
        printf("      UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("     VOLUME => %.0f  MU => %.4E  MU*VOL => %.6f\n", current_volume, mu, mu*current_volume);
        printf("LIMIT: %5i DEADLIMIT: %5i D_KL: %.4E ENRGY: %.4E \n", workingLimit, deadLimit, currentKL, current_energy);
        cout << "*******************                                        *******************" << endl;

        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }
//
        // test bead_indices and vector are in sync
//        for(int i=0; i<workingLimit; i++){
//            if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << endl;
//                exit(0);
//            }
//        }
//
//        for(std::set<int>::iterator bit = beads_in_use_tree.begin(); bit != beads_in_use_tree.end(); ++bit){
//            if(std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, *bit) == bead_indices.begin()+workingLimit){
//                cout << "!!!!! " << " NOT FOUND in BEAD_INDICES " << *bit << endl;
//                exit(0);
//            }
//        }
//
//        if(beads_in_use_tree.size() != workingLimit){
//            cout << "!!!!! " << " Incorrect size WORKINGLIMIT "  << endl;
//            exit(0);
//        }

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowest_energy = current_energy;
            lowestWorkingLimit = workingLimit;
            lowestDeadLimit = deadLimit;
            lowest_volume = current_volume;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            averaging[averagingIndex] = workingLimit;
            averagingIndex++;
            if (averagingIndex == capacity){
                averagingIndex =0;
            }
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }

    int sum_x_squared=0, sum_x=0;
    for(int i=0; i<capacity; i++){
       sum_x_squared += averaging[i]*averaging[i];
       sum_x += averaging[i];
    }

    double average_x = sum_x/(double)capacity;
    double stdev = sqrt(sum_x_squared/(double)capacity - average_x*average_x);

    cout << "AVERAGE " << average_x << " " << stdev << " CNTR " << counter << endl;

    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);

    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setStartingDeadLimit(lowestDeadLimit);

    //pModel->setBeadAverageAndStdev(pModel->getStartingWorkingLimit(), 0.17*pModel->getStartingWorkingLimit());
    pModel->setBeadAverageAndStdev(average_x, ceil(stdev));

    cout << "*******************                                        *******************" << endl;
    cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << endl;
    printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    printf("      TRK AVERAGE => %.2f (%i) \n", (sumIt/counter), (int)counter);
    cout << "*******************                                        *******************" << endl;
    cout << "*******************                 VOLUME                 *******************" << endl;
    printf("          INITIAL => %.0f \n", initialVolume);
    printf("           LOWEST => %.0f \n", lowest_volume);
    printf("            FINAL => %.0f\n", current_volume);
    cout << "*******************                                        *******************" << endl;

    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << endl;
    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << endl;
        cout << "INCREASE highTempRounds, g" << endl;
        return false;
    }

}


std::string Anneal::refineHomogenousBodyASACVX(Model *pModel, Data *pData, std::string outputname){
    cout << "########################<<<<<>>>>>>############################## " << endl;
    cout << "#                                                               # " << endl;
    cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << endl;
    cout << "#                                                               # " << endl;
    cout << "########################<<<<<>>>>>>############################## " << endl;
    //float oldN = pModel->getVolumeAverage();     // need to reset this for modeling via rounds
    //float oldStdev = pModel->getVolumeStdev();
    float runningAverage = pModel->getVolumeAverage();
    float runningVariance = pModel->getVolumeStdev();
    int lowerLimit = runningAverage - runningVariance;

    double lowTempStop = (double)highTempStartForCooling;
    int swap1;

    int alterMe;
    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));

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
    //pModel->centerLatticeModel(workingLimit, bead_indices);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // set deadLimit of the selected set
    int deadLimit;
    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    std::shuffle(bead_indices.begin() + workingLimit, bead_indices.begin() + deadLimit, gen); // randomize search space

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
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;
    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float tempAverageContacts;

    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;

    float step_limit = (updateCount < 10000) ? 27000.0 : (float)updateCount;
    int deadUpdate = std::ceil(updateCount*0.091);


    std::vector<float> acceptanceRateDuringRun((int)step_limit);
    std::vector<double> tempDuringRun((int)step_limit);
    std::vector<float> divergenceDuringRun((int)step_limit);
    std::vector<int> workingLimitDuringRun((int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float sum_x_squared=0, sum_x=0, divideBy=0, average_x, stdev, acceptRate = 0.5, inv500 = 1.0/500.0;

    int counter=1;
    double inv_kb_temp = 1.0/lowTempStop;
    int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    float this_energy, lowestKL = currentKL;
    //char addRemoveText[75];

    double runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, pModel);

    // slowly increase weight of total Contact energy over D_kl
    double base = 1.65125;
    double finalbase = 2.27452; // emphasize packing as temp gets lower

//    double base = 3.65125;
//    double finalbase = 1.4712;

    //double etaConstant = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) + 2 );
    double etaConstant = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit))  )*base;
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);

    double etaFactor = std::pow(finalbase/base, 1.0/floor(step_limit/(float)deadUpdate));
    //double baseFactor = etaConstant;
    double baseFactor = base;

    //double etaFactor = 1.0/std::pow(100, 1.0/(step_limit/(float)deadUpdate));

    float low_temp_limit = step_limit*0.91;
    //float half_limit = step_limit-1.2*deadUpdate;
    float half_limit = 1.27*deadUpdate;
    // want D_KL to settle < 10^-5 to 10^-6
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    std::normal_distribution<float> volumeGen(runningAverage, runningVariance);

    std::clock_t startTime;
    int attempts=0, failures=0;
    double newSum;
    double runtime;

    int numberOfCoolingTempSteps;

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
        beginIt = bead_indices.begin();
        endIt = bead_indices.end();
        std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased

//        if (numberOfCoolingTempSteps < half_limit && distribution(gen) < percentAddRemove  ) { //add or remove bead within working Set (exclude deadzone)
        if (distribution(gen) < percentAddRemove  ) { //add or remove bead within working Set (exclude deadzone)
            //std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               ADD?REMOVE               *******************" << endl;
            alterMe = (int) volumeGen(gen);
            // build a list of indices within defined region of convex hull
            //if (distribution(gen) < 0.5){ // ADD BEAD?
            //if (alterMe > workingLimit || workingLimit < lowerLimit){ // ADD BEAD?
                if (alterMe > workingLimit ){ // ADD BEAD?
                cout << "*******************                  ADD                   *******************" << endl;
                startTime = std::clock();
                double afterAdding;

                original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while(original == -1){
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

                itIndex = std::find(bead_indices.begin(), bead_indices.end(), original);
                // select first element from randomized active_set
                // check if available neighbor can be added
                newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);

                afterAdding = etaConstant*newSum/(double)(workingLimit+1);
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);
                addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = pData->calculateKLDivergence(binCount);

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                this_energy = testKL + lambda*connectivityPotential(1);

                beads_in_use_tree.insert(original);
                tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    //rePopulateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, &deadLimit, pModel, original);
//                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    //rePopulateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, &deadLimit, pModel, original);
//                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(original);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }


            } else { // REMOVE BEADS?
                cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
//                sprintf(addRemoveText, "     REMOVE => %i", 1);

                double afterRemoving;
                startTime = std::clock();
                int randomSpot = randomIndex(gen);
                original = bead_indices[randomSpot];

                tempNumberOfComponents = eulerTour.removeNode(original);

                if (tempNumberOfComponents == 1){
                    // grab from randomized active_indices list

                    newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                    afterRemoving = etaConstant*newSum/(double)(workingLimit-1);
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                    // still need to sort, swap changes the order
                    testKL = pData->calculateKLDivergence(binCount);
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
                        //removeFromdDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel, original);
                    } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                        currentKL = testKL;
                        current_energy = this_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalContactEnergy = afterRemoving;
                        runningContactsSum = newSum;
                        isUpdated = true;
                        //removeFromdDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel, original);
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
            bool isSwapped;
            attempts +=1;
            // shuffling should not change the location of the iterator
            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               POSITIONAL               *******************" << endl;
            cout << "*******************                                        *******************" << endl;
            int originalSwap2Value;
            startTime = std::clock();
            //int wldl = (numberOfCoolingTempSteps > low_temp_limit) ? (workingLimit + (int)((deadLimit-workingLimit)*0.37)) : (workingLimit + (int)((deadLimit-workingLimit)*0.13));
            int wldl = (numberOfCoolingTempSteps > low_temp_limit) ? (int)((diffContacts*workingLimit)*0.13) :  (int)(diffContacts*workingLimit)*0.0157;

            // randomly select an index to move
            int position = randomIndex(gen);
            swap1 = bead_indices[ position ];

            // test bead_indices and vector are in sync
//        for(int i=0; i<workingLimit; i++){
//            if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << endl;
//                exit(0);
//            }
//        }
//
//        for(std::set<int>::iterator bit = beads_in_use_tree.begin(); bit != beads_in_use_tree.end(); ++bit){
//            if(std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, *bit) == bead_indices.begin()+workingLimit){
//                cout << "!!!!! " << " NOT FOUND in BEAD_INDICES " << *bit << endl;
//                exit(0);
//            }
//        }
//            std::cout << " FINISHED CHECKING " << endl;

            if (eulerTour.removeNode(swap1) == 1){
                // prefer to move index if it has too few contacts
                // remove contribution of swap1
                double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);

                isSwapped = false; // swap
                itIndex = bead_indices.begin() + position;

//                for(int i=0; i<workingLimit; i++){
//                    if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                        cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << endl;
//                        exit(0);
//                    }
//                }


                // remove selected index from P(r)
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

                beads_in_use_tree.erase(swap1); // remove from in_use tree

                double newPotential, invWorkingLimit = etaConstant/(double)(workingLimit);
                // possible multi-thread the loop
                // break on first success
                for (int bb = 0; bb < wldl; bb++) { // go through and find a better position within unused beads

                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while (originalSwap2Value == -1 || originalSwap2Value == swap1){
                        originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }

                    // get available neighbor position
                    pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), originalSwap2Value);
                    // make the swap, sort and update P(r)
                    // randomly pick a point in beads_in_use
                    // look at neighborhood for point.
                 //   if (numberOfContactsFromSet(&beads_in_use_tree, pModel, originalSwap2Value) > 0){
                        newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
                        newPotential = newSum*invWorkingLimit;

                        beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                        std::iter_swap(itIndex, pSwap2);
                        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted
//                        sortMe(&bead_indices, itIndex, pSwap2, &workingLimit);

                        addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                        testKL = pData->calculateKLDivergence(workingBinCount);

                        this_energy = testKL + lambda * connectivityPotential(1);
                        tempTotalContactEnergy = newPotential-totalContactEnergy;

                        if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = 1;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
//                            sprintf(addRemoveText, "     SWAPPED => %i to %i ( %i )", swap1, originalSwap2Value, wldl);
                            break;
                        } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) >
                                   distribution(gen)) {

                            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                            currentKL = testKL;
                            current_energy = this_energy;
                            currentNumberOfComponents = 1;
                            totalContactEnergy = newPotential;
                            runningContactsSum = newSum;
                            isUpdated = true;
                            isSwapped = true;
                            eulerTour.addNode(originalSwap2Value, pModel);
//                            sprintf(addRemoveText, "  SA SWAPPED => %i to %i ( %i )", swap1, originalSwap2Value, wldl);
                            break;
                        }
                        // reverse changes; find index in sorted list and replace
                        // binCount is P(r) distribution without swap1 value
                        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                        std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                        itIndex = bead_indices.begin() + position;
                        //itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                        beads_in_use_tree.erase(originalSwap2Value);
                   // }
                }

                // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
                if (!isSwapped){ // itIndex already reverses on exit from loop
                    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
//                    sprintf(addRemoveText, "      FAILED => %i", swap1);
                    eulerTour.addNode(swap1, pModel);
                }
            } else {
//                sprintf(addRemoveText, "   EU FAILED => %i", swap1);
                eulerTour.addNode(swap1, pModel);
            }
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("       TEMP => %-.4E \n     ACCEPT => %.5f  FAILURES => %i\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop, acceptRate, failures, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
        printf("  UPDATECNT => %8d TIME => %.5f (SECONDS) \n", deadUpdate, runtime);
        printf("      GRAPH => %3d \n", currentNumberOfComponents);
//        printf(" LATTCE AVG => %.0f        STDEV => %.3f  %s\n", runningAverage, runningVariance, addRemoveText);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.4f BASE => %.2f\n", totalContactEnergy, etaConstant, runningContactsSum, baseFactor);
        printf("LIMIT: %5i DEADLIMIT : %5i D_KL : %.4E ENRGY : %.4E\n", workingLimit, deadLimit, currentKL, (current_energy+totalContactEnergy));

        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = etaConstant;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1.0);
            isUpdated = false;
            if (currentKL < lowestKL ){
                lowestKL = currentKL;
            }
            // I am only moving/removing one bead at a time, should this be updated every 10 steps?
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
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
            //Race Condition
            //runningAverage = 0.5*(average_x + runningAverage);
            runningAverage = (average_x + 2.0*runningAverage)*0.333333;
            runningVariance = stdev;
            //pModel->setBeadAverageAndStdev(0.5*(average_x + pModel->getVolumeAverage()), stdev);
            //recalculate eta, increase last term to decrease weight of eta
            // if too much eta, shape is determined by contact potential
            // need right balance between KL and contactPotential
            //eta = pow(10, ceil(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)) - 0.911);
            //runningContactsSum = calculateTotalContactSum( &beads_in_use_tree, pModel);
            //etaConstant *= etaFactor;
            //std::string name = "model_" + std::to_string(numberOfCoolingTempSteps)+"_" +  std::to_string(baseFactor);
            //pModel->writeModelToFile(workingLimit, bead_indices, name, numberOfCoolingTempSteps);

            baseFactor *= etaFactor;
            etaConstant = pow(10, floor(log10(currentKL) - log10(runningContactsSum /(double)workingLimit)))*baseFactor;

            totalContactEnergy = etaConstant*(runningContactsSum / (double)workingLimit);
            //std::normal_distribution<float> volumeGen(runningAverage, runningVariance);
        }
        counter++;
    } // end of steps
//
//    // At end of each temp, update a probability model for volume?  Use this to select
//    // string tempName = "rough";
//    // pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, tempName, this, pData);
//
//    pModel->writeModelToFile(workingLimit, bead_indices, "before_final", numberOfCoolingTempSteps);
//    cout << "------------------------------------------------------------------------------" << endl;
//    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
//    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
//    cout << "------------------------------------------------------------------------------" << endl;
//
//    tempAverageContacts=0.0;
//    for (int i=0; i<workingLimit; i++){
//        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
//        tempAverageContacts += temp;
//    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << " POSITIONAL SUCCESS/FAILURE " << ((double)failures/(double)attempts)*100 << endl;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    int totalLessThanAverage=0, numberContacts;
    for (int i=0; i<workingLimit; i++){
        numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        if (numberContacts < average_number_of_contacts ){
            totalLessThanAverage +=1;
        }
    }
//
//    int finalRounds = (totalLessThanAverage*std::log((double)totalLessThanAverage) + 0.5772156649*totalLessThanAverage + 0.5);
//    cout << " TOTAL LATTICE POINTS LESS THAN AVERAGE " << totalLessThanAverage << endl;
//    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);
//
//    int swappedCounter =0;
//    int round=0;
//    float average_number_of_contacts_before = tempAverageContacts/(float)workingLimit;
//    tempAverageContacts=0.0;
//    for (int i=0; i<workingLimit; i++){
//        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
//    }
//    average_number_of_contacts = tempAverageContacts/(float)workingLimit;
//
    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " FINAL LOW TEMP REFINEMENT" << endl;
//    printf("    SWAPPED => %i             \n", swappedCounter);
//    printf("     ROUNDS => %i (%i)            \n", round, finalRounds);
    cout << "              KL DIVERGENCE   ENERGY " << endl;
    //printf("     BEFORE =>    %.4E           %.4E  \n", beforeFinalKL, beforeFinalEnergy);
    printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    cout << "    AVERAGE" << endl;
    //printf("     BEFORE =>    %.2f  \n", average_number_of_contacts_before);
    printf("      AFTER =>    %.2f  \n", average_number_of_contacts);
    cout << "------------------------------------------------------------------------------" << endl;
//
//    //this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);
//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF
    int numpoints = workingLimit;
    coordT points[3*numpoints];
    char flags[]="qhull FA";
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(cvx);
    pData->printKLDivergence(binCount);

    //float average_number_of_contacts = 4.1;
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
    std::cout << " FINISHED TEMP " << std::endl;
    return nameOfModel;
}
