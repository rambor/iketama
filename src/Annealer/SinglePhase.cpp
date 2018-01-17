//
// Created by Robert Rambo on 13/01/2017.
//
#include <Keeper.h>
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

    maxbin += 1;  // maximum bin index, so vector size must be +1

    std::cout << "MAX BIN : " << maxbin << std::endl;

    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    //std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50


    std::cout << "    TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "    MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "              BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    //std::vector<int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<int> backUpState(totalBeadsInSphere);
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
    unsigned int lowerN = round(lowerV*invBeadVolume);
    unsigned int upperN = round(upperV*invBeadVolume);

    // connectivity potential
    // volume potential
    // pick random number between lowerN and upperN
    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN);
    unsigned int workingLimit = number_of_beads_to_use(gen);
    pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "universe", 0);
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::cout << " POPULATING DEAD LIMIT" << std::endl;
//    int deadLimit; // as bead indices are discarded, set upper limit of vector
//    populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    //populateLayeredDeadlimit(bead_indices.begin(), workingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
//    std::cout << "    DEADLIMIT SET => " << deadLimit << std::endl;
    std::cout << " WORKINGLIMIT SET => " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";
    //int numpoints = 3*totalBeadsInSphere;
    coordT points[3*(upperN+(int)(0.30*upperN))];

    float lowest_volume, test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float initialVolume = current_volume;
//
//    // area of the model P(r) distribution
    std::cout << "       CVX VOLUME => " << current_volume  << std::endl;
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);

    int tempNumberOfComponents, currentNumberOfComponents  = eulerTour.getNumberOfComponents();
    //int currentNumberOfComponents;
    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, currentNumberOfComponents);
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    std::vector<unsigned int>::iterator beginBinCount = binCount.begin();
    std::vector<int>::iterator itIndex;
    std::vector<unsigned int>::iterator endBinCount = binCount.end();

    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
//    float lowest_volume, test_volume, current_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);//calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    float initialVolume = current_volume;
    // area of the model P(r) distribution
    std::cout << "       CVX VOLUME => " << current_volume  << " WL => " << workingLimit << std::endl;
    float hlambda = 0.001;
    float muConstant = mu*currentKL/((1.0-mu)*current_volume);
    float temp_volume_energy, current_volume_energy = muConstant*current_volume;
    float testKL, test_energy, current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_volume_energy;//(float)workingLimit;

    std::cout << " STARTING KL => " << currentKL << std::endl;
    float lowest_energy = current_energy;

    int testIndex;
    //float inv_kb_temp = 1.0/(float)highT;
    float inv_kb_temp = 1.0/0.0001;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    unsigned int lowestWorkingLimit;
    lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    double sumIt = 0.0, counter=0.0;

    bool updateCVX =false;

    //int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();  seems on average there are less than 12 available to look at
    std::cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< endl;
    int high;
   // std::clock_t startTime;

    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
        std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased

        if (distribution(gen) < 0.371){

            //alterMe = number_of_beads_to_use(gen);  // uniform distribution

            if ((distribution(gen) < 0.5) && workingLimit < upperN){
                std::cout << "*******************                  ADD                   ******************* "  << endl;
                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while ( testIndex == -1) { // find a new position
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

                itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
                //make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);

                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->calculateKLDivergence(binCount);

                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                //startTime = std::clock();
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                //std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
                temp_volume_energy = muConstant*test_volume;
                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    current_volume_energy = temp_volume_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX = true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    current_volume_energy = temp_volume_energy;
                    currentNumberOfComponents = tempNumberOfComponents;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX = true;
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                }

            } else {

                std::cout << "*******************                 REMOVE                 ******************* "  << endl;

                if (workingLimit > lowerN){

                    // grab from randomized active_indices list
                    //int randomSpot = randomIndex(gen);
                    testIndex = bead_indices[randomIndex(gen)];

                    std::vector<int>::iterator beginIt = bead_indices.begin();
                    removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = pData->calculateKLDivergence(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                    //test_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);
                    //startTime = std::clock();
                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    //std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
                    //isConnectedComponent(&bead_indices, workingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                    temp_volume_energy = muConstant*test_volume;
                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                    if (test_energy < current_energy ) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        current_volume_energy = temp_volume_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else if (exp(-(test_energy - current_energy)*inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        current_volume_energy = temp_volume_energy;
                        currentNumberOfComponents = tempNumberOfComponents;
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        beginIt = bead_indices.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    }
                }
            }

        } else {

            std::cout << "*******************               POSITIONAL               *******************" << endl;

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
            for (unsigned int i = 0; i < workingLimit; i++) {
                beadToPoint(&points[i*3], pModel->getBead(bead_indices[i]));
                active_indices[i] = bead_indices[i];
            }

            // needs to be optimized
            qh_new_qhull(3, workingLimit, points, 0, flags, NULL, NULL);
            vertexT * vertices = qh vertex_list;
            int totalV = qh num_vertices;
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

            //std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            int neighbor;

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
                //int randind;

                for (bb = 0; bb < 1; bb++) { // go through and find a better position within unused beads
                    // find a bead in use, and see if it has an empty neighbor to use
                    //randind = randomIndex(gen);
                    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while (neighbor == -1 || neighbor == swap1){
                        //randind = randomIndex(gen);
                        neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }
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

                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

                    temp_volume_energy = muConstant*test_volume;
                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                    if (test_energy < current_energy) { //
                        beads_in_use_tree.insert(neighbor);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                        std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        current_volume_energy = temp_volume_energy;
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
                        current_volume_energy = temp_volume_energy;
                        current_energy = test_energy;
                        isSwapped = true;
                        updateCVX = true;
                        break;
                    }  // reverse changes; find index in sorted list and replace

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

        std::cout << "*******************                                        *******************" << endl;
        printf("           MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        printf("              GRAPH => %3i LIMIT => %5i\n", currentNumberOfComponents, workingLimit);
        printf("              UPPER => %i LOWER => %i \n", upperN, lowerN);
        printf("             VOLUME => %.0f  MU => %.4E \n", current_volume, muConstant);
        printf("               D_KL : %.4E E_VOL : %.4E ENRGY : %.4E \n", currentKL, current_volume_energy, current_energy);
        std::cout << "*******************                                        *******************" << endl;

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

            lowestWorkingLimit = workingLimit;
            lowest_volume = current_volume;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

            // update D_kl
            muConstant = mu*currentKL/((1.0-mu)*current_volume);
            //current_energy = currentKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + muConstant*current_volume;
            current_energy = currentKL + muConstant*current_volume;
            lowest_energy = current_energy;
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }


    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);
    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
//    pModel->setStartingDeadLimit(lowestDeadLimit);
    //pModel->setBeadAverageAndStdev(pModel->getStartingWorkingLimit(), 0.17*pModel->getStartingWorkingLimit());
    pModel->setBeadAverageAndStdev(lowestWorkingLimit, lowestWorkingLimit*0.2);

    std::cout << "*******************                                        *******************" << endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << endl;
    printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    printf("      TRK AVERAGE => %.2f (%i) \n", (sumIt/counter), (int)counter);
    std::cout << "*******************                                        *******************" << endl;
    std::cout << "*******************                 VOLUME                 *******************" << endl;
    printf("          INITIAL => %.0f \n", initialVolume);
    printf("           LOWEST => %.0f \n", lowest_volume);
    printf("            FINAL => %.0f\n", current_volume);
    std::cout << "*******************                                        *******************" << endl;


    std::cout << "     MAX MODEL N_S BINS : " << maxbin << endl;
    std::cout << "   SHANNON BINS IN DATA : " << pData->getShannonBins() << endl;
    std::cout << "            ZERO BIN AT : " << pData->getZeroBin() << endl;

    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << endl;
    float    tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << endl;
//    exit(0);

    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << endl;
        cout << "INCREASE highTempRounds, g" << endl;
        return false;
    }
}


std::string Anneal::refineHomogenousBodyASAHybrid(Model *pModel, Data *pData, std::string outputname){
    cout << "########################<<<<<>>>>>>############################## " << endl;
    cout << "#                                                               # " << endl;
    cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << endl;
    cout << "#                                                               # " << endl;
    cout << "########################<<<<<>>>>>>############################## " << endl;
    //float oldN = pModel->getVolumeAverage();     // need to reset this for modeling via rounds
    //float oldStdev = pModel->getVolumeStdev();

    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.0001 : (double)highTempStartForCooling;

    int swap1, totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

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

    unsigned int workingLimit = pModel->getStartingWorkingLimit(), deadLimit;

    // reset iterators to internal bead_indices
    //pModel->centerLatticeModel(workingLimit, bead_indices);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    this->fillPrBinsAndAssignTotalBin(pBin,  pDistance,  totalDistancesInSphere,  pData);
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

    std::cout << "    TOTAL EXP N_S BINS: " << totalBins << endl;
    std::cout << "    MAX MODEL N_S BINS: " << maxbin << endl;
    std::cout << "              BINWIDTH: " << pData->getBinWidth() << endl; // lattice should be limited by binwidth

    auto beginBinCount = binCount.begin();
    auto endBinCount = binCount.end();

    float testKL, currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    float lowestKL = currentKL;
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

    //
    deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
    enlargeDeadLimit(bead_indices, &deadLimit, pModel);
    //pModel->writeModelToFile(deadLimit, bead_indices, "dead", 0);

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;
    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);

    int updateCount = ccmultiple*coupons;
    int updatedCount = 0;

    //float step_limit = (updateCount < 10000) ? 27000.0 : (float)updateCount;
    float step_limit = updateCount/0.1;

   // int deadUpdate = std::ceil(step_limit/updateSteps);
    // if renormalization is too frequent, then a arrangemetn that may be too out of norm will be reset to low
    // so, need to span few steps?
    //
    // if 0.2, essentially divides step_limit into 5 steps for recalibrating weight
    // if 0.1 => 10 steps
    // if 0.05 => 20 steps
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
    float acceptRate = 0.5, inv500 = 1.0/500.0;
    float inv500slash499 = 499.0/500.0;

    int counter=1;
    double inv_kb_temp = 1.0/lowTempStop;
    //int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    float this_energy;
    char addRemoveText[75];

    double runningContactsSum = calculateTotalContactSumPotential(&beads_in_use_tree, pModel);

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    double tempKLDivContacts, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    // slowly increase weight of total Contact energy over D_kl
    double percentDiff=0;

    //float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + currentVolE;
    double cDWConstant = 0;
    double cotactDistributionWeight = currentKL*cDWConstant/(currentKLDivContacts - cDWConstant*currentKLDivContacts);
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);

    //double etaConstant = (eta == 0) ? 0 : eta;
    double invWL = 1.0/(double)workingLimit;

    double etaConstant = current_energy*eta/((1.0 - eta)*runningContactsSum*invWL);
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*runningContactsSum*invWL;
    float lowest_energy = current_energy + totalContactEnergy;

    std::clock_t startTime;
    int attempts=0, failures=0, animateCount=0, updated=0, successess=0;
    //int upperN = round(upperV/pModel->getBeadVolume());
    int upperN = workingLimit;

    double newSum;
    double runtime;

    std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<int> randomDeadIndex (workingLimit, deadLimit-1);

    std::vector<int> selections(totalBeadsInSphere);
    for(unsigned long int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    int numberOfCoolingTempSteps=0;
    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());

        if (distribution(gen) < percentAddRemove  ) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << endl;
            // build a list of indices within defined region of convex hull
            if (distribution(gen) < 0.50 ){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << endl;
                startTime = std::clock();
                double afterAdding;

                if (distribution(gen) > acceptRate ){
                    // get lattice point not in use but touching a position in use, pick random bead to add to
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while(original == -1){
                        original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }
                } else {
                    original = bead_indices[randomDeadIndex(gen)];
                }

                itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

                // check if available neighbor can be added
                newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, original, runningContactsSum);
                afterAdding = etaConstant*newSum/((double)workingLimit+1.0);

//                addToContactsDistribution(original, contactsDistributionOfModel, &beads_in_use_tree, pModel);
//                tempKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);  // alters backUpState
                addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = pData->calculateKLDivergence(binCount);

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                //this_energy = testKL + lambda*connectivityPotential(1);
                tempNumberOfComponents = eulerTour.addNode(original, pModel);
                this_energy = testKL +  lambda*connectivityPotential(tempNumberOfComponents);
//                this_energy = testKL +  contactDistributionWeight*tempKLDivContacts + lambda*connectivityPotential(tempNumberOfComponents);
                tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                    beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
//                    currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    isUpdated = true;
//                    eulerTour.addNode(original, pModel);
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                    beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
//                    currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    isUpdated = true;
//                    eulerTour.addNode(original, pModel);
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    eulerTour.removeNode(original);
//                    beads_in_use_tree.insert(original);
//                    removeFromContactsDistribution(original, contactsDistributionOfModel, &beads_in_use_tree, pModel);
//                    beads_in_use_tree.erase(original);
                    sprintf(addRemoveText, "     ADD => %i", 0);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
                sprintf(addRemoveText, "     REMOVE => %i", 1);
                startTime = std::clock();
                original = bead_indices[randomIndex(gen)];

//                bool tourtest = true;
//                while(tourtest){
//                    if (eulerTour.removeNode(original) == 1){
//                        tourtest = false;
//                    } else {
//                        eulerTour.addNode(original, pModel);
//                        original = bead_indices[ randomIndex(gen) ]; // potential to select the same thing twice
//                    }
//                }
// why would remove be favored over add?
                newSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, original, runningContactsSum);
                double afterRemoving = etaConstant*newSum/((double)workingLimit-1.0);

//                removeFromContactsDistribution(original, contactsDistributionOfModel, &beads_in_use_tree, pModel);
//                tempKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

                beginIt = bead_indices.begin();
                removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                // still need to sort, swap changes the order
                testKL = pData->calculateKLDivergence(binCount);

                tempNumberOfComponents = eulerTour.removeNode(original);
                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);
//                this_energy = testKL + contactDistributionWeight*tempKLDivContacts + lambda*connectivityPotential(tempNumberOfComponents);

                tempTotalContactEnergy = (afterRemoving - totalContactEnergy);

                if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
                    beads_in_use_tree.erase(original);
                    currentKL = testKL;
                    current_energy = this_energy;
//                    currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    isUpdated = true;

                } else if ((testKL > 0 ) && exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                    beads_in_use_tree.erase(original);
                    currentKL = testKL;
                    current_energy = this_energy;
//                    currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    isUpdated = true;

                } else { // undo changes and move to next bead (rejecting)
//                    beads_in_use_tree.erase(original);
//                    addToContactsDistribution(original, contactsDistributionOfModel, &beads_in_use_tree, pModel);
//                    beads_in_use_tree.insert(original);
                    beginIt = bead_indices.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);

                    eulerTour.addNode(original, pModel);
                }
            }

//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            startTime = std::clock();

            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            // randomly select an index to move
            int position = selections[0];
            swap1 = bead_indices[position];
            for(unsigned int i=1;i<workingLimit; i++){
                if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                    break;
                } else { // if last node, will break out
                    eulerTour.addNode(swap1, pModel);
                    position = selections[i];
                    swap1 = bead_indices[position]; // what happens if I fail all the way?
                }
            }

            //std::cout << "RUN TIME 2 " << ((std::clock() - startTime)/(double) CLOCKS_PER_SEC) << std::endl;
            // remove contribution of swap1
            double oldSum = recalculateContactsPotentialSumRemove(&beads_in_use_tree, pModel, swap1, runningContactsSum);
            itIndex = bead_indices.begin() + position;

            // remove selected index from P(r)
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

            // remove contribution to Contacts Distribution
//            for(int c=0; c<13; c++){ // make BackUp Copy
//                tempContactsDistributionOfModel[c] = contactsDistributionOfModel[c];
//            }
//            removeFromContactsDistribution(swap1, contactsDistributionOfModel, &beads_in_use_tree, pModel);

            beads_in_use_tree.erase(swap1); // remove from in_use tree
            // possible multi-thread the loop
            // break on first success
            int originalSwap2Value;
            for(unsigned int i=0;i<workingLimit; i++){
                int newPosition = bead_indices[selections[i]];
                if (newPosition != swap1){
                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, newPosition);
                    if (originalSwap2Value != -1 && originalSwap2Value != swap1){
                        break;
                    }
                }
            }

            // make the swap, sort and update P(r)
            // randomly pick a point in beads_in_use
            // look at neighborhood for point.
            newSum = recalculateContactsPotentialSumAdd(&beads_in_use_tree, pModel, originalSwap2Value, oldSum);
            double newPotential = etaConstant*newSum*invWL;

//            addToContactsDistribution(originalSwap2Value, contactsDistributionOfModel, &beads_in_use_tree, pModel);
//            tempKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

            // get available neighbor position
            pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->calculateKLDivergence(workingBinCount);

            if (currentNumberOfComponents == 1){

                this_energy = testKL;
                tempTotalContactEnergy = newPotential-totalContactEnergy;

                if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    totalContactEnergy = newPotential;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
                    sprintf(addRemoveText, "     SWAPPED => %i to %i ", swap1, originalSwap2Value);
                } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) >
                           distribution(gen)) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    totalContactEnergy = newPotential;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    eulerTour.addNode(originalSwap2Value, pModel);
                    sprintf(addRemoveText, "  SA SWAPPED => %i to %i", swap1, originalSwap2Value);
                } else {
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                    sprintf(addRemoveText, "      FAILED => %i", swap1);
                    eulerTour.addNode(swap1, pModel);
                }
            } else {
                tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);

//            this_energy = testKL + contactDistributionWeight*tempKLDivContacts + lambda*connectivityPotential(tempNumberOfComponents);
                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);
                tempTotalContactEnergy = newPotential-totalContactEnergy;

                if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                    currentKL = testKL;
                    current_energy = this_energy;
//                currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = newPotential;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    //eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
                    sprintf(addRemoveText, "     SWAPPED => %i to %i ", swap1, originalSwap2Value);
                } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) >
                           distribution(gen)) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                    currentKL = testKL;
                    current_energy = this_energy;
//                currentKLDivContacts = tempKLDivContacts;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalContactEnergy = newPotential;
                    runningContactsSum = newSum;
                    isUpdated = true;
                    //eulerTour.addNode(originalSwap2Value, pModel);
                    sprintf(addRemoveText, "  SA SWAPPED => %i to %i", swap1, originalSwap2Value);
                } else {
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    beads_in_use_tree.insert(swap1); // add back swap one to tree
                    sprintf(addRemoveText, "      FAILED => %i", swap1);
                    eulerTour.removeNode(originalSwap2Value);
                    eulerTour.addNode(swap1, pModel);

//                for(int c=0; c<13; c++){ // restore from BackUp Copy
//                    contactsDistributionOfModel[c] = tempContactsDistributionOfModel[c];
//                }
                }
            }



        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("      TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        printf("      TIME : %.5f (SECONDS)  %s\n", runtime, addRemoveText);
        printf("     LIMIT : %5i  DEADLIMIT => %5i     COMP => %3i\n", workingLimit, upperN, currentNumberOfComponents);
        printf(" CONTACTE => %-5.4E ETA => %.4E AVG => %.2f N => %i\n", totalContactEnergy, etaConstant, runningContactsSum, updatedCount);
        printf("     D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, lowestKL, (current_energy+totalContactEnergy));

        // random access
        pTempDuringRun[numberOfCoolingTempSteps] = runningContactsSum;
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = (current_energy + totalContactEnergy);
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            invWL = 1.0/(double)workingLimit;
            updated++;

            if (updated % coupons == 0 && numberOfCoolingTempSteps > 1){
                deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
                enlargeDeadLimit(bead_indices, &deadLimit, pModel);
                //pModel->writeModelToFile(deadLimit, bead_indices, "dead", updated);
                //pModel->writeModelToFile(workingLimit, bead_indices, "cooling", updated);
            }

            randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
            randomDeadIndex = std::uniform_int_distribution<int>(workingLimit, deadLimit-1);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            //selections.resize(workingLimit);
            int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned long int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

//        populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
//        // check distribution
//        for(int c=0; c<13; c++){
//            if (tempContactsDistributionOfModel[c] != contactsDistributionOfModel[c]){
//                std::cout << c << " contacts " << tempContactsDistributionOfModel[c] << " != " << contactsDistributionOfModel[c] << std::endl;
//                exit(0);
//            }
//        }
//        checkSetAndVector(workingLimit, &bead_indices, &beads_in_use_tree);

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        //updateASAModTemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        // update for running average
        // correlate the updating based on D_kl, which means only update if we get lower values
        // if we correlated with overall energy, then improvements in volume or contacts will cause constants to be re-updated
        if ( (current_energy + totalContactEnergy) < lowest_energy){ // if counter too small, add/remove may not sample sufficiently
            // pModel->setBeadAverageAndStdev(0.5*(average_x + pModel->getVolumeAverage()), stdev);
            // renormalization
            // percentDiff = (lowest_energy - (current_energy + totalContactEnergy))/(0.5*((current_energy + totalContactEnergy)+lowest_energy));
//            percentDiff = (contactDistributionWeight*currentKLDivContacts)/current_energy;
//            if ( percentDiff > cDWConstant  ) {
//                contactDistributionWeight = currentKL*cDWConstant/(currentKLDivContacts - cDWConstant*currentKLDivContacts);
//                current_energy = currentKL + contactDistributionWeight*currentKLDivContacts + lambda*connectivityPotential(currentNumberOfComponents);;
//                diffCount++;
//            }
//
            // if outside window, renormalize
           // double test = totalContactEnergy/current_energy - eta;
           // if ( std::abs(test) > 0.05 && acceptRate > 0.31) {
           if ( totalContactEnergy/current_energy < eta) {
                //invWL = 1.0/(double)workingLimit;
                etaConstant = current_energy*eta/((1.0 - eta)*runningContactsSum*invWL);
                totalContactEnergy = etaConstant*runningContactsSum*invWL;
               updatedCount++;
            }

            lowestKL = currentKL;
            lowest_energy = current_energy + totalContactEnergy;
        }

        counter++;
    } // end of steps
//
    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    std::string name = "ar_" + outputname;
//    pModel->writeModelToFile2(
//            currentKL,
//            workingLimit,
//            bead_indices,
//            binCount,
//            name,
//            this,
//            pData,
//            numberOfCoolingTempSteps,
//            10000,
//            0);

    double invWorkingLimit = 1.0/(double)workingLimit;

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);
    std::cout << " currentKLDivContacts " << currentKLDivContacts << std::endl;

    int totalCounts=0;
    for (int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (int c=1; c<13; c++){
        int totalContactsAt = 0;
        for (int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::cout << c << " " << totalContactsAt/(double)totalCounts << std::endl;
    }

    std::cout << "KL DIVERGENCE " << std::endl;
    pData->printKLDivergence(binCount);
//    contactDistributionWeight = 0.0001;
//    current_energy = currentKL + contactDistributionWeight*currentKLDivContacts;

    // relax the connectivity constraint ?
//    if ((current_energy + totalContactEnergy) < lowest_energy){
//        etaConstant *= 0.1;
//        etaConstant = current_energy*eta/(runningContactsSum - eta*runningContactsSum);
//        totalContactEnergy = etaConstant*runningContactsSum;
//        lowestKL = currentKL;
//    }

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(current_volume);
    // reset temp for low temp annealing (strictly positional)
    float tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    cout << " POSITIONAL SUCCESS/FAILURE " << ((double)failures/(double)attempts)*100 << endl;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    name = "annealed_" + outputname;
    pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            name,
            this,
            pData,
            numberOfCoolingTempSteps,
            10000,
            0);

    // CONSTANT TEMP REFINEMENT
    // make Contact Potential a percentage of final energy
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




    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << " FINAL LOW TEMP REFINEMENT" << std::endl;
    std::cout << "              KL DIVERGENCE   ENERGY " << std::endl;
    std::printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    std::cout << "    CONTACTS" << std::endl;
    std::printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    std::printf("    BIN [ONE] => %d AVG => %.3f \n", binCount[0], binCount[0]/(float)workingLimit);
    std::printf("    VOLUME  START => %.4E FINAL => %.4E \n", starting_volume, current_volume);
    std::printf("E_CONTACTS  FINAL => %.4E \n", runningContactsSum/(float)workingLimit);
    std::cout << "------------------------------------------------------------------------------" << std::endl;
//
   this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);

//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF


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
            current_volume,
            average_number_of_contacts);

    return nameOfModel;
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
    int swap1, totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

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

    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    // reset iterators to internal bead_indices
    //pModel->centerLatticeModel(workingLimit, bead_indices);
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // set deadLimit of the selected set
    unsigned int deadLimit;
    //populateLayeredDeadlimitUsingSet(&bead_indices, &beads_in_use_tree, workingLimit, &deadLimit, pModel);
    std::shuffle(bead_indices.begin() + workingLimit, bead_indices.begin() + deadLimit, gen); // randomize search space

    // convert distances in Search Sphere to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    this->fillPrBinsAndAssignTotalBin(pBin,  pDistance,  totalDistancesInSphere,  pData);
    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
//    std::vector<int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS: " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS: " << maxbin << endl;
    cout << "              BINWIDTH: " << pData->getBinWidth() << endl; // lattice should be limited by binwidth

    std::vector<unsigned int>::iterator beginBinCount = binCount.begin();
    std::vector<unsigned int>::iterator endBinCount = binCount.end();

    float testKL, currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    char flags[] = "qhull FA";
    float tempVolume=0, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    float lowestKL = currentKL;
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << endl;
    int original;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float updateSteps = 7.0;
    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;

    float step_limit = (updateCount < 10000) ? 27000.0 : (float)updateCount;
    int deadUpdate = std::ceil(step_limit/updateSteps);
    // if renormalization is too frequent, then a arrangemetn that may be too out of norm will be reset to low
    // so, need to span few steps?
    //
    // if 0.2, essentially divides step_limit into 5 steps for recalibrating weight
    // if 0.1 => 10 steps
    // if 0.05 => 20 steps

    std::vector<float> acceptanceRateDuringRun((int)step_limit);
    std::vector<double> tempDuringRun((int)step_limit);
    std::vector<float> divergenceDuringRun((int)step_limit);
    std::vector<int> workingLimitDuringRun((int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float sum_x_squared=0, sum_x=0, divideBy=0, average_x, stdev, acceptRate = 0.5, inv500 = 1.0/500.0;
    float inv500slash499 = 499.0/500.0;

    int counter=1;
    double inv_kb_temp = 1.0/lowTempStop;
    //int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    float this_energy;
    char addRemoveText[75];

    // slowly increase weight of total Contact energy over D_kl
    // maximize binCount[0] ?
    // binCount[0] total number of pairwise contacts
    // maximize 2*binCount[0]/workingLimit => maximizes connectivity
    // minimize Dkl while keeping contacts within a specified value
    // constant eta is acheived by setting:
    // baseFactor => eta
    // etaFactor = 1;
    //double baseFactor = 0.1;
    //double finalMu = 0.2714;
    //double muFactor = std::pow(finalMu/mu, 1.0/updateSteps);
    double muFactor = 1;
    double muAt = mu;

    // exit(0);
    // want D_KL to settle < 10^-5 to 10^-6
    float muConstant = muAt*currentKL/((1.0-muAt)*current_volume);
    float tempVolE=0, currentVolE = muConstant*current_volume;
    //float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + currentVolE;
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + currentVolE;
    float lowest_energy = current_energy;

    std::clock_t startTime;
    int attempts=0, failures=0, animateCount=0, updated=0;
    double runtime;

    int numberOfCoolingTempSteps;
    std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
        beginIt = bead_indices.begin();
        endIt = bead_indices.end();

        if (distribution(gen) < percentAddRemove  ) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            cout << "______________________________________________________________________________" << endl;
            cout << "*******************               ADD?REMOVE               *******************" << endl;
            // build a list of indices within defined region of convex hull
            if (distribution(gen) < 0.5){ // ADD BEAD?
                cout << "*******************                  ADD                   *******************" << endl;
                startTime = std::clock();
                double afterAdding;
                // get lattice point not in use but touching a position in use
                // pick random bead to add to
                original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while(original == -1){
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }
                itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);
                //std::vector<int>::iterator lit = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);
                // select first element from randomized active_set
                // check if available neighbor can be added
                //afterAdding = etaConstant*averageContactsPotential(newSum/(double)(workingLimit+1));
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);
                addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = pData->calculateKLDivergence(binCount);

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                //this_energy = testKL + lambda*connectivityPotential(1);
                //tempNumberOfComponents = eulerTour.addNode(original, pModel);
                tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                tempVolE = muConstant * tempVolume;///(float)workingLimit;
                this_energy = testKL + tempVolE;// + lambda*connectivityPotential(tempNumberOfComponents);

                if ( this_energy < current_energy ) {
                    beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents = 1; // tempNumberOfComponents;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp(-(this_energy - current_energy) * inv_kb_temp) > distribution(gen) ) {
                    beads_in_use_tree.insert(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents = 1; // tempNumberOfComponents;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                cout << "*******************                 REMOVE                 *******************" << endl;
                // test for deletion
                sprintf(addRemoveText, "     REMOVE => %i", 1);
                startTime = std::clock();
                original = bead_indices[randomIndex(gen)];

                bool tourtest = true;
                while(tourtest){
                    if (eulerTour.removeNode(original) == 1){
                        tourtest = false;
                    } else {
                        eulerTour.addNode(original, pModel);
                        original = bead_indices[ randomIndex(gen) ]; // potential to select the same thing twice
                    }
                }

                // double afterRemoving = etaConstant*averageContactsPotential(newSum/(double)(workingLimit-1));
                removeLatticePositionToModel(&beginIt, bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);
                // still need to sort, swap changes the order
                testKL = pData->calculateKLDivergence(binCount);
                tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                tempVolE = muConstant * tempVolume;
                this_energy = testKL + tempVolE; // + lambda*connectivityPotential(tempNumberOfComponents);
                //tempNumberOfComponents = eulerTour.removeNode(original);
                //this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);

                if (this_energy < current_energy ) {
                    beads_in_use_tree.erase(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents = 1; //  tempNumberOfComponents;
                    isUpdated = true;
                    randomIndex = std::uniform_int_distribution<int> (0,workingLimit-1); // guaranteed unbiased
                } else if ((testKL > 0 ) && exp(-(this_energy - current_energy)*inv_kb_temp) > distribution(gen) ){
                    beads_in_use_tree.erase(original);
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents = 1; // tempNumberOfComponents;
                    isUpdated = true;
                    randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
                } else { // undo changes and move to next bead (rejecting)
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.addNode(original, pModel);
                }
            }

//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }

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
            int wldl = 1;

            // randomly select an index to move
            int position = randomIndex(gen);
            swap1 = bead_indices[ position ];

            // select only node I can move?
            bool tourtest = true;
            while(tourtest){
                if (eulerTour.removeNode(swap1) == 1){
                    tourtest = false;
                } else {
                    eulerTour.addNode(swap1, pModel);
                    position = randomIndex(gen);
                    swap1 = bead_indices[ position ];
                }
            }

            // prefer to move index if it has too few contacts
            // remove contribution of swap1
            isSwapped = false; // swap
            itIndex = bead_indices.begin() + position;

            // remove selected index from P(r)
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

            beads_in_use_tree.erase(swap1); // remove from in_use tree

            //double newPotential, invWorkingLimit = 1.0/(double)(workingLimit);
            // possible multi-thread the loop
            // break on first success
            for (int bb = 0; bb < wldl; bb++) { // go through and find a better position within unused beads
                int tempRandom = bead_indices[randomIndex(gen)];
                originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, tempRandom);
//                    originalSwap2Value = getUseableNeighborFromSetDepth(&beads_in_use_tree, pModel, swap1);

                while (originalSwap2Value == -1 || tempRandom == swap1 || originalSwap2Value == swap1){
//                    while (originalSwap2Value == -1 || originalSwap2Value == swap1){
                    tempRandom = bead_indices[randomIndex(gen)];
                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, tempRandom);
//                        originalSwap2Value = getUseableNeighborFromSetDepth(&beads_in_use_tree, pModel, swap1);
                }

                // get available neighbor position
                pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);

                // make the swap, sort and update P(r)
                // randomly pick a point in beads_in_use
                // look at neighborhood for point.
                std::iter_swap(itIndex, pSwap2);
                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted
//                        sortMe(&bead_indices, itIndex, pSwap2, &workingLimit);
                addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
                // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                testKL = pData->calculateKLDivergence(workingBinCount);
                tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                tempVolE = muConstant*tempVolume;
                this_energy = testKL + tempVolE; // + lambda*connectivityPotential(tempNumberOfComponents);
//                    tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
//                    this_energy = testKL;// if connectivity is one, potential is 0

                if (this_energy < current_energy) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents = 1;// tempNumberOfComponents;
                    isUpdated = true;
                    isSwapped = true;
                    eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
                    sprintf(addRemoveText, "     SWAPPED => %i to %i ( %i )", swap1, originalSwap2Value, wldl);
                    break;
                } else if (exp(-(this_energy - current_energy) * inv_kb_temp) >
                           distribution(gen)) {
                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                    std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
                    currentKL = testKL;
                    current_energy = this_energy;
                    current_volume = tempVolume;
                    currentVolE = tempVolE;
                    currentNumberOfComponents =  1;// tempNumberOfComponents;
                    isUpdated = true;
                    isSwapped = true;
                    eulerTour.addNode(originalSwap2Value, pModel);
                    sprintf(addRemoveText, "  SA SWAPPED => %i to %i ( %i )", swap1, originalSwap2Value, wldl);
                    break;
                }
                // reverse changes; find index in sorted list and replace
                // binCount is P(r) distribution without swap1 value
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                itIndex = bead_indices.begin() + position;
                //eulerTour.removeNode(originalSwap2Value);
                //itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                //beads_in_use_tree.erase(originalSwap2Value);
            }

            // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
            if (!isSwapped){ // itIndex already reverses on exit from loop
                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // may not be necessary
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.addNode(swap1, pModel);
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("      TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        printf("      TIME : %.5f (SECONDS)  %s\n", runtime, addRemoveText);
        printf("     LIMIT : %5i  LATTCE AVG => %.0f     STDEV => %.3f\n", workingLimit, runningAverage, runningVariance);
        printf("     VOLE => %-5.4E ( %.0f ) AVG => %.3f \n", currentVolE, current_volume, binCount[0]/(float)workingLimit);
        printf("     D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, lowestKL, current_energy);

        // random access
        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        //pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = totalContactEnergy;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            // I am only moving/removing one bead at a time, should this be updated every 10 steps?
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            failures=0;

            sum_x_squared += workingLimit*workingLimit;
            sum_x += workingLimit;
            divideBy += 1.0;
            updated++;
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

//        checkSetAndVector(workingLimit, &bead_indices, &beads_in_use_tree);
        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        // update for running average
        // correlate the updating based on D_kl, which means only update if we get lower values
        // if we correlated with overall energy, then improvements in volume or contacts will cause constants to be re-updated
        if (acceptRate > 0.40 && currentKL < lowestKL ){ // if counter too small, add/remove may not sample sufficiently
            float inv_divideBy = 1.0/divideBy;
            average_x = sum_x*inv_divideBy;

            float sum_x_square_inv_divideBy = sum_x_squared*inv_divideBy;
            stdev = sqrt(sum_x_square_inv_divideBy - average_x*average_x);

            sum_x_squared = 0.0;
            sum_x = 0.0;
            divideBy= 0.0;

            runningAverage = average_x;
            runningVariance = stdev;
            // pModel->setBeadAverageAndStdev(0.5*(average_x + pModel->getVolumeAverage()), stdev);
            // recalculate eta, increase last term to decrease weight of eta
            // renormalization
            lowestKL = currentKL;
            // reset eta and mu based on lowestKL
            current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
            muAt *= muFactor;
            muConstant = muAt *currentKL/((1.0-muAt)*current_volume);
            current_energy = currentKL;// + muConstant * current_volume;
            lowest_energy = current_energy;

            pModel->writeModelToFile2(
                    currentKL,
                    workingLimit,
                    bead_indices,
                    binCount,
                    "cooling_" + to_string(counter),
                    this,
                    pData,
                    numberOfCoolingTempSteps,
                    10000,
                    0);
        }
        counter++;
    } // end of steps
//

    pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            "before",
            this,
            pData,
            numberOfCoolingTempSteps,
            10000,
            0);

    double invWorkingLimit = 1.0/(double)workingLimit;
    // keep the particle volume near the current_volume?
    double old_volume = current_volume;
    muConstant = mu*currentKL/((1.0-mu) * current_volume);// * invWorkingLimit);  // volume per bead
    current_energy = currentKL ;//+ muConstant * current_volume;// * invWorkingLimit;

    randomIndex = std::uniform_int_distribution<int>(0,workingLimit-1); // guaranteed unbiased
    acceptRate = 0.5;
    updated = 0;
    isUpdated=false;

    // reset temp for low temp annealing (strictly positional)
    lowTempStop *= 100;
    //lowTempStop = (double)highTempStartForCooling;
    inv_kb_temp = 1.0/lowTempStop;
    // coupon collector problem
    step_limit= 7*(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
    beginBinCount = binCount.begin();

    std::vector<int> selections(workingLimit);
    for(int i=0;i<workingLimit; i++){
        selections[i] = i;
    }


    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){
        // shuffling should not change the location of the iterator
        std::shuffle(selections.begin(), selections.end(), gen);
        cout << "______________________________________________________________________________" << endl;
        cout << "*******************                FIXED WL                *******************" << endl;
        cout << "*******************                                        *******************" << endl;
        startTime = std::clock();
        // randomly select an index to move
        int position = selections[0];
        swap1 = bead_indices[ position ];
        for(int i=1;i<workingLimit; i++){
            if (eulerTour.removeNode(swap1) > 1){
                eulerTour.addNode(swap1, pModel);
                position = selections[i];
                swap1 = bead_indices[position];
            } else {
                break;
            }
        }

        //runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        // remove contribution of swap1
        itIndex = bead_indices.begin() + position;
        // remove selected index from P(r)
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
        removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

        beads_in_use_tree.erase(swap1); // remove from in_use tree

        int originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, swap1);
        while (originalSwap2Value == -1 || originalSwap2Value == swap1){
            originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, swap1);
        }
        // make the swap, sort and update P(r)
        pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), originalSwap2Value);
        std::iter_swap(itIndex, pSwap2);
        std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

        addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);
        // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
        testKL = pData->calculateKLDivergence(workingBinCount);
        tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
        this_energy = testKL + muConstant * tempVolume; //*invWorkingLimit;// + lambda*connectivityPotential(tempNumberOfComponents);

        if (this_energy < current_energy) {
            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
            currentKL = testKL;
            current_energy = this_energy;
            current_volume = tempVolume;
            isUpdated = true;
            eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
            sprintf(addRemoveText, "     SWAPPED => %i to %i ", swap1, originalSwap2Value);
        } else if (exp(-(this_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount);
            currentKL = testKL;
            current_energy = this_energy;
            current_volume = tempVolume;
            isUpdated = true;
            eulerTour.addNode(originalSwap2Value, pModel);
            sprintf(addRemoveText, "  SA SWAPPED => %i to %i ", swap1, originalSwap2Value);
        } else {
            std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
            std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
            beads_in_use_tree.insert(swap1); // add back swap one to tree
            sprintf(addRemoveText, "      FAILED => %i", swap1);
            eulerTour.addNode(swap1, pModel);
        }

        //runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("       TEMP => %-.4E \n     ACCEPT => %.5f  FAILURES => %i\n      INVKB => %.3E\n   MAXSTEPS => %.0f (%4i) \n", lowTempStop, acceptRate, failures, inv_kb_temp, step_limit, numberOfCoolingTempSteps);
        printf("    UPDATED => %8d TIME => %.5f (SECONDS)\n", updated, runtime);
        printf("     VOLUME => %.0f   %s\n", current_volume, addRemoveText);
        printf("  BIN [ONE] => %6d AVG => %.3f \n", binCount[0], binCount[0]*invWorkingLimit);
        printf("LIMIT: %5i D_KL : %.4E ( %.4E ) ENRGY : %.4E\n", workingLimit, currentKL, lowestKL, current_energy);

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            failures=0;
            updated+=1;
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
    }


    pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            "low_temp_annealed",
            this,
            pData,
            numberOfCoolingTempSteps,
            10000,
            0);


    int lowTempSuccessCount=0;
    int finalCoupons = 7*(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    for(int round=0; round < finalCoupons; round++){

        bool isSwapped=false;
        attempts +=1;
        // shuffling should not change the location of the iterator
        cout << "______________________________________________________________________________" << endl;
        cout << "*******************       CONSTANT TEMP POSITIONAL         *******************" << endl;
        cout << "*******************                                        *******************" << endl;
        // randomly select an index to move
        int originalSwap2Value, position = randomIndex(gen);
        swap1 = bead_indices[ position ];

        // test bead_indices and vector are in sync
        // select only node I can move?
        bool tourtest = true;
        while(tourtest){
            if (eulerTour.removeNode(swap1) == 1){
                tourtest = false;
            } else {
                eulerTour.addNode(swap1, pModel);
                position = randomIndex(gen);
                swap1 = bead_indices[ position ];
            }
        }

        // remove contribution of swap1
        itIndex = bead_indices.begin() + position;

        // remove selected index from P(r)
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
        removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
        std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

        beads_in_use_tree.erase(swap1); // remove from in_use tree

        // possible multi-thread the loop
        // break on first success
        for (int bb = 0; bb < 1; bb++) { // go through and find a better position within unused beads
            int tempRandom = bead_indices[randomIndex(gen)];
            originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, tempRandom);
//                    originalSwap2Value = getUseableNeighborFromSetDepth(&beads_in_use_tree, pModel, swap1);

            while (originalSwap2Value == -1 || tempRandom == swap1 || originalSwap2Value == swap1){
//                    while (originalSwap2Value == -1 || originalSwap2Value == swap1){
                tempRandom = bead_indices[randomIndex(gen)];
                originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, tempRandom);
//                        originalSwap2Value = getUseableNeighborFromSetDepth(&beads_in_use_tree, pModel, swap1);
            }

            // get available neighbor position
            pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), originalSwap2Value);
            // make the swap, sort and update P(r)
            // randomly pick a point in beads_in_use
            // look at neighborhood for point.
            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, workingBinCount);

            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->calculateKLDivergence(workingBinCount);
            tempVolume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
            tempVolE = muConstant * tempVolume;//*invWorkingLimit;
            this_energy = testKL + tempVolE;// + lambda*connectivityPotential(tempNumberOfComponents);
            //this_energy = testKL;// if connectivity is one, potential is 0

            if (this_energy < current_energy) {
                beads_in_use_tree.insert(originalSwap2Value); // add new lattice
                std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                currentKL = testKL;
                current_energy = this_energy;
                current_volume = tempVolume;
                currentVolE = tempVolE;
                isSwapped = true;
                eulerTour.addNode(originalSwap2Value, pModel); // faster than removing
                lowTempSuccessCount++;
                sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
                break;
            }
            // reverse changes; find index in sorted list and replace
            // binCount is P(r) distribution without swap1 value
            std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
            std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
            itIndex = bead_indices.begin() + position;
        }

        // if no suitable location is found, return swap1 value to P(r) by copying binCountBackUp;
        if (!isSwapped){ // itIndex already reverses on exit from loop
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // may not be necessary
            std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
            beads_in_use_tree.insert(swap1); // add back swap one to tree
            sprintf(addRemoveText, "      FAILED => %i", swap1);
            eulerTour.addNode(swap1, pModel);
        }

        printf("   MAXSTEPS => %8i (%4i) \n", round, finalCoupons);
        printf("  SUCCESSES => %8d  %s\n", lowTempSuccessCount, addRemoveText);
        printf("  BIN [ONE] => %d VOL => %.0f\n", binCount[0], current_volume);
        printf("LIMIT: %5i D_KL : %.4E ( %.4E ) ENRGY : %.4E\n", workingLimit, currentKL, lowestKL, current_energy);
    }


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
    float    tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
//    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << " POSITIONAL SUCCESS/FAILURE " << ((double)failures/(double)attempts)*100 << endl;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

//    int totalLessThanAverage=0, numberContacts;
//    for (int i=0; i<workingLimit; i++){
//        numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
//        if (numberContacts < average_number_of_contacts ){
//            totalLessThanAverage +=1;
//        }
//    }

    // LOW TEMP ADD?REMOVE REFINEMENT

//
//    int finalRounds = (totalLessThanAverage*std::log((double)totalLessThanAverage) + 0.5772156649*totalLessThanAverage + 0.5);
//    cout << " TOTAL LATTICE POINTS LESS THAN AVERAGE " << totalLessThanAverage << endl;
    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);
//
//    int swappedCounter =0;
//    int round=0;
//    float average_number_of_contacts_before = tempAverageContacts/(float)workingLimit;
//    tempAverageContacts=0.0;
//    for (int i=0; i<workingLimit; i++){
//        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
//    }
//    average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    double runningContactsSum = calculateTotalContactSumPotential(&beads_in_use_tree, pModel);
//
    cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << endl;
    cout << "------------------------------------------------------------------------------" << endl;
    cout << " FINAL LOW TEMP REFINEMENT" << endl;
    cout << "              KL DIVERGENCE   ENERGY " << endl;
    printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    cout << "    CONTACTS" << endl;
    printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    printf("    BIN [ONE] => %d AVG => %.3f \n", binCount[0], binCount[0]/(float)workingLimit);
    printf("    VOLUME    OLD => %.0f FINAL => %.0f \n", old_volume, current_volume);
    printf("E_CONTACTS  FINAL => %.4E \n", runningContactsSum/(float)workingLimit);
    cout << "------------------------------------------------------------------------------" << endl;
//
    this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);
//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF
    int numpoints = workingLimit;
    coordT points[3*numpoints];
//    char flags[]="qhull FA";
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



bool Anneal::initializeModelToRefine(Model *pModel, Data *pData, std::string name, std::string PDBFilename) {

    srand(time(0));
    contactCutOff = interconnectivityCutOff;

    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * pBin = pModel->getPointerToBins(); // initialized as empty in Model class

    // convert distances within the large search space to bins based on input P(R)-DATA file
    this->fillPrBinsAndAssignTotalBin( pBin,  pDistance,  totalDistancesInSphere,  pData);
    std::cout << "MAX BIN : " << maxbin << std::endl;

    // initialize Universe and fill indices
    const int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    std::vector<int> bead_indices(totalBeadsInSphere); // large vector ~1000's

    int * ptr = (totalBeadsInSphere != 0) ? &bead_indices.front() : NULL;
    for(int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // returns only seed, needs to be mapped to entire Universe
    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    int workingLimit = pModel->getTotalInSeed();

    // sort reassigned bead model into working universe
    int locale=0;
    for(std::vector<int>::const_iterator it = pModel->getSeedBegin(); it != pModel->getSeedEnd(); ++it) {
        std::vector<int>::iterator fit = std::find(bead_indices.begin(), bead_indices.end(), *it); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+locale, fit);
        locale++;
    }

    // create working observed probability distribution that encompasses search sphere
    pData->createWorkingDistribution(maxbin);
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50

    std::cout << "    TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "    MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "              BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    //
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::cout << " WORKINGLIMIT SET => " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int currentNumberOfComponents  = eulerTour.getNumberOfComponents();

//    float lowest_volume, test_volume, current_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);//calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    float initialVolume = current_volume;

    std::cout << "    STARTING KL => " << currentKL << std::endl;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< std::endl;
    int high;

    pModel->writeModelToFile(workingLimit, bead_indices, name, high);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setBeadAverageAndStdev(workingLimit, workingLimit*0.1);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    std::printf("            FINAL => %.0f\n", current_volume);
    std::cout << "*******************                                        *******************" << std::endl;

    std::cout << "     MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "   SHANNON BINS IN DATA : " << pData->getShannonBins() << std::endl;
    std::cout << "            ZERO BIN AT : " << pData->getZeroBin() << std::endl;

    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    float    tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;

    return true;
}