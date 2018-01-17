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
    violation_limit = 1.51*pModel->getBeadRadius();

    //float inv_kb_temp = 1.0/highT;
    float inv_kb_temp = 1.0/0.0001;

    // convert distances to ShannonBin membership
    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
    float * pDistance = pModel->getPointerToDistance();
    int * const pBin = pModel->getPointerToBins();

    maxbin=0;
    int violations, temp_violations;
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
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
//    int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();
    // as bead indices are discarded, set upper limit of vector
    std::vector<int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<int> test_indices(totalBeadsInSphere);        // large vector ~1000's
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

    const float invBeadVolume = 1.0/pModel->getBeadVolume();

    //float inv_kb = totalBeadsInSphere*beadVolume/(pData->getShannonBins());
    int lowerN = round(lowerV*invBeadVolume/(float)pModel->getNumberOfSubUnits());
    int upperN = round(upperV*invBeadVolume/(float)pModel->getNumberOfSubUnits());

    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN);   // number of beads in ASU

    unsigned int subUnitWorkingLimit = number_of_beads_to_use(gen);
    pModel->writeModelToFile(totalBeadsInSphere, subUnit_indices, "universe", 0);

    std::cout << "Bead Search Limited to: " << lowerN << " <= N <= " << upperN << std::endl;
    std::cout << "      INITIAL MODEL WL: " << subUnitWorkingLimit << std::endl;
    std::cout << "              SYMMETRY: " << pModel->getSymmetryString() << std::endl;
    std::cout << "        TOTAL SUBUNITS: " << pModel->getNumberOfSubUnits() << std::endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    std::cout << "        CREATING INITIAL RANDOM MODEL " << std::endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());

    std::set<int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::uniform_int_distribution<int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    char flags[25];
    std::sprintf(flags, "qhull s FA");
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
    //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
    //populateLayeredDeadlimit(subUnit_indices.begin(), subUnitWorkingLimit, &deadLimit, pModel, totalBeadsInSphere); // resets deadLimit
    unsigned int deadLimit = recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, pModel, totalBeadsInSphere);
    enlargeDeadLimit(subUnit_indices, &deadLimit, pModel);
    std::uniform_int_distribution<int> randomDeadIndex (subUnitWorkingLimit, deadLimit-1);
    // calculate starting energy

    std::clock_t start = std::clock();

    // fill binCount for first time
    float currentKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
    float hlambda = 0.001;
    float muConstant = mu*currentKL/((1.0-mu)*current_volume);

    std::cout << " THREAD TIME "  << (std::clock() - start)/(double) CLOCKS_PER_SEC << " seconds " << std::endl;
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");

    float testKL, test_energy; // sets alpha as a constant during high temp eq

    int lowestWorkingLimit = subUnitWorkingLimit, lowestDeadLimit;
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
    auto beginBinCount = binCount.begin();
    auto endBinCount = binCount.end();

    char addRemoveText[50];
    std::vector<int> averageV(highTempRounds);
    float sum_x_squared=0;
    int volumeCount=0, testIndex;
    float volumeSum=0, workingLimitSum=0;

    int totalViolations = violations;//(double)subUnitWorkingLimit;
    beta = 0.01;
    current_energy += beta * totalViolations;

    bool updateCVX = false;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup cop

    int high;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) < 0.36551){
            // create points from workSet to determine HULL
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

            std::shuffle(indices_to_check.begin(), indices_to_check.end(), gen);
            int swapCount = 0, swap1;

            for (int v = 0; v < 1; v++) { // try n times

                swap1 = indices_to_check[v];
                isSwapped = false;
                // find bead to swap in active set
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);
                // remove selected index from P(r)
                std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
                int tempValue = totalViolations - removeFromPrSym(swap1, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
                std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // make copy of altered P(r)

                beads_in_use_tree.erase(swap1);
                eulerTour.removeNode(swap1);

                for (int bb = 0; bb < 3; bb++) {
                    //int neighbor = subUnit_indices[subUnitWorkingLimit+bb];
                    int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                    while (neighbor == -1 || neighbor == swap1){
                        neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                    }

                    // make the swap, sort and update P(r)
                    pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);
                    std::iter_swap(itIndex, pSwap2);
                    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit); // bead_indices needs to be sorted

                    temp_violations = tempValue + addToPrSym(neighbor, subUnit_indices, subUnitWorkingLimit, workingBinCount, pModel, pData);
                    testKL = pData->calculateKLDivergence(workingBinCount); // 100x faster than calculateKLEnergySymmetry
                    //testKL = calculateKLEnergySymmetry(&subUnit_indices, &workingBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );

                    tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

                    test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + test_volume*muConstant + beta*temp_violations;

                    if (test_energy < current_energy) {
                        beads_in_use_tree.insert(neighbor);
                        //move swapped position to the deadLimit
                       // std::iter_swap(pSwap2, subUnit_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                        //std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        totalViolations = temp_violations;
                        current_energy = test_energy;
                        isSwapped = true;
                        updateCVX = true;
                        swapCount++;
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        beads_in_use_tree.insert(neighbor);
                        //std::iter_swap(pSwap2, subUnit_indices.begin()+deadLimit-1-swapCount);
                        std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin());
                       // std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
                        currentNumberOfComponents = tempNumberOfComponents;
                        currentKL = testKL;
                        current_volume = test_volume;
                        totalViolations = temp_violations;
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
            // alterMe = number_of_beads_to_use(gen); // pick number between lower and upper bounds
            auto beginIt = subUnit_indices.begin();

            if (subUnitWorkingLimit > lowerN) { // REMOVE beads from sorted list into useable range < deadLimit
                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                sprintf(addRemoveText, "  REMOVE  ");

                testIndex = subUnit_indices[randomIndex(gen)];

                temp_violations = totalViolations - removeLatticePositionToModelSym(&beginIt, subUnit_indices, binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);

//                testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);

                test_volume = calculateCVXHULLVolume(flags,
                                                     &subUnit_indices,
                                                     subUnitWorkingLimit,
                                                     pModel);

                tempNumberOfComponents = eulerTour.removeNode(testIndex);

                test_energy = testKL + hlambda * (tempNumberOfComponents - 1) * (tempNumberOfComponents - 1) +
                              muConstant * test_volume + beta*temp_violations;

                if (test_energy < current_energy) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = temp_violations;
                    beads_in_use_tree.erase(testIndex);
                    updateCVX=true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    currentKL = testKL;
                    current_energy = test_energy;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = temp_violations;
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

//                testIndex = subUnit_indices[randomDeadIndex(gen)];

                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                while ( testIndex == -1) { // find a new position
                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                }
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.end(), testIndex);

                addLatticPositionToModel(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &itIndex);
                temp_violations = totalViolations + addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);
                //testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);

                test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1)
                              + muConstant*test_volume + beta*temp_violations;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = temp_violations;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentKL = testKL;
                    current_volume = test_volume;
                    currentNumberOfComponents = tempNumberOfComponents;
                    totalViolations = temp_violations;
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else { // undo changes (rejecting)
                    restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                    currentNumberOfComponents = eulerTour.removeNode(testIndex);
                }
            }
        } // END OF ADD?REMOVE

        if (updateCVX){
            //std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
            deadLimit = recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, pModel, totalBeadsInSphere);
            enlargeDeadLimit(subUnit_indices, &deadLimit, pModel);

            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<int>(0,subUnitWorkingLimit-1); // guaranteed unbiased
            randomDeadIndex = std::uniform_int_distribution<int> (subUnitWorkingLimit, deadLimit-1);
            updateCVX=false;
        }

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts
        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("         MAXSTEPS => %i (%4i) \n", highTempRounds, high);
        std::printf("            GRAPH => %3i\n", currentNumberOfComponents);
        std::printf("            UPPER => %i LOWER => %i \n", upperN, lowerN);
        std::printf("           VOLUME => %.0f  MU => %.4E  MU*VOL => %.3E\n", current_volume, muConstant, muConstant*current_volume);
        std::printf("       VIOLATIONS => %d  BETA => %.4E  \n", totalViolations, beta);
        std::printf("        LIMIT : %5i DEADLIMIT : %5i D_KL : %.4E ENRGY: %.4E \n", subUnitWorkingLimit, deadLimit, currentKL, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            //std::string name = "initial_sym_search_" + std::to_string(high);
            //pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, name, high);
            //pModel->writeModelToFile(deadLimit, subUnit_indices, name+"_skin", high);
            //pModel->writeSymModelToFile(subUnitWorkingLimit, subUnit_indices, name+"sym_");
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;

            muConstant = mu*currentKL/((1.0-mu)*current_volume);

            current_energy = currentKL + muConstant*current_volume + beta*totalViolations;

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
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            cout << high  << "           Failed " << pData->calculateKLDivergence(binCount) << " != " << testKL << endl;
//            cout << high  << " Failed currentKL " << currentKL << " != " << testKL << endl;
//            printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//
//            int bins = testBinCount.size();
//            std::cout << " BIN : " << std::endl;
//            for(int i=0; i<bins; i++){
//                if (testBinCount[i] != binCount[i]){
//                    std::cout << i << " : " << testBinCount[i] << " == " << binCount[i] << std::endl;
//                }
//            }
//
//            exit(0);
//        }

    } // end of HIGH TEMP EQUILIBRATION

    pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage;

    // remove points close to hull
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, "initial_CVX_subunit_" + filenameprefix, high);
    pModel->writeSymModelToFile(currentKL, lowestWorkingLimit, lowest_subUnit_indices, binCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, current_volume, 4.1);

    pModel->setStartingSet(lowest_subUnit_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);

    std::cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << std::endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::printf("    LOWEST => %d \n", lowestWorkingLimit);
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
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY" << std::endl;

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

    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    // reset iterators to internal bead_indices
    beginIt = bead_indices.begin();
    endIt = bead_indices.end();

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    //pModel->writeModelToFile(workingLimit, bead_indices, "restart" + filenameprefix, 0);

    unsigned int deadLimit = recalculateDeadLimit(workingLimit, bead_indices, pModel, totalBeadsInSphere);
    enlargeDeadLimit(bead_indices, &deadLimit, pModel);
    // set deadLimit of the selected set
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
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    std::cout << "   TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "   MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "             BINWIDTH : " << pData->getBinWidth() << std::endl; // lattice should be limited by binwidth

    auto beginBinCount = binCount.begin();
    auto endBinCount = binCount.end();

    int violations;
    float testKL, currentKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
    char flags[] = "qhull FA"; // CVX HULL STUFF
    float temp_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    //std::shuffle(bead_indices.begin()+workingLimit, bead_indices.end(), gen);
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float tempAverageContacts;//, updateSteps = 7.0;;
    // coupon collector's problem
    int coupons = (workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5);
    int updateCount = ccmultiple*coupons;
    float step_limit = (updateCount < 30000) ? 30000 : (float)updateCount;

    std::vector<float> tempDuringRun(step_limit);
    std::vector<float> divergenceDuringRun(step_limit);
    std::vector<int> workingLimitDuringRun(step_limit);

    float * const pTempDuringRun = &tempDuringRun.front();
    float * const pDivergenceDuringRun = &divergenceDuringRun.front();
    int * const pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int tempNumberOfComponents, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float acceptRate = 0.5, inv500 = 1.0/500.0;
    float inv500slash499 = 499.0/500.0;
    int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    float this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    double newSum, runningContactsSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);
    //double invWL = 1.0/(double)workingLimit;
    double invWL = 1.0;

    double etaConstant = eta*currentKL/((1.0-eta)*runningContactsSum*invWL);
    double tempTotalContactEnergy, totalContactEnergy = etaConstant*(runningContactsSum*invWL);

    // want D_KL to settle < 10^-5 to 10^-6
    // int span = (int)(std::floor(log10f(currentKL)) - std::floor(log10f(totalContactEnergy/(float)workingLimit))) + 3;
    //mu = 0;
    float muConstant = 1;// mu*currentKL/((1.0-mu)*current_volume);
    float tempVolE=0, currentVolE = muConstant*current_volume;

    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + currentVolE;

    beta=0.001;
    int totalViolations = violations, temp_violations;
    current_energy += beta*totalViolations;

    float lowest_energy = current_energy + totalContactEnergy, afterAdding, afterRemoving;;

    std::cout << "           VIOLATIONS : " << violations << std::endl; // lattice should be limited by binwidth
    std::cout << "        CURRENT ENRGY : " << current_energy << std::endl; // lattice should be limited by binwidth

    std::clock_t startTime;
    std::uniform_int_distribution<int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<int> randomDeadIndex (workingLimit, deadLimit-1);

    std::vector<int> selections(workingLimit);
    for(unsigned long int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    int numberOfCoolingTempSteps=0, updatedETA=0;
    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
        beginBinCount = binCount.begin();

        if (distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // additional points to expand deadlimit will occur via enlarging CVX Hull

            startTime = std::clock();
            if (distribution(gen) < 0.5){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;

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
                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &backUpState, &workingLimit, &itIndex);

                beads_in_use_tree.insert(original);
                newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);
                afterAdding = etaConstant*newSum;
                //afterAdding = etaConstant*newSum/(double)(workingLimit);

                temp_violations = totalViolations + addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
                testKL = pData->calculateKLDivergence(binCount);
//                testKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );

                tempNumberOfComponents = eulerTour.addNode(original, pModel);
//                temp_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                tempVolE = muConstant*temp_volume;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + tempVolE;
                tempTotalContactEnergy = (afterAdding - totalContactEnergy);

                if ( (this_energy + afterAdding) < (current_energy + totalContactEnergy) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    totalViolations = temp_violations;
//                    current_volume = temp_volume;
//                    currentVolE = tempVolE;
                    isUpdated = true;
                    // eulerTour.addNode(original, pModel);
                    currentNumberOfComponents = tempNumberOfComponents;
                    sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterAdding;
                    runningContactsSum = newSum;
                    totalViolations = temp_violations;
//                    current_volume = temp_volume;
//                    currentVolE = tempVolE;
                    isUpdated = true;
                    // eulerTour.addNode(original, pModel);
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
                temp_violations = totalViolations - removeLatticePositionToModelSym(&beginIt, bead_indices, binCount, &workingLimit, &original, pModel, pData);

                beads_in_use_tree.erase(original);
                newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);

                //afterRemoving = etaConstant*newSum/(double)(workingLimit);
                afterRemoving = etaConstant*newSum;
                testKL = pData->calculateKLDivergence(binCount);

//                testKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );

                tempNumberOfComponents = eulerTour.removeNode(original);
//                temp_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                tempVolE = muConstant*temp_volume;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + tempVolE;

                tempTotalContactEnergy = (afterRemoving - totalContactEnergy);

                if ((this_energy + afterRemoving) < (current_energy + totalContactEnergy) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    totalViolations = temp_violations;
//                    current_volume = temp_volume;
//                    currentVolE = tempVolE;
                    isUpdated = true;
                    currentNumberOfComponents = tempNumberOfComponents;

                } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy)*inv_kb_temp) > distribution(gen) ){
                    currentKL = testKL;
                    current_energy = this_energy;
                    totalContactEnergy = afterRemoving;
                    runningContactsSum = newSum;
                    totalViolations = temp_violations;
//                    current_volume = temp_volume;
//                    currentVolE = tempVolE;
                    isUpdated = true;
                    currentNumberOfComponents = tempNumberOfComponents;

                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
                    beginIt = bead_indices.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.addNode(original, pModel);
                }
            }

        } else { // positional refinement

            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            //float localPotentialAtOldPosition, localPotentialAtOldPositionWithOutBead, localPotentialAtNewPosition, localPotentialAtNewPositionWithBead;
            // find bead to swap in active set
            // only refine a single position with 2 or less contacts

            startTime = std::clock();

            // randomly select an index to move
            // select only node I can move?
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            int position = selections[0];
            swap1 = bead_indices[ position ];

            for(int i=1; i<workingLimit; i++){
                if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                    break;
                } else { // if last node, will break out
                    eulerTour.addNode(swap1, pModel);
                    position = selections[i];
                    swap1 = bead_indices[ position ]; // what happens if I fail all the way?
                }
            }

            // swap
            itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
            int tempValue = totalViolations - removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
            std::copy(binCount.begin(), binCount.end(), workingBinCount.begin()); // make copy of altered P(r)

            beads_in_use_tree.erase(swap1);

            // find better position
            int originalSwap2Value;
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

            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            temp_violations = tempValue + addToPrSym(originalSwap2Value, bead_indices, workingLimit, workingBinCount, pModel, pData);

            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->calculateKLDivergence(workingBinCount);

            beads_in_use_tree.insert(originalSwap2Value);
            newSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);
            double newPotential = etaConstant*newSum*invWL;

            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
//            temp_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//            tempVolE = muConstant*temp_volume;

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + tempVolE;
            tempTotalContactEnergy = newPotential - totalContactEnergy;

            if ((this_energy + newPotential) < (current_energy + totalContactEnergy)) {
//                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
//                    eulerTour.addNode(originalSwap2Value, pModel);
                std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin()); // final copy so next round make backup
                currentKL = testKL;
                current_energy = this_energy;
                currentNumberOfComponents = tempNumberOfComponents;
                totalContactEnergy = newPotential;
                runningContactsSum = newSum;
                isUpdated = true;
                totalViolations = temp_violations;
//                current_volume = temp_volume;
//                currentVolE = tempVolE;
                sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else if (exp(-(this_energy - current_energy + tempTotalContactEnergy) * inv_kb_temp) > distribution(gen)) {
//                    beads_in_use_tree.insert(originalSwap2Value); // add new lattice
//                    eulerTour.addNode(originalSwap2Value, pModel);
                std::copy(workingBinCount.begin(), workingBinCount.end(), binCount.begin());
                currentKL = testKL;
                current_energy = this_energy;
                currentNumberOfComponents = tempNumberOfComponents;
                totalContactEnergy = newPotential;
                runningContactsSum = newSum;
                isUpdated = true;
                totalViolations = temp_violations;
//                current_volume = temp_volume;
//                currentVolE = tempVolE;
                sprintf(addRemoveText, "  SA SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.erase(originalSwap2Value);

                beads_in_use_tree.insert(swap1); // add back swap one to tree
                sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()


        printf("       TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("     ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        printf("       TIME : %.5f (SECONDS)  %s\n", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC), addRemoveText);
        printf("      LIMIT : %5i  VOL => %-4.3E BIN[1] AVG => %.3f \n", workingLimit, current_volume, binCount[0]/(float)workingLimit);
        printf(" VIOLATIONS => %d  BETA => %.3E COMP => %d \n", totalViolations, beta, currentNumberOfComponents);
        printf("   CONTACTE => %-5.4E ETA => %.4E AVG => %.2f N => %d \n", totalContactEnergy, etaConstant, runningContactsSum, updatedETA);
        printf("       D_KL => %-5.4E ( %.4E ) ENRGY : %.4E\n", currentKL, lowest_energy, (current_energy+totalContactEnergy));

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;
//        float testKL1 = calculateKLEnergySymmetry(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
//        if (currentKL != testKL1){// || checkForRepeats(bead_indices)){
//            testKL = pData->calculateKLDivergence(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            return "stopped";
//        }

        // update run time parameters
        *(pTempDuringRun+numberOfCoolingTempSteps) = lowTempStop;
        *(pDivergenceDuringRun+numberOfCoolingTempSteps) = currentKL;
        *(pWorkingLimitDuringRun+numberOfCoolingTempSteps) = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            //invWL = 1.0/(double)workingLimit;
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
                updatedETA++;
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

    tempAverageContacts = tempAverageContacts*invWL;
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }

    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    // At end of each temp, update a probability model for volume?  Use this to select
    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
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
inline float Anneal::calculateKLEnergySymmetry(std::vector<int> *subUnit_indices, std::vector<unsigned int> *binCount, const int indicesWorkingLimit, const int totalBeadsInSphere, int &violation, Model *pModel, Data *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);
    violation=0;

    const int subUnits = pModel->getNumberOfSubUnits();
    const int totalCoordinates = subUnits*indicesWorkingLimit;
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
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    for (int i=0; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                violation++;
            }
        }
    }

    return pData->calculateKLDivergence(*binCount);
}


/**
 * The new position to Add MUST BE WITHIN THE WORKING LIMIT
 * beadsInUse must be sorted up to workingLimit
 *
 * @param addMeSubUnitIndex
 * @param beadsInUse
 * @param workingLimit
 * @param prBins
 * @param pModel
 * @param pData
 * @return
 */
inline int Anneal::addToPrSym(int addMeSubUnitIndex, std::vector<int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){
    // each bead has a sym related partner in pModel
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
//    const int totalSubUnits = pModel->getNumberOfSubUnits();
//    vector3 * refVec, * targetVec;
//    const unsigned long int  indexInMasterBeadArray = pModel->mapSubUnitIndex(addMeSubUnitIndex);
    int findIt=0;
    const int * ptr = &beadsInUse.front();
    for(int i=0; i<workingLimit; i++){
        if (ptr[i] == addMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const int stopAt = findIt;
    //
    // add interdomain distances
    //
//    const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);
//    // for selected bead (given by removeMe), calculate distance to all other selected beads
//    int next = 0;
//    int bin;
//
//    while(next < stopAt && (next < workingLimit)){ // stops when addMeSubUnitIndex is reached
//        // for each index in subunit
//        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[next]) );
//
//        prBins[ pData->convertToBin((*prefVec - *targetVec).length()) ] += totalSubUnits;
//        //totalAdded += totalSubUnits;
//        //for(int ss=0; ss<totalSubUnits; ss++){
//        //    prBins[ bin ]++;
//        //}
//        next++;
//    }
//
//    next++; // move past the addMeSubUnitIndex by incrementing
//    for(int b=next; b<workingLimit; b++){
//
//        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]));
//
//        prBins[ pData->convertToBin((*prefVec - *targetVec).length()) ] += totalSubUnits;
//        //totalAdded += totalSubUnits;
//        //for(int ss=0; ss<totalSubUnits; ss++){
//        //    prBins[ bin ]++;
//        //}
//    }
//
//    //
//    // add intra domain distances (between subunits)
//    //
//    int subUnitLimit = totalSubUnits-1;
//    int outer;
//    int violations=0;
//    float distance_to;
//
//
//    for(int s=0; s < subUnitLimit; s++) {
//
//        // iterate over each position in outer subunit
//        outer=0;
//        for(; outer < stopAt; outer++){ // go over all the atoms in the subunit
//
//            const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
//            for(int ss=(s+1); ss < totalSubUnits; ss++) {
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
//                distance_to = (*prefVec - *targetVec).length();
//                prBins[ pData->convertToBin(distance_to) ]++;
//
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                //  totalAdded++;
//            }
//        }
//
//        // outer is pointing to new position in outer subunit
//        const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[outer]) + s );
//        for(int ss=(s+1); ss < totalSubUnits; ss++) {
//
//            for(int subIndex=0; subIndex < workingLimit; subIndex++){
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex(beadsInUse[subIndex]));
//                distance_to = (*prefVec - *targetVec).length();
//                prBins[ pData->convertToBin(distance_to) ]++;
//
//
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                //  totalAdded++;
//            }
//        }
//
//        outer++;
//        for(; outer < workingLimit; outer++){ // go over all the atoms in the subunit
//
//            const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[outer]) + s);
//            for(int ss=(s+1); ss < totalSubUnits; ss++) {
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + indexInMasterBeadArray);
//                distance_to = (*prefVec - *targetVec).length();
//                prBins[ pData->convertToBin(distance_to) ]++;
//
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                // totalAdded++;
//            }
//        }
//    }


    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin

    int violations=0;
    float distance_to;

    const int subUnits = pModel->getNumberOfSubUnits();
    const int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;

    // create first subunit from selected indices and populate coordinates
    for (int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(beadsInUse[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    int count = workingLimit;
    for (int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    // adjust Pr
    for (int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                ++violations;
            }
        }

        for (int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                ++violations;
            }
        }
    }


    for (int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
            if (distance_to < violation_limit){
                --violations;
            }
        }
    }

    return violations;
}






/**
 * beadsInUse does not have to be sorted
 *
 * prBins is the model P(r) distribution
 */
inline int Anneal::removeFromPrSym(int const &removeMeSubUnitIndex, std::vector<int> & beadsInUse, int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){

    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    //const int totalSubUnits = pModel->getNumberOfSubUnits();
    //vector3 * targetVec;
    //const unsigned long int indexInMasterBeadArray = pModel->mapSubUnitIndex(removeMeSubUnitIndex);
    int findIt=0;
    //int totalCalculations = 0;
    // expected => n*intradomain + interdomain
    //int expected = totalSubUnits*(workingLimit-1)

    // find index of bead to remove
    int * ptr = &beadsInUse.front();
    for(int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const int stopAt = findIt;
    //
    // remove interdomain distances
    //
//    const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray);
//
//    Bead * tempBead = pModel->getBead(removeMeSubUnitIndex);
//    if ((tempBead->getVec() - *prefVec).length() != 0){
//        std::cout << " Not the same 1" << std::endl;
//        exit(0);
//    }
//
//    // for selected bead (given by removeMe), calculate distance to all other selected beads
//    unsigned int next = 0;
//    for(; next<stopAt; next++){ // stops when addMeSubUnitIndex is reached
//        // for each index in subunit
//        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse( pModel->mapSubUnitIndex(beadsInUse[next]) );
//
//        tempBead = pModel->getBead(beadsInUse[next]);
//        if ((tempBead->getVec() - *targetVec).length() != 0){
//            std::cout << " Not the same 2" << std::endl;
//            exit(0);
//        }
//
//        prBins[ pData->convertToBin((*prefVec - *targetVec).length()) ] -= totalSubUnits;
//        totalCalculations++;
//    }
//
//    next++; // move past the addMeSubUnitIndex by incrementing
//    for(; next<workingLimit; next++){
//        targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[next]));
//
//        tempBead = pModel->getBead(beadsInUse[next]);
//        if ((tempBead->getVec() - *targetVec).length() != 0){
//            std::cout << " Not the same 3" << std::endl;
//            exit(0);
//        }
//
//        prBins[ pData->convertToBin((*prefVec - *targetVec).length()) ] -= totalSubUnits;
//        totalCalculations++;
//    }
//
//    //
//    // remove intra domain distances (between subunits)
//    //
//    //int subUnitLimit = totalSubUnits-1;
    int violations=0;
    float distance_to;
//    // remove distances of refVec to subUnits
//    for(unsigned int s=0; s < totalSubUnits; s++) {
//        const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray+s);
//
//        for(int ss=(s+1); ss < totalSubUnits; ss++) {
//
//            for(int subIndex=0; subIndex < workingLimit; subIndex++){
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex(beadsInUse[subIndex]));
//                distance_to = (*prefVec - *targetVec).length();
//                //prBins[ pData->convertToBin(distance_to) ]--;
//                --prBins[ pData->convertToBin(distance_to) ];
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                totalCalculations++;
//            }
//        }
//    }
//
//    // for each remaining bead, remove contributions to the removed points
//    for(unsigned int s=1; s < totalSubUnits; s++) {
//
//        const vector3 * prefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray+s);
//
//        for(int ss=0; ss < s; ss++) {
//
//            int b=0;
//            for(; b<stopAt; b++){
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]) + ss);
//                distance_to = (*prefVec - *targetVec).length();
//                --prBins[ pData->convertToBin(distance_to) ];
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                totalCalculations++;
//            }
//
//            b++;
//            for(; b<workingLimit; b++){
//                targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(pModel->mapSubUnitIndex(beadsInUse[b]) + ss);
//                distance_to = (*prefVec - *targetVec).length();
//                --prBins[ pData->convertToBin(distance_to) ];
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//                totalCalculations++;
//            }
//        }
//    }
//
//    std::cout << "TOTAL CALCULATONS on REMOVE " << totalCalculations << std::endl;


    const int subUnits = pModel->getNumberOfSubUnits();
    const int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;

    // create first subunit from selected indices and populate coordinates
    for (int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(beadsInUse[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    int count = workingLimit;
    for (int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    // adjust Pr
    for (int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                ++violations;
            }
        }

        for (int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                ++violations;
            }
        }
    }

    /**
     *
     * we double remove the symmetry related vectors, need to add back
     */
    for (int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                --violations;
            }
        }
    }

    return violations;
}

/**
 *
 *
 */
inline int Anneal::removeLatticePositionToModelSym(std::vector<int>::iterator * pBeginIt,
                                                    std::vector<int> & bead_indices,
                                                    std::vector<unsigned int> & modelPrBins,
                                                    unsigned int * pWorkingLimit,
                                                    const int * pLatticePointToRemove, Model * pModel, Data *pData){

    std::vector<int>::iterator itIndex = std::find(*pBeginIt, *pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    int violations = removeFromPrSym(*pLatticePointToRemove, bead_indices, *pWorkingLimit, modelPrBins, pModel, pData);
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
    return violations;
}




//need to determine, for each bead, number of contacts it makes

/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSumPotentialSym(std::set<int> *beads_in_use, std::vector<int> * beads_indices, int workingLimit, Model *pModel){

    //calculate contacts per bead for selectedIndex
    double sum=0;
    const int totalSubUnits = pModel->getNumberOfSubUnits();
    vector3 * targetVec, * tempVec1;

    //int * ptr = &(beads_indices->front());
    int sizeOfSubUnitNeighbors = pModel->getSizeOfSubUnitNeighbors();
    std::vector<int> goodNeighbors(totalSubUnits, 0);

    // create coordinate space
    const int totalCoordinates = totalSubUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    // create first subunit from selected indices and populate coordinates
    const int * ptr = &(beads_indices->front());
    for (int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(ptr[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }


    // multithread this part
    // can be made thread safe since no indices will overlap in writing to coordinates
    int count = workingLimit;
    for (int s=1; s < totalSubUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    // determine contacts
    float distance_to;

    for (std::set<int>::iterator it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {

        int subTotal = numberOfContactsFromSet(beads_in_use, pModel, *it);
        // now include contacts from symmetry
        // int indexInMasterBeadArray = pModel->mapSubUnitIndex(*it);
        // const vector3 * pRefVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(indexInMasterBeadArray); // calculate distances from ref to other
        const vector3 * pRefVec = &pModel->getBead(*it)->getVec();
        // check first neighbor =1
        // Need to chehck that at least one lattice point is in contact with neighboring subunit
        for(int ss=1; ss < totalSubUnits; ss++) {

            if (pModel->isSubUnitNeighbor(ss)){
//                bool tempTest = true;
                int * pGood = &goodNeighbors[ss];

                for(int subIndex=0; subIndex < workingLimit; subIndex++){

                    targetVec = &coordinates[subIndex + workingLimit*ss];
                    //targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex( ptr[subIndex] ));
                    distance_to = (*pRefVec - *targetVec).length();

                    if (distance_to < contactCutOff && distance_to > violation_limit){ /*!< distance between lattice points to be counted as a neighbor */
//                        subTotal++;
                        // if no distance can be considered a contact between neighbors must consider not touching
                        (*pGood) = 1;
                    }
                }

            }
            //else {

//                for(int subIndex=0; subIndex < workingLimit; subIndex++){
//                    targetVec = &coordinates[subIndex + workingLimit*ss];
//                    //targetVec = pModel->getVectorOfBeadCoordinateInSymBeadUniverse(ss + pModel->mapSubUnitIndex( ptr[subIndex] ));
//                    distance_to = (*pRefVec - *targetVec).length();
//
//                    if (distance_to < contactCutOff && distance_to > violation_limit){ /*!< distance between lattice points to be counted as a neighbor */
//                        subTotal++;
//                        // std::cout << *it << " subunit : " << ss << " index " << subIndex << " = " << partial << " radius " << radius << " =: " << distance_to << std::endl;
//                    }
//                }
//            }
        }


        //if (subTotal > 12){
        //    sum += totalContactsPotential(12);
        //} else {
            sum += totalContactsPotential(subTotal);
        //}
        //std::cout << *it << " CNT : " << subTotal << " partial " << partial << " = " << (subTotal - partial) << std::endl;
    }

    // Determine penalty if base subunit is in contact with its neighbors
    int sumIt = 0;
    for(int i=0; i<totalSubUnits; i++){
        sumIt += goodNeighbors[i];
    }

    int diff = sizeOfSubUnitNeighbors - sumIt;
    if (diff > 0){
        sum += 2*totalContactsPotential(0);
    } else if (sumIt < 0){
        std::cout << "less than zero " << sumIt ;
        exit(0);
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

    return (sum + totalContactsPotential(newContactsFromSelectedIndex));
    //return (sum + newContactsFromSelectedIndex);
}


bool Anneal::initializeModelToRefineSym(Model *pModel, Data *pData, std::string name, std::string PDBFilename) {

    srand(time(0));
    contactCutOff = interconnectivityCutOff;
    violation_limit = 1.51*pModel->getBeadRadius();

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
    pModel->createSeedFromPDBSym(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    unsigned int workingLimit = pModel->getTotalInSeed();

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

    //float currentKL = calculateKLEnergy(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    int violations;
    float currentKL = calculateKLEnergySymmetry(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, violations, pModel, pData );
    //
    std::set<int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::cout << " WORKINGLIMIT SET => " << workingLimit << std::endl;
    std::cout << "       VIOLATIONS => " << violations << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // area of the model P(r) distribution
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    int currentNumberOfComponents  = eulerTour.getNumberOfComponents();

//    float lowest_volume, test_volume, current_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);//calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    float initialVolume = current_volume;
    std::cout << "      STARTING KL => " << currentKL << std::endl;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    //int high;
    //pModel->writeModelToFile(workingLimit, bead_indices, name, high);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setBeadAverageAndStdev(workingLimit, 0);

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
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

    double runningContactsSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;


    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);

    return true;
}