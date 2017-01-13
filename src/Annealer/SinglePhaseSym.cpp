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
void Anneal::createInitialModelSymmetry(Model *pModel, Data *pData) {

    srand(time(0));

    //interconnectivityCutOff = pModel->getBeadRadius()*2.0001;
    //interconnectivityCutOff = pModel->getBeadRadius()*2.8285;
    contactCutOff = interconnectivityCutOff;

    float inv_kb_temp = 1.0/highT;

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
    std::vector<int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> workingBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<int> binCountBackUp(maxbin);        // smallish vector, typically < 50

    cout << "    TOTAL EXP N_S BINS : " << totalBins << endl;
    cout << "    MAX MODEL N_S BINS : " << maxbin << endl;
    cout << "              BINWIDTH : " << pData->getBinWidth() << endl;
    cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << endl;

    int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

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

    std::uniform_int_distribution<> number_of_beads_to_use (lowerN, upperN); // number of beads in ASU

    int subUnitWorkingLimit = number_of_beads_to_use(gen);
    pModel->writeModelToFile(totalBeadsInSphere, subUnit_indices, "universe");

    cout << "Bead Search Limited to: " << lowerN << " <= N <= " << upperN << endl;
    cout << "      INITIAL MODEL WL: " << subUnitWorkingLimit << endl;
    cout << "              SYMMETRY: " << pModel->getSymmetryString() << endl;
    cout << "        TOTAL SUBUNITS: " << pModel->getNumberOfSubUnits() << endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    std::vector<int>::iterator beginIt = subUnit_indices.begin();
    std::vector<int>::iterator endIt = subUnit_indices.end();

    cout << "        CREATING INITIAL RANDOM MODEL " << endl;
    std::shuffle(beginIt, endIt, gen);
    std::sort(beginIt, beginIt+subUnitWorkingLimit);
    std::set<int> beads_in_use_tree(beginIt, beginIt + subUnitWorkingLimit);

    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    int numpoints = 3*totalBeadsInSphere;
    coordT points[numpoints];
    vertexT * vertices;  // not sure if this needs to be explicitly deleted

    //int indexOfHullpt, count, potentialContacts, index;
    char flags[25];
    sprintf(flags, "qhull s FA");
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, points, pModel);
    refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
    // calculate starting energy

    std::clock_t start;
    start = std::clock();

    // fill binCount for first time
    float currentKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
    lambda = 0.001;
    mu = 0.00001; // scale this to volume of search space

    cout << " THREAD TIME "  << (std::clock() - start)/(double) CLOCKS_PER_SEC << " seconds " << endl;
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");
    // BIN calculation is  0.000865
    // HULL calculation is 0.0015 seconds

    float testKL, test_energy; // sets alpha as a constant during high temp eq

    int lowestWorkingLimit, lowestDeadLimit;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    int tempNumberOfComponents, currentNumberOfComponents;
    EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
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
    vector<int> averageV(highTempRounds);
    std::vector<int>::iterator pSwap2;
    int volumeCount=0, testIndex, originalSwap2Value;
    float volumeSum=0, workingLimitSum=0;

    float tempViolations, totalViolations=0;
    beta = 0.0001;
    for(int i=0; i<subUnitWorkingLimit; i++){
        totalViolations += symViolationsPotential(subUnit_indices[i], &subUnit_indices, subUnitWorkingLimit, pModel);
    }
    current_energy += beta * totalViolations;

    cout << "Total Violations " << totalViolations << endl;


    for (int high=0; high < highTempRounds; high++){ // iterations during the high temp search

        std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
        std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy

        //cout << high << " currentKL => " << currentKL << " calc => " << calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData ) << endl;

        if (distribution(gen) < 0){
            // create points from workSet to determine HULL
            float invSubUnitworkingLimit = mu/(float)subUnitWorkingLimit;
            bool isSwapped;

            sprintf(addRemoveText, "POSITIONAL");
            std::copy(beginIt, beginIt+subUnitWorkingLimit, active_indices.begin());   // make copy
            std::shuffle(active_indices.begin(), active_indices.begin() + subUnitWorkingLimit, gen);

            for (int v=0; v < 1; v++){
                isSwapped = false;
                testIndex = active_indices[v];
                itIndex = std::find(beginIt, beginIt + subUnitWorkingLimit, testIndex);

                std::copy(beginIt, endIt, backUpState.begin());
                std::copy(beginBinCount, endBinCount, binCountBackUp.begin());  // unaltered P(r)
                // remove contribution from P(r)-distribution, subUnitWorkingLimit stays same
                removeFromPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
                std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // make copy of altered P(r)

                beads_in_use_tree.erase(testIndex);
                // find better position
                // eulerTour.removeNode(testIndex);

                for (int bb = subUnitWorkingLimit; bb < deadLimit; bb++) {

                    pSwap2 = (beginIt + bb); //&subUnit_indices[bb]; // located in workingZone
                    originalSwap2Value = *pSwap2;
                    beads_in_use_tree.insert(*pSwap2);

                    std::iter_swap(itIndex, pSwap2); // iterators will point to same spot relative to beginIt
                    std::sort(beginIt, beginIt + subUnitWorkingLimit); // bead_indices needs to be sorted

                    addToPrSym(originalSwap2Value, subUnit_indices, subUnitWorkingLimit, workingBinCount, pModel, pData);

                    testKL = pData->calculateKLDivergence(workingBinCount);
                    //testKL = calculateKLEnergySymmetry(&subUnit_indices, &workingBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
                    //tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
                    isConnectedComponent(&subUnit_indices, subUnitWorkingLimit, pDistance, totalBeadsInSphere, tempNumberOfComponents);
                    test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, points, pModel);
                    test_energy = testKL + lambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + test_volume*invSubUnitworkingLimit;



                    if (test_energy < current_energy){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                        //std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
                        cout << "violations  " << symViolationsPotential(originalSwap2Value, &subUnit_indices, subUnitWorkingLimit, pModel) << endl;
                        isSwapped = true;
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)){
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        std::copy(workingBinCount.begin(), workingBinCount.end(), beginBinCount); // final copy so next round make backup
                        cout << "violations  " << symViolationsPotential(originalSwap2Value, &subUnit_indices, subUnitWorkingLimit, pModel) << endl;
                        //std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
                        isSwapped = true;
                        break;
                    }
                    // undo and test next position within deadLimit
                    std::copy(beginBinCount, endBinCount, workingBinCount.begin()); // removes swapped contribution from Pr in workingBinCount
                    std::copy(backUpState.begin(), backUpState.end(), beginIt);
                    itIndex = std::find(beginIt, beginIt + subUnitWorkingLimit, testIndex);
                    // original values of pSwap2 and itIndex are restored
                    beads_in_use_tree.erase(originalSwap2Value);
                    // eulerTour.removeNode(originalSwap2Value);
                    // itIndex should be pointing back to old vertex and bb is restored
                }

                if (!isSwapped) { // return to original
                    std::sort(beginIt, beginIt + subUnitWorkingLimit);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), beginBinCount); // add back swap 1 from backup
                    // eulerTour.addNode(testIndex, pModel);
                    beads_in_use_tree.insert(testIndex); // add back swap one to tree
                }
            }

            testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData );
            if (currentKL != testKL){
                cout << "POSITIONAL ERROR " << isSwapped << endl;
                // if this fails, update of P(r) is wrong or subUnit_indices is corrupted
                cout << high  << "           Failed " << pData->calculateKLDivergence(binCount) << " != " << testKL << endl;
                cout << high  << " Failed currentKL " << currentKL << " != " << testKL << endl;
                return;
            }

        } else { // add or remove

            // dead limit is set by convex hull
            alterMe = number_of_beads_to_use(gen); // pick number between lower and upper bounds

            if (subUnitWorkingLimit > alterMe && alterMe > lowerN) { // REMOVE beads from sorted list into useable range < deadLimit

                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                std::copy(beginIt, beginIt+subUnitWorkingLimit, active_indices.begin());
                std::shuffle(active_indices.begin(), active_indices.begin() + subUnitWorkingLimit, gen);
                sprintf(addRemoveText, "  REMOVE  ");

                // go through entire list
                for (int i = 0; i < subUnitWorkingLimit; i++) { // removing bead always decreases total volume?
                    // grab from randomized active_indices list
                    testIndex = active_indices[i];

                    tempViolations = (totalViolations - symViolationsPotential(testIndex, &subUnit_indices, subUnitWorkingLimit, pModel));

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
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        currentKL = testKL;
                        current_energy = test_energy;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalViolations = tempViolations;
                        beads_in_use_tree.erase(testIndex);
                        break;
                    } else { // undo changes and move to next bead (rejecting)
                        // std::copy(backUpState.begin(), backUpState.begin()+workingLimit+1, beginIt);
                        //currentNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                        restoreRemovingLatticePointFromBackUp(&beginIt, &subUnitWorkingLimit, &binCountBackUp,
                                                              &beginBinCount);
                    }

                } // end of for loop, break if successful


            } else { // ADD beads

                std::shuffle(beginIt+subUnitWorkingLimit, beginIt+deadLimit, gen);
                std::copy(beginIt, endIt, backUpState.begin());   // make backup copy
                sprintf(addRemoveText, "   ADD    ");

                for(int i=subUnitWorkingLimit; i<deadLimit; i++){

                    //testIndex = active_indices[i]; // remove from active_indices if used
                    //itIndex = find(beginIt + subUnitWorkingLimit, beginIt + deadLimit, testIndex);
                    itIndex = beginIt + i;
                    testIndex = *itIndex;
                    //cout << i << " " << active_indices[i] << " == " << subUnit_indices[i] << endl;
                    addLatticPositionToModel(&beginIt, &endIt, &backUpState, &subUnitWorkingLimit, &itIndex);
                    addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);

                    //testKL = calculateKLEnergySymmetry(&subUnit_indices, &binCount, subUnitWorkingLimit, totalBeadsInSphere, violations, pModel, pData);
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
                        //adjustDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                        beads_in_use_tree.insert(testIndex);
                        break;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        current_energy = test_energy;
                        currentKL = testKL;
                        current_volume = test_volume;
                        currentNumberOfComponents = tempNumberOfComponents;
                        totalViolations = tempViolations;
                        //adjustDeadLimitPerBead(beginIt, workingLimit, &deadLimit, pModel, totalBeadsInSphere, testIndex);
                        beads_in_use_tree.insert(testIndex);
                        break;
                    } else { // undo changes (rejecting)
                        restoreAddingFromBackUp(&beginIt, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                        //currentNumberOfComponents = eulerTour.removeNode(testIndex);
                    }
                }
            }

        } // END OF ADD?REMOVE

        //refineCVXHull(subUnit_indices, active_indices, totalBeadsInSphere, subUnitWorkingLimit, &deadLimit, pModel);
        populateLayeredDeadlimitUsingSet(&subUnit_indices, &beads_in_use_tree, subUnitWorkingLimit, &deadLimit, pModel);

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

        if (current_energy < lowest_energy){
            string name = "initial_sym_search_" + std::to_string(high);
            pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, name);
            pModel->writeModelToFile(deadLimit, subUnit_indices, name+"_skin");

            //pModel->writeSymModelToFile(subUnitWorkingLimit, subUnit_indices, name+"sym_");
            //counter++;
            std::copy(beginIt, endIt, lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;
            lowest_energy = current_energy;

            workingLimitSum += subUnitWorkingLimit;
            volumeSum += current_volume;
            volumeCount++;
        }

    } // end of HIGH TEMP EQUILIBRATION

    //pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, "best_HT_sym_" + filenameprefix);
    //std::sort(lowest_subUnit_indices.begin(), lowest_subUnit_indices.begin()+lowestWorkingLimit);

    pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage, volumeStdev;

    volumeAverage = workingLimitSum/(float)volumeCount;
    volumeStdev = 0.17*volumeAverage;

    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);
    //std::normal_distribution<float> volumeGen(beadAverage, beadStDev);
    cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << endl;
    cout << "Average R " << calculateAverageDistance(pDistance, &lowestWorkingLimit, &lowest_subUnit_indices, pModel) << endl;

    // remove points close to hull
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, "best_HT_sym_" + filenameprefix);

    pModel->setStartingSet(lowest_subUnit_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
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
inline void Anneal::addToPrSym(int addMeSubUnitIndex, vector<int> & beadsInUse, int workingLimit, vector<int> & prBins, Model *pModel, Data *pData){
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


// for a specified bead position, determine number of violations (within contactCutoff)
/**
 * for the specified lattice index, count number of violations when symmetry mates are created
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
            if ((*refVec - *targetVec).length() < cutOff){
                violations++;
            }
        }
    }

    return violations/(float)workingLimit;
}