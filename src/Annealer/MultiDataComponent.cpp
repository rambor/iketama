//
// Created by Robert Rambo on 07/08/2017.
//


#include "MultiDataComponent.h"


// have to set the data with the highest resolution as the base universe
// each dataset gets its own universe



/**
 * read segmentation file
 * @param setup_filename
 */
MultiDataComponent::MultiDataComponent(
        std::string segmentation_file,
        int highTempRounds,
        float contacts,
        int totalSASteps,
        float eta,
        float lambda,
        float alpha,
        float mu,
        float percentAddRemove,
        float contactsPerBead):
        setupFile(segmentation_file),
        highTempRounds(highTempRounds),
        contacts(contacts),
        totalSASteps(totalSASteps),
        eta(eta),
        lambda(lambda),
        alpha(alpha),
        mu(mu),
percentAddRemove(percentAddRemove){

    checkPhaseFormat();

    boost::regex prDataFormat("(prdata)|(PRDATA)");
    boost::regex nameFormat("(name)|(NAME)");
    boost::regex phaseFormat("(phase)|(PHASE)");
    boost::regex volumeFormat("(volume)|(VOLUME)|(VOL)|(vol)");
    boost::regex contrastFormat("(contrast)|(CONTRAST)");
    boost::regex sigmaFormat("(sigma)|(SIGMA)");
    boost::regex belongsToFormat("(BELONGS_TO)|(belongs_to)|(belongsto)|(BELONGSTO)");
    boost::regex multiFormat("(mult)|(MULT)");
    boost::regex linkedFormat("(linked_to)|(LINKED_TO)");

    // parse the file and add create Universes (based on PRDATA)
    /*
     * read setup file
     * Steps:
     * 1. parse the file for PRDATA and create Universe for each PRDATA
     * 2. parse the file and make a phase for each
     * delimiter is =>
     */
    std::string line;
    std::vector<std::string> tempLine;
    std::vector<std::string> fileLines;
    link_potential.resize(4);
    link_potential[0]=100;
    link_potential[1]=10;
    link_potential[2]=3;
    link_potential[3]=1;

    std::ifstream segmentationFileStream (this->setupFile);
    if (segmentationFileStream.is_open()) {
        while(!segmentationFileStream.eof()) {
            std::getline(segmentationFileStream, line);
            fileLines.push_back(line);
        }
        segmentationFileStream.close();
    }

    int totalLines = fileLines.size();

    for(int i=0; i<totalLines; i++){
        line = fileLines[i];
        // if line contains PRDATA
        if (boost::regex_search(line, prDataFormat) && boost::regex_search(line, nameFormat)) {
            // split the line and make new model object
            std::string tempFilename, tempName;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
            // where in vector is model?
            std::vector<std::string>::iterator it = std::find(tempLine.begin(), tempLine.end(), "PRDATA");
            int index = std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempFilename = tempLine[index+2];
            }

            it = std::find(tempLine.begin(), tempLine.end(), "NAME");
            index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempName = tempLine[index+2];
            }

            this->addDataset(tempFilename, tempName);
        }
    }

    if (datasets.size() < 1){
        throw std::invalid_argument( "**** ERROR => NO PRDATA FILE FOUND CHECK INPUT FILE: \n");
    }

    // create Universes, need to sort out which PRDATA had largest dmax
    universalDmax=0.0d;
    universalBinWidth=111.0d;
    for(it_type mit = datasets.begin(); mit != datasets.end(); ++mit){
        Data * tempData = &mit->second;
        universalBinWidth = (tempData->getBinWidth() < universalBinWidth) ? tempData->getBinWidth() : universalBinWidth;
        universalDmax = (tempData->getDmax() > universalDmax) ? tempData->getDmax() : universalDmax;
    }

//    char tempText[80];
//    std::sprintf(tempText,"%.2f", universalDmax);
//    FUNCTIONS_RPR::printInfo("UNIVERSAL DMAX FOR ALL UNIVERSES : " + std::string(tempText));
    // make Universes, each will have same dmax and binWidth and number of beads
    for(it_type mit = datasets.begin(); mit != datasets.end(); ++mit){
        std::string name = mit->first;
        FUNCTIONS_RPR::printInfo("CREATING UNIVERSE : " + name);
        universes.insert(std::pair<std::string, Universe> (name, Universe(mit->first, &(mit->second), this->universalDmax, this->universalBinWidth, contactsPerBead)));
    }

    inv_total_universes = 1.0/universes.size();

    //
    // now, need to make the Segments within each universe
    // grab line that starts with PHASE
    // split
    // set volume
    // set contrast -> used to calculate D_KL
    FUNCTIONS_RPR::printInfo("INITIALIZING SEGMENTATION : ");
    boost::regex arrowFormat("=>");
    for(int i=0; i<totalLines; i++){

        line = fileLines[i];
        // if line contains PHASE
        if (boost::regex_search(line, phaseFormat) && boost::regex_search(line, volumeFormat) && boost::regex_search(line, belongsToFormat)) {
            // split the line and make new model object
            float tempSigma = 0.1, contrast = 1.0f; // default
            std::string tempPhase, tempBelongsTo;
            int tempVol;

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            // where in vector is model? need to correct for case?
            std::vector<std::string>::iterator it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("PHASE") );
            int index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempPhase = tempLine[index+2];
                boost::to_upper(tempPhase);
                FUNCTIONS_RPR::printInfo("PROCESSING SEGMENT: " + tempPhase);
            }

            it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("VOL")); // case insensitive find in vector?
            index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempVol = std::stoi(tempLine[index+2]); // limit on volume > 2*beadVolume?
                FUNCTIONS_RPR::printInfo("    VOLUME : " + tempLine[index+2]);
            }

            it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("SIGMA")); // case insensitive find in vector?
            index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempSigma = std::stof(tempLine[index+2]);
                FUNCTIONS_RPR::printInfo("    SIGMA : " + tempLine[index+2]);
            }

            it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("CONTRAST")); // case insensitive find in vector?
            index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                 contrast = std::stof(tempLine[index+2]);
                FUNCTIONS_RPR::printInfo("    CONTRAST : " + tempLine[index+2]);
            }

            // check the belongs_to
            std::vector<std::string> tempStringVector;
            boost::split(tempStringVector, line, boost::is_any_of("\t  "), boost::token_compress_on);
            int totalInSplit = tempStringVector.size();

            bool notFound = false;

            for(int ii=0; ii<totalInSplit; ii++){

                if (boost::regex_match(tempStringVector[ii], belongsToFormat) && boost::regex_match(tempStringVector[ii+1], arrowFormat) ){ // get file name
                    std::string tempBelongsToUni = tempLine[ii+2]; // belongs_to name
                    boost::to_upper(tempBelongsToUni); // must be uppercase

                    // make a segment in the Universe with name => tempBelongsToUni
                    // multiplicity means make the same thing multiple times
                    // should the updating of model be linked, like add/remove be coordinated?
                    Universe * pUniverse = &(universes.find(tempBelongsToUni)->second);
                    notFound = pUniverse->createPhase(tempPhase, tempVol, 1, contrast, tempSigma);
                }
            }

            if (notFound){
                phases.insert(tempPhase);
            } else {
                FUNCTIONS_RPR::printError("PHASE MAPPING ERROR CHECK SEG FILE");
                exit(0);
            }
        }
    }

    // if two phases are explicitly linked/connected, then the two sets must have at last one neighbor in common
    for(int i=0; i<totalLines; i++) {

        line = fileLines[i];
        if (boost::regex_search(line, phaseFormat) &&  boost::regex_search(line, linkedFormat)) {
            std::string tempLink, tempPhase;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            std::vector<std::string>::iterator it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("PHASE") );
            int index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempPhase = tempLine[index+2];
                boost::to_upper(tempPhase);
                FUNCTIONS_RPR::printInfo("PROCESSING SEGMENT: " + tempPhase);
            }

            // where in vector is model? need to correct for case?
            it = std::find_if(tempLine.begin(), tempLine.end(), find_by_string("LINKED_TO") );
            index =  std::distance(tempLine.begin(),it);
            if (it != tempLine.end()){
                tempLink = tempLine[index+2];
                boost::to_upper(tempLink);
                FUNCTIONS_RPR::printInfo("PROCESSING SEGMENT: " + tempLink);
                links.insert(std::pair<std::string, std::string>(tempPhase, tempLink));
            }

        }
    }

    std::cout << "**************************-----------------------**************************" << std::endl;
    //FUNCTIONS_RPR::printInfo("PROCESSING SEGMENT: " + tempLink);


    /*
     * create baseUniverses for each phase
     * 1. Iterate through each phase
     * 2. check all the universe to see if phase is in the universe
     * 3. if phase is in unverse, check the resolution (bead radius)
     * 4. select universe with smallest bead_radius => set as base
     */
    for(std::set<std::string>::iterator itPhaseName=phases.begin(); itPhaseName != phases.end(); ++itPhaseName){
        float tempWidth = 1020203; // arbitrary number

        std::string baseUniverseForPhase=universes.begin()->first; // default to first entry
        // find the universes that contain smallest radius, this will be the base
        for(auto &entry : universes) {
            auto &key = entry.first;
            auto &uni = entry.second;

            if (uni.containsPhase(*itPhaseName) && (uni.getBaseWidth() < tempWidth)){
                tempWidth = uni.getBaseWidth();
                baseUniverseForPhase=key;
//                FUNCTIONS_RPR::printInfo("    Phase : " + *itPhaseName + " UNI : " + uni.getNameOfUniverse() + " binwidth : " + std::to_string(tempWidth) );
//                std::cout << " Phase => " << *itPhaseName << " UNI : " << uni.getNameOfUniverse() << " biWidth " << tempWidth << std::endl;
            }
        }

        // base universe could be base for multiple phases
        // or base could be for one phase but not for another
        baseUniversesForEachPhase.insert(std::pair<std::string, Universe *> (*itPhaseName, &(universes.find(baseUniverseForPhase)->second)));
        /*
         * each baseUniverse holds pointers to non-base Universes
         */
        for(auto &entry : universes) {
            auto &key = entry.first;
            auto &uni = entry.second;

            if (uni.containsPhase(*itPhaseName) && (key != baseUniverseForPhase)  ) {
                FUNCTIONS_RPR::printInfo("PARALLEL UNIVERSE : " + key + " | FOR BASE : " + baseUniverseForPhase );
                (universes.find(baseUniverseForPhase)->second).addParallelUniverse(&uni);
                //universes[baseUniverseForPhase].addParallelUniverse(&uni);  // this will map base onto parallel universes
            }
        }
    }

    totalInvPhases = 1.0/(double)baseUniversesForEachPhase.size();

    // set weights for each Segment
    // get largest volume
    double maxVolume = 0;
    for(auto &entry : baseUniversesForEachPhase) {
        auto const &phase = entry.first;
        auto &baseUniverse = entry.second;
        // make initial model of the phase
        if (baseUniverse->containsPhase(phase)){
            if (baseUniverse->getSegments(phase)->getExpectedVolume() > maxVolume){
                maxVolume = baseUniverse->getSegments(phase)->getExpectedVolume();
            }
        } else {
            std::cout << " ERROR : MISSING PHASE => " << phase << std::endl;
            exit(0);
        }
    }

    // set weights
    for(auto &entry : baseUniversesForEachPhase) {
        auto const &phase = entry.first;
        auto &baseUniverse = entry.second;
        // make initial model of the phase
        if (baseUniverse->containsPhase(phase)){
            // calculate weight for phase
            double temp = maxVolume/baseUniverse->getSegments(phase)->getExpectedVolume();
            std::cout << " SETTING WEIGHT " << phase << " => " << temp << std::endl;
            for(auto &entry : universes) { // go through each universe and see if segment exists
                auto &uni = entry.second;
                if (uni.containsPhase(phase)){
                    uni.getSegments(phase)->setWeight(temp);
                }
            }
        }
    }
}

/**
 * PRDATA => file1.dat NAME => A
 * PRDATA => file2.dat NAME => AB
 * PRDATA => file3.dat NAME => AC
 * @param filename
 * @param name
 */
void MultiDataComponent::addDataset(std::string filename, std::string name) {

    // open file
    // create IofQ object

    boost::to_upper(name); // change to uppercase
    FUNCTIONS_RPR::printInfo("ADDING DATASET " + filename + " NAME => " + name );

    datasets.insert(std::pair<std::string, Data> (name, Data()) );
    datasets.find(name)->second.addPofRData(filename);

    // each phase needs to be assigned a base universe
//    if (datasets.find(name)->second.getBinWidth() < smallestBinWidth){
//        this->smallestBinWidth = datasets.find(name)->second.getBinWidth();
//        baseUniverseKey = name;
//    }

    // for now, everything has the same weight
    weights.insert(std::pair<std::string, float> (name, 1.0f));
}



bool MultiDataComponent::createInitialModel() {

    Universe * pUniverseToModify;
    // for each round
    //    randomize phases list
    //    for through each phase
    //         1. add or remove
    //         2. positional
    //         3. swap

    // for each phase in each universe, make randomize selection limited by upperlimit
    // * phase is not aware of the other phases
    // * universe has list of beads_in_use
    srand(time(0));

    // energy terms
    // cvxhull
    // pr distribution

    // have to make model for each Universe
    // starting with the highest resolution as base
    // if next Universe contains component from prior Universe in common, map the selected beads into the next
    // and create new component if specified

    // go through each phase and make initial model in respective bases;
    // each base universe will propagate changes to parallel universes

    FUNCTIONS_RPR::printInfo("CREATING INITIAL DISTRIBUTION OF POINTS FOR EACH PHASE");

    for(auto &entry : baseUniversesForEachPhase) {
        auto const &phase = entry.first;
        auto &baseUniverse = entry.second;
        // make initial model of the phase
        FUNCTIONS_RPR::printInfo(" PHASE " + phase + " BASE UNIVERSE : " + baseUniverse->getNameOfUniverse());
        if (baseUniverse->containsPhase(phase)){
            baseUniverse->createInitialModelOfPhase(phase);
        } else {
            std::cout << " ERROR : MISSING PHASE => " << phase << std::endl;
            exit(0);
        }
    }

    for(auto &entry : universes) {
        auto &uni = entry.second;
        // make initial model of the phase
        uni.initializeContactsPotential();
    }

    // at this point all segments/phases should be constructed in each of the respective universes
    // indices should be sorted for each universe
    // we need to calculate P(r)-distribution for each but need to determine contrast weights

//    for (auto &entry : universes){
//        auto const &name = entry.first;
//        auto &uni = entry.second;
//        std::cout << " name " << name << " " << uni.getTotalSegments() << std::endl;
//        //uni.validateIndices();
//    }

    /**
     *
     * Energy is :
     * D_kl
     * Volume
     * Components (Euler Tour)
     *
     */
    double inv_kb_temp = 1.0d/0.01d;
    double tempTotalKLEnergy, currentTotalKLEnergy = calculateTotalKLEnergy();

    // create eulur tour
    double hlambda = 0.1d;

    this->initializeEulerTours();
    double tempConnectivityEnergy, currentConnectivityEnergy=calculateTotalComponentEnergy(hlambda);
    /*
     * CONVEX HULL
     * Can be calculated over
     * 1. entire set of selected indices in Universe
     * 2. selected indices within Segment
     */
    double tempVolumeSum, currentVolumeSum = this->calculateVolumeEnergy();
    double muConstant = mu*currentTotalKLEnergy/((1.0-mu)*currentVolumeSum);
    double temp_volume_energy, current_volume_energy = muConstant*currentVolumeSum;

    double temp_link, current_link = this->calculateLinkPotential(1.0);

    double betaConstant = (current_link == 0) ? 1 : (0.2*current_link/((1.0-0.2)*current_link));
    double current_link_energy = betaConstant*current_link;

    double test_energy = 0, current_energy = currentTotalKLEnergy + current_volume_energy + currentConnectivityEnergy + current_link_energy;
    double lowest_energy = current_energy;

    int totalPhases = phases.size();
    std::string phaseToModify;
    std::vector<std::string> phases_vector(totalPhases);
    std::copy(phases.begin(), phases.end(), phases_vector.begin());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    int high, upperbound, lowerbound, cwl;
    char addRemoveText[75];
    std::clock_t startTime;
    double runtime;
    int lowCounter = 0;

//    exit(0);

    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search
        // randomly phase vector and iterate through
        std::shuffle(phases_vector.begin(), phases_vector.end(), gen);
        startTime = std::clock();

        for (int ph = 0; ph < totalPhases; ph++) {

            phaseToModify = phases_vector[ph];
            pUniverseToModify = baseUniversesForEachPhase.find(phaseToModify)->second;
            std::cout << "*******************----------------------------------------*******************" << std::endl;

            if (distribution(gen) < 0.51) { // add or remove from phase

                upperbound = pUniverseToModify->getSegments(phaseToModify)->getUpperNumberOfBeads();
                lowerbound = pUniverseToModify->getSegments(phaseToModify)->getLowerNumberOfBeads();
                cwl = pUniverseToModify->getSegments(phaseToModify)->getWorkingLimit();

                if ( (distribution(gen) > 0.5 && cwl < upperbound) ) {
                    std::sprintf(addRemoveText, "     %s => %s [ %s ] ", "  ADD ", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());
                    pUniverseToModify->addLatticePointToPhase(phaseToModify);
                    // get energy of new state
                    tempTotalKLEnergy = calculateTotalKLEnergy();
                    // calculate connectivity potential
                    tempConnectivityEnergy = calculateTotalComponentEnergy(hlambda);
                    // calculate new volume
                    tempVolumeSum = calculateVolumeEnergy();
                    temp_volume_energy = muConstant * tempVolumeSum;
                    temp_link = calculateLinkPotential(1.0);

                    test_energy = tempTotalKLEnergy + temp_volume_energy + tempConnectivityEnergy + betaConstant*temp_link;

                    if (test_energy < current_energy) {
                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else { // undo changes (rejecting)
                        pUniverseToModify->restoreFromBackupAdd(phaseToModify);
                    }


                } else if (cwl > lowerbound) {
                    sprintf(addRemoveText, "     %s => %s [ %s ]", "REMOVE", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());
                    pUniverseToModify->removeLatticePointFromPhase(phaseToModify);
                    tempTotalKLEnergy = calculateTotalKLEnergy();
                    // calculate connectivity potential
                    tempConnectivityEnergy = calculateTotalComponentEnergy(hlambda);
                    // calculate new volume
                    tempVolumeSum = calculateVolumeEnergy();
                    temp_volume_energy = muConstant * tempVolumeSum;
                    temp_link = calculateLinkPotential(1.0);

                    test_energy = tempTotalKLEnergy + temp_volume_energy + tempConnectivityEnergy + betaConstant*temp_link;

                    if (test_energy < current_energy) {
                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else { // undo changes (rejecting)
                        pUniverseToModify->restoreFromBackRemove(phaseToModify);
                    }

                }
            } else { // positional
                    sprintf(addRemoveText, "     %s => %s [ %s ]", " REPOS", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());
                    pUniverseToModify->rePositionLatticePointInPhase(phaseToModify);

                    tempTotalKLEnergy = calculateTotalKLEnergy();
                    // calculate connectivity potential
                    tempConnectivityEnergy = calculateTotalComponentEnergy(hlambda);
                    // calculate new volume
                    tempVolumeSum = calculateVolumeEnergy();
                    temp_volume_energy = muConstant * tempVolumeSum;
                    temp_link = calculateLinkPotential(1.0);

                    test_energy = tempTotalKLEnergy + temp_volume_energy + tempConnectivityEnergy + betaConstant*temp_link;

                    if (test_energy < current_energy) {

                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                        current_energy = test_energy;
                        currentTotalKLEnergy = tempTotalKLEnergy;
                        currentVolumeSum = tempVolumeSum;
                        current_link = temp_link;
                        currentConnectivityEnergy = tempConnectivityEnergy;
                        current_volume_energy = temp_volume_energy;
                    } else { // undo changes (rejecting)
                        pUniverseToModify->undoRePositionLatticePointInPhase(phaseToModify);
                    }

            }

            runtime = (std::clock() - startTime) / (double) CLOCKS_PER_SEC;
            std::cout << "*******************                                        *******************" << std::endl;
            printf("      TIME : %.5f (SECONDS) %s\n", runtime, addRemoveText);
            printf("             INVKBT => %.5f \n", inv_kb_temp);
            printf("           MAXSTEPS => %i (%4i)       [ %i ] \n", highTempRounds, high, lowCounter);
            printf("            E_GRAPH => %.4E E_LINK %.1f\n", currentConnectivityEnergy, current_link);
            printf("             VOLUME => %.4E  MU => %.4E \n", currentVolumeSum, muConstant);
            printf("              E_VOL => %.4E  \n", current_volume_energy);
            printf("               D_KL : %.4E TEST : %.4E ENRGY : %.4E \n", currentTotalKLEnergy, test_energy, current_energy);
            std::cout << "*******************                                        *******************" << std::endl;

            //pUniverseToModify->checkIndices();

            if (currentConnectivityEnergy == 0 && current_energy < lowest_energy) {
                // update D_kl
                muConstant = mu * currentTotalKLEnergy / ((1.0 - mu) * currentVolumeSum);

                current_link = this->calculateLinkPotential(1.0);

                betaConstant = (current_link == 0) ? 1 : (0.2*current_link/((1.0-0.2)*current_link));
                current_link_energy = betaConstant*current_link;

                // update current energy with new terms
                current_energy = currentTotalKLEnergy + muConstant * currentVolumeSum  + current_link_energy;

                lowest_energy = current_energy;

                this->makeCopyOfCurrentState();

                lowCounter++;
                if (lowCounter > 31){
                    goto finished;
                }
//                double temp = 1.0/inv_kb_temp*0.91;
//                inv_kb_temp = 1.0/temp;
//                inv_kb_temp *= 100/97.0;
            }
        }
    }

finished:
    // Housekeeping
    for (auto &entry : universes){
        auto const &name = entry.first;
        auto &uni = entry.second;
        uni.printSegments("highTemp");
        uni.deletePointers();
    }

}


void MultiDataComponent::anneal(int ccmultiple, float highTempStart){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::clock_t startTime;
    char addRemoveText[75];
    double runtime;

    bool isUpdated = false;
    float divideBy=0, acceptRate = 0.5, inv500 = 1.0/500.0;
    float inv500slash499 = 499.0/500.0;
    int failures=0, updated=0;

    Universe * pUniverseToModify;

    // re-initialize models in universe using lowest store state from createInitialModel
    for(auto &entry : universes) {
        entry.second.initializeFromStoredState();
    }

    int totalPhases = phases.size();
    std::string phaseToModify;
    std::vector<std::string> phases_vector(totalPhases);
    std::copy(phases.begin(), phases.end(), phases_vector.begin());
    std::uniform_int_distribution<int> phasesSet(0,totalPhases-1);

    // find the universe with the largest workingLimit
    int max_selected_indices = 0, lowerbound, upperbound, cwl;;

    double currentTemp = highTempStart;
    double inv_kb_temp = 1.0/highTempStart;

    for(auto &entry : universes) {
        auto &uni = entry.second;
        max_selected_indices = (uni.getTotalWorkingLimit() > max_selected_indices) ? uni.getTotalWorkingLimit() : max_selected_indices;
        uni.initializeContactsPotential();
    }

    float updateSteps = 7.0;
    int coupons = (max_selected_indices*std::log((double)max_selected_indices) + 0.5772156649*max_selected_indices + 0.5);
    int updateCount = ccmultiple*coupons;

    float step_limit = (updateCount < 10000) ? 27000.0 : (float)updateCount;
    int deadUpdate = std::ceil(step_limit/updateSteps);

    double tempTotalKLEnergy, currentTotalKLEnergy = calculateTotalKLEnergy();
    /*
     * Contacts Potential
     */
    double runningContactsSum = calculateContactPotential();
    double baseFactor = eta;
    double etaFactor = 1;
    double etaConstant = baseFactor*currentTotalKLEnergy/((1.0-baseFactor)*runningContactsSum);
    double tempTotalContactEnergy, current_totalContactEnergy = etaConstant*runningContactsSum;
    /*
     * Volume Potential (default mu : 0.57)
     * specific to the segment/phase
     */
    double tempVolumeSum, currentVolumeSum = this->calculateVolumeEnergy();
    double muConstant = mu*currentTotalKLEnergy/((1.0-mu)*currentVolumeSum);
    double temp_volume_energy, current_volume_energy = muConstant*currentVolumeSum;

    /*
     * Link Potential
     */
    double temp_link, current_link = this->calculateLinkPotential(1.0);
    double betaConstant = (current_link == 0) ? 1 : (0.2*current_link/((1.0-0.2)*current_link));
    double current_link_energy = betaConstant*current_link;

    double test_energy, current_energy = currentTotalKLEnergy + current_volume_energy + current_link_energy + current_totalContactEnergy;
    double lowest_energy = current_energy;

    int counter=1;
    int numberOfCoolingTempSteps, workingLimit;
    Segment * pSegment;

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::shuffle(phases_vector.begin(), phases_vector.end(), gen);
        startTime = std::clock();

        phaseToModify = phases_vector[phasesSet(gen)];
        pUniverseToModify = baseUniversesForEachPhase.find(phaseToModify)->second;
        pSegment = pUniverseToModify->getSegments(phaseToModify);

        std::cout << "*******************----------------------------------------*******************" << std::endl;

        if (distribution(gen) < percentAddRemove  ) { //add or remove bead within working Set (exclude deadzone)

            upperbound = pSegment->getUpperNumberOfBeads();
            lowerbound = pSegment->getLowerNumberOfBeads();
            cwl = pSegment->getWorkingLimit();

            if ( (distribution(gen) < 0.5 && (cwl < upperbound))){ // ADD BEAD?
//            if ( (distribution(gen) < 0.5)){ // ADD BEAD?
                std::sprintf(addRemoveText, "     %s => %s [ %s ] ", "  ADD ", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());
                pUniverseToModify->addLatticePointToPhase(phaseToModify);
                // get energy of new state
                tempTotalKLEnergy = calculateTotalKLEnergy();
                // calculate new volume
                tempVolumeSum = calculateVolumeEnergy();
                temp_volume_energy = muConstant * tempVolumeSum;
                temp_link = calculateLinkPotential(1.0);

                tempTotalContactEnergy = etaConstant*calculateContactPotential();

                test_energy = tempTotalKLEnergy + temp_volume_energy + betaConstant*temp_link + tempTotalContactEnergy;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentTotalKLEnergy = tempTotalKLEnergy;
                    currentVolumeSum = tempVolumeSum;
                    current_link = temp_link;
                    current_volume_energy = temp_volume_energy;
                    current_totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentTotalKLEnergy = tempTotalKLEnergy;
                    currentVolumeSum = tempVolumeSum;
                    current_link = temp_link;
                    current_volume_energy = temp_volume_energy;
                    current_totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else { // undo changes (rejecting)
                    pUniverseToModify->restoreFromBackupAdd(phaseToModify);
                }

            } else  { // remove

                std::sprintf(addRemoveText, "     %s => %s [ %s ] ", "REMOVE", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());
                pUniverseToModify->removeLatticePointFromPhaseConstrained(phaseToModify);
                tempTotalKLEnergy = calculateTotalKLEnergy();
                // calculate new volume
                tempVolumeSum = calculateVolumeEnergy();
                temp_volume_energy = muConstant * tempVolumeSum;
                temp_link = calculateLinkPotential(1.0);

                tempTotalContactEnergy = etaConstant*calculateContactPotential();

                test_energy = tempTotalKLEnergy + temp_volume_energy + betaConstant*temp_link + tempTotalContactEnergy;

                if (test_energy < current_energy) {
                    current_energy = test_energy;
                    currentTotalKLEnergy = tempTotalKLEnergy;
                    currentVolumeSum = tempVolumeSum;
                    current_link = temp_link;
                    current_volume_energy = temp_volume_energy;
                    current_totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                    current_energy = test_energy;
                    currentTotalKLEnergy = tempTotalKLEnergy;
                    currentVolumeSum = tempVolumeSum;
                    current_link = temp_link;
                    current_volume_energy = temp_volume_energy;
                    current_totalContactEnergy = tempTotalContactEnergy;
                    isUpdated = true;
                } else { // undo changes (rejecting)
                    pUniverseToModify->restoreFromBackRemove(phaseToModify);
                }
            }

        } else { // positional
            std::sprintf(addRemoveText, "     %s => %s [ %s ]", " REPOS", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());

            pUniverseToModify->rePositionAnyLatticePointInPhase(phaseToModify);

            tempTotalKLEnergy = calculateTotalKLEnergy();
            // calculate new volume
            tempVolumeSum = calculateVolumeEnergy();
            temp_volume_energy = muConstant * tempVolumeSum;
            temp_link = calculateLinkPotential(1.0);

            tempTotalContactEnergy = etaConstant*calculateContactPotential();

            test_energy = tempTotalKLEnergy + temp_volume_energy + betaConstant*temp_link + tempTotalContactEnergy;

            if (test_energy < current_energy) {
                current_energy = test_energy;
                currentTotalKLEnergy = tempTotalKLEnergy;
                currentVolumeSum = tempVolumeSum;
                current_link = temp_link;
                current_volume_energy = temp_volume_energy;
                current_totalContactEnergy = tempTotalContactEnergy;
                isUpdated = true;
            } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                current_energy = test_energy;
                currentTotalKLEnergy = tempTotalKLEnergy;
                currentVolumeSum = tempVolumeSum;
                current_link = temp_link;
                current_volume_energy = temp_volume_energy;
                current_totalContactEnergy = tempTotalContactEnergy;
                isUpdated = true;
            } else { // undo changes (rejecting)
                pUniverseToModify->undoRePositionLatticePointInPhase(phaseToModify);
            }
        }

       // pUniverseToModify->testContactsUpdating();
        runtime = (std::clock() - startTime) / (double) CLOCKS_PER_SEC;
        std::cout << "*******************                                        *******************" << std::endl;
        printf("            TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", currentTemp, step_limit, numberOfCoolingTempSteps);
        printf("          ACCEPT : %.5f  FAILURES => %i  \n", acceptRate, failures);
        printf("            TIME : %.5f (SECONDS) \n", runtime);
        printf("              WL : %s \n", addRemoveText);
        printf("             E_LINK => %.1f\n", current_link);
        printf("          E_CONTACT => %.4E  ETA => %.4E \n", current_totalContactEnergy, etaConstant);
        printf("           E_VOLUME => %.4E  ( %.4E ) MU => %.4E \n", current_volume_energy, currentVolumeSum, muConstant);
        printf("            D_KL : %.4E TEST : %.4E ENRGY : %.4E \n", currentTotalKLEnergy, test_energy, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // Adaptive simulated annealing part
        if (isUpdated){
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            divideBy += 1.0;
            updated++;
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, currentTemp, inv_kb_temp);

        if (counter % deadUpdate == 0 && acceptRate < 0.44){ // if counter too small, add/remove may not sample sufficiently
            // reset eta and mu based on lowestKL
            //baseFactor *= etaFactor;
            double tempCon = calculateContactPotential();
            etaConstant = baseFactor*currentTotalKLEnergy/((1.0-baseFactor)*tempCon);

            current_totalContactEnergy = etaConstant*tempCon;

            muConstant = mu*currentTotalKLEnergy/((1.0-mu)*currentVolumeSum);
            current_volume_energy = muConstant*currentVolumeSum;

            current_energy = currentTotalKLEnergy + current_volume_energy + betaConstant*current_link_energy + current_totalContactEnergy;
//            for (auto &entry : universes){
//                auto const &name = entry.first;
//                auto &uni = entry.second;
//                uni.printSegments("annealed_rd_"+std::to_string(numberOfCoolingTempSteps));
//            }
        }
        counter++;

    } // end of steps

/*
 *
 * constant positional annealing
 *
 */
    // reset eta and mu based on lowestKL
    //baseFactor *= 0.333332;
    double tempCon = calculateContactPotential();
    etaConstant = baseFactor*currentTotalKLEnergy/((1.0-baseFactor)*tempCon);
    current_totalContactEnergy = etaConstant*tempCon;

    muConstant = mu*currentTotalKLEnergy/((1.0-mu)*currentVolumeSum);
    current_volume_energy = muConstant*currentVolumeSum;

    current_energy = currentTotalKLEnergy + current_volume_energy + current_link_energy + current_totalContactEnergy;

    for (auto &entry : universes){
        //auto const &name = entry.first;
        auto &uni = entry.second;
        uni.printSegments("annealed_HT");
    }

    for(auto &entry : universes) {
        auto &uni = entry.second;
        max_selected_indices = (uni.getTotalWorkingLimit() > max_selected_indices) ? uni.getTotalWorkingLimit() : max_selected_indices;
    }
    coupons = (max_selected_indices*std::log((double)max_selected_indices) + 0.5772156649*max_selected_indices + 0.5);
    updateCount = 11*coupons;

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < updateCount; numberOfCoolingTempSteps++){

        std::shuffle(phases_vector.begin(), phases_vector.end(), gen);
        startTime = std::clock();

        phaseToModify = phases_vector[phasesSet(gen)];
        pUniverseToModify = baseUniversesForEachPhase.find(phaseToModify)->second;
        pSegment = pUniverseToModify->getSegments(phaseToModify);
        workingLimit = pSegment->getWorkingLimit();

        std::cout << "*******************-------------CONSTANT TEMP--------------*******************" << std::endl;


            std::sprintf(addRemoveText, "     %s => %s [ %s ]", " REPOS", phaseToModify.c_str(), pUniverseToModify->getNameOfUniverse().c_str());

            pUniverseToModify->rePositionAnyLatticePointInPhase(phaseToModify);

            tempTotalKLEnergy = calculateTotalKLEnergy();
            // calculate new volume
            tempVolumeSum = calculateVolumeEnergy();
            temp_volume_energy = muConstant * tempVolumeSum;
            temp_link = calculateLinkPotential(1.0);

            tempTotalContactEnergy = etaConstant*calculateContactPotential();

            test_energy = tempTotalKLEnergy + temp_volume_energy + betaConstant*temp_link + tempTotalContactEnergy;

            if (test_energy < current_energy) {
                current_energy = test_energy;
                currentTotalKLEnergy = tempTotalKLEnergy;
                currentVolumeSum = tempVolumeSum;
                current_link = temp_link;
                current_volume_energy = temp_volume_energy;
                current_totalContactEnergy = tempTotalContactEnergy;
                isUpdated = true;
            } else if (exp(-(test_energy - current_energy) * inv_kb_temp) > distribution(gen)) {
                current_energy = test_energy;
                currentTotalKLEnergy = tempTotalKLEnergy;
                currentVolumeSum = tempVolumeSum;
                current_link = temp_link;
                current_volume_energy = temp_volume_energy;
                current_totalContactEnergy = tempTotalContactEnergy;
                isUpdated = true;
            } else { // undo changes (rejecting)
                pUniverseToModify->undoRePositionLatticePointInPhase(phaseToModify);
            }

        runtime = (std::clock() - startTime) / (double) CLOCKS_PER_SEC;
        std::cout << "*******************                                        *******************" << std::endl;
        printf("            TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", currentTemp, numberOfCoolingTempSteps, updateCount);
        printf("            TIME : %.5f (SECONDS) \n", runtime);
        printf("              WL : %-6i  %s \n", workingLimit, addRemoveText);
        printf("             E_LINK => %.1f\n", current_link);
        printf("          E_CONTACT => %.4E  ETA => %.4E \n", current_totalContactEnergy, etaConstant);
        printf("           E_VOLUME => %.4E  ( %.4E ) MU => %.4E \n", current_volume_energy, currentVolumeSum, muConstant);
        printf("            D_KL : %.4E TEST : %.4E ENRGY : %.4E \n", currentTotalKLEnergy, test_energy, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;
    }

    this->printSummary();

    for (auto &entry : universes){
       // auto const &name = entry.first;
        auto &uni = entry.second;
        uni.printSegments("annealed");
        uni.deletePointers();
    }
}

bool MultiDataComponent::fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

/**
 * Check that all the belongs_to specified are represented in the Universes
 * @param file
 * @return
 */
bool MultiDataComponent::checkPhaseFormat() {

    std::ifstream phaseFileStream (setupFile);

    boost::regex prDataFormat("(prdata)|(PRDATA)");
    boost::regex arrowFormat("=>");

    boost::regex phaseFormat("(phase)|(PHASE)");
    boost::regex belongsToFormat("(BELONGS_TO)|(belongs_to)|(belongsto)|(BELONGSTO)");

    std::string line;
    std::vector<std::string> tempLine;

    std::vector<std::string> fileLines;
    if (phaseFileStream.is_open()) {
        while(!phaseFileStream.eof()) {
            std::getline(phaseFileStream, line);
            fileLines.push_back(line);
        }
        phaseFileStream.close();
    }

    for(unsigned int f=0; f<fileLines.size(); f++) {
        line = fileLines[f];

        // Does line contain PRDATA and does file exist and is NAME found in PHASE LINE?
        if (boost::regex_search(line, prDataFormat)) {
            int nameCount = 0;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
            int totalInSplit = tempLine.size();

            for(int i=0; i<totalInSplit; i++){

                if (boost::regex_match(tempLine[i], prDataFormat) && boost::regex_match(tempLine[i+1], arrowFormat) ){ // get file name

                    std::string tempFile = tempLine[i+2];

                    if (fileExists(tempFile)){
                        // get name and make sure Name is present in phase lines
                        // should only have one occurence of name
                        for(unsigned int j=0; j<tempLine.size(); j++){

                            if (tempLine[j] == "NAME" || tempLine[j] == "name"){
                                std::string tempName = tempLine[j+2];
                                bool notFound = true;
                                // make sure name is in phases
                                for(unsigned int ff=0; ff<fileLines.size(); ff++) {
                                    std::string tempPhaseline = fileLines[ff];

                                    if (boost::regex_search(tempPhaseline, phaseFormat)) { // if phase line then see if Name is present in belongs to
                                        std::vector<std::string> vectorOfWords;
                                        boost::split(vectorOfWords, tempPhaseline, boost::is_any_of("\t  "), boost::token_compress_on);
                                        for(unsigned int fff=0; fff<vectorOfWords.size(); fff++) { //
                                            if ((vectorOfWords[fff] == "BELONGS_TO" || vectorOfWords[fff] == "belongs_to") && vectorOfWords[fff+1] =="=>"){  // get name
                                                std::string belongs_to = vectorOfWords[fff+2];
                                                if (belongs_to == tempName){
                                                    notFound = false;
                                                    break;
                                                }

                                            }
                                        }
                                    }
                                }

                                if (notFound){
                                    throw std::invalid_argument( "**** ERROR => NAME NOT FOUND IN PHASES: \n\t" + line  + " \n");
                                }

                                nameCount++;
                            }
                        }

                        if (nameCount != 1){
                            throw std::invalid_argument( "**** ERROR => PRDATA CAN ONLY HAVE ONE NAME : \n\t" + line  + " \n");
                        }

                    } else {
                        throw std::invalid_argument( "**** ERROR => FILE CANNOT BE FOUND : \n\t" + tempFile  + " \n");
                    }
                }
            }
        }
    }

    // check that each belongs_to is found
    for(unsigned int f=0; f<fileLines.size(); f++) {

        line = fileLines[f];

        // Does line contain PHASE and belongs_to target found in PRDATA LINE?
        if (boost::regex_search(line, phaseFormat)) {

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
            int totalInSplit = tempLine.size();

            for(int i=0; i<totalInSplit; i++){

                if (boost::regex_match(tempLine[i], belongsToFormat) && boost::regex_match(tempLine[i+1], arrowFormat) ){ // get file name

                    std::string tempBelongsTo = tempLine[i+2]; // belongs_to name
                    bool notFound = true;

                    for(unsigned int j=0; j<fileLines.size(); j++){
                        std::string prLine = fileLines[j];

                        if (boost::regex_search(prLine, prDataFormat) ){
                            // get name
                            std::vector<std::string> vectorOfWords;
                            boost::split(vectorOfWords, prLine, boost::is_any_of("\t  "), boost::token_compress_on);
                            for(unsigned int fff=0; fff<vectorOfWords.size(); fff++) { //

                                if ((vectorOfWords[fff] == "NAME" || vectorOfWords[fff] == "name") && vectorOfWords[fff+1] =="=>"){  // get name
                                    std::string pr_belongs_to = vectorOfWords[fff+2];
                                    if (pr_belongs_to == tempBelongsTo){
                                        notFound = false;
                                        goto outLoop;
                                    }
                                }
                            }
                        }
                    }

                    outLoop:
                    if (notFound){
                        throw std::invalid_argument( "**** ERROR => BELONGS_TO NOT MATCHED : \n\t" + tempBelongsTo  + " \n");
                    }
                }
            }
        }
    }

    return true;
}

void MultiDataComponent::printSummary() {

    for(auto &entry : baseUniversesForEachPhase) {
        //const auto &phaseName = entry.first;
        auto &uni = entry.second;
        std::cout << "*******************----------------------------------------*******************" << std::endl;

        //totalVolumeEnergy += uni->getSegments(phaseName)->calculateCVXHULLVolume();
        FUNCTIONS_RPR::printInfo("  PRDATA : " + uni->getNameOfUniverse());
        char text[100];
        std::sprintf(text, "    D_KL : %.5f", uni->getKLEnergy());
        FUNCTIONS_RPR::printInfo(text);
        std::sprintf(text, " TOTAL N : %i", uni->getTotalWorkingLimit());
        FUNCTIONS_RPR::printInfo(text);
        uni->getNormalizedCVXVolume();
        std::sprintf(text, " CVX VOL : %.1f", uni->calculateCVXHULLVolume());
        FUNCTIONS_RPR::printInfo(text);
        std::sprintf(text, "  PER PT : %.2f", uni->getNormalizedVolume());
        FUNCTIONS_RPR::printInfo(text);


        std::map<std::string, Segment> * pSegs = uni->getPointerToSegments();
        for(auto &entry : *pSegs) {
            //auto &phaseName = entry.first;
            auto &segment = entry.second;
            std::sprintf(text, "       PHASE : %s", segment.getID().c_str());
            FUNCTIONS_RPR::printInfo(text);
            std::sprintf(text, "     TOTAL N : %i", segment.getWorkingLimit());
            FUNCTIONS_RPR::printInfo(text);
            std::sprintf(text, "   N CVX VOL : %.5f", segment.getNormalizedVolume());
            FUNCTIONS_RPR::printInfo(text);
            std::sprintf(text, " EULER TOURS : %i", segment.getEulerComponents());
            FUNCTIONS_RPR::printInfo(text);
            std::cout << "                ----------------------------------------" << std::endl;
        }
    }
}
