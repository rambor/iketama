//
// Created by Robert Rambo on 13/01/2017.
//

#include "Anneal.h"
#include "PDBModel.h"

//using namespace std;

Anneal::Anneal(float highT,
               float percent,
               int lowerV,
               int upperV,
               int highTempRounds,
               float contactsPerBead,
               std::string fileprefix,
               int totalSteps,
               float eta,
               float lambda,
               float alpha,
               float mu,
               int multiple,
               float accRate) {

    this->highT = highT;
    this->percentAddRemove = percent;
    this->lowerV = lowerV;
    this->upperV = upperV;
    this->highTempRounds = highTempRounds;
    this->stepsPerTemp = 5;
    this->contactsPerBead = contactsPerBead;
    this->totalTempStop = totalSteps;

    //this->highTempStartForCooling = 0.00001; //0.00001
    this->highTempStartForCooling = 0.0000016; //0.00001
    this->filenameprefix = fileprefix;

//    this->numberOfCoolingTempSteps = (int) ceil(log(lowTempStop/highTempStartForCooling)/(log(expSlowCoolConstant)));
    this->eta = eta;
    this->lambda = lambda;
    this->alpha = alpha;
    this->mu = mu;
    this->ccmultiple = multiple;

    // expansionSlope = (stepsPerTemp - 5)/numberOfCoolingTempSteps;
    // printf("  SETTING EXPANSION SLOPE => %.2f (%.0f -> %.0f)\n", expansionSlope, highT, lowTempStop);
    contactsDistribution.resize(13);
    contactsDistribution[0] = 0.0d;
    contactsDistribution[1] = 0.266030833333333;
    contactsDistribution[2] = 0.384529;
    contactsDistribution[3] = 0.2497085;
    contactsDistribution[4] = 0.05954545;
    contactsDistribution[5] = 0.02736574;
    contactsDistribution[6] = 0.00854701666666667;
    contactsDistribution[7] = 0.0042735;
    contactsDistribution[8] = 0.0000001;
    contactsDistribution[9] = 0.0d;
    contactsDistribution[10] = 0.0d;
    contactsDistribution[11] = 0.0d;
    contactsDistribution[12] = 0.0d;
    double totaltemp =0;
    for(int i=0; i<13; i++){
        totaltemp += contactsDistribution[i];
    }
    for(int i=1; i<9; i++){ // normalize
        contactsDistribution[i] = contactsDistribution[i]/totaltemp;
    }

    this->asaAcceptanceRate = accRate;
    complementASAAcceptanceRate = 1.0 - accRate;
    intASAAcceptanceRate = (int)1000*accRate;
    intComplementASAAcceptanceRate = (int)1000*complementASAAcceptanceRate;
}




float Anneal::calculateAverageDistance(float * pDistance, int *stopAt, std::vector<int> *bead_indices, Model * pModel){

    unsigned long int totalBeads = (unsigned long int)pModel->getTotalNumberOfBeadsInUniverse();

    int row, count=0;

    unsigned long int row2;
    float distanceSum=0.0;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    for(int m=0; m < *stopAt; m++){
        row = (*bead_indices)[m];
        row2 = row * (totalBeads) - row*(row+1)*0.5 - row - 1;

        // parallel
        for(int n=m+1; n < *stopAt; n++){
            distanceSum +=  *(pDistance + row2 + (*bead_indices)[n]);
            count++;
        }
    }
    return distanceSum/count;
}

//float Anneal::calculateCVXHULLVolume(char *flags, std::vector<int> *bead_indices, int upTo, Model *pModel) {
float Anneal::calculateCVXHULLVolume(char *flags, std::vector<int> *bead_indices, int upTo, Model *pModel) {

    // points should already be set to 3*upTo
//    int dim = 3;
//    char flagCVX[25];
//    sprintf(flagCVX, "qhull s FA");
    //coordT pointSet[3*upTo];
    int numpoints = 3*upTo;
    coordT points[numpoints];

    for (int i=0; i<upTo; i++){
        beadToPoint(&points[i*3], pModel->getBead((*bead_indices)[i]));
        //  beadToPoint(&(pointSet[i*3]), pModel->getBead((*bead_indices)[i]));
    }

    qh_new_qhull (3, upTo, points, 0, flags, NULL, NULL);
    int volume_test = qh totvol;
    //qh_freebuffers(qh);
    qh_freeqhull(true);

    return volume_test;
}


void Anneal::populateContactsDistribution(std::vector<double> & distribution, std::set<int> *beads_in_use, Model *pModel){
    std::set<int>::iterator it;
    for(int i=0; i<13; i++){
        distribution[i] =0;
    }

    for (it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {
        distribution[ numberOfContactsFromSet(beads_in_use, pModel, *it) ] += 1;
    }
}


/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSumPotential(std::set<int> *beads_in_use, Model *pModel){

    //calculate contacts per bead for selectedIndex
    double sum=0;

    std::set<int>::iterator it;
    for (it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {
        //int test = numberOfContactsFromSet(beads_in_use, pModel, *it);
        sum += totalContactsPotential(numberOfContactsFromSet(beads_in_use, pModel, *it));
    }

    return sum;
}


/**
 * go through each lattice point within working limit and determine total contact potential
 *
 */
double Anneal::calculateTotalContactSum(std::set<int> *beads_in_use, Model *pModel){

    //calculate contacts per bead for selectedIndex
    double sum=0;

    std::set<int>::iterator it;
    for (it = beads_in_use->begin(); it != beads_in_use->end(); ++it) {
        sum += (double)numberOfContactsFromSet(beads_in_use, pModel, *it);
    }

    return sum;
}


// need to return the selected index
void Anneal::createPlacesToCheck(int workingLimit,
                                          int average_number_of_contacts,
                                          int swap1,
                                          std::set<int> * beads_in_use,
                                          std::set<int> * returnMe,
                                          Model * pModel){

    std::vector<int>::iterator primaryNeighborhood, secondaryNeighborhood;
    int secondaryN;
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    // select first element from randomized active_set
    // check if available neighbor can be added
    primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);

    for (int n=0; n < totalNeighbors; n++){ // run into problem of check
        int neighbor = *(primaryNeighborhood + n);

        if (neighbor > -1){
            std::set<int>::iterator inSet = beads_in_use->find(neighbor);

            if (inSet == beads_in_use->end()){ // not in use, so its in deadlimit
                // if number of contacts at new position (less current) is 0, skip
                returnMe->insert(neighbor);

            } else { // its already a neighbor within WorkSet so check its neighborhood
                secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);

                for (int s=0; s < totalNeighbors; s++){
                    // adding to any of these already insures at least 1 contact
                    // check if secondary neighbor is in beads_in_use
                    secondaryN = *(secondaryNeighborhood + s);
                    inSet = beads_in_use->find(secondaryN);
                    if (secondaryN > -1 && inSet == beads_in_use->end()){ // if .end(), means not in use
                        returnMe->insert(secondaryN);
                    } else if (secondaryN == -1){
                        break;
                    }
                }
            }
        } else if (neighbor == -1) {
            break;
        }
    }


//    for (int i=0; i < workingLimit; i++) {
//
//        swap1 = *(activeIt + i);
//        //swap1 = *active_indices[i];
//
//        numberContacts = numberOfContactsFromSet(beads_in_use, pModel, swap1);
//        if (numberContacts < average_number_of_contacts){
//            // select first element from randomized active_set
//            // check if available neighbor can be added
//            primaryNeighborhood = pModel->getPointerToNeighborhood(swap1);
//
//            for (int n=0; n < totalNeighbors; n++){ // run into problem of check
//                int neighbor = *(primaryNeighborhood + n);
//
//                if (neighbor > -1){
//                    set<int>::iterator inSet = beads_in_use->find(neighbor);
//
//                    if (inSet == beads_in_use->end()){ // not in use, so its in deadlimit
//                        // if number of contacts at new position (less current) is 0, skip
//                        placesToCheck.insert(neighbor);
//
//                    } else { // its already a neighbor within WorkSet so check its neighborhood
//                        secondaryNeighborhood = pModel->getPointerToNeighborhood(neighbor);
//
//                        for (int s=0; s < totalNeighbors; s++){
//                            // adding to any of these already insures at least 1 contact
//                            // check if secondary neighbor is in beads_in_use
//                            secondaryN = *(secondaryNeighborhood + s);
//                            inSet = beads_in_use->find(secondaryN);
//                            if (secondaryN > -1 && inSet == beads_in_use->end()){ // if .end(), means not in use
//                                placesToCheck.insert(secondaryN);
//                            } else if (secondaryN == -1){
//                                break;
//                            }
//                        }
//                    }
//                } else if (neighbor == -1) {
//                    break;
//                }
//            }
//            break;
//        }
//    }
}

bool Anneal::checkForRepeats(std::vector<int> beads) {
    bool state = false;
    int beadSize = beads.size();

    std::set<int> testSet(beads.begin(), beads.end());
    std::cout << "______________________________________________________________________________" << std::endl;
    std::cout << "*******************                 TEST                   *******************" << std::endl;
    std::cout << "*******************              -----------               *******************" << std::endl;
    std::cout << " TEST SET " << testSet.size() << " vector set " << beadSize << std::endl;
    if (testSet.size() != beadSize){
        for(int i=0; i<10; i++){
            std::cout << "                 !!!!!!!!!!!!!!!!DUPLICATE ENTRIES FOUND! " << testSet.size() << " " << beadSize <<  std::endl;
        }
        state = true;
    }
    return state;
}

/*!
 *
 */
float Anneal::connectivityPotential(int numberOfComponents){

    switch(numberOfComponents) {
        case 1:
            return 0.f;
        case 2:
            return 10.f;
        case 3:
            return 1000.f;
        case 4:
            return 10000.f;
        case 5:
            return 100000.f;
        case 6:
            return 1000000.f;
        default:
            return 1000000.f*numberOfComponents;
    }
}

/*
 * go through each vertex of hull and move one unit further away
 *
 */
void Anneal::enlargeDeadLimit(std::vector<int> &bead_indices,
                              unsigned int *pDeadLimit,
                              Model *pModel) {

    // for each point upto Deadlimit, add its neighbor if not withinDeadLimit

    std::set<int> beads_in_use(bead_indices.begin(), bead_indices.begin() + *pDeadLimit);

    unsigned int total = *pDeadLimit;
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<int>::iterator it;

    for (int i=0; i < total; i++){
        // go through each bead, add neigbor if not found
        it = pModel->getPointerToNeighborhood(bead_indices[i]);

        for (int n=0; n< totalNeighbors; n++){
            int neighbor = *(it+n);
            if (beads_in_use.find(neighbor) == beads_in_use.end() && (neighbor > -1)){ // positions not in use should not be found
                // swap to deadlimit and increment
                std::vector<int>::iterator foundIt = std::find(bead_indices.begin() + *pDeadLimit, bead_indices.end(), neighbor);
                std::iter_swap(bead_indices.begin() + *pDeadLimit, foundIt);
                beads_in_use.insert(neighbor);
                (*pDeadLimit)++;
            } else if (neighbor == -1) {
                break;
            }
        }
    }
}

/**
 * has to be sorted
 * &nonSeed, nonSeedCount, pDistances, totalBeads, numberOfComponentsToBeModeled
 */
bool Anneal::isConnectedComponent(std::vector<int> *activeIndices, int availableWorkingLimit, float *pDistances,
                                  int totalBeads, int &numberOfComponents) {

//    if (!std::is_sorted(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit)){
//        cout << " NOT SORTED FROM isConnectedComponent " << endl;
//        return false;
//    }

    bool test=false;
    std::vector<int> tempIndices(availableWorkingLimit);

    std::copy(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit, tempIndices.begin());
//    std::sort(tempIndices.begin(), tempIndices.end());

    // int lastBeadIndex =tempIndices[availableWorkingLimit-1];
    //Graph graph(lastBeadIndex+1);  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    Graph graph;  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    boost::graph_traits<Graph>::edge_descriptor edge;
    bool flag;
    // float distance;
    // boost::tie(edge, flag) = add_edge(0, 1, graph);  // tie is a reference based tuple
    // ds.union_set(0,1);
    // cout << " SIZE OF : " << availableWorkingLimit << endl;

    //int count=0;

    // create edges within cutoff
    std::vector< Vertex> vertices(availableWorkingLimit); ////typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    for(int i=0; i < availableWorkingLimit; i++){
        vertices[i] = boost::add_vertex(graph); // adding vertices to the graph, function returns vertex descriptor
    }

    Vertex * pOne, * pTwo;

    // create edges between points that are within cutoff
    for(int i=0; i < availableWorkingLimit; i++){

        int firstBead = tempIndices[i];
        pOne = &vertices[i];
        int next = i+1;

        for(int j=next; j<availableWorkingLimit; j++){
            int secondBead = tempIndices[j];
            // is secondBead a neighbor of firstBead
            // distance = *(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead);
            if ((*(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead)) <= interconnectivityCutOff){
                pTwo = &vertices[j];
                //cout << " j " << j << " " << secondBead << " " << distance << endl;
                //boost::tie(edge, flag) = boost::add_edge(firstBead, secondBead, graph);
                boost::tie(edge, flag) = boost::add_edge(*pOne, *pTwo, graph); // returns pointer to new edge
                //cout << " j " << j << " EDGE ADDED " << endl;
                //ds.union_set(firstBead, secondBead);
                //count++;
            }

        }
    }

    std::vector<VertexIndex> rank(num_vertices(graph));
    std::vector<Vertex> parent(num_vertices(graph));

    boost::disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]); // Rank and Parent are pointers to Vertex and VertexIndex

    initialize_incremental_components(graph, ds); // graph of edges and empty ds
    incremental_components(graph, ds);

    boost::component_index<VertexIndex> components(parent.begin(), parent.end());

    // one single component should be
    numberOfComponents = components.size();
/*
    BOOST_FOREACH(VertexIndex current_index, components) {
                    std::cout << "component " << current_index << " contains: ";
                    // std::cout << boost::num_vertices(components[current_index].first) << endl;
                    // Iterate through the child vertex indices for [current_index]
                    BOOST_FOREACH(VertexIndex child_index,
                                  components[current_index]) {
                                    std::cout << child_index << " ";
                                }

                    std::cout << std::endl;
                }
*/
    // map index of vertex back to active indices via sorted temp indices
    // go through each compnent and redistribute to largest to make connected
    if (numberOfComponents == 1){
        test = true;
    }

    return test;
}

/**
 * contacts calculation must be performed on a sorted list and includes the beadIndex in the list
 */
int Anneal::numberOfContacts(int &beadIndex, std::vector<int> *bead_indices, unsigned int const &workingLimit, Model *pModel, float * pDistance){

    int count=0;
    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();
    // beadsInUse must be sorted
    int row;
    unsigned long int row2;

    int i=0;
    // add down column (column is beadIndex
    while ((*bead_indices)[i] < beadIndex && i < workingLimit){
        row = (*bead_indices)[i];
        row2 = row*totalBeads - (row*(row+1)*0.5) - row - 1;
        if (*(pDistance + row2 + beadIndex) < contactCutOff) {
            count++;
        }
        i++;
    }
    // out of first loop, assume bead_indices[i] == beadIndex at this point
    i++;
    // Add across row (constant row)
    row2 = beadIndex*totalBeads - beadIndex*(beadIndex+1)*0.5 - beadIndex - 1;
    while (i < workingLimit){
        if (*(pDistance + row2 + (*bead_indices)[i]) < contactCutOff) {
            count++;
        }
        i++;
    }
    return count;
}



int Anneal::getRandomNeighbor(int &locale, std::set<int> *beads_in_use, Model * pModel){

    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(locale);
    int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<int> temp(totalNeighbors);
    std::set<int>::iterator endOfSet = beads_in_use->end();

    int count=0; // possible beads to use
    for (int i=0; i< totalNeighbors; i++){
        int neighbor = *(it+i);
        if (beads_in_use->find(neighbor) == endOfSet && (neighbor > -1)){ // positions not in use should not be found
            temp[count] = neighbor;
            count++;
        } else if (neighbor == -1) {
            break;
        }
    }
    // if count == 0, no neighbors to use
    if (count == 0){
        return 0;
    }else {
        return (temp[rand()%count]);
    }
}




/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within deadLimit
 */
void Anneal::populateLayeredDeadlimit(std::vector<int>::iterator iteratorBeadIndices, const int workingLimit,
                                      int *pDeadLimit, Model *pModel, const int totalBeads) {
    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex;
    int distance;
    int neighbor;

    for (int i = 0; i<workingLimit; i++){

        it = pModel->getPointerToNeighborhood(*(iteratorBeadIndices + i));

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);

            itIndex = std::find(iteratorBeadIndices+(*pDeadLimit), iteratorBeadIndices + totalBeads, neighbor);
            distance = (itIndex - iteratorBeadIndices); // distance from beginning of vector
            // if not found, itIndex will report last

            if ((neighbor > -1) && (distance >= *pDeadLimit) && (distance < totalBeads)){
                std::iter_swap(iteratorBeadIndices + (*pDeadLimit), itIndex);
                (*pDeadLimit)++;
            } else if (neighbor == -1){
                break;
            }

//            if ( distance >= *pDeadLimit && (neighbor != -1) && (distance < totalBeads)) {
//                //cout << j << " SELECTED INDEX: " << *(iteratorBeadIndices + i) << " NEIGHHBOR: "  << neighbor << " " << distance << " DL " << *pDeadLimit << endl;
//                std::iter_swap(iteratorBeadIndices + (*pDeadLimit), itIndex);
//                (*pDeadLimit)++;
//            } else if (neighbor == -1) {
//                break;
//            }
        }
    }
}

/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within deadLimit
 */
void Anneal::populateLayeredDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, const int workingLimit,
                                              int * pDeadLimit, Model * pModel) {
    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex;
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int neighbor;

    for (int i = 0; i < workingLimit; i++){ // excludes true model positions

        it = pModel->getPointerToNeighborhood((*bead_indices)[i]);

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);

            if ((neighbor > -1) && (beads_in_use->find(neighbor) == endOfSet)){ // if true, add neighbor to search space

                itIndex = std::find( (bead_indices->begin() + (*pDeadLimit)), bead_indices->end(), neighbor);

                // if neighbor was already added from previous point, it will not be found
                if (itIndex != bead_indices->end()){
                    std::iter_swap( (bead_indices->begin() + (*pDeadLimit)), itIndex);
                    (*pDeadLimit)++;
                }
            } else if (neighbor == -1) {
                break;
            }
        }
    }
}


void Anneal::removeFromdDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, const int workingLimit,
                                                int * pDeadLimit, Model * pModel, int indexOfRemovedPosition) {

    std::vector<int>::iterator it, itIndex;
    int neighbor;

    it = pModel->getPointerToNeighborhood(indexOfRemovedPosition);

    for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
        // if neighbor is inside workinglimit, don't add
        neighbor = *(it+j);
        // check for each neighbor that it is within contact of workingLimit
        if (neighbor > -1 && numberOfContactsFromSet(beads_in_use, pModel, neighbor) == 0){

            itIndex = std::find( (bead_indices->begin() + workingLimit), bead_indices->end(), neighbor);

            if (std::distance(bead_indices->begin(), itIndex) < *pDeadLimit){
                *pDeadLimit = (*pDeadLimit) - 1;
                std::iter_swap( itIndex, (bead_indices->begin() + (*pDeadLimit)));
            }

        } else if (neighbor == -1) {
            break;
        }
    }

}

unsigned int Anneal::recalculateDeadLimit(unsigned int workingLimit, std::vector<int> &bead_indices, Model * pModel, int totalBeadsInSphere ){

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[] = "qhull FA";

    coordT hullPoints2[3*workingLimit];

    // calculate CVX Hull from selected indices
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints2[i*3], pModel->getBead(bead_indices[i]));
    }

    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, NULL, NULL);
//    vertexT * vertices = qh vertex_list;
//    int totalV = qh num_vertices;

    // UPDATE DEADZONE : move points not selected that are outside of hull to deadZone
    std::vector<int>::iterator foundIt, beginIt = bead_indices.begin();
    int * beadPosition;

    unsigned int deadLimit = workingLimit;
    for(int i = workingLimit; i < totalBeadsInSphere; i++){ // for beads no selected, are they inside or outside hull
        beadPosition = &bead_indices[i];
        beadToPoint(testPoint, pModel->getBead(*beadPosition));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (!isoutside){
            foundIt = std::find(beginIt + deadLimit, bead_indices.end(), *beadPosition);
            std::iter_swap(beginIt + deadLimit, foundIt);
            deadLimit++;
        }
    }

    qh_freeqhull(true);

    // repartition
//    std::copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
//    deadLimit = workingLimit + insideCount;
//    std::copy(outside.begin(), outside.begin()+outsideCount, beginIt + deadLimit);

    return deadLimit;
}


/*
 * only want non-used positions that are near component positions
 */
void Anneal::populatedDeadLimitExcludingSet(std::vector<int> * bead_indices,
                                            std::set<int> * beads_in_use,
                                            std::vector<int> * beads_in_component,
                                            int totalInComponent,
                                            const int workingLimit,
                                            int * pDeadLimit,
                                            Model * pModel){

    *pDeadLimit = workingLimit;
    std::vector<int>::iterator it, itIndex, endIt = bead_indices->end();
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int neighbor;

    for (int i = 0; i < totalInComponent; i++){ // excludes true model positions

        it = pModel->getPointerToNeighborhood((*beads_in_component)[i]);

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);
            // need to check if neighbor is found in use, if not add it to possible
            if ((neighbor > -1) && (beads_in_use->find(neighbor) == endOfSet)){ // if true, add neighbor to search space

                itIndex = std::find( (bead_indices->begin() + (*pDeadLimit)), endIt, neighbor);
                // if neighbor was already added from previous point in beads_in_component, it will not be found
                if (itIndex != endIt){
                    std::iter_swap( (bead_indices->begin() + (*pDeadLimit)), itIndex);
                    (*pDeadLimit)++;
                }
            } else if (neighbor == -1) {
                break;
            }
        }
    }
}

void Anneal::refineCVXHull(std::vector<int> &bead_indices,
                           std::vector<int> &active_indices,
                           int totalBeadsInSphere,
                           int workingLimit,
                           int *pDeadLimit,
                           Model *pModel){

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[] = "qhull FA";

    coordT hullPoints[3*workingLimit];

    // create set to determine CVX Hull
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints[i*3], pModel->getBead(bead_indices[i]));
        active_indices[i] = bead_indices[i]; // copy of the indices before any sorting up to workingLimit
    }

    // needs to be optimized
    qh_new_qhull(3, workingLimit, hullPoints, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;
    pModel->setCVXHullVolume(qh totvol);

    // needs to be optimized
    // UPDATE DEADZONE : move points not selected that are outside of hull to deadZone
    *pDeadLimit = totalBeadsInSphere;

    std::vector<int>::iterator beginIt;
    std::vector<int> inside(totalBeadsInSphere-workingLimit);
    std::vector<int> outside(totalBeadsInSphere-workingLimit);
    int insideCount=0, outsideCount=0;
    int * beadPosition;

    // sort unselected beads into those inside and outside CVX Hull
    for(int i=workingLimit; i < *pDeadLimit; i++){
        beadPosition = &bead_indices[i];
        beadToPoint(testPoint, pModel->getBead(*beadPosition));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (isoutside){
            outside[outsideCount] = *beadPosition;
            outsideCount++;
        } else {
            inside[insideCount] = *beadPosition;
            insideCount++;
        }
    }

    beginIt = bead_indices.begin();
    // copy sorted indices back to bead_indices
    std::copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
    *pDeadLimit = workingLimit + insideCount;
    std::copy(outside.begin(), outside.begin()+outsideCount, beginIt + *pDeadLimit);

    // make copy of points from vertices and then free qh_hull
//    std::vector<int> vertexIndices(totalV);
//    //int vertexCount=0;
//    for (int v = 0; v < totalV; v++) { // SWAP VERTICES OF HULL TO BETTER SPOT
//        vertexIndices[v] = active_indices[qh_pointid(vertices->point)];
//        vertices = vertices->next;
//    }
    qh_freeqhull(true);
    // return vertexIndices to see if they are at the edges
    //cout <<"TOTAL POINTS IN HULL : " << totalV  << endl;
}

/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASATemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=asaAcceptanceRate;

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate = asaAcceptanceRate+complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.15 && stepEval < 0.65){
        lamRate = asaAcceptanceRate;
    } else if (stepEval >= 0.65){
        lamRate = asaAcceptanceRate*pow(intASAAcceptanceRate, -(stepEval - 0.65)*2.857142857);
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
        changed=true;
    } else {
        temp = temp*1.001001001001;
        changed=true;
    }

    if (changed){
        inv_temp = 1.0/temp;
    }
}

/**
 * change acceptance rate to 0.1
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASALowTemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=0.01;

    if (stepEval < 0.15) { // if less than 15% into run
        lamRate = 0.01 + 0.990*pow(990, -stepEval*6.666667);
    } else if (stepEval >= 0.15 && stepEval < 0.65){
        lamRate = 0.01;
    } else if (stepEval >= 0.65){
        lamRate = 0.01*pow(10, -(stepEval - 0.65)*2.857142857);
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
        changed=true;
    } else {
        temp = temp*1.001001001001;
        changed=true;
    }

    if (changed){
        inv_temp = 1.0/temp;
    }
}

void Anneal::updateASAModTemp(int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

//    bool changed = false;
    double stepEval = (double)index/evalMax;
    double lamRate=0.27;
    double finalEval = 0.65;

    if (stepEval < 0.10) { // 80% acceptance rate
        lamRate = 0.80+0.20*pow(200, -stepEval*6.666667); // exponential increase to get to
    } else if (stepEval >= 0.10 && stepEval < 0.20){
        lamRate = 0.70; // 70% acceptance rate
    } else if (stepEval >= 0.20 && stepEval < 0.30){
        lamRate = 0.60; // 60% acceptance rate
    } else if (stepEval >= 0.30 && stepEval < 0.40){
        lamRate = 0.50; // 50% acceptance rate
    } else if (stepEval >= 0.40 && stepEval < 0.55){
        lamRate = 0.44; // 44% acceptance rate
    } else if (stepEval >= 0.55 && stepEval < finalEval){
        lamRate = 0.27; // 27% acceptance rate
    } else if (stepEval >= finalEval){ // exponential decay from 0.27
        lamRate = 0.27*pow(270, -(stepEval - finalEval)*2.857142857);
    }


    if (acceptRate > lamRate){
        temp = 0.999*temp;
        inv_temp = 1.0/temp;
//        changed=true;
    } else {
        temp = temp*1.001001001001;
        inv_temp = 1.0/temp;
//        changed=true;
    }

//    if (changed){
//        inv_temp = 1.0/temp;
//    }
}


/**
 * Diagnostic for calculating potential around each lattice in the model
 */
void Anneal::printContactList(std::vector<int> &bead_indices, std::set<int> * beads_in_use_tree, int workingLimit, Model * pModel){
    std::cout << "CONTACT LIST" << std::endl;
    for(int i=0; i<workingLimit; i++){
        int index = bead_indices[i];
        float value = calculateLocalContactPotentialPerBead(beads_in_use_tree, pModel, index);
        int cc = numberOfContactsFromSet(beads_in_use_tree, pModel, index);
        std::cout<< i << " " << index << " => " << cc <<  "  POTENTIAL => " << value << " <=> " << totalContactsPotential(cc) <<  std::endl;
    }
}


/**
 * Treat the connectivity potential as a look up table.
 *
 */
void Anneal::populatePotential(int totalNeighbors){

    connectivityPotentialTable.resize(totalNeighbors+1);

    for(int i=0; i<(totalNeighbors+1); i++){
        double diff = 1.0-exp(-alpha*(i-contactsPerBead));
//        double diff = exp(abs(i-contactsPerBead));
//        double diff = (i-contactsPerBead);
//         connectivityPotentialTable[i]=diff;
         connectivityPotentialTable[i]=1.0*(diff*diff);
//        connectivityPotentialTable[i] = pow(abs(i-contactsPerBead),3);
    }

     connectivityPotentialTable[0] *= 100.0;
     connectivityPotentialTable[1] *= 10;
//     connectivityPotentialTable[2] *= 7;
//     connectivityPotentialTable[3] *= 2;
//    connectivityPotentialTable[4] = 0;
//    connectivityPotentialTable[5] = 0;
//    connectivityPotentialTable[6]  *= 3;
//    connectivityPotentialTable[7]  *= 10;
//    connectivityPotentialTable[8]  *= 10;
//    connectivityPotentialTable[9]  *= 1000;//
//    connectivityPotentialTable[10] *= 1000;
//    connectivityPotentialTable[11] *= 1000;
//    connectivityPotentialTable[12] *= 1000;

//    connectivityPotentialTable[0] = 100.0;
//    connectivityPotentialTable[1] = 400;
//    connectivityPotentialTable[2] = 20;
//    connectivityPotentialTable[3] = 0;
//    connectivityPotentialTable[4] = 0;
//    connectivityPotentialTable[5] = 0;
//    connectivityPotentialTable[6]  = 100;
//    connectivityPotentialTable[7]  = 400;
//    connectivityPotentialTable[8]  = 1000;
//    connectivityPotentialTable[9]  = 1000;
//    connectivityPotentialTable[10] = 1000;
//    connectivityPotentialTable[11] = 1000;
//    connectivityPotentialTable[12] = 1000;
}



void Anneal::printContactsFromSet(std::vector<int> &bead_indices, int workingLimit, std::set<int> *beads_in_use,
                                           Model *pModel,
                                           int const selectedIndex){

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> tempVec(workingLimit);
    std::copy(bead_indices.begin(), bead_indices.begin()+workingLimit, tempVec.begin());
    std::shuffle(tempVec.begin(), tempVec.end(), gen);


    std::vector<int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    int neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int totalNeighbors = pModel->getSizeOfNeighborhood();

    for (int i=0; i< totalNeighbors; i++){

        int neighbor = *(it+i);

        if (beads_in_use->find(neighbor) != endOfSet){
            neighborContacts += 1;
        } else if (neighbor == -1) {
            break;
        }
    }

}



/**
 * modelPR and targetPR are the same size
 * targetPR is derived from PDB
 */
float Anneal::calculateKLDivergenceAgainstPDBPR(std::vector<unsigned int> &modelPR, std::vector<double> &targetPR){

    float totalCounts = 0.0;
    float kl=0.0;
    double prob, * value;

    int totalm = modelPR.size();
    std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    //last nonzero bin
    int last=0;
    for (int i=0; i<totalm; i++){
        if (targetPR[last] <= 0){
            break;
        }
        last++;
    }

    for (int i=0; i<totalm; i++){
        totalCounts += modelPR_float[i];
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    for (int i=0; i < last; i++){
        prob = targetPR[i];  // get experimental bin value
        value = &modelPR_float[i]; // get model value
        if (prob > 0 && *value > 0){
            kl += prob * log(prob/(*value) * totalCounts);
        } else if (prob > 0 && *value <= 0){ // severely penalize any model bin that is zero
            kl += 10000430;
        }
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    // severely penalize choices that make a bin zero for values < dmax
    /*
    int last = totalm-1;
    for (int i=0; i < last; i++){

        if (modelPR_float[i] <= 0){
            kl += 10000000000000000000;
        } else {
            prob = targetPR[i];  //
            kl += prob * log(prob/(modelPR_float[i])*totalCounts);
        }
    }

    // assume if last bin in dataset is 0, then 0*log0 = 0
    if (probability_per_bin[last] > 0){
        prob = targetPR[last];
        kl += prob * log(prob/(modelPR_float[last])*totalCounts);
        //kl += prob * log(prob/(modelPR_float[last])*totalCounts);
    }
     */

    return kl;  // returns value per bin
}

// convert PDB model into lattice model
void Anneal::readPDB(Model *pModel, std::vector<int> * keptBeads, std::string filename){

    PDBModel pdbModel(filename, true, true, pModel->getBeadRadius());

    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    Bead * currentBead;
    float diffx, diffy, diffz, bx, by, bz;
    int totalAtoms = pdbModel.getTotalAtoms();
    float beadradius = pModel->getBeadRadius();
    float b2 = beadradius*beadradius;

    for(int i=0; i<totalBeads; i++){

        currentBead = pModel->getBead(i);
        bx = currentBead->getX();
        by = currentBead->getY();
        bz = currentBead->getZ();

        for (int t=0; t<totalAtoms; t++){
            diffx = bx - *(pdbModel.getCenteredX()+t);
            diffy = by - *(pdbModel.getCenteredY()+t);
            diffz = bz - *(pdbModel.getCenteredZ()+t);

            if ((diffx*diffx + diffy*diffy + diffz*diffz) < b2){
                keptBeads->push_back(i);
                break;
            }
        }
    }

    int totalKept = keptBeads->size();
    //pModel->printSelectedBeads(0, totalKept, *keptBeads);
    pModel->writeSubModelToFile(0, totalKept, *keptBeads, "ideal");
}











bool Anneal::checkSetAndVector(int workingLimit, std::vector<int> * indices, std::set<int> * beads_in_use){

    // check all beadss in set
    for(std::set<int>::iterator it = beads_in_use->begin(); it != beads_in_use->end(); ++it){

        std::vector<int>::iterator vit = std::find(indices->begin(), indices->end(), *it);

        if (vit == indices->end()){
            std::cout << "************** ==> SET INDEX exceeds VECTOR " << *it << std::endl;
            exit(0);
            return true;
        }else{
            int dis = std::distance(indices->begin(), vit);

            if (dis > workingLimit){
                std::cout << "************** ==> SET INDEX exceeds workinglimit " << *it << std::endl;
                exit(0);
                return true;
            }
        }
    }


    for(int i=0;i<workingLimit; i++){
        int index = (*indices)[i];
        std::set<int>::iterator vit = beads_in_use->find(index);
        if (vit == beads_in_use->end()){
            std::cout << "************** ==> VECTOR INDEX NOT IN SET " << index << std::endl;
            exit(0);
            return true;
        }
    }

    return false;
}

void Anneal::printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<int> * wl){

    int total = accept->size();

    FILE * pFile;
    const char *outputFileName;
    std::string nameOf = "run_parameters.txt";
    outputFileName = nameOf.c_str() ;
    pFile = std::fopen(outputFileName, "w");

    std::string index;

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );

    for (int i=0; i < total; i++){
        //residue_index = std::to_string(selectedBeads[i]);
        //index = std::to_string(i + 1);
        fprintf(pFile, "%i %.5E %0.9f %0.9f %i\n", (i+1), (*accept)[i], (*temp)[i], (*divergence)[i], (*wl)[i] );
    }

    fclose(pFile);

}

