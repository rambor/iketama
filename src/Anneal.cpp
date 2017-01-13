//
// Created by Robert Rambo on 13/01/2017.
//

#include "Anneal.h"
#include "Model.h"

using namespace std;

Anneal::Anneal(float highT, float percent, int lowerV, int upperV, int highTempRounds, int contactsPerBead, string fileprefix, int totalSteps, float stepFactor, float eta, float lambda) {

    this->highT = highT;
    this->percentAddRemove = percent;
    this->lowerV = lowerV;
    this->upperV = upperV;
    this->highTempRounds = highTempRounds;
    this->stepsPerTemp = 5;
    this->contactsPerBead = contactsPerBead;
    this->totalTempStop = totalSteps;
    this->stepFactor = stepFactor;
    this->highTempStartForCooling = 0.00001; //0.00001
    this->filenameprefix = fileprefix;

    this->numberOfCoolingTempSteps = (int) ceil(log(lowTempStop/highTempStartForCooling)/(log(expSlowCoolConstant)));
    this->eta = eta;
    this->lambda = lambda;

    // expansionSlope = (stepsPerTemp - 5)/numberOfCoolingTempSteps;
    // printf("  SETTING EXPANSION SLOPE => %.2f (%.0f -> %.0f)\n", expansionSlope, highT, lowTempStop);
}




float Anneal::calculateAverageDistance(float * pDistance, int *stopAt, vector<int> *bead_indices, Model * pModel){

    int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    int row, row2, count = 0;
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


float Anneal::calculateCVXHULLVolume(char *flags, vector<int> *bead_indices, int upTo, coordT * points, Model *pModel) {

    // points should already be set to 3*upTo

//    int dim = 3;
//    char flagCVX[25];
//    sprintf(flagCVX, "qhull s FA");

    //coordT pointSet[3*upTo];

    for (int i=0; i<upTo; i++){
        beadToPoint(&points[i*3], pModel->getBead((*bead_indices)[i]));
        //  beadToPoint(&(pointSet[i*3]), pModel->getBead((*bead_indices)[i]));
    }

    qh_new_qhull(3, upTo, points, 0, flags, NULL, NULL);

    int volume_test = qh totvol;
    qh_freeqhull(true);
    return volume_test;
}


/*
 * go through each vertex of hull and move one unit further away
 */
void Anneal::enlargeDeadLimit(std::vector<int> &vertexIndices,
                              int totalV,
                              std::vector<int> &outsidePoints,
                              int outsideCounts,
                              std::vector<int> &bead_indices,
                              int workingLimit,
                              int *deadLimit,
                              int totalBeadsInSphere,
                              Model *pModel) {

    int * pIndexOfHullpt;
    vector<int>::iterator beginIt = bead_indices.begin(), it;
    vector<int>::iterator endIt = bead_indices.end();

    float sum_x=0.0, sum_y=0.0, sum_z=0.0, countf=0.0;
    Bead * currentBead;
    // calculate center of the hull
    for (int v = 0; v < totalV; v++) {
        currentBead = pModel->getBead(vertexIndices[v]);
        sum_x += currentBead->getX();
        sum_y += currentBead->getY();
        sum_z += currentBead->getZ();
        countf += 1.0;
    }

    float ave_x, ave_y, ave_z, distance_from_center, temp_distance, diff_x, diff_y, diff_z;
    // calculate center of HULL based on vertices
    float inv_count = 1.0/countf;
    ave_x = sum_x*inv_count;
    ave_y = sum_y*inv_count;
    ave_z = sum_z*inv_count;
    // end calculate center of hull

    vector<int> neighbors(60);
    int totalNeighbors, replaceWith, distanceTo;
    int modifiedWorkingLImit = workingLimit;

    // Search region between workingLimit and DeadLimit for a better point
    for (int v = 0; v < totalV; v++) { // SWAP VERTICES OF HULL TO BETTER SPOT

        pIndexOfHullpt = &vertexIndices[v];
        // find all contact neighboring beads
        // for through each one and accept first one that increase hull volume
        pModel->getNeighbors(*pIndexOfHullpt, neighbors, totalNeighbors);

        currentBead = pModel->getBead(*pIndexOfHullpt);

        diff_x = currentBead->getX()-ave_x;
        diff_y = currentBead->getY()-ave_y;
        diff_z = currentBead->getZ()-ave_z;
        distance_from_center = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z; // current distance from center of HULL

        replaceWith = *pIndexOfHullpt; // if nothing greater use same vertex

        // PARALLELIZE
        for (int i=0; i < totalNeighbors; i++){ // not all beads are
            currentBead = pModel->getBead(neighbors[i]);
            diff_x = currentBead->getX()-ave_x;
            diff_y = currentBead->getY()-ave_y;
            diff_z = currentBead->getZ()-ave_z;
            temp_distance = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z;

            if (temp_distance > distance_from_center){
                // add bead to workinglist
                distance_from_center = temp_distance;
                replaceWith = neighbors[i];
            }
        }

        it = find(beginIt, endIt, replaceWith);
        distanceTo = distance(beginIt, it);

        if ( distanceTo >= *deadLimit){ // outside feasible region
            std::iter_swap(beginIt + *deadLimit, it);
            std::iter_swap(beginIt + *deadLimit, beginIt + modifiedWorkingLImit);
            (*deadLimit)++;
            modifiedWorkingLImit++;
        } else if (distanceTo >= workingLimit) { // point is within Feasible region
            std::iter_swap(it, beginIt + modifiedWorkingLImit);
            modifiedWorkingLImit++;
        }
    }

    // add these points to region within workinglimit and deadlimit (working zone)
    // add points, increase working limit
    *deadLimit = recalculateDeadLimit(modifiedWorkingLImit, bead_indices, pModel, totalBeadsInSphere);
    // calculate convex HULL with this set of points
}

/**
 * has to be sorted
 * &nonSeed, nonSeedCount, pDistances, totalBeads, numberOfComponentsToBeModeled
 */
bool Anneal::isConnectedComponent(std::vector<int> *activeIndices, int availableWorkingLimit, float *pDistances,
                                  int totalBeads, int &numberOfComponents) {

    if (!std::is_sorted(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit)){
        cout << " NOT SORTED FROM isConnectedComponent " << endl;
        return false;
    }

    bool test=false;
    std::vector<int> tempIndices(availableWorkingLimit);

    std::copy(activeIndices->begin(), activeIndices->begin()+availableWorkingLimit, tempIndices.begin());
    std::sort(tempIndices.begin(), tempIndices.end());

    // int lastBeadIndex =tempIndices[availableWorkingLimit-1];
    //Graph graph(lastBeadIndex+1);  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    Graph graph;  // typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    boost::graph_traits<Graph>::edge_descriptor edge;
    bool flag;
    float distance;
    // boost::tie(edge, flag) = add_edge(0, 1, graph);  // tie is a reference based tuple
    // ds.union_set(0,1);
    // cout << " SIZE OF : " << availableWorkingLimit << endl;

    int count=0;

    // create edges within cutoff
    std::vector< Vertex> vertices(availableWorkingLimit); ////typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    for(int i=0; i < availableWorkingLimit; i++){
        vertices[i] = boost::add_vertex(graph);
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
            //distance = *(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead);
            if ((*(pDistances + (int)(firstBead*totalBeads - 0.5*firstBead*(firstBead+1)) - firstBead - 1 + secondBead)) <= interconnectivityCutOff){
                pTwo = &vertices[j];
                //cout << " j " << j << " " << secondBead << " " << distance << endl;
                //boost::tie(edge, flag) = boost::add_edge(firstBead, secondBead, graph);
                boost::tie(edge, flag) = boost::add_edge(*pOne, *pTwo, graph);
                //cout << " j " << j << " EDGE ADDED " << endl;
                //ds.union_set(firstBead, secondBead);
                count++;
            }

        }
    }

    std::vector<VertexIndex> rank(num_vertices(graph));
    std::vector<Vertex> parent(num_vertices(graph));

    boost::disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);

    initialize_incremental_components(graph, ds);
    incremental_components(graph, ds);

    boost::component_index<VertexIndex> components(parent.begin(), parent.end());

    // one single component should be
    // parent.size() - availableWorkingLimit + 1
    //cout << " Components Size " << components.size() << " Parent size " << parent.size() << " Working limit " << availableWorkingLimit << endl;
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
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within deadLimit
 */
void Anneal::populateLayeredDeadlimitUsingSet(std::vector<int> * bead_indices, std::set<int> * beads_in_use, const int workingLimit,
                                              int * pDeadLimit, Model * pModel) {
    *pDeadLimit = workingLimit;

    std::vector<int>::iterator it, itIndex, iteratorBeadIndices = bead_indices->begin(), endIt = bead_indices->end();
    std::set<int>::iterator endOfSet = beads_in_use->end();
    int neighbor;

    for (int i = 0; i < workingLimit; i++){ // excludes true model positions

        it = pModel->getPointerToNeighborhood((*bead_indices)[i]);

        for (int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is inside workinglimit, don't add
            neighbor = *(it+j);

            if ((neighbor > -1) && (beads_in_use->find(neighbor) == endOfSet)){ // if true, add neighbor to search space

                itIndex = std::find( (bead_indices->begin() + (*pDeadLimit)), endIt, neighbor);

                // if neighbor was already added from previous point, it will not be found
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


int Anneal::recalculateDeadLimit(int workingLimit, vector<int> &bead_indices, Model * pModel, int totalBeadsInSphere ){

    int deadLimit = totalBeadsInSphere;

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[25];
    sprintf(flags, "qhull s FA");

    coordT hullPoints2[3*workingLimit];

    // can be threaded
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints2[i*3], pModel->getBead(bead_indices[i]));
    }

    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;

    // UPDATE DEADZONE : move points not selected that are outside of hull to deadZone
    std::vector<int>::iterator beginIt;
    beginIt = bead_indices.begin();

    std::vector<int> inside(totalBeadsInSphere-workingLimit);
    std::vector<int> outside(totalBeadsInSphere-workingLimit);
    int insideCount=0, outsideCount=0;
    int * beadPosition;

    for(int i=workingLimit; i < deadLimit; i++){
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

    qh_freeqhull(true);

    copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
    deadLimit = workingLimit + insideCount;
    copy(outside.begin(), outside.begin()+outsideCount, beginIt + deadLimit);

    return deadLimit;
}


void Anneal::refineCVXHull(vector<int> &bead_indices,
                           vector<int> &active_indices,
                           int totalBeadsInSphere,
                           int workingLimit,
                           int *pDeadLimit,
                           Model *pModel){

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[25];
    sprintf(flags, "qhull s FA");

    coordT hullPoints[3*workingLimit];

    // create set to determine CVX Hull
    for (int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints[i*3], pModel->getBead(bead_indices[i]));
        active_indices[i] = bead_indices[i];
    }

    // needs to be optimized
    qh_new_qhull(3, workingLimit, hullPoints, 0, flags, NULL, NULL);
    vertexT * vertices = qh vertex_list;
    int totalV = qh num_vertices;

    pModel->setCVXHullVolume(qh totvol);

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
    copy(inside.begin(), inside.begin()+insideCount, beginIt + workingLimit);
    *pDeadLimit = workingLimit + insideCount;
    copy(outside.begin(), outside.begin()+outsideCount, beginIt + *pDeadLimit);

    // make copy of points from vertices and then free qh_hull
    vector<int> vertexIndices(totalV);

    for (int v = 0; v < totalV; v++) { // SWAP VERTICES OF HULL TO BETTER SPOT
        vertexIndices[v] = active_indices[qh_pointid(vertices->point)];
        vertices = vertices->next;
    }

    qh_freeqhull(true);
    enlargeDeadLimit(vertexIndices, totalV, outside, outsideCount, bead_indices, workingLimit, pDeadLimit, totalBeadsInSphere, pModel);
    //cout <<"TOTAL POINTS IN HULL : " << totalV  << endl;
}




