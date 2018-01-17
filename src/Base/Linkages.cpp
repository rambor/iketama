//
// Created by xos81802 on 26/09/2017.
//

#include "Linkages.h"

//Linkages::Linkages() {}

Linkages::Linkages(std::string baseUniverse, std::string parallelUniverse, Model * pParallelUniverseModel) : baseUniverse(baseUniverse), parallelUniverse(parallelUniverse), pParallelUniverseModel(pParallelUniverseModel) {
    totalBeadsInParallelUniverse = pParallelUniverseModel->getTotalNumberOfBeadsInUniverse();
    for(int i=0; i<totalBeadsInParallelUniverse; i++){
        parallelUniverseToBaseMappings.insert( std::pair<int, Links >( i, Links(i) ));
    }
}

/**
 * Parameters that are passed in are with respect to base universe
 * @param beadIndexInBase
 * @param beadradius
 * @param beadInBase
 */
void Linkages::createMappingOfBaseBead(int beadIndexInBase, float beadradius, Bead * beadInBase){

    baseToParallelUniverseMappings.insert( std::pair<int, Links >( beadIndexInBase, Links(beadIndexInBase) ));
    Links * pBaseLink = &baseToParallelUniverseMappings.find(beadIndexInBase)->second;
    vector3 diffvec;

    float limit = beadradius + pParallelUniverseModel->getBeadRadius();

    for(int i=0; i<totalBeadsInParallelUniverse; i++){ // find all beads in ParallelUniverse tha are within limit
        Bead * bead = pParallelUniverseModel->getBead(i);

        diffvec = beadInBase->getVec() - bead->getVec();

        if (diffvec.length() < limit){
            pBaseLink->addLinks(i);
            // create internal mapping
            parallelUniverseToBaseMappings.find(i)->second.addLinks(beadIndexInBase);
        }
    }
}

