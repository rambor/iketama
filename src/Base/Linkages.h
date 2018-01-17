//
// Created by xos81802 on 26/09/2017.
//

#ifndef IKETAMA_LINKAGES_H
#define IKETAMA_LINKAGES_H

#include <string>
#include <vector>
#include <algorithm>
#include <random>
//#include <set>
#include <utility>
#include <map>
#include "Bead.h"
#include "../Model.h"
#include "Links.h"


class Linkages {

    const std::string baseUniverse;
    const std::string parallelUniverse;

    int totalBeadsInParallelUniverse;

    Model * pParallelUniverseModel;
    std::map<int, Links> parallelUniverseToBaseMappings; // from internal indicee to base
    std::map<int, Links> baseToParallelUniverseMappings; // base index to internal indices of universe

public:

    Linkages(std::string baseUniverse, std::string parallelUniverse, Model * pParallelUniverseModel);

    void createMappingOfBaseBead(int beadIndexInBase, float beadradius, Bead * beadInBase);

    int getTotalLinksForBaseBead(int index){ return baseToParallelUniverseMappings.find(index)->second.getTotalLinks();}
    int getTotalLinksForParallelBead(int index){ return parallelUniverseToBaseMappings.find(index)->second.getTotalLinks();}

    //std::vector<int>::iterator getStartOfLinkVectorInBase(int beadIndex){ baseToParallelUniverseMappings[beadIndex].;}
    int * getPointerToLinkVectorInBase(int index){ return baseToParallelUniverseMappings.find(index)->second.getPointerToMapping(); }

    //
    //Links * getPointerToLinkInBase(int beadIndex){ return &baseToParallelUniverseMappings[beadIndex]; }
    //std::vector<int>::iterator getEndOfLinkVectorInBase(int beadIndex){ baseToParallelUniverseMappings[beadIndex].mappings.end();}

};


#endif //IKETAMA_LINKAGES_H
