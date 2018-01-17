//
// Created by xos81802 on 27/09/2017.
//

#ifndef IKETAMA_LINKS_H
#define IKETAMA_LINKS_H


#include <vector>

class Links {
    const int parentBeadIndex;
    int totalLinks;

    std::vector<int> mappings;

public:

    Links(int parentIndex) : parentBeadIndex(parentIndex) {};

    void addLinks(int index){
        mappings.push_back(index);
        totalLinks = mappings.size();
        mappings.resize(totalLinks);
    }

    int * getPointerToMapping(){return &mappings[0];}

    int getIndexByIndex(int index){ return mappings[index];}
    //void addLinks(int index);
    int getTotalLinks(){return totalLinks;}

};


#endif //IKETAMA_LINKS_H
