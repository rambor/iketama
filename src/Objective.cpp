//
// Created by Robert Rambo on 12/01/2017.
//

#include "Objective.h"

using namespace std;
using namespace boost::filesystem;

Objective::Objective(){

}

/**
 * returns a pointer to Data Object
 * holds the memory address to element in datasets vector
 *
 */
Data * Objective::getMainDataset() {

    Data * preturnMe;

    int phaseCount = 0;

    if (datasets.size() == 1){
        // get address
        preturnMe = &(datasets.begin()->second);
    } else {
        // iterate over each dataset
        // determine total number of phases in each dataset
        // identify dataset with largest number of phases
        for(auto iterator = datasets.begin(); iterator != datasets.end(); iterator++) {
            // getData object
            Data temp = iterator->second;
            if (temp.getTotalPhases() > phaseCount){
                preturnMe = &(iterator->second);
                phaseCount = temp.getTotalPhases();
            }
        }

    }

    return preturnMe;
}


/**
 *
 */
//void Objective::createWorkingSets(Model & model) {
//
//    int totalDatasets = datasets.size();
//
//    // generates different sets of q-values per dataset
//    for(int i=0; i<totalDatasets; i++){
//        this->getDataObjectByIndex(i)->creatingWorkingSet(model);
//    }
//
//}

void Objective::addDataObject(std::string iofqfile, std::string pofrfile) {
    // load file and create Data object
    // get name of iofqfile, strip away

    string line;
    path p = iofqfile;
    boost::regex slashes("/|\\\\");

    string identifierKey;
    vector<string> tempLine;

    try {
        if (exists(p)){
            //split name if forward or backslash is present
            if(boost::regex_search(iofqfile, slashes)){
                boost::split(tempLine, iofqfile, boost::is_any_of("/|\\\\"), boost::token_compress_on);
                int last = tempLine.size()-1;
                identifierKey = tempLine[last];
            } else {
                identifierKey = iofqfile;
            }

            // create IofQ object
            datasets.insert(pair<string, Data> (identifierKey, Data(iofqfile)));
            keys.push_back(identifierKey);

            // Add PofR Dataset
            Data * pdata;
            pdata = &datasets[identifierKey];
            (*pdata).addPofRData(pofrfile);
        }

    } catch (const filesystem_error& ex) {
        cout << ex.what() << '\n';
    }

    /*
    cout  <<  "  root_name()----------: " << p.root_name() << '\n';
    cout  <<  "  root_directory()-----: " << p.root_directory() << '\n';
    cout  <<  "  root_path()----------: " << p.root_path() << '\n';
    cout  <<  "  relative_path()------: " << p.relative_path() << '\n';
    cout  <<  "  parent_path()--------: " << p.parent_path() << '\n';
    cout  <<  "  filename()-----------: " << p.filename() << '\n';
    */
}


/**
 * return last dataset
 */
Data * Objective::getLastDataset(){
    return &datasets[keys[keys.size()-1]];
}