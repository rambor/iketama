#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <regex>
#include "Base/Phase.h"
#include "Base/Data.h"
#include "Model.h"
#include "Objective.h"
#include <thread>
#include "ThreadPool.h"
#include <iostream>
#include "Anneal.h"


using namespace std;
namespace po = boost::program_options;

namespace {
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}


bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

Phase & findPhaseByID(string id, vector<Phase *> phases){

    int size = phases.size();
    for(int i=0; i<size; i++){
        if ((phases[i]->getID()).compare(id) == 0){
            return *phases[i];
        }
    }
}


int main(int argc, char** argv) {

    unsigned int hw=std::thread::hardware_concurrency();
    unsigned int hwConcurr= (hw != 0)? hw : 4;

    std::string datFile;
    std::string multiphaseFile;
    std::string chemicalType;
    std::vector<std::string> pdbFiles;
    std::vector<std::string> datFiles;

    int upperV, lowerV;
    float bead_radius;

    int volume, totalSteps = 10000;
    float sigma;

    float contactsPerBead;
    float highT, lowT, alpha;
    float percentAddRemove = 0.07, lambda = 0.01, decayRate = 100000;
    float eta = 0.000001; // 10^-4 to 10^-5 seems to be acceptable
    int highTempRounds;
    int multiple = 51;

    int totalModels, totalPhasesForSeeded=1;

    std::string mode, prefix, ref, seedFile, anchorFile;
    bool mirror, superimpose = false, isSeeded = false, unconstrainedVolume, latticeResolution=false, fast=true, refine=false;

    po::options_description desc("\n  Usage: bubbles myRefinedScatterData.dat  \n To increase steps per temp, increase exponentialSlowCoolConstant to 0.90\n");

    desc.add_options()
            ("help,h", "Print help messages")
            ("dat", po::value<std::vector<std::string> >(&datFiles), "ScAtter *.dat files from P(r) refinement (_sx and _pr)")
            ("contacts,c", po::value<float>(&contactsPerBead)->default_value(4), " 0 < contacts < 6")
            ("ccmultiple,u", po::value<int>(&multiple)->default_value(51), "Multiple of the Coupon Collector Sampling")
            ("totalModels,d", po::value<int>(&totalModels)->default_value(1))
            ("highTempForSearch", po::value<float>(&highT)->default_value(0.0005)) //0.00001
            ("lowTempForRefine", po::value<float>(&lowT)->default_value(0.0000002))
            ("totalCoolingSteps", po::value<int>(&totalSteps), "Default is 10000 steps")

            ("highTempRounds,g", po::value<int>(&highTempRounds)->default_value(6719))
            ("percentAddRemove", po::value<float>(&percentAddRemove), "Sets probability of Add/Remove versus positional refinement")
            ("alpha", po::value<float>(&alpha)->default_value(0.0173), "Morse potential depth")
            ("type,t", po::value<string>(&chemicalType)->default_value("protein"), "Protein or nucleic or both")
            ("vol,v", po::value<int>(&volume), "Volume of particle, default is read from dat file for single phase")
            ("sigma,s", po::value<float>(&sigma)->default_value(0.59), "Sigma for volume of protein or nucleic or both")
            ("components,a", po::value<std::string>(&multiphaseFile), "File describing mapping of components with dat files")
            ("sym,o", po::value<std::string>(&mode)->default_value("C1"), "Sets symmetry operator for single model run")
            ("prefix,x", po::value<std::string>(&prefix)->default_value("run"), "Name to tag output files")
            ("fast", po::value<bool>(&fast)->default_value(true), "Default radius is half binwidth, slow mode is radius 1/(2*root3)*binwidth")
            ("seed", po::value<std::string>(&seedFile), "Use specified input PDB as seed")
            ("totalPhasesForSeeded", po::value<int>(&totalPhasesForSeeded), "Total number of unique, interconnected phases to be modeled")
            ("refine", po::value<bool>(&refine), "Refine input PDB model, file is specified using seed flag")
            //("eta,e", po::value<float>(&eta)->default_value(eta), "compactness weight, default is 10^-6")
            ("decayRate,k", po::value<float>(&eta)->default_value(eta), "compactness weight, default is 10^-6")
            ("lambda,l", po::value<float>(&lambda), "connectivity weight, default is 10^-5")
            ("anchor", po::value<string>(&anchorFile), "anchor ")
            ;

    // specify anneal parameter file
    po::positional_options_description positionalOptions;
    positionalOptions.add("dat", 2);
    //positionalOptions.add("pdbs", -1); // -1 refers to unlimited

    po::variables_map vm;
    std::vector<Phase *> phases;

    try {

        po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(), vm);

        // checking options list
        po::notify(vm);

        Objective minFunction;

        // Must be between 1 and 10
        if (contactsPerBead > 40 || contactsPerBead < 2){
            cout << "Incorrect contactsPerBead : " << contactsPerBead << "\nRe-Setting contacts per bead to 4" << endl;
            contactsPerBead = 4;
        }

        /*
        cout << "Before test: " << vm.count("align") << endl;
        if ( vm.count("align") ){

            superimpose = true;

            pdbFiles = vm["pdbs"].as<vector<string> >();

            if (pdbFiles.size() == 0){
             cout << "Must specify files using -f option such as run_*.pdb" << endl;
             } else {
                cout << "Aligning files: " << endl;
               for (int i=0; i<pdbFiles.size(); i++){
                   cout << pdbFiles[i] << endl;
               }
            }
        }
       */
        // single phase system (only one iofq and pofr dat file)
        if (vm.count("dat") == 1) {

            datFiles = vm["dat"].as<vector<string> >();

            string pofrfile;
            string iofqfile;
            pofrfile = datFiles[1];
            iofqfile = datFiles[0];

            if (Data::checkPofRFile(datFiles[0])) {
                pofrfile = datFiles[0];
                iofqfile = datFiles[1];
            } else if (Data::checkPofRFile(datFiles[1])) {
                pofrfile = datFiles[1];
                iofqfile = datFiles[0];
            }

            if (vm.count("seed")){
                isSeeded = true;
                cout << "Using : " << seedFile << " as input PDB model " << endl;
            }


            cout << "Using : " << iofqfile << " as Intensity file " << endl;
            cout << "Using : " << pofrfile << " as PofR file " << endl;

            // fast mode, neeed to use different function for adding dataObject
            minFunction.addDataObject(iofqfile, pofrfile);

            Data *myDataObject = minFunction.getDataObjectByIndex(0); // store address of the DataObject

            // single phase system
            float contrast;
            if (chemicalType.compare("protein") == 0) {
                contrast = 1;
            } else if (chemicalType.compare("nucleic") == 0) {
                contrast = 1.31;
            } else {
                contrast = 1.5;
            }

            // get volume from data file
            volume = myDataObject->getVolume();
            // reset volume to command line value if specified
            if (vm.count("vol")) {
                volume = vm["vol"].as<float>();
            }

            if (vm.count("totalPhasesForSeeded") < 1){
                totalPhasesForSeeded == 1;
            }

            // create the object on the stack/heap?
            phases.push_back(new Phase(volume, sigma, contrast));

            lowerV = (int)(volume - volume*sigma);
            upperV = (int)(volume + volume*0.05);

            // make association of dataset to phase
            myDataObject->addPhase(*phases[0]);

        } else if (vm.count("components") == 1) {
            // read in file

            if (fileExists(multiphaseFile)){
                ifstream phaseFileStream (multiphaseFile);
                std::string line;
                std::vector<std::string> tempLine;
                Data * tempData;

                boost::regex modelFormat("(model)|(MODEL)");
                boost::regex volumeFormat("(volume)|(VOLUME)");
                boost::regex contrastFormat("(contrast)|(CONTRAST)");
                boost::regex sigmaFormat("(sigma)|(SIGMA)");
                boost::regex prfileFormat("(POFRFILE)|(pofrfile)");
                boost::regex belongsToFormat("(belongs_to)|(BELONGS_TO)");

                if (phaseFileStream.is_open()) {
                    while(!phaseFileStream.eof()) {
                        getline(phaseFileStream, line); //this function grabs a line and moves to next line in file

                        // if line contains MODEL, CONTRAST VOLUME then make new model object
                        // boost::to_upper(line);
                        std::vector<std::string>::iterator it, it2;
                        int tempVolume, index, index2;
                        float tempContrast, tempSigma;
                        std::string tempId, tempprfile, tempiofqfile;

                        if (boost::regex_search(line, modelFormat) && boost::regex_search(line, volumeFormat) && boost::regex_search(line, contrastFormat) && boost::regex_search(line, sigmaFormat)) {
                            // split the line and make new model object
                            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                            // where in vector is model?
                            it = std::find(tempLine.begin(), tempLine.end(), "MODEL");
                            index =  distance(tempLine.begin(),it);
                            if (it != tempLine.end()){
                                tempId = tempLine[index+1];
                            }

                            it = std::find(tempLine.begin(), tempLine.end(), "VOLUME");
                            index =  distance(tempLine.begin(),it);
                            if (it != tempLine.end()){
                                tempVolume = stoi(tempLine[index+1]);
                            }

                            it = std::find(tempLine.begin(), tempLine.end(), "CONTRAST");
                            index =  distance(tempLine.begin(),it);
                            if (it != tempLine.end()){
                                tempContrast =  stof(tempLine[index+1]);
                            }

                            it = std::find(tempLine.begin(), tempLine.end(), "SIGMA");
                            index =  distance(tempLine.begin(),it);
                            if (it != tempLine.end()){
                                tempSigma =  stof(tempLine[index+1]);
                            }

                            phases.push_back(new Phase(tempVolume, tempSigma, tempContrast, tempId));

                        } else if (boost::regex_search(line, prfileFormat) && boost::regex_search(line, belongsToFormat)) { // associate PRFILE with models
                            // boost::split(tempLine, line, boost::is_any_of("\t\s=>"), boost::token_compress_on);
                            boost::trim_if(line, boost::is_any_of("\t "));
                            boost::split(tempLine, line, boost::is_any_of("\t =>"), boost::token_compress_on);

                            // add DataObject to minFunction
                            it = std::find(tempLine.begin(), tempLine.end(), "POFRFILE");
                            it2 = std::find(tempLine.begin(), tempLine.end(), "IOFQFILE");
                            index =  distance(tempLine.begin(),it);
                            index2 =  distance(tempLine.begin(),it2);

                            if (it != tempLine.end() && it2 != tempLine.end()){
                                tempprfile = tempLine[index+1];
                                tempiofqfile = tempLine[index2+1];
                                cout << " READING DATA FILE " << endl;
                                cout << " IOFQFILE : " << tempiofqfile << endl;
                                cout << " POFRFILE : " << tempiofqfile << endl;
                                minFunction.addDataObject(tempiofqfile, tempprfile);
                                cout << "\n"<< endl;
                            }

                            // associate phases with DataObject
                            tempData = minFunction.getLastDataset();

                            it = std::find(tempLine.begin(), tempLine.end(), "WEIGHT");
                            index =  distance(tempLine.begin(),it);
                            float testValue = std::stof(tempLine[index+1]);

                            if (it != tempLine.end() && (testValue > 0)){
                                tempData->setWeight(testValue);
                            } else {
                                cout << "  \nCHECK COMPONENT PHASE FILE: WEIGHT MUST BE GREATER THAN ZERO\n" << endl;
                                return ERROR_IN_COMMAND_LINE;
                            }

                            // for each element of BELONGS_TO
                            it = std::find(tempLine.begin(), tempLine.end(), "BELONGS_TO");
                            index =  distance(tempLine.begin(),it);
                            if (it != tempLine.end()){
                                for(int i=(index+1); i<tempLine.size(); i++){
                                    // as long as tempLine[i] != AND
                                    if ( !(boost::iequals(tempLine[i], "AND")) && (tempLine[i].compare("!") != 0)){
                                        boost::trim_if(tempLine[i], boost::is_any_of("\t "));
                                        tempData->addPhase(findPhaseByID(tempLine[i], phases));
                                    } else if ((tempLine[i].compare("!") == 0)){ //  ignore comments marked by !
                                        break;
                                    }
                                }
                            }
                        } // end of line parsing if statements
                    }
                }

                for(int p=0; p<phases.size(); p++){
                    phases[p]->setNumberOfBeads(bead_radius);
                }

                cout << "\n Setting P(r) Dataset with most components as main dataset " << endl;
                cout << "\n Main P(r) Dataset for multi-component modeling => " << minFunction.getMainDataset()->getPofRFilename() << endl;

                volume = minFunction.getMainDataset()->getVolume();
                cout << " VOLUME : " << volume << endl;
                lowerV = (int)(volume - volume*sigma); // SIGMA IS FROM command line
                upperV = (int)(volume + volume*0.05);
                cout << " Initial VOLUME search constrainted by " << lowerV << " => " << upperV << endl;

            }

            // create model for each entry and then associate with PRFILE
            // need to determine volume of the entire object for initial model

        } else if (vm.count("help")) {
            std::cout << desc << "\n";
            return SUCCESS;
        }

        // determine a suitable bead_radius for multicomponent
        std::cout << " TOTAL COMPONENTS : " << phases.size() << endl;

        std::cout << "\n" << endl;
        float min_bin_width = 1000.0;
        for(int t=0; t < minFunction.getTotalDatasets(); t++){
            Data * pData = minFunction.getDataObjectByIndex(t);
            int totalPhasesTemp = pData->getTotalPhases();

            if (pData->getBinWidth() < min_bin_width) {
                min_bin_width = pData->getBinWidth();
            }

            // cout << "DATASET " <<  (t+1) << " " << pData->getVolume() << "\n TOTAL PHASES => " << pData->getTotalPhases() << " " << endl;
            for (int tt=0; tt<totalPhasesTemp; tt++){
                cout << "\t" << (tt+1) << " PHASE ID : " << pData->getPhase(tt)->getID() << " VOLUME => " <<  pData->getPhase(tt)->getVolume() << endl;
            }
        }

        Data *mainDataset = minFunction.getMainDataset();

        float interconnectivityCutOff;

        if (fast){
            bead_radius = 0.499999999999995*(mainDataset->getBinWidth());
            interconnectivityCutOff = bead_radius*2.001;
        } else {

            bead_radius = 0.5*1.0/(sqrt(2))*(mainDataset->getBinWidth());
            interconnectivityCutOff = bead_radius*2.8285;
            alpha = 0.07;
            contactsPerBead = 4.5;
            //bead_radius = 0.5*1.0/(sqrt(3))*(mainDataset->getBinWidth());
            //interconnectivityCutOff = bead_radius*3.464*1.001;
            //bead_radius = 0.5*1.0/(sqrt(2))*(mainDataset->getBinWidth());
        }

        cout << "       BEAD RADIUS : " << bead_radius << endl;
        cout << "         BIN WIDTH : " << mainDataset->getBinWidth() << " HALF => " << 0.5*mainDataset->getBinWidth() << endl;

        // bead_radius must be adjusted in the multi-phase refinement (so maybe regrid the model)
        // Model model(minFunction.getMainDataset()->getDmax(), bead_radius, mode);
        int searchSpace;

        if (isSeeded && !refine){
            // we have to use an enlarged search space since we don't know how the object is arranged relative to dmax
            searchSpace = minFunction.getMainDataset()->getDmax()*1.5;
        } else {
            searchSpace = minFunction.getMainDataset()->getDmax();
        }

        cout << "              DMAX : " << searchSpace << endl;
        Model model(searchSpace, bead_radius, fast, mode);

        if (fileExists(anchorFile) && fileExists(seedFile)){
            // check that anchor points exists in seedFile
            model.setAnchorPoints(anchorFile, seedFile);
            //return SUCCESS;
        }

        // create Instance of Annealer SYMMETRY or Not
        Anneal mainAnneal(highT,
                          percentAddRemove,
                          lowerV,
                          upperV,
                          highTempRounds,
                          contactsPerBead,
                          prefix,
                          totalSteps,
                          eta,
                          lambda,
                          alpha,
                          multiple);

        mainAnneal.setInterconnectivityCutoff(interconnectivityCutOff);
        // GENERATE INITIAL MODEL
        // Determine dataset that contains all the phases
        // A dataset object contains pointers to all the associated phases
        // minFunction contains all datasets in a vector
        //
        //
        //
        //
        if (std::regex_match(mode, std::regex("(C|D)[0-9]+")) && mode.compare("C1") != 0 ) { // if sym is set

            cout << "      SYMMETRY SET : " << mode << endl;
            if (!mainAnneal.createInitialModelSymmetry(&model, mainDataset)){
                return ERROR_UNHANDLED_EXCEPTION;
            }

        } else if (isSeeded && !refine) {
//
//            // make bead model from PDB
            cout << "*** CREATING INITIAL MODEL FROM SEED ***" << endl;
            mainAnneal.createSeedFromPDB(&model, mainDataset, "reduced_seed", seedFile, totalPhasesForSeeded);
            return SUCCESS;
//
        } else if (refine && isSeeded){
//            // create lattice model from input PDB, can be proper PDB or bead model from another program such as DAMMIN
//            mainAnneal.reAssignLatticeModel(seedFile, &model, mainDataset);
//
        } else { // no symmetry
            //mainAnneal.createInitialModelCVXHull(&model, mainDataset, "initialCVX");
            if (!mainAnneal.createInitialModelCVXHull(&model, mainDataset, "initialCVX")){
                return ERROR_UNHANDLED_EXCEPTION;
            }
        }


        // REFINE Initial model
        if (phases.size() > 1 && !refine){ // multiphase refinement
//            //
//            // refine the homogenous body before partitioning into phases
//            //
//            // mainAnneal.refineHomogenousBodyASA(&model, mainDataset, 1);
//            // partition in phases
//
//            for(int i=0; i<totalModels; i++){ // use resulting model as input for several simulated annealing runs
//                //void Anneal::refineMultiphase(Model * pModel, Objective * pObjective, vector<Phase *> * pPhases)
//                mainAnneal.refineMultiphase(&model, &minFunction, &phases);
//            }
//
        } else { // if no phase file, just refine the initial model

//            // PARALLELize the loops when CVX hull is rewritten
//            // refinement will be ran in a parallel loop for independent refinements
//            // #pragma omp parallel for
//            // for (int i=0; i<rounds; i++){
//            //     string nameOfOutfile = "fast_"+static_cast<ostringstream*>( &(ostringstream() << i) )->str();
//            //     mainAnneal.createInitialModel(&model, mainDataset, nameOfOutfile);
//            // }

//            // symmetry model
            if (std::regex_match(mode, std::regex("(C|D)[0-9]+")) && mode.compare("C1") != 0){

                mainAnneal.populatePotential(model.getSizeOfNeighborhood());

                if (totalModels == 1){
                    mainAnneal.refineSymModel(&model, mainDataset, prefix);
                } else { // distribute multiple models over many threads
                    ThreadPool pool(hw);
                    std::vector< std::future<std::string> > results;
                    for(int i = 0; i < totalModels; ++i) {
                        std::string nameTo = prefix + "_" + std::to_string(i+1);
                        std::packaged_task<std::string(Model *, Data *, std::string)> task(std::bind(&Anneal::refineSymModel, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));

                        results.emplace_back(pool.enqueue(
                                &Anneal::refineSymModel, &mainAnneal, &model, mainDataset, nameTo
                        )
                        );
                    }

                    for(auto && result: results)
                        std::cout << "     FINISHED => " << result.get() << ' ' << endl;
                    std::cout << std::endl;
                }

            } else {
//                //cout << "isSeeded " << isSeeded << endl;
//                if (isSeeded && !refine){
//                    // build back missing part if seeded, seeded is fixed set of lattice positions that are used throughout search
//                    for(int i=0; i < totalModels; i++) { // use resulting model as input
//                        mainAnneal.initialHighTempSearchSeeded(&model, mainDataset, "initial_seeded_model");
//                        mainAnneal.refineHomogoenousBodyFromPDBSeed(&model, mainDataset, i+1);
//                    }
//                } else {
//                    // if refining from input PDB, the input PDB sets initial search space, all lattice positions are refineable
//

                // define the tasks
//                std::packaged_task<std::string(Model *, Data *, std::string)> task(std::bind(&Anneal::refineHomogenousBodyASACVX, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
//                std::packaged_task<std::string(Model *, Data *, std::string)> task2(std::bind(&Anneal::refineHomogenousBodyASACVX, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
//                std::packaged_task<std::string(Model *, Data *, std::string)> task3(std::bind(&Anneal::refineHomogenousBodyASACVX, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
//                std::thread t1(std::move(task), &model, mainDataset, "one");
//                std::thread t2(std::move(task2), &model, mainDataset, "two");
//                std::thread t3(std::move(task3), &model, mainDataset, "three");
//                t1.join();
//                t2.join();
//                t3.join();

                mainAnneal.populatePotential(model.getSizeOfNeighborhood());
                // launch annealing in separate thread(s)
                if (totalModels == 1){
                        mainAnneal.refineHomogenousBodyASACVX(&model, mainDataset, prefix);
                } else {

                    ThreadPool pool(hw);
                    std::vector< std::future<std::string> > results;
                    for(int i = 0; i < totalModels; ++i) {
                        std::string nameTo = prefix + "_" + std::to_string(i+1);
                        std::packaged_task<std::string(Model *, Data *, std::string)> task(std::bind(&Anneal::refineHomogenousBodyASACVX, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));

                        results.emplace_back(pool.enqueue(
                                &Anneal::refineHomogenousBodyASACVX, &mainAnneal, &model, mainDataset, nameTo
                        )
                        );
                    }

                    for(auto && result: results)
                        std::cout << "     FINISHED => " << result.get() << ' ' << endl;
                    std::cout << std::endl;

                }

//                }
            }
        }


    } catch (boost::program_options::required_option& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (boost::program_options::error& e){

        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    }


    return 0;
}