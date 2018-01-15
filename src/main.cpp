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
//#include <iostream>
#include "Anneal.h"
#include "Annealer/MultiDataComponent.h"


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

//Phase & findPhaseByID(string id, vector<Phase *> phases){
//
//    int size = phases.size();
//    for(int i=0; i<size; i++){
//        if ((phases[i]->getID()).compare(id) == 0){
//            return *phases[i];
//        }
//    }
//}


int main(int argc, char** argv) {

    unsigned int hw=std::thread::hardware_concurrency();
    unsigned int hwConcurr= (hw != 0)? hw : 4;

    std::string datFile;
    std::string multiphaseFile;
    std::string chemicalType;
    std::vector<std::string> pdbFiles;
    std::vector<std::string> datFiles;

    int upperV=0, lowerV=0;
    float bead_radius;

    int volume, totalSteps = 10000;
    float sigma;

    float contactsPerBead, xaxis, yaxis, zaxis;
    float highT, alpha;
    float percentAddRemove = 0.1, lambda = 0.001;
    //float eta = 0.46, mu = 0.57; // 10^-4 to 10^-5 seems to be acceptable
    float acceptanceRate = 0.44, eta = 0.47, mu = 0.001; // 10^-4 to 10^-5 seems to be acceptable
    int highTempRounds;
    int multiple = 1;

    int totalModels, totalPhasesForSeeded=1;

    std::string sym, prefix, ref, seedFile, phaseFile, anchorFile;
    bool mirror, superimpose = false, isSeeded = false, unconstrainedVolume, latticeResolution=false, fast=true, refine=false;
    std::string descText =
                "\n          USAGE : iketama myRefinedScatterData_sx.dat myRefinedScatterData_pr.dat";
    descText += "\n FOR REFINEMENT : --refine 1 --seed infile.pdb -r 0.001 -u 3 -p 0.1 -x rfd";
    descText += "\n   FOR SYMMETRY : -o D2 -u 70 --mu 0.001\n";
    descText += "\n ADAPTIVE SIMULATED ANNEALING (ASA)";
    descText += "\n ====>  length of ASA is set using ccmultiple";
    descText += "\n ====>  acceptanceRate controls the starting temperature of the ASA run\n";
    descText += "\n CONSTRAINTS";
    descText += "\n ====>  contacts should be between 2 and 7, default is 4.1";
    descText += "\n ====>  strength of contacts constraint is set using eta";
    descText += "\n ====>  strength of volume constraint set by mu";
    descText += "\n ====>  high mu or eta flattens model\n\n";
    descText += "\n SEQUENTIAL RUNS IN BASH SHELL";
    descText += "\n ====>  for i in {0..10} ; do iketama file_sx.dat file_pr.dat -x run_${i} ; done\n\n";
    descText += "\n MODELS WITH SYMMETRY";
    descText += "\n ====>  increase mu to 0.1 using --mu 0.1 and set -u to a higher value";
    descText += "\n";
    po::options_description desc(descText);

    desc.add_options()
            ("help,h", "Print help messages")
            ("dat", po::value<std::vector<std::string> >(&datFiles), "ScAtter *.dat files from P(r) refinement (_sx and _pr)")
            ("contacts,c", po::value<float>(&contactsPerBead)->default_value(4.1), " 0 < contacts < 12")
            ("ccmultiple,u", po::value<int>(&multiple)->default_value(30), "Multiple of the Coupon Collector Sampling, increase u to increase number of steps in ASA")
            ("totalModels,d", po::value<int>(&totalModels)->default_value(1))
            ("highTempForSearch", po::value<float>(&highT)->default_value(0.0005)) //0.00001
            ("totalCoolingSteps", po::value<int>(&totalSteps), "Default is 10000 steps")
            ("highTempRounds,g", po::value<int>(&highTempRounds)->default_value(21357))
            ("percentAddRemove,p", po::value<float>(&percentAddRemove), "Sets probability of Add/Remove versus positional refinement")
            ("alpha", po::value<float>(&alpha)->default_value(0.031), "Morse Potential depth/curvature")
            ("type,t", po::value<string>(&chemicalType)->default_value("protein"), "Protein or nucleic or both")
            ("vol,v", po::value<int>(&volume), "Volume of particle, default is read from dat file for single phase")
            ("sigma,s", po::value<float>(&sigma)->default_value(0.59), "Sigma for volume of protein or nucleic or both")
            ("components,a", po::value<std::string>(&multiphaseFile), "File describing mapping of components with dat files")
            ("sym,o", po::value<std::string>(&sym)->default_value("C1"), "Sets symmetry operator for single model run")
            ("prefix,x", po::value<std::string>(&prefix)->default_value("run"), "Name to tag output files")
            ("seed", po::value<std::string>(&seedFile), "Use specified input PDB as seed")
            ("totalPhasesForSeeded", po::value<int>(&totalPhasesForSeeded), "Total number of unique, interconnected phases to be modeled")
            ("phaseFile", po::value<std::string>(&phaseFile), "File that maps data to components")
            ("refine,f", po::value<bool>(&refine), "Refine input PDB model, file is specified using seed flag")
            ("eta,e", po::value<float>(&eta)->default_value(eta), "compactness weight")
            ("mu,m", po::value<float>(&mu)->default_value(mu), "convex hull weight, percentage of KL divergence")
            ("lambda,l", po::value<float>(&lambda), "connectivity weight, default is 10^-2")
            ("anchor", po::value<string>(&anchorFile), "anchor ")
            ("acceptanceRate,r", po::value<float>(&acceptanceRate), "ASA acceptance rate, use a rate < 0.1 for refinement")
            ("xaxis", po::value<float>(&acceptanceRate), "length of x-axis, Angstroms")
            ("yaxis", po::value<float>(&acceptanceRate), "length of y-axis, Angstroms")
            ("zaxis", po::value<float>(&acceptanceRate), "length of z-axis, Angstroms")
            ;

    // specify anneal parameter file
    po::positional_options_description positionalOptions;
    positionalOptions.add("dat", 2);
    //positionalOptions.add("pdbs", -1); // -1 refers to unlimited

    po::variables_map vm;
    //std::vector<Phase *> phases;

    try {

        po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(), vm);

        // checking options list
        po::notify(vm);

        Objective minFunction;

        // Must be between 1 and 10
        if (contactsPerBead > 10 || contactsPerBead < 2){
            std::cout << "Incorrect contactsPerBead : " << contactsPerBead << "\nRe-Setting contacts per bead to 4" << std::endl;
            contactsPerBead = 4.1;
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
        if (vm.count("dat") == 0 ){

        }

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

            Data * myDataObject = minFunction.getDataObjectByIndex(0); // store address of the DataObject

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
            //phases.push_back(new Phase(volume, sigma, contrast));

            lowerV = (int)(volume - volume*sigma);
            upperV = (int)(volume + volume*0.05);

            // make association of dataset to phase
            //myDataObject->addPhase(*phases[0]);

        } else if (vm.count("components") == 1) {

            // read in file

            // validate multiphase file format
            try {
                if (fileExists(multiphaseFile)){
                    // check all the names and that they match belongs_to field and likewise all belongs_to have a defined dataset
                    // each component must have a belongs_to
                    MultiDataComponent segmentationMap = MultiDataComponent(
                            multiphaseFile,
                            highTempRounds,
                            contactsPerBead,
                            totalSteps,
                            eta,
                            lambda,
                            alpha,
                            mu,
                            percentAddRemove,
                            contactsPerBead
                    );

                    segmentationMap.createInitialModel();
                    segmentationMap.anneal(multiple, highT);

                } else { // throw excepction
                    throw std::invalid_argument( "**** ERROR => FILE CANNOT BE FOUND or READ : \n\t" + multiphaseFile  + " \n");
                }
            } catch (std::exception &err) {
                std::cerr<<"Caught "<<err.what()<< std::endl;
                std::cerr<<"Type "<<typeid(err).name()<< std::endl;
                exit(0);
            }

            return 1;

        } else if (vm.count("help")) {
            std::cout << desc << "\n";
            return SUCCESS;
        }

        // determine a suitable bead_radius for multicomponent
        // std::cout << " TOTAL COMPONENTS : " << phases.size() << endl;

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
        } else { // not in use, only fast mode available
            bead_radius = 0.5*1.0/(sqrt(2))*(mainDataset->getBinWidth());
            interconnectivityCutOff = bead_radius*2.8285;
            //alpha = 0.07;
            //contactsPerBead = 4.5;
            //bead_radius = 0.5*1.0/(sqrt(3))*(mainDataset->getBinWidth());
            //interconnectivityCutOff = bead_radius*3.464*1.001;
            //bead_radius = 0.5*1.0/(sqrt(2))*(mainDataset->getBinWidth());
        }

        std::cout << "       BEAD RADIUS : " << bead_radius << std::endl;
        std::cout << "         BIN WIDTH : " << mainDataset->getBinWidth() << " HALF => " << 0.5*mainDataset->getBinWidth() << std::endl;

        // bead_radius must be adjusted in the multi-phase refinement (so maybe regrid the model)
        // Model model(minFunction.getMainDataset()->getDmax(), bead_radius, mode);
        int searchSpace;

        if (isSeeded && !refine){
            // we have to use an enlarged search space since we don't know how the object is arranged relative to dmax
            searchSpace = minFunction.getMainDataset()->getDmax()*1.5;
        } else {
            searchSpace = minFunction.getMainDataset()->getDmax()*1.2;
        }

        std::cout << "              DMAX : " << searchSpace << std::endl;

        //Model model(searchSpace, bead_radius, fast, mode);
        Model model;

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
                          mu,
                          multiple,
                          acceptanceRate);

        mainAnneal.setInterconnectivityCutoff(interconnectivityCutOff);


        /**
         * Generate Initial Model
         * Typically a constant temp search looking for a single continuous set of lattice points
         *
         */
        if (std::regex_match(sym, std::regex("(C|D)[0-9]+")) && sym.compare("C1") != 0 ) { // SYMMETRY CONSTRAINT

            model = Model(searchSpace, bead_radius, fast, sym);

            std::cout << "      SYMMETRY SET : " << sym << std::endl;

            if ((refine && isSeeded)){
                /*
                 * align all models using a program like damaver -a *.pdb
                 * realign damfilt or damstart to reference model if necessary
                 * Easy check is to open damfilt against the reference, should superimpose
                 * use this damfilt/damstart model as input for refined.
                 */
                mainAnneal.populatePotential(model.getSizeOfNeighborhood());
                mainAnneal.initializeModelToRefineSym(&model, mainDataset, prefix, seedFile);
            } else {
                if (!mainAnneal.createInitialModelSymmetry(&model, mainDataset)){
                    return ERROR_UNHANDLED_EXCEPTION;
                }
            }


        } else if (isSeeded && !refine) { // input seed model, will build off of it

            model = Model(searchSpace, bead_radius, fast);

            if (fileExists(anchorFile)){
                // check that anchor points exists in seedFile
                std::cout << "*** READING ANCHOR FILE ***" << std::endl;
                mainAnneal.setAnchorPoints(anchorFile, seedFile, &model);
            }

            // make bead model from PDB
            // to thread this for multiple runs, the bead model must be made each time independently
            /*
             * initial model is created from the input PDB before searching
             */
            std::cout << "*** CREATING INITIAL MODEL FROM SEED ***" << std::endl;
            mainAnneal.populatePotential(model.getSizeOfNeighborhood());
            mainAnneal.createSeedFromPDB(&model, mainDataset, "reduced_seed", seedFile, totalPhasesForSeeded);

            bool redo = false;
            while (!redo){
                redo = mainAnneal.createInitialModelCVXHullSeeded(&model, mainDataset, "initialCVXSeeded");
            }

            exit(0);
            //mainAnneal.refineHomogenousBodyASACVXSeeded(&model, mainDataset, "refined");

            // create multiple models off of same refined SEED
            ThreadPool pool(hw);
            std::vector< std::future<std::string> > results;
            for(int i = 0; i < totalModels; ++i) {
                std::string nameTo = prefix + "_" + std::to_string(i+1);
                std::packaged_task<std::string(Model *, Data *, std::string)> task(std::bind(&Anneal::refineHomogenousBodyASACVXSeeded, &mainAnneal, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));

                results.emplace_back(pool.enqueue(
                        &Anneal::refineHomogenousBodyASACVXSeeded,
                        &mainAnneal,
                        &model,
                        mainDataset,
                        nameTo
                )
                );
            }

            for(auto && result: results)
                std::cout << "     FINISHED => " << result.get() << ' ' << endl;

            std::cout << std::endl;

            return SUCCESS;
        } else if (refine && isSeeded){ // refine input model with no symmetry
            /*
             * create lattice model from input PDB, can be proper PDB or bead model from another program such as DAMMIN
             * low temp refine that includes add/remove and positional
             */
            model = Model(searchSpace, bead_radius, fast);

            mainAnneal.populatePotential(model.getSizeOfNeighborhood());
            mainAnneal.initializeModelToRefine(&model, mainDataset, prefix, seedFile);
//            mainAnneal.reAssignLatticeModel(seedFile, &model, mainDataset);
        } else { // no symmetry

            model = Model(searchSpace, bead_radius, fast); // set model

            if (!mainAnneal.createInitialModelCVXHull(&model, mainDataset, "initialCVX")){
                return ERROR_UNHANDLED_EXCEPTION;
            }

//            mainAnneal.populatePotential(model.getSizeOfNeighborhood());
//            mainAnneal.refineHomogenousBodyMaxEntropyCVX(&model, mainDataset, prefix);
        }


        /**
         * REFINING INITIAL MODEL
         */
        std::cout << " REFINING STARTED  " << endl;

//        if (phases.size() > 1 && !refine){ // multiphase refinement
          if (0 > 1 && !refine){ // multiphase refinement
//            //
//            // refine the homogenous body before partitioning into phases
//            //
            //mainAnneal.createInitialModelCVXHullSeeded(&model, mainDataset, "initialCVXSeeded");
//            // mainAnneal.refineHomogenousBodyASA(&model, mainDataset, 1);
//            // partition in phases
//
//            for(int i=0; i<totalModels; i++){ // use resulting model as input for several simulated annealing runs
//                //void Anneal::refineMultiphase(Model * pModel, Objective * pObjective, vector<Phase *> * pPhases)
//                mainAnneal.refineMultiphase(&model, &minFunction, &phases);
//            }
//
        } else { // if no phase file, just refine the initial model

            /**
             * Refine the Initial Model with Symmetry contraints
             */
            if (std::regex_match(sym, std::regex("(C|D)[0-9]+")) && sym.compare("C1") != 0){

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
                                &Anneal::refineSymModel,
                                &mainAnneal,
                                &model,
                                mainDataset,
                                nameTo
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
                std::cout << " POPULATING POTENTIAL " << std::endl;
                // launch annealing in separate thread(s)
                if (totalModels == 1){
                    mainAnneal.refineHomogenousBodyASAHybrid(&model, mainDataset, prefix);
                   // mainAnneal.refineHomogenousBodyMaxEntropyCVX(&model, mainDataset, prefix);
                        //mainAnneal.refineHomogenousBodyASACVX(&model, mainDataset, prefix);
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

//        for (std::vector< Phase * >::iterator it = phases.begin() ; it != phases.end(); ++it) {
//            delete (*it);
//        }
//        phases.clear();

    } catch (boost::program_options::required_option& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (boost::program_options::error& e){

        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    }


    return 0;
}