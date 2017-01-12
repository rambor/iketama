#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <regex>

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


int main() {
    string datFile;
    string multiphaseFile;
    string chemicalType;
    std::vector<string> pdbFiles;
    std::vector<string> datFiles;

    int upperV, lowerV;
    float bead_radius;

    int volume, totalSteps = 10000;
    float sigma;
    float stepFactor;

    int contactsPerBead;
    float highT, lowT;
    float percentAddRemove = 0.05, lambda = 0.0001;
    float eta = 0.00001; //
    int highTempRounds;
    int oversampling=31;
    int totalModels, totalPhasesForSeeded=1;

    string mode, prefix, ref, seedFile, anchorFile;
    bool mirror, superimpose = false, isSeeded = false, unconstrainedVolume, latticeResolution=false, fast=true, refine=false;

    po::options_description desc("\n  Usage: bubbles myRefinedScatterData.dat  \n To increase steps per temp, increase exponentialSlowCoolConstant to 0.90\n");

    desc.add_options()
            ("help,h", "Print help messages")
            ("dat", po::value<std::vector<std::string> >(&datFiles), "ScAtter Dat file from P(r) refinement")
            ("contacts,c", po::value<int>(&contactsPerBead)->default_value(5), " 0 < contacts < 6")
            ("totalModels,d", po::value<int>(&totalModels)->default_value(1))
            ("highTempForSearch", po::value<float>(&highT)->default_value(0.0005)) //0.00001
            ("lowTempForRefine", po::value<float>(&lowT)->default_value(0.0000002))
            ("totalCoolingSteps", po::value<int>(&totalSteps), "Default is 10000 steps")
            ("overSampling,q",  po::value<int>(&oversampling), "Default is 23")
            ("stepFactor,m", po::value<float>(&stepFactor)->default_value(1.045))
            ("highTempRounds,g", po::value<int>(&highTempRounds)->default_value(4700))
            ("percentAddRemove", po::value<float>(&percentAddRemove), "Sets probability of Add/Remove versus positional refinement")
            ("type,t", po::value<string>(&chemicalType)->default_value("protein"), "Protein or nucleic or both")
            ("vol,v", po::value<int>(&volume), "Volume of particle, default is read from dat file for single phase")
            ("sigma,s", po::value<float>(&sigma)->default_value(0.59), "Sigma for volume of protein or nucleic or both")
            ("components,a", po::value<string>(&multiphaseFile), "File describing mapping of components with dat files")
            ("sym,o", po::value<string>(&mode)->default_value("C1"), "Sets symmetry operator for single model run")
            ("prefix,x", po::value<string>(&prefix)->default_value("run"), "Name to tag output files")
            ("fast", po::value<bool>(&fast)->default_value(true), "Default radius is half binwidth, slow mode is radius 1/(2*root3)*binwidth")
            ("seed", po::value<string>(&seedFile), "Use specified input PDB as seed")
            ("totalPhasesForSeeded", po::value<int>(&totalPhasesForSeeded), "Total number of unique, interconnected phases to be modeled")
            ("refine", po::value<bool>(&refine), "Refine input PDB model, file is specified using seed flag")
            ("eta", po::value<float>(&eta)->default_value(eta), "compactness weight, default is 10^-6")
            ("lambda", po::value<float>(&lambda), "connectivity weight, default is 10^-4")
            ("anchor", po::value<string>(&anchorFile), "anchor ")
            ;

    // specify anneal parameter file
    po::positional_options_description positionalOptions;
    positionalOptions.add("dat", 2);
    //positionalOptions.add("pdbs", -1); // -1 refers to unlimited

    po::variables_map vm;

    cout << "Hello, World!" << endl;
    return 0;
}