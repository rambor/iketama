//
// Created by Robert Rambo on 12/01/2017.
//

#include "PDBModel.h"

using namespace std;

// constructor
PDBModel::PDBModel(string file, bool rna, bool discardWaters, float lower) {

    filename = file;
    ifRNA = rna;
    ifWaters = false;
    watersPerResidue = 30;

    string line, tempResi;
    // open file and read in contents
    ifstream scxFile (filename.c_str());
    int fileLength = 0;
    volume=0.0;
    cutoff = lower*lower;

    atomType.reserve(5000);
    chainID.reserve(5000);
    resID.reserve(5000);

    waterLines.reserve(500*watersPerResidue);

    x.reserve(5000);
    y.reserve(5000);
    z.reserve(5000);

    boost::regex pdbStart("ATOM");
    boost::regex hetatm("HETATM");
    boost::regex wat("HOH");
    boost::regex pdbX("-*[0-9]+.[0-9]+");
    boost::regex ifHydrogen("H[A-Z]+");

    cout << "\tReading PDB File : " << filename << endl;

    if (scxFile.is_open()) {
        while(!scxFile.eof()) {
            getline(scxFile, line); //this function grabs a line and moves to next line in file
            // string::substring(position,length)
            // Check if line starts with ATOM or HETATM and exclude hydrogens

            if ((line.length() > 0 && (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)) && !boost::regex_search(line.substr(17,3), wat)) && boost::regex_search(line.substr(31,8),pdbX) && !boost::regex_search(line.substr(12,4), ifHydrogen)) {

                x.push_back(strtod((char*)line.substr(30,8).c_str(), NULL));
                y.push_back(strtod((char*)line.substr(38,8).c_str(), NULL));
                z.push_back(strtod((char*)line.substr(46,8).c_str(), NULL));

                //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType.push_back(line.substr(12,4));//

                // residue name, three letter for protein, two for nucleic acids
                tempResi = line.substr(17,3);

                // reassign residue abbreviations for RNA
                if (ifRNA){  // A => ALA, G => GLY, C => CYS
                    if ((tempResi == "  A") || (tempResi == "ADE") || (tempResi == " rA") || (tempResi == "A  ")){
                        tempResi = " rA";
                    } else if ((tempResi == "  G") || (tempResi == "GUA") || (tempResi == " rG") || (tempResi == "G  ")) {
                        tempResi = " rG";
                    } else if ((tempResi == "  U") || (tempResi == "URI") || (tempResi == " rU") || (tempResi == "U  ")) {
                        tempResi = " rU";
                    } else if ((tempResi == "  C") || (tempResi == "CYT") || (tempResi == " rC") || (tempResi == "C  ")) {
                        tempResi = " rC";
                    }
                }

                resi.push_back(tempResi);             // residue name Protein (ALA, GLY, ...), RNA (rA, rG, rU, rC), DNA (dA, dG, dU, rC)
                resID.push_back( atoi( line.substr(22,4).c_str()) ); // residue sequence number
                chainID.push_back(line.substr(21,1));

                atomVolume.push_back(FUNCTIONS_RPR::residueToVolume( atomType.back(), tempResi ));

                volume += atomVolume.back(); //FUNCTIONS_RPR::residueToVolume( line.substr(17,3), atomType.back() );
                fileLength++;

            } else if (!discardWaters && line.length() > 20 && (line.substr(17,3) == "HOH")) {
                waterLines.push_back(line);
            }
            // keep HETATM flag - say you have a heme?
            // WATERS r in lines containing HETATM
            // if include waters is set, must
        }

        totalAtoms = fileLength;
        cout << "\tTotal atoms: " << totalAtoms << endl;
    }

    scxFile.close();
    scxFile.clear();

    trimmedAtomType = atomType;
    trimmedResi = resi;

    // center coordinates and convert to r, theta, phi
    rThetaPhiAtomType = new float[totalAtoms*5];

    atomicNumbers = new int[totalAtoms];
    atomicExpVolume = new float[totalAtoms];

    centeredX = new float[totalAtoms];
    centeredY = new float[totalAtoms];
    centeredZ = new float[totalAtoms];

    // send in x,y,z coordintes return centered coordinates
    FUNCTIONS_RPR::dmaxFromPDB(x,y,z, &dmax, centeredX, centeredY, centeredZ, totalAtoms);

    int locale;
    int atomicNumber;

    for (int n=0; n < totalAtoms; n++) {
        /*
         Convert to spherical coordinates
         xyz_to_rtp(x[],y[],z[])
         pass in &address to the array elements and dereference in the function to get value
         */
        // float * xyz_to_rtp(const float &tempx, const float &tempy, const float &tempz);


        pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp((const float &) centeredX[n], (const float &) centeredY[n], (const float &) centeredZ[n]);

        atomicNumber = convertToAtomicNumber(atomType[n], resi[n]);

        locale = n*5;
        rThetaPhiAtomType[locale] = *pconvertXYZ;               // [0] r
        rThetaPhiAtomType[locale+1] = cosf(*(pconvertXYZ+1));   // [1] cos(theta)
        rThetaPhiAtomType[locale+2] = *(pconvertXYZ+2);         // [2] phi
        rThetaPhiAtomType[locale+3] = atomicNumber;             // [3] atomic number
        rThetaPhiAtomType[locale+4] = 1.0;                      // [4]

        atomList.insert( atomicNumber );                        // Set of unique atomic numbers for calculating f(q)
        atomicNumbers[n] = atomicNumber;

        atomicExpVolume[n] = cbrt( atomVolume[n]*atomVolume[n] )*invFourPI; // volume^(2/3)*1/(4*PI)

        // call trim on resi names and atom types
        boost::algorithm::trim (trimmedAtomType[n]);
        boost::algorithm::trim (trimmedResi[n]);
    }
    // atomList is by atomic number, but water will be a molecule, give it a special number of 99
    atomList.insert(99); // Add unique for water.
    totalResidues = resID.size();
    // maximum theoretical amount
    keptWaters = new Coords[watersPerResidue*totalResidues];
    waterCount = 0;

    if (!discardWaters && (waterLines.size() > 0)){
        waterCount = waterLines.size();
        float invwaterCount = 1.0/waterCount;
        float waterX=0, waterY=0, waterZ=0, aveX, aveY, aveZ;

        for(int w=0; w < waterCount; w++){
            keptWaters[w].x = strtod((char*)waterLines[w].substr(30,8).c_str(), NULL);
            keptWaters[w].y = strtod((char*)waterLines[w].substr(38,8).c_str(), NULL);
            keptWaters[w].z = strtod((char*)waterLines[w].substr(46,8).c_str(), NULL);
            keptWaters[w].type = "O";
            keptWaters[w].occ = 1.0;
            waterX += keptWaters[w].x;
            waterY += keptWaters[w].y;
            waterZ += keptWaters[w].z;
        }

        aveX= waterX*invwaterCount;
        aveY= waterY*invwaterCount;
        aveZ= waterZ*invwaterCount;

        for(int w=0; w < waterCount; w++){
            keptWaters[w].x += -aveX;
            keptWaters[w].y += -aveY;
            keptWaters[w].z += -aveZ;
        }

        Coords * waterCoordinate;
        rWater = new float[waterCount];
        phiWater = new float[waterCount];
        waterCosTheta = new float[waterCount];

        for(int rt=0; rt<waterCount; rt++){
            waterCoordinate = &keptWaters[rt];
            //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp((const float &) centeredX[n], (const float &) centeredY[n], (const float &) centeredZ[n]);
            //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp(&(*waterCoordinate).x, &(*waterCoordinate).y, &(*waterCoordinate).z);
            pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp((const float &) (*waterCoordinate).x, (const float &) (*waterCoordinate).y, (const float &) (*waterCoordinate).z);
            rWater[rt] = *pconvertXYZ;
            waterCosTheta[rt] = cosf(*(pconvertXYZ+1));
            phiWater[rt] = *(pconvertXYZ+2);
        }

        ifWaters = true;
        cout << "\tUsing => " << waterCount << " waters from PDB" << endl;
    }
} // end of constructor

void PDBModel::writeKeptWaters() const {

    int count = 1;
    for(int i=0; i< waterCount; i++){
        if (keptWaters[i].type.size() == 0){
            break;
        }

        printf("%-3s%7i%4s%5s%2s%4i     %7.3f %7.3f %7.3f  %4.2f %4.2f\n", "ATOM", count, "O", "HOH", "A", count, keptWaters[i].x, keptWaters[i].y, keptWaters[i].z, keptWaters[i].occ, 100.0 );
        count++;
    }
}

void PDBModel::writeCenteredCoordinatesToFile(string name) const {

    string residue_index;
    FILE * pFile;

    const char *outputFileName;
    name = name + ".pdb";
    outputFileName = name.c_str() ;
    const char *originalPDBFilename;
    originalPDBFilename = this->filename.c_str();
    pFile = fopen(outputFileName, "w");

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", originalPDBFilename);
    for (int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<string>(resID[n]);
        fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}


void PDBModel::writeCenteredCoordinates() const {

    string residue_index;

    for (int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<string>(resID[n]);
        printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
    }
}
// for each residue, calculate waters

//void PDBModel::addWaters(Waters & waterModel) {
//
//    int resid, j, atomsCount;
//    string residue;
//
//    vector<string> atomsInResidue;
//    //string * atomsInResidue = NULL;
//    vector<float> res_x;
//    vector<float> res_y;
//    vector<float> res_z;
//
//    float neigh_x [1200]; // this is a guess at the number of atoms that will be within 15 angstroms of a water
//    float neigh_y [1200]; // if too small, it will cause a SEGFAULT
//    float neigh_z [1200];
//
//    bool keepIt;
//    int neigh_count, water_count=1;
//    float diff_x, diff_y, diff_z, d2, temp_x, temp_y, temp_z;
//    Coords * waterCoordinate;
//    Coords tempCoords;
//
//    for(int i=0; i<totalAtoms; i++) {
//
//        atomsInResidue.clear();
//        res_x.clear();
//        res_y.clear();
//        res_z.clear();
//
//        atomsInResidue.reserve(watersPerResidue);
//        res_x.reserve(watersPerResidue);
//        res_y.reserve(watersPerResidue);
//        res_z.reserve(watersPerResidue);
//
//        //residue = resi[i];
//        residue = trimmedResi[i];
//        resid = resID[i];
//        atomsCount = 1;
//
//        atomsInResidue.push_back(trimmedAtomType[i]);
//        res_x.push_back(centeredX[i]);
//        res_y.push_back(centeredY[i]);
//        res_z.push_back(centeredZ[i]);
//
//        // assemble atoms of residue, includes backbone and sidechains
//        j = i + 1;
//        while ((j < totalAtoms) && (resid == resID[j])) {
//            atomsInResidue.push_back(trimmedAtomType[j]);
//            res_x.push_back(centeredX[j]);
//            res_y.push_back(centeredY[j]);
//            res_z.push_back(centeredZ[j]);
//            i++;
//            j++;
//            atomsCount++;
//        }
//
//        // create water model
//        // residue will have to be centered
//        // hydrate residue (residue, atomType, x, y, z, alignedWaters)
//        waterModel.hydrateResidue(residue, atomsInResidue, atomsCount, res_x, res_y, res_z);
//        // waterModel.writeWaterCoords(residue, water_count);
//        // check waters, if acceptable keep in model
//        int totalWatersAllocated = waterModel.getTotalTempWaterArraySize();
//
//        int n = 0;
//        neigh_count = 0;
//        keepIt = true;
//
//        tempCoords = waterModel.getWaterCoords(0);
//
//        if (tempCoords.type.size() > 0) {
//
//        // go through first water and make neighbor list
//        if (n == 0) {
//            tempCoords = waterModel.getWaterCoords(n);
//
//            n++;
//            for (int m = 0; m < totalAtoms; m++) {
//                temp_x = centeredX[m];
//                temp_y = centeredY[m];
//                temp_z = centeredZ[m];
//                diff_x = tempCoords.x - temp_x;
//                diff_y = tempCoords.y - temp_y;
//                diff_z = tempCoords.z - temp_z;
//                d2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
//
//                // if structure atom is within 15 Angstrom of water, add to neighbor list
//                if (d2 <= 225) {
//                    neigh_x[neigh_count] = temp_x;
//                    neigh_y[neigh_count] = temp_y;
//                    neigh_z[neigh_count] = temp_z;
//                    //cout << "neigh_count " << neigh_count << endl;
//                    neigh_count++;
//                }
//
//                // reject water
//                if (d2 <= cutoff) {
//                    keepIt = false;
//                    waterModel.setOccupancyToZero(0);
//                }
//            }
//
//            if (keepIt) {
//                waterCoordinate = &tempCoords;
//                keptWaters[waterCount].x = (*waterCoordinate).x;
//                keptWaters[waterCount].y = (*waterCoordinate).y;
//                keptWaters[waterCount].z = (*waterCoordinate).z;
//                keptWaters[waterCount].type = (*waterCoordinate).type;
//                keptWaters[waterCount].occ = (*waterCoordinate).occ;
//                waterCount++;
//            }
//        }
//
//        while ((n < totalWatersAllocated) && (waterModel.getWaterCoords(n).type.size() > 0)) {
//            tempCoords = waterModel.getWaterCoords(n);
//            keepIt = true;
//
//            for (int m = 0; m < neigh_count; m++) {
//                diff_x = tempCoords.x - neigh_x[m];
//                diff_y = tempCoords.y - neigh_y[m];
//                diff_z = tempCoords.z - neigh_z[m];
//                d2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
//                // reject water
//                if (d2 <= cutoff) {
//                    keepIt = false;
//                    waterModel.setOccupancyToZero(n);
//                    break;
//                }
//            }
//
//            if (keepIt) {
//                waterCoordinate = &tempCoords;
//                keptWaters[waterCount].x = (*waterCoordinate).x;
//                keptWaters[waterCount].y = (*waterCoordinate).y;
//                keptWaters[waterCount].z = (*waterCoordinate).z;
//                keptWaters[waterCount].type = (*waterCoordinate).type;
//                keptWaters[waterCount].occ = (*waterCoordinate).occ;
//                waterCount++;
//            }
//            n++;
//        }
//    }
//    }
//
//    rWater = new float[waterCount];
//    phiWater = new float[waterCount];
//    waterCosTheta = new float[waterCount];
//
//    for(int rt=0; rt<waterCount; rt++){
//        waterCoordinate = &keptWaters[rt];
//        //(const float &)
//        //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp(&(*waterCoordinate).x, &(*waterCoordinate).y, &(*waterCoordinate).z);
//        // pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp((const float &)(*waterCoordinate).x, (const float &)(*waterCoordinate).y, (const float &)(*waterCoordinate).z);
//        rWater[rt] = *pconvertXYZ;
//        waterCosTheta[rt] = cosf(*(pconvertXYZ+1));
//        phiWater[rt] = *(pconvertXYZ+2);
//    }
//
//    ifWaters = true;
//}
