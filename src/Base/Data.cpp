//
// Created by Robert Rambo on 12/01/2017.
//
#include "Data.h"

//using namespace std;
// empty constructor
Data::Data(){}

Data::Data(std::string filename) {
    this->iofqfilename = filename;
    // read in file
    // determine qmin and qmax
    std::ifstream data (this->iofqfilename.c_str());
    values.reserve(1800);
    std::string line;
    std::vector<std::string> tempLine;

    float qvalue, iofq, sigma;

    boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");

    boost::regex volFormat("(volume)|(VOLUME)");

    totalDataPoints = 0;
    if (data.is_open()) {
        while(!data.eof()) // if qmax is not present (qmaxPresent), do all the values in file
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) && boost::regex_search(tempLine.at(1), dataFormat)  ) {

                qvalue = std::stof(tempLine[0].c_str() );
                iofq = std::stof(tempLine[1].c_str());
                sigma = std::stof(tempLine[2].c_str());

                values.push_back(Datum(qvalue, iofq, sigma));
                totalDataPoints++;

            } else if (boost::regex_search(line, volFormat)){

                std::vector<std::string> volLine;
                boost::split(volLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                int index = findInArray("VOLUME", &volLine); // need to make case insensitive
                if (index > 0){
                    this->volume = std::stof(volLine[index + 2]);
                    std::cout << this->volume << std::endl;
                }
                //cout << " find in vector: => " << findInArray("VOLUME", &volLine) << endl;
                //this->volume = stof(volLine[4]);
            }
        }
    }

    qmax = values[totalDataPoints-1].getQ();
}



float Data::convertToBinProbability(float distance){

    float prob = 0.0;

    int binlocale = this->convertToBin(distance);

    int length = probability_per_bin.size();

    if (binlocale >= length){
        //cout << "length " << length << " " << binlocale << " => distance " << distance << endl;
        prob = 0.0;
    } else {
        prob = probability_per_bin[binlocale];
    }

    return prob;
}

/**
 * convert distance to Shannon bin of experimental P(r) distribution
int Data::convertToBin(float distance){

    int binlocale = (int) floor(distance/bin_width);

    if (binlocale >= round(shannon_bins)){
        binlocale = shannon_bins + 1;
    }

    return binlocale;
}
*/

/**
 * convert Shannon bin to distance
 */
float Data::convertBinToDistance(int bin){

    return (bin*bin_width + 0.5*bin_width);
}

float Data::calculateVolumeFromDistribution(int totalBins, std::vector<int> * binCount){

    float tempvol = 0;

    for(unsigned long int i=0; i < totalBins; i++){
        tempvol += (*binCount)[i];
    }

    tempvol *= bin_width*0.5;
    return tempvol;
}

/**
 * calculate PofR using Moore coefficients and normalize to 1
 */
void Data::normalizeMoorePofR() {

    int n;
    float norm=0;
    for (int i = 0; i < moore_coefficients.size(); i++) {
        n = i+1;
        norm += moore_coefficients[i]/(float)n*pow(-1,(n+1));
    }

    norm *= dmax*dmax/M_PI;
    // can I down-sample distribution ?
    probability_per_bin.reserve(shannon_bins);
    // bin_width = dmax/total_bins;
    // round up when calculating shannon number and use the d_max from the round up.

    // for each bin, calculate area
    float lower, upper, sum=0, value;
    for(int i=0; i<shannon_bins; i++){
        lower = i*bin_width;     //r-value
        upper = (i+1)*bin_width; //r-value
        // integrate between lower and upper
        value = integrateMoore(lower, upper);
        //value = integrateMooreSingle(bin_width*(0.5+i), upper);
        sum += value;
        probability_per_bin.push_back(value);
        //probability_per_bin.push_back(calculatePofRUsingMoore(bin_width*(0.5+i))/norm);
        std::cout << bin_width*(0.5+i) << " " << probability_per_bin[i] << std::endl;
    }

    // Normalilize
    //std::cout << "NORM " << norm << " " << sum << std::endl;
    normalize(1.0/sum);
}


/**
 * calculate PofR using Moore coefficients and normalize to 1
 */
void Data::normalizePrBins() {

    float norm=0;
    for (int i = 0; i < bin_coefficients.size(); i++) {
        norm += bin_coefficients[i];
    }

    // can I down-sample distribution ?
    probability_per_bin.reserve(bin_coefficients.size());
    // bin_width = dmax/total_bins;
    // round up when calculating shannon number and use the d_max from the round up.

    // for each bin, calculate area
    float value;
    std::cout << "  => NORMALIZED BINNED VALUES " << std::endl;
    for(int i=0; i<bin_coefficients.size(); i++){
        // integrate between lower and upper
        value = bin_coefficients[i]/norm;
        probability_per_bin.push_back(value);

        std::cout << printf("    RVALUE %7.3f : %.6f", (bin_width*(0.5+i)), (value)) << std::endl;
    }
}

void Data::normalize(float norm){

    FILE * pFile;

    const char *outputFileName;
    std::string nameOf = "normalized_expPr.txt";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    // Normalilize
    std::cout << " NORMALIZATION CONSTANT " << norm << std::endl;
    std::cout << " NORMALIZED P(r) : " << std::endl;
    fprintf(pFile, "# NORMALIZED P(r)\n");
    fprintf(pFile, "# REMARK  BINWIDTH => %.3f \n", bin_width);
    fprintf(pFile, "# REMARK      DMAX => %d \n", dmax);
    fprintf(pFile, "# REMARK   NS DMAX => %.0f \n", ns_dmax);
    fprintf(pFile, "  %.3f %.8f\n", 0.0, 0.0);
    for(int i=0; i<probability_per_bin.size(); i++){
        probability_per_bin[i] = probability_per_bin[i]*norm;
        //printf("  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]);
        fprintf(pFile, "  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]); // should match bin -> P(r) mapping in input file
    }

    fclose(pFile);
}

/**
 * Integrate sections of the P(r) distribution bounded by lower and upper
 */
float Data::integrateMoore(float lower, float upper){
    float lowerSum=0, upperSum=0, inva;
    for(int i=0; i<moore_coefficients.size(); i++){
        int n = i + 1;
        inva = dmax/(M_PI*n);
        // integrate between lower and upper
        lowerSum+= moore_coefficients[i]*(-inva*(lower*cos(lower*n*M_PI/dmax) - inva*sin(lower*n*M_PI/dmax)));
        upperSum+= moore_coefficients[i]*(-inva*(upper*cos(upper*n*M_PI/dmax) - inva*sin(upper*n*M_PI/dmax)));
    }
    return (upperSum - lowerSum);
}


/**
 * Integrate sections of the P(r) distribution bounded by lower and upper
 */
float Data::integrateMooreSingle(float midpointrvalue, float upper){
    float upperSum=0;
    for(int i=0; i<moore_coefficients.size(); i++){
        int n = i + 1;
        // integrate between lower and upper
        upperSum+= moore_coefficients[i]*(sin(upper*n*M_PI/dmax));
    }
    return 0.5/(double)dmax*midpointrvalue*upperSum;
}



float Data::calculatePofRUsingMoore(float rvalue){

    float sum=0;
    for(int i=0; i<shannon_bins; i++){
        int n = i+1;
        sum += moore_coefficients[i]*sin(M_PI*rvalue*n/dmax);
    }

    return rvalue*sum;
}


/**
 * creates Shannon limited binning of input P(r) distribution
 * use trapezoid rule to integerate for normalization
 * not sure if I should just interpolate the middle point of the bin using interpolation theory rather than average
 * param int count is the number of r-values in P(r) file
 */
void Data::normalizePofR(int count) {
    // float sum = pofr[0].pr + pofr[count-1].pr; first and last are 0 by default
    float sum;
    // integrate per bin
    probability_per_bin.reserve(shannon_bins);
    // round up when calculating shannon number and use the d_max from the round up.
    float r_value;
    float totalSum=0.0, lower, upper;

    for(int i=0; i<shannon_bins; i++){ // integrate between lower and upper
        lower = i*bin_width;     //q-value
        upper = (i+1)*bin_width; //q-value
        // cc = 0; // counts ticks
        // integrate bin (area of bin)
        // interpolate P(r) values for lower and upper
        int lowerIndex=0, upperIndex=count;
        float lowerRValue, upperRValue;

        for(int j=0; j<count; j++){
            if(pofr[j].r > lower){ // equal when r is zero
                lowerIndex = j-1;  // index of first value greater than lower bound on bin
                break;
            }
        }

        for(int j=0; j<count; j++){
            upperIndex = j;
            if(pofr[j].r >= upper){ // equal when r is dmax
                //upperIndex = j; // index of first value greater than upper bound on bin
                break;
            }
        } // if upper exceeds dmax, upperIndex is last element

        // perform linear interpolation of missing value
        lowerRValue = 0;
        if (lowerIndex > 0){
            lowerRValue = pofr[lowerIndex].pr + (pofr[lowerIndex+1].pr - pofr[lowerIndex].pr)*(lower - pofr[lowerIndex].r)/(pofr[lowerIndex+1].r- pofr[lowerIndex].r);
        }

        upperRValue = 0;
        if (upperIndex < count){
            upperRValue = pofr[upperIndex-1].pr + (pofr[upperIndex].pr - pofr[upperIndex-1].pr)*(upper - pofr[upperIndex-1].r)/(pofr[upperIndex].r- pofr[upperIndex-1].r);
        }

        // do integration of bin
        sum = 0.0;
        for (int j=lowerIndex; j<count; j++){
            r_value = pofr[j].r;
            int lastIndex;
            if (r_value > lower && r_value < upper){
                float height = std::min(lowerRValue, pofr[j].pr);
                sum += height*(r_value-lower) + abs(lowerRValue - pofr[j].pr)*(r_value-lower)*0.5;

                lower = r_value;
                lowerRValue = pofr[j].pr;
                lastIndex = j;
            }

            if (r_value >= upper) { // add last trapezoid to area
                // area of rectable + area of triangle
                // sum += (upper - lower)*min(upperRValue, pofr[lastIndex].pr) + abs(upperRValue - pofr[lastIndex].pr)*(r_value-lower)*0.5;
                sum += (upper - lower)*std::min(upperRValue, pofr[lastIndex].pr) + abs(upperRValue - pofr[lastIndex].pr)*(upper-lower)*0.5;
                break;
            }
        }

        probability_per_bin.push_back(sum);
        totalSum += sum;
    }

    // trapezoid rule (integrate Pr-distribution to calculate normalization constant)
    // assume last r-value is 0 at d_max
    float partialSum=0.0f;
    for(int i=0; i< count; i++){
        partialSum += 2*pofr[i].pr;
    }

    //float invPartialSum = (2.0*(float)(count-1))/(dmax*partialSum);
    float invPartialSum = 1.0/totalSum;
    //cout << "Partial SUM DIRECT => " << dmax/(2.0*(float)(count-1))*partialSum << " " << totalSum << endl;
/*
    // for each Shannon point, find first that is larger than it
    // then do weighted average based on distance from mid-point to upper and lower bounds of bin-width
    int nextIndex = 0;

    for(int i=0; i<total_bins; i++) {
        float mid = i*bin_width + 0.5*bin_width;
        int lowerIndex;
        float upperRValue;
        float upperPrValue;

        for(int j=nextIndex; j<count; j++){
            r_value = pofr[j].r;
            if(r_value > mid){
                upperRValue = r_value;
                upperPrValue = pofr[j].pr;
                lowerIndex = j-1;
                nextIndex = j;
                break;
            }
        }

        float diffLow = (mid - pofr[lowerIndex].r);
        float diffUp = (upperRValue - mid);
        float invWidth = 1.0/(diffUp + diffLow);
        // set normalized P(r) by multiplying normalization constant
        probability_per_bin.push_back(invPartialSum*(diffUp*upperPrValue + diffLow*pofr[lowerIndex].pr)*invWidth);
    }
*/
    // Normalilize
    //int stopLimit = probability_per_bin.size();
//    for(int i=0; i<probability_per_bin.size(); i++){
//        probability_per_bin[i] = probability_per_bin[i]*invPartialSum;
//    }
    normalize(invPartialSum);
}

void Data::calculateRatioPr(std::vector<float> &modelPR){

    float rv;
    int total = modelPR.size();
    float totalSum=0.0;

    for (int i=0; i<total; i++){
        rv = this->getBinRValue(i);
        totalSum += probability_per_bin[i]/modelPR[i];
        std::cout << rv << "  " << probability_per_bin[i]/modelPR[i] << std::endl;
    }

    for (int i=0; i<total; i++){
        rv = this->getBinRValue(i);
        std::cout << rv << "  " << probability_per_bin[i]/modelPR[i]/totalSum << std::endl;
    }

}

/**
 * model Pr can have more bins than probability_per_bin
 * if model exceeds the experimental bin, Divergence is zero for that bin
 * since zero x log [zero/number] => zero
 */
float Data::calculateKLDivergence(std::vector<int> &modelPR){

    double totalCounts = 0.0;
    float kl=0.0;
    double prob, *value;
    int totalm = modelPR.size();
    std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    // trapezoid rule for normalizing model PR
    // modelPR does not include the end points, r=0, r = dmax
    //float lastValueModelPR = modelPR_float[totalm-1];
    //float firstPart = 0.5*(0.5*bin_width)*modelPR_float[0];
    //float lastPart  = 0.5*(0.5*bin_width)*lastValueModelPR;
    for(int i=0; i< totalm; i++){
        totalCounts += modelPR_float[i];//*bin_width;
    }

    //totalCounts *= 2.0; //multiple by 2
    //totalCounts -= (modelPR_float[0] + lastValueModelPR); //
    //totalCounts = ((totalm*bin_width-bin_width)*totalCounts)/(2.0*(totalm-1)) + firstPart + lastPart;
    //totalCounts = (totalm*bin_width)/(2.0*totalm)*2.0*totalCounts;
    //totalCounts = bin_width*totalCounts;
    // end trapezoid rule
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    //int last = totalObsInPr-1;

    for (int i=0; i < zeroBin; i++){
        // i know every value in working_probability up to zeroBin is nonzero
        value = &modelPR_float[i];
        if (*value > 0){
            prob = working_probability_per_bin[i];  // bounded by experimental Shannon Number
            kl += prob * log(prob/(*value) * totalCounts);
        } else { // severely penalize any model bin that is zero
            //kl += 33554430;
            kl += 10000000;
        }
    }


    /*
    for (int i=0; i < zeroBin; i++){

        if (modelPR_float[i] <= 0){
            kl += 10000000000000000000;
        } else {
            prob = working_probability_per_bin[i];  //
            //kl += prob * log(prob/(bin_width*modelPR_float[i])*totalCounts);
            kl += prob * log(prob/(modelPR_float[i])*totalCounts);
        }
    }

    // assume if last bin in dataset is 0, then 0*log0 = 0
    //if (probability_per_bin[last] > 0){
    if (working_probability_per_bin[last] > 0){
        prob = working_probability_per_bin[last];
        //kl += prob * log(prob/(bin_width*modelPR_float[last])*totalCounts);
        kl += prob * log(prob/(modelPR_float[last])*totalCounts);
        //kl += prob * log(prob/(modelPR_float[last])*totalCounts);
    }
*/
    return kl*1.0/(double)shannon_bins;  // returns value per bin
    //return kl*1.0/(double)totalm;  // returns value per bin
}

/**
 * model Pr can have more bins than probability_per_bin
 * if model exceeds the experimental bin, Divergence is zero for that bin
 * since zero x log [zero/number] => zero
 */
float Data::calculateKLDivergenceMultiComponent(std::vector<float> &modelPR) {

    double totalCounts = 0.0;
    float kl=0.0;
    double prob;
    float *value;
    int totalm = modelPR.size();
    //std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    // trapezoid rule for normalizing model PR
    // modelPR does not include the end points, r=0, r = dmax

    for(int i=0; i< totalm; i++){
        totalCounts += modelPR[i];//*bin_width;
    }

    //totalCounts *= 2.0; //multiple by 2
    //totalCounts -= (modelPR_float[0] + lastValueModelPR); //
    //totalCounts = ((totalm*bin_width-bin_width)*totalCounts)/(2.0*(totalm-1)) + firstPart + lastPart;
    //totalCounts = (totalm*bin_width)/(2.0*totalm)*2.0*totalCounts;
    //totalCounts = bin_width*totalCounts;
    // end trapezoid rule
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0

    for (int i=0; i < zeroBin; i++){
        // i know every value in working_probability up to zeroBin is nonzero
        //value = &modelPR_float[i];
        value = &modelPR[i];
        if (*value > 0){
            prob = working_probability_per_bin[i];  // bounded by experimental Shannon Number
            kl += prob * log(prob/(*value) * totalCounts);
        } else { // severely penalize any model bin that is zero
            //kl += 33554430;
            kl += 1000000;
        }
    }

    return kl*1.0/(float)shannon_bins;  // returns value per bin
    //return kl*1.0/(double)totalm;  // returns value per bin
}

void Data::printKLDivergence(std::vector<int> &modelPR){

    float totalCounts = 0.0;
    float prob;
    int totalm = modelPR.size();

    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    for (int i=0; i<totalm; i++){
        totalCounts += modelPR[i];
    }

    float tempPR, r;
    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    std::cout << "FINAL MODEL " << std::endl;
    std::cout << "       r       MODEL      EXP        D_KL" << std::endl;
    for (int i=0; i < shannon_bins; i++){
        prob = probability_per_bin[i];  // bounded by experimental Shannon Number
        tempPR = modelPR[i];
        r = bin_width*i+0.5*bin_width;
        printf(" %10.4f  %.6f  %.6f  % .6f\n", r, prob, tempPR/totalCounts, (prob * log(prob/tempPR*totalCounts)));
        //cout << r << " " << prob << " " << tempPR/totalCounts << " " << prob * log(prob/tempPR*totalCounts) << endl;
    }
}

/**
 * bins can be 1 or more times the number of shannon bins
 */
void Data::setDataBinSize(int bins) {

    lowerBounds.clear();
    data_bins = round(bins*shannon_bins);

    float binwidth = qmax/(float)data_bins;

    int count, index = 0;
    float sumSN, averageNoise, invCount, deltaQ;

    // determine points per bin
    float lower, upper;

    for (int i=0; i < data_bins; i++){
        lower = binwidth*i;
        upper = binwidth*(i+1);
        float shannonHartley;
        count = 0;
        sumSN = 0;
        // go through data and count and determine average noise in bin

        while( (values[index].getQ() >= lower && values[index].getQ() < upper) && index < totalDataPoints ){
            Datum temp = values[index];
            sumSN += temp.getI()/temp.getSigma();
            count++; // number of points within range
            index++;
        }

        invCount = 1.0/(float)count;
        deltaQ = (upper-lower)*invCount;

        averageNoise = sumSN*invCount;

        shannonHartley = 2.0*M_PI/(float)dmax*log(1 + averageNoise);

        if (count > 0){
            if (deltaQ > shannonHartley){

                if (count < 5){
                    pointsPerBin.push_back(count); // use all points in bin
                } else {
                    pointsPerBin.push_back(5);
                }
            } else {
                pointsPerBin.push_back(2);
            }

            lowerBounds.push_back(lower);
        }
    }
}

void Data::addPofRData(std::string filename) {

    this->pofrfilename = filename;
    // read in file

    std::ifstream data (this->pofrfilename.c_str());
    pofr.reserve(200);
    std::string line;
    std::vector<std::string> tempLine;

    float rvalue, pvalue, sigma;

    boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
    boost::regex dmaxFormat("(dmax)|(DMAX)");
    boost::regex qmaxFormat("(qmax)|(QMAX)");
    boost::regex volFormat("(volume)|(VOLUME)");
    boost::regex raveFormat("<r>");
    boost::regex rgFormat("REAL Rg :");

    // Read in Pr data - if sascif, this will need to change for now, three column format
    int tempCount = 0;
    if (data.is_open()) {

        // if qmax is not present (qmaxPresent), do all the values in file
        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) && boost::regex_search(tempLine.at(1), dataFormat)  ) {

                rvalue = strtof(tempLine[0].c_str(), NULL);
                pvalue = strtof(tempLine[1].c_str(), NULL);
                sigma = strtof(tempLine[2].c_str(), NULL);

                pofr.push_back(RealSpace(rvalue, pvalue, sigma));
                tempCount++;

            } else if (boost::regex_search(line, dmaxFormat)){
                std::vector<std::string> dmaxLine;
                boost::split(dmaxLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->dmax = std::stof(dmaxLine[4]);
            } else if (boost::regex_search(line, qmaxFormat)){
                std::vector<std::string> qmaxLine;
                boost::split(qmaxLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->qmax = stof(qmaxLine[6]);
            } else if (boost::regex_search(line, volFormat)){
                std::vector<std::string> volLine;
                boost::split(volLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->volume = stof(volLine[4]);
            } else if (boost::regex_search(line, raveFormat)){
                std::vector<std::string> raveLine;
                boost::split(raveLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->rave = stof(raveLine[5]);
            } else if (boost::regex_search(line, rgFormat)){
                std::vector<std::string> rgLine;
                boost::split(rgLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
                this->rg = stof(rgLine[5]);
            }

        }
        std::cout << "  => PRDATA : " << this->pofrfilename << " QMAX => " << this->qmax << " DMAX => " << this->dmax << " VOL " << this->volume << std::endl;
    }

    data.close();

    if (dmax <= 0 || qmax < 0.1){
        // throw exception
        throw "QMAX/DMAX NOT FOUND OR EQUAL TO ZERO";
    }

    //
    // shannon_bins = round(dmax*qmax/M_PI); // based on the data
    // qmax is specified in pr.dat file
    // shannon bins can be <= Moore Coefficents, e.g. downsampling the distribution
    //

    shannon_bins = ceil(dmax*qmax/M_PI);  // shannon bins <= Moore Coefficients
    ns_dmax = shannon_bins*M_PI/qmax;     // dmax corresponding to Shannon Bin

    //bin_width = M_PI/qmax;

    if (shannon_bins < 3){// throw exception
        throw "Too few shannon bins, USE ANOTHER PROGRAM";
    }

    // use experimental Shannon number before scaling in setBinSize

    this->parseBins();

    //this->parseMooreCoefficients();

    if (bin_coefficients.size() > 2){

        if (shannon_bins > bin_coefficients.size()){
            std::cout << "-----------------------     WARNING     -----------------------" << std::endl;
            std::cout << "-- qmax NOT SUPPORTED by Shannon Number and bin size in file          --" << std::endl;
            std::cout << "-- Validate input files                                         --" << std::endl;
            std::cout << "-- Shannon_BINS : " << shannon_bins << " != " << bin_coefficients.size() << std::endl;
            exit(0);
        }

        bin_width = dmax/(double)bin_coefficients.size();
        //normalizeMoorePofR();
        normalizePrBins();
    } else {
        this->normalizePofR(tempCount);
    }

    //this->setDataBinSize(2);
}

void Data::parseBins() {

    std::ifstream data (this->pofrfilename.c_str());
   // boost::regex binValueLine("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
    boost::regex coefficient("BIN_[0-9]+");
    boost::regex remarkFormat("REMARK");

    std::string line;
    std::vector<std::string> tempLine;

    int binCount=0;
    std::cout << "  => READING POFR FILE " << this->pofrfilename << std::endl;

    if (data.is_open()) {

        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
            //&& !boost::regex_search(line, background)


            if(boost::regex_search(tempLine.at(0), remarkFormat) && boost::regex_search(line, coefficient) ){
                std::vector<std::string> mline;
                boost::trim(line);
                boost::split(mline, line, boost::is_any_of("\t  "), boost::token_compress_on);
                int total = mline.size();

                try {
                    float value = std::abs(stof(mline.back()));

                    if(std::strcmp(mline[total-2].c_str(), ":")==0 && value > 0){
                        bin_coefficients.push_back(stof(mline.back()));
                        binCount++;
                    } else {
                        throw std::invalid_argument( "Improper Bin Value at : \n\t" + line  + " \n Can not be zero or negative");
                    }

                } catch (std::exception &err) {
                    std::cerr<<"Caught "<<err.what()<<std::endl;
                    std::cerr<<"Type "<<typeid(err).name()<<std::endl;
                    exit(0);
                };
            }
        }
    }

    data.close();

}


void Data::parseMooreCoefficients(){
    // read in file
    std::ifstream data (this->pofrfilename.c_str());
    boost::regex mooreLine("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
    boost::regex coefficient("m_\\(\\s?[0-9]+\\)");
    boost::regex background("m_\\(0\\)");
    boost::regex remarkFormat("REMARK");

    std::string line;
    std::vector<std::string> tempLine;

    int mooreCount=0;
    std::cout << "  => READING POFR FILE " << this->pofrfilename << std::endl;

    if (data.is_open()) {

        // if qmax is not present (qmaxPresent), do all the values in file
        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);
            //&& !boost::regex_search(line, background)

            if(boost::regex_search(tempLine.at(0), remarkFormat) && boost::regex_search(line, coefficient) ){
                std::vector<std::string> mline;
                boost::trim(line);
                boost::split(mline, line, boost::is_any_of("\t  "), boost::token_compress_on);
                int total = mline.size();
                std::cout << "USING LINE: " << line << std::endl;

                try {
                    float value = abs(stof(mline.back()));
                    if(std::strcmp(mline[total-2].c_str(), ":")==0 && value != 0){
                        moore_coefficients.push_back(stof(mline.back()));
                        mooreCount++;
                    } else {
                        throw std::invalid_argument( "Improper Moore Value at : \n\t" + line  + " \n Can not be zero");
                    }
                } catch (std::exception &err) {
                    std::cerr<<"Caught "<< err.what()<<std::endl;
                    std::cerr<<"Type "<<typeid(err).name()<<std::endl;
                    std::exit(0);
                };
            }
        }
    }

    data.close();
}


/**
 * i indexes from 0...(probability_per_bin.size()-1)
 * lower => i*bin_width
 * upper => (i+1)*bin_width
 * r_value in bin => 0.5*bin_width + i*bin_width
 */
float Data::getBinRValue(int binIndex){
    float rvalue = 0.5*bin_width + binIndex*bin_width;
    return rvalue;
}

/**
 * add memory address of the phases
 */
void Data::addPhase(Phase & phase){
    phases.push_back(&phase);
    totalPhases = phases.size();
}

bool Data::checkPofRFile(std::string file) {

    // read in file
    std::ifstream data (file.c_str());
    bool returnMe = false;
    std::string line;

    boost::regex prFormat("P\\(r|R\\)");

    if (data.is_open()) {

        std::getline(data, line);

        if(boost::regex_search(line, prFormat)){
            returnMe = true;
        }
    }

    data.close();
    return returnMe;
}

bool Data::isPhasePresent(std::string id){
    bool test = false;

    for(int i=0;i<getTotalPhases(); i++){
        if ((phases[i]->getID()).compare(id) == 0){
            test = true;
            break;
        }
    }

    return test;
}

int Data::findInArray(std::string value, std::vector<std::string> * strings) {

    int index = -1;
    std::vector<std::string>::iterator beginIt = strings->begin();
    std::vector<std::string>::iterator it = std::find(beginIt, strings->end(), value);
    if (it != strings->end()){
        index = std::distance(beginIt, it);
    }
    return index;
}
