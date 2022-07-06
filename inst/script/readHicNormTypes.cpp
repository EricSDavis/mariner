#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <streambuf>
#include "zlib.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// returns whether or not this is valid HiC file
bool readMagicString(ifstream& fin) {
    string str;
    getline(fin, str, '\0' );
    return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

// [[Rcpp::export]]
Rcpp::CharacterVector readHicNormTypes(std::string fname) {
    ifstream fin(fname, ios::in | ios::binary);
    if (!fin) {
        fin.close();
        stop("file %s cannot be opened for reading.", fname);
    }

    // hicMagicString
    if (!readMagicString(fin)) {
        fin.close();
        stop("Hi-C magic string is missing, does not appear to be a hic file.");
    }

    // Verison
    int version;
    fin.read((char*)&version, sizeof(int));
    if (version < 6) {
        fin.close();
        stop("Version %d no longer supported.", version);
    }

    // footerPosition
    long master;
    fin.read((char*)&master, sizeof(long));

    // Skip to footer
    fin.seekg(master, ios::beg);

    // Number of bytes for the "version 5" footer
    int32_t nBytes;
    fin.read((char*)&nBytes, sizeof(int32_t));

    // Initialize variable to store norm types
    Rcpp::CharacterVector normTypes;

    // Read through footer section to get to norm types
    int32_t nEntries;
    fin.read((char*)&nEntries, sizeof(int32_t));
    for (int i=0; i<nEntries; i++) {
        string str;
        getline(fin, str, '\0');
        int64_t fpos;
        fin.read((char*)&fpos, sizeof(int64_t));
        int32_t sizeinbytes;
        fin.read((char*)&sizeinbytes, sizeof(int32_t));
    }

    int32_t nExpectedValues;
    fin.read((char*)&nExpectedValues, sizeof(int32_t));
    for (int i=0; i<nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0'); //unit
        int32_t binSize;
        fin.read((char*)&binSize, sizeof(int32_t));

        int32_t nValues;
        fin.read((char*)&nValues, sizeof(int32_t));
        for (int j=0; j<nValues; j++) {
            double v;
            fin.read((char*)&v, sizeof(double));
        }

        int32_t nNormalizationFactors;
        fin.read((char*)&nNormalizationFactors, sizeof(int32_t));
        for (int j=0; j<nNormalizationFactors; j++) {
            int32_t chrIdx;
            fin.read((char*)&chrIdx, sizeof(int32_t));
            double v;
            fin.read((char*)&v, sizeof(double));
        }
    }

    fin.read((char*)&nExpectedValues, sizeof(int32_t));
    for (int i=0; i<nExpectedValues; i++) {

        // Record available norm types
        string type;
        getline(fin, type, '\0'); //typeString
        normTypes.push_back(type);

        string unit0;
        getline(fin, unit0, '\0'); //unit
        int32_t binSize;
        fin.read((char*)&binSize, sizeof(int32_t));

        int32_t nValues;
        fin.read((char*)&nValues, sizeof(int32_t));
        for (int j=0; j<nValues; j++) {
            double v;
            fin.read((char*)&v, sizeof(double));
        }
        int32_t nNormalizationFactors;
        fin.read((char*)&nNormalizationFactors, sizeof(int32_t));
        for (int j=0; j<nNormalizationFactors; j++) {
            int32_t chrIdx;
            fin.read((char*)&chrIdx, sizeof(int32_t));
            double v;
            fin.read((char*)&v, sizeof(double));
        }
    }

    // Close .hic file
    fin.close();

    // Replace empty string with NONE
    for(int i=0; i<normTypes.length(); i++){
        if(normTypes[i] == ""){
           normTypes[i] = "NONE";
        }
    }

    // Get unique
    normTypes = unique(normTypes);

    return normTypes;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
hicFile <- "/Users/eric/Phanstiel Lab Dropbox/Eric Davis/projects/mariner/mariner/inst/extdata/test.hic"
# hicFile <- "/Users/eric/Phanstiel Lab Dropbox/Eric Davis/projects/CHON/data/hic/CHON_HiC_C28_WT_NA_0_S_inter_30.hic"
# hicFile <- "/Users/eric/Phanstiel Lab Dropbox/Eric Davis/largeData/hic/LEUK/condition/HEK_HiC_NUP_IDR_FS_A9_megaMap_inter_30.hic"
readHicNormTypes(fname = hicFile)
*/
