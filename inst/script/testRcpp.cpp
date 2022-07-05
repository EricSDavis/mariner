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

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// returns whether or not this is valid HiC file
bool readMagicString(ifstream& fin) {
    string str;
    getline(fin, str, '\0' );
    return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

// [[Rcpp::export]]
void readHicNorms(std::string fname) {
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
    Rcout << "Version: " << version << std::endl;

    // footerPosition
    long footerPosition;
    fin.read((char*)&footerPosition, sizeof(long));
    Rcout << "Footer pos: " << footerPosition << std::endl;

    // Genome
    string genomeId;
    getline(fin, genomeId, '\0' );
    Rcout << "genomeId: " << genomeId << std::endl;

    // Skip to footer
    istream& seekg(streampos footerPosition);

    // Number of bytes for the "version 5" footer
    int nBytesV5;
    fin.read((char*)&nBytesV5, sizeof(int));
    Rcout << "nBytesV5: " <<  nBytesV5 << endl;

    // Master index
    int nEntries;
    fin.read((char*)&nEntries, sizeof(int));
    Rcout << "nEntries: " << nEntries << endl;
    // reading and ignoring attribute-value dictionary
    // for (int i=0; i<nEntries; i++) {
    //     string key, position, size;
    //     getline(fin, key, '\0');
    //     getline(fin, position, '\0');
    //     getline(fin, position, '\0');
    // }

    // Skip master index entries
    istream& seekg(streampos nEntries);

    // Expected value vectors
    int nExpectedValueVectors;
    fin.read((char*)&nExpectedValueVectors, sizeof(int));
    Rcout << "nExpectedValueVectors: " << nExpectedValueVectors << endl;

    // // nattributes
    // int nattributes;
    // fin.read((char*)&nattributes, sizeof(int));
    // Rcout << nattributes << std::endl;
    //
    // // reading and ignoring attribute-value dictionary
    // for (int i=0; i<1; i++) {
    //     string key, value;
    //     getline(fin, key, '\0');
    //     getline(fin, value, '\0');
    // }
    //
    // int nChrs;
    // fin.read((char*)&nChrs, sizeof(int));
    // Rcout << nChrs << std::endl;
    //
    // StringVector chrom_names(nChrs);
    // NumericVector chrom_lengths(nChrs);
    // for (int i=0; i<nChrs; i++) {
    //     string name;
    //     int length;
    //     getline(fin, name, '\0');
    //     fin.read((char*)&length, sizeof(int));
    //     chrom_names[i] = name;
    //     chrom_lengths[i] = length;
    // }

    fin.close();

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
hicFile <- "/Users/eric/Phanstiel Lab Dropbox/Eric Davis/projects/mariner/mariner/inst/extdata/test.hic"
readHicNorms(fname = hicFile)
*/
