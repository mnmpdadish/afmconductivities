//
//  Utilities.h
//  oneBuddy
//

#pragma once
#include <cassert>
//#include <cctype>
//#include <cstdlib>
#include <cstring>
#include <iostream>
//#include <iosfwd>
#include <fstream>
#include <string>
//#include <sstream>
#include <limits>
#include <math.h>
#include <complex>
#include <sys/stat.h>


#define MAX_LENGHT std::numeric_limits<std::streamsize>::max()
using namespace std;


// NUMBER UTILITIES
bool dbEqual(const double &a, const double &b) {
    return abs(a - b) < 0.000000001;
}

// FILES UTILITIES
inline bool exists(const string &name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

void nextLine(istream &flux) {
    flux.ignore(MAX_LENGHT,'\n');
}

void reset(istream &flux){
    flux.clear();
    flux.seekg(0);
}

bool reach(istream &flux, const string &search) {
    reset(flux);
    string buff;
    while (flux.peek() != EOF)
    {
        flux >> buff;
        if(buff.compare(search)==0){
            return true;
        }
        nextLine(flux);
    }
    return false;
}

template<typename TYPE>
bool readNumber(istream &flux, const string &search, TYPE &value) {
    if(reach(flux, search)){
        string line;
        getline(flux,line);
        istringstream buffer(line);
        string word;
        
        if(buffer.peek() == EOF || buffer.peek() == '\n') {
            cout << search << " = " << value << '\n';
            return false;
        }    
        buffer >> value;
        cout << search << " := " << value << '\n';
        return true;
    }
    cout << search << " = " << value << '\n';
    return false;
}        



