#include "Fibre.h"
#include <stdlib.h>
#include <iostream>
#include <vector>

// Add to fibre vertices placed in the 'forwards' direction
void fibre::addPosF(double xxF, double yyF, double zzF){
    fibre::xF.push_back(xxF);
    fibre::yF.push_back(yyF);
    fibre::zF.push_back(zzF);
}

// Add to fibre vertices placed in the 'backwards' direction
void fibre::addPosB(double xxB, double yyB, double zzB){
    fibre::xB.push_back(xxB);
    fibre::yB.push_back(yyB);
    fibre::zB.push_back(zzB);
}

// Returns fibre vertices placed in the 'forwards' direction
std::vector<double> fibre::getxF(){return fibre::xF;}
std::vector<double> fibre::getyF(){return fibre::yF;}
std::vector<double> fibre::getzF(){return fibre::zF;}

// Returns fibre vertices placed in the 'backwards' direction
std::vector<double> fibre::getxB(){return fibre::xB;}
std::vector<double> fibre::getyB(){return fibre::yB;}
std::vector<double> fibre::getzB(){return fibre::zB;}

// Add to the fibre length measurement
void fibre::addLength(double length){fibre::fibreLength += length;}

// Returns fibre length
double fibre::getFibreLength(){return fibre::fibreLength;}

// Clear fibre data
void fibre::clearFibre(){
    fibre::xB.clear();
    fibre::yB.clear();
    fibre::zB.clear();
    fibre::xF.clear();
    fibre::yF.clear();
    fibre::zF.clear();
	fibre::fibreLength = 0;
}

fibre::fibre(){fibre::fibreLength = 0;}
