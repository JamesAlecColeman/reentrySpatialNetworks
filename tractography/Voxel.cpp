#include "Voxel.h"
#include <stdlib.h>
#include <iostream>
#include <vector>

std::vector<double> voxel::getxs(){return voxel::xs;}
std::vector<double> voxel::getys(){return voxel::ys;}
std::vector<double> voxel::getzs(){return voxel::zs;}
std::vector<int> voxel::getFibIndex(){return voxel::fibIndex;}
std::vector<int> voxel::getCellIndex(){return voxel::cellIndex;}

void voxel::addVertex(float xVertex, float yVertex, float zVertex, int fibreIndex, int cellNumber){
    voxel::xs.push_back(xVertex);
    voxel::ys.push_back(yVertex);
    voxel::zs.push_back(zVertex);
    voxel::fibIndex.push_back(fibreIndex);
    voxel::cellIndex.push_back(cellNumber);
}

voxel::voxel(){}
