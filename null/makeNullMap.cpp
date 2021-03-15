#include <iostream>
#include <cmath>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include "numpy.hpp"
#include <math.h>
#include <fstream>
#include "Cell.h"
#include <algorithm>

using namespace std;

 // Coarse graining factor
 const double g = 3;

 // Human fibre orientation dataset dimensions
 //const int Lx = (int)ceil(303/g);
 //const int Ly = (int)ceil(275/g);
 //const int Lz = (int)ceil(244/g);

 // Sheep1 fibre orientation dataset dimensions
 const int Lx = (int)ceil(325/g);
 const int Ly = (int)ceil(949/g);
 const int Lz = (int)ceil(611/g);
 
 int cellIndex = 0;
 
 vector<cell> cells;
 
 // Stores flattened final fibre map	   
 vector<double> cellsData;

// Stores fibre orientation dataset as 4D array
 vector<vector<vector<vector<double> > > > data(3, vector<vector<vector<double> > >
                                           (Lx, vector<vector<double> >
                                           (Ly, vector<double>
                                           (Lz))));
										   
 vector<vector<vector<int> > > cellIndices(Lx, vector<vector<int> >
                                           (Ly, vector<int>
                                           (Lz)));


// Returns pseudorandom number between 0 and 1
float rnd(){return static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * 1;}

// Returns distance^2 between vectors (x, y, z) and (u, v, w)
float squaredDistance(float x, float y, float z, float u, float v, float w){
    return (pow(x - u, 2) + pow(y - v, 2) + pow(z - w, 2));
}

// Returns whether point resides in a full voxel (part of atrial volume) or not (blank space)
bool full(float i, float j, float k){

    i = floor(i);
    j = floor(j);
    k = floor(k);

    if(data[0][i][j][k] != 0 || data[1][i][j][k] != 0 || data[2][i][j][k] != 0){return 1;}
    else{return 0;}

}

// Returns whether point resides within the dataset boundaries
bool inBoundaries(float i, float j, float k){

    i = floor(i);
    j = floor(j);
    k = floor(k);

    if(i >= Lx || j >= Ly || k >= Lz || i < 0 || j < 0 || k < 0){return 0;}
    else{return 1;}
}

// Loads and unflattens the flat fibre orientation dataset of form (i, j, k, v1, v2, v3, i, j, k, v1, v2, v3, ...)
void loadData(){ 

    vector<double> flatData;
    vector<int> flatMeta;
	
    aoba::LoadArrayFromNumpy("PATH TO .NPY FILE HERE", flatMeta, flatData);

    long long int counter = 0;
    int looper = 0;
    int i;
    int j;
    int k;
    double v1;
    double v2;
    double v3;
	
    while(counter < flatData.size()){
		
        if(looper == 0){i = flatData[counter]; looper++;}
        else if(looper == 1){j = flatData[counter]; looper++;}
        else if(looper == 2){k = flatData[counter]; looper++;}
        else if(looper == 3){v1 = flatData[counter]; looper++;}
        else if(looper == 4){v2 = flatData[counter]; looper++;}
        else if(looper == 5){
            v3 = flatData[counter];
            data[0][i][j][k] = v1;
            data[1][i][j][k] = v2;
            data[2][i][j][k] = v3;
            looper = 0;
        }

        counter++;
    }
	

}

// Calculates distances between all nodes and their nearby nodes (<2 voxel units away)
void doProximal(){

	float x;
	float y;
	float z;
	
	float xVoxel;
	float yVoxel;
	float zVoxel;
	
	float xOther;
	float yOther;
	float zOther;
	
	float distance;
    
	for(long long int i = 0; i < cells.size(); i++){
		
		x = cells[i].getx();
		y = cells[i].gety();
		z = cells[i].getz();
		
		xVoxel = floor(x);
		yVoxel = floor(y);
		zVoxel = floor(z);
		
		for(int p = -2; p < 3; p++){
		    for(int q = -2; q < 3; q++){
				for(int r = -2; r < 3; r++){
					
					if(inBoundaries(xVoxel + p, yVoxel + q, zVoxel + r)){
						
						xOther = cells[cellIndices[xVoxel + p][yVoxel + q][zVoxel + r]].getx();
						yOther = cells[cellIndices[xVoxel + p][yVoxel + q][zVoxel + r]].gety();
						zOther = cells[cellIndices[xVoxel + p][yVoxel + q][zVoxel + r]].getz();
						
						distance = pow(pow(x - xOther, 2) + pow(y - yOther, 2) + pow(z - zOther, 2), 0.5);		
									
						if(distance < 2){
							
							cells[i].addProximalCell(cellIndices[xVoxel + p][yVoxel + q][zVoxel + r], distance);
						
						}
					}
				}
			}
		}
	}
		
}


int main (int argc, char **argv) {
	
    srand (time(NULL));
	
	loadData();
	
	cellIndex=0;
    
	for(int i = 0; i < Lx; i++){
        for(int j = 0; j < Ly; j++){
            for(int k = 0; k < Lz; k++){
				if (full(i,j,k)){
					cells.push_back(cell((double) i + rnd(), (double) j + rnd(), (double) k + rnd(), 0));
					cellIndices[i][j][k]=cellIndex;
					cellIndex++;
				}
			}
		}
	}

    doProximal();
	
	// Flatten and save cell objects
    for(int i = 0; i < cells.size(); i++){
        cellsData.push_back(cells[i].getx());
        cellsData.push_back(cells[i].gety());
        cellsData.push_back(cells[i].getz());
        cellsData.push_back(cells[i].getFib());
        cellsData.push_back(0);
        cellsData.push_back(0);
        cellsData.push_back(cells[i].getcons().size());
		
        for(int c = 0; c < cells[i].getcons().size(); c++){
            cellsData.push_back(cells[i].getcons()[c]);
        }

        cellsData.push_back(cells[i].getCellDistances().size());
        for(int c = 0; c < cells[i].getCellDistances().size(); c++){
            cellsData.push_back(cells[i].getCellDistances()[c]);
        }

        cellsData.push_back(cells[i].getProximalCells().size());
        for(int c = 0; c < cells[i].getProximalCells().size(); c++){
            cellsData.push_back(cells[i].getProximalCells()[c]);
        }
    }

    aoba::SaveArrayAsNumpy("SAVE PATH FOR FLAT .NPY NULL MAP HERE", cellsData);
	
    return 0;
}



