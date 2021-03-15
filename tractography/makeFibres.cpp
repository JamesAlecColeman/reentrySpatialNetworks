/*
Modified version of Evenly Spaced Streamlines fibre tractography (Merhof et al. 2005)
Adapted to increase fibre length and coverage.
*/

#include <iostream>
#include <cmath>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <fstream>

#include "Fibre.h"
#include "Voxel.h"
#include "Cell.h"
#include "numpy.hpp"

using namespace std;

 // Coarse graining factor
 const double g = 1;

 // Human fibre orientation dataset dimensions
 //const int Lx = (int)ceil(303/g);
 //const int Ly = (int)ceil(275/g);
 //const int Lz = (int)ceil(244/g);

 // Sheep1 fibre orientation dataset dimensions
 //const int Lx = (int)ceil(325/g);
 //const int Ly = (int)ceil(949/g);
 //const int Lz = (int)ceil(611/g);

 // Sheep2 fibre orientation dataset dimensions
 const int Lx = (int)ceil(89/g);
 const int Ly = (int)ceil(345/g);
 const int Lz = (int)ceil(208/g);
 
 // Tractography parameters
 float dStep = 1; // Fixed step size tracked along voxel orientation
 float dSep = 0.7; // Generation of fibre terminated if it gets closer than dSep to a different fibre
 int threshold = 30; // Minimum number of vertices a fibre must have to be added to fibre map
 int seedsFactor = 1; // Generates (number of full voxels)*seedsFactor seeds per round
 int increment = 5;
 
 // Fibres generated from seed points randomly placed over atria
 vector<double> xSeeds;
 vector<double> ySeeds;
 vector<double> zSeeds;

 // Stores fibres and nodes as objects
 vector<fibre> fibres;
 vector<cell> cells;
 
 // Fibre map metadata
 int fibreCount;
 int fullCount;
 double nOnes;
 double nZeros;
 int totalLength = 0;
 
 // Temporary fibre class and its vertices
 fibre tempFibre;
 vector<float> xVertices;
 vector<float> yVertices;
 vector<float> zVertices;

// Stores fibre orientation dataset as 4D array
 vector<vector<vector<vector<double> > > > data(3, vector<vector<vector<double> > >
                                           (Lx, vector<vector<double> >
                                           (Ly, vector<double>
                                           (Lz))));

// Stores fibre vertices corresponding to each voxel on the dataset
 vector<vector<vector<voxel> > > voxels(Lx, vector<vector<voxel> >
                                           (Ly, vector<voxel>
                                           (Lz)));

// Stores (equally spaced 1 voxel apart) nodes corresponding to each voxel on the dataset
 vector<vector<vector<voxel> > > rVoxels(Lx, vector<vector<voxel> >
                                           	(Ly, vector<voxel>
                                           	(Lz)));

 // Stores coverage of fibres over the atria
 vector<vector<vector<bool> > > complete(Lx, vector<vector<bool> >
                                           (Ly, vector<bool>
                                           (Lz)));
 // Stores flattened final fibre map	   
 vector<double> cellsData;

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

// Loads and unflattens the flat fibre orientation dataset of form (i, j, k, v1, v2, v3, i, j, k, v1, v2, v3, ...) stored as .npy
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
							
						for(int n = 0; n < rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getxs().size(); n++){
								
							xOther = rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getxs()[n];
							yOther = rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getys()[n];
							zOther = rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getzs()[n];
									
							distance = pow(pow(x - xOther, 2) + pow(y - yOther, 2) + pow(z - zOther, 2), 0.5);
									
									
							if(rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getFibIndex()[n] != cells[i].getFib() && distance < 2){
										
								cells[i].addProximalCell(rVoxels[xVoxel + p][yVoxel + q][zVoxel + r].getCellIndex()[n], distance);
										
										
							}
						}
					}
				}
		    }
		}	
	}
}

// Returns whether a fibre vertex is too close (<dSep) away from a vertex on a different fibre
bool tooClose(float i, float j, float k){

    int iVoxel = floor(i);
    int jVoxel = floor(j);
    int kVoxel = floor(k);

    float xOther;
    float yOther;
    float zOther;

    for(int p = -1; p < 2; p++){
        for(int q = -1; q < 2; q++){
            for(int r = -1; r < 2; r++){
                if(inBoundaries(iVoxel + p, jVoxel + q, kVoxel + r)){
                    for(int n = 0; n < voxels[iVoxel + p][jVoxel + q][kVoxel + r].getxs().size(); n++){

                        xOther = voxels[iVoxel + p][jVoxel + q][kVoxel + r].getxs()[n];
                        yOther = voxels[iVoxel + p][jVoxel + q][kVoxel + r].getys()[n];
                        zOther = voxels[iVoxel + p][jVoxel + q][kVoxel + r].getzs()[n];

                        if(pow(xOther - i, 2) + pow(yOther - j, 2) + pow(zOther - k, 2) < pow(dSep, 2)){return 1;}
                    }
                }
            }
        }
    }
    return 0;
}

// Generates a seed randomly placed on the atrial volume
void generateSeed(){

    bool seeded = 0;
    float xInitSeed;
    float yInitSeed;
    float zInitSeed;
    while(seeded == 0){
        xInitSeed = rnd()*(Lx-1);
        yInitSeed = rnd()*(Ly-1);
        zInitSeed = rnd()*(Lz-1);
        if(full(xInitSeed, yInitSeed, zInitSeed) && (voxels[floor(xInitSeed)][floor(yInitSeed)][floor(zInitSeed)].getxs().size() == 0)){seeded = 1;}
    }

    xSeeds.push_back(xInitSeed);
    ySeeds.push_back(yInitSeed);
    zSeeds.push_back(zInitSeed);

}

// Generates the forwards or backwards part of a fibre
void makeFibre(float x, float y, float z, bool forwards){

    int xVoxel = floor(x);
    int yVoxel = floor(y);
    int zVoxel = floor(z);

    int currentDirection;

    vector<float> prevData;

    bool running = 1;

    if(forwards){
        tempFibre.clearFibre();
        xVertices.clear();
        yVertices.clear();
        zVertices.clear();
    }

    if(!inBoundaries(x, y, z)){running = 0;}
    else if(!full(x, y, z)){running = 0;}
    else if(tooClose(x, y, z)){running = 0;}

    if(running){

        xVertices.push_back(x);
        yVertices.push_back(y);
        zVertices.push_back(z);

        if(forwards){
            tempFibre.addPosF(x, y, z);
            prevData.push_back(data[0][xVoxel][yVoxel][zVoxel]);
            prevData.push_back(data[1][xVoxel][yVoxel][zVoxel]);
            prevData.push_back(data[2][xVoxel][yVoxel][zVoxel]);
        }
        else{

            tempFibre.addPosB(x, y, z);
            prevData.push_back(-data[0][xVoxel][yVoxel][zVoxel]);
            prevData.push_back(-data[1][xVoxel][yVoxel][zVoxel]);
            prevData.push_back(-data[2][xVoxel][yVoxel][zVoxel]);
        }
    }

    while(running){


        if (prevData[0]*data[0][xVoxel][yVoxel][zVoxel]+prevData[1]*data[1][xVoxel][yVoxel][zVoxel]+prevData[2]*data[2][xVoxel][yVoxel][zVoxel] > 0){
            currentDirection = 1;
        }
        else {
            currentDirection = -1;
        }


        x = x + dStep * currentDirection * data[0][xVoxel][yVoxel][zVoxel];
        y = y + dStep * currentDirection * data[1][xVoxel][yVoxel][zVoxel];
        z = z + dStep * currentDirection * data[2][xVoxel][yVoxel][zVoxel];

        xVoxel = floor(x);
        yVoxel = floor(y);
        zVoxel = floor(z);

        if(!inBoundaries(x, y, z)){running = 0;}
        else if(!full(x, y, z)){running = 0;}
		else if(tooClose(x, y, z)){running = 0;}
        else if(prevData[0]*currentDirection*data[0][xVoxel][yVoxel][zVoxel]+prevData[1]*currentDirection*data[1][xVoxel][yVoxel][zVoxel]+prevData[2]*currentDirection*data[2][xVoxel][yVoxel][zVoxel] <= 0.766){running = 0;}
        


        if(running){

            xVertices.push_back(x);
            yVertices.push_back(y);
            zVertices.push_back(z);

            if(forwards){tempFibre.addPosF(x, y, z);}
            else{tempFibre.addPosB(x, y, z);}


            prevData[0] = currentDirection * data[0][xVoxel][yVoxel][zVoxel];
            prevData[1] = currentDirection * data[1][xVoxel][yVoxel][zVoxel];
            prevData[2] = currentDirection * data[2][xVoxel][yVoxel][zVoxel];
        }
    }

    if(!forwards){

        if(tempFibre.getxF().size() > threshold || tempFibre.getxB().size() > threshold){

            totalLength += tempFibre.getxF().size() + tempFibre.getxB().size() - 1;

            for(int n = 0; n < xVertices.size(); n++){

                voxels[floor(xVertices[n])][floor(yVertices[n])][floor(zVertices[n])].addVertex(xVertices[n], yVertices[n], zVertices[n], fibres.size(), 0);

            }

            fibres.push_back(tempFibre);
        }
    }
}


int main (int argc, char **argv) {

    srand (time(NULL));

    loadData();

    fibreCount = 0;
    fullCount = 0;

    for(int i = 0; i < Lx; i++){ // Count full voxels
        for(int j = 0; j < Ly; j++){
            for(int k = 0; k < Lz; k++){
                complete[i][j][k] = 0;
                if(full(i, j, k)){fullCount += 1;}
            }
        }
    }

    cout << "There are " << fullCount << " full voxels." << endl;

    for(int p = 0; p < fullCount * seedsFactor; p++){generateSeed();}

    cout << "Seeds generated" << endl;

    while(xSeeds.size() > 0 && threshold >= 0){ // Fibre creation

        if(xSeeds.size() % 100000 == 0){cout << "nFibres = " << fibres.size() << ", nSeeds = " << xSeeds.size() << endl;}

        if(xSeeds.size() > 0 && threshold >= 0){ // Make new fibre as long as seeds exist
            makeFibre(xSeeds[0], ySeeds[0], zSeeds[0], 1);
            makeFibre(xSeeds[0], ySeeds[0], zSeeds[0], 0);

            xSeeds.erase(xSeeds.begin());
            ySeeds.erase(ySeeds.begin());
            zSeeds.erase(zSeeds.begin());

        }

        if(xSeeds.size() == 0 && threshold > 0){
			
			// Changes fibre length threshold

            for(int p = 0; p < fullCount * seedsFactor; p++){generateSeed();}

            threshold-=increment;

        }

		else if(xSeeds.size() == 0 && threshold == 0){

			// Tests coverage

			nOnes = 0;
			nZeros = 0;

			for(int n = 0; n < fibres.size(); n++){

				for(int l = 0; l < fibres[n].getxF().size();l++){

					complete[floor(fibres[n].getxF()[l])][floor(fibres[n].getyF()[l])][floor(fibres[n].getzF()[l])] = 1;

				}

				for(int v = 0; v < fibres[n].getxB().size(); v++){

					complete[floor(fibres[n].getxB()[v])][floor(fibres[n].getyB()[v])][floor(fibres[n].getzB()[v])] = 1;
				}
			}

			for(int i = 0; i < Lx; i++){
				for(int j = 0; j < Ly; j++){
					for(int k = 0; k < Lz; k++){
						if(data[0][i][j][k]!=0 ||data[1][i][j][k]!=0 ||data[2][i][j][k]!=0){
							if(complete[i][j][k]){nOnes++;}
							else{nZeros++;}
						}
					}
				}
			}
			
			if(nOnes/(nOnes + nZeros) < 0.99){
				
				cout << "Coverage still " << nOnes/(nOnes + nZeros) << endl;
				
				for(int p = 0; p < fullCount; p++){generateSeed();}
			}
			
			
			
		}
    }

    int firstCell;

	// Convert fibre vertices (dStep units apart) to nodes (1 voxel apart) and store in cell class
	
    for(int i = 0; i < fibres.size(); i++){ 

        firstCell = cells.size();

        for(int f = 0; f < fibres[i].getxF().size(); f+=1){

            cells.push_back(cell(fibres[i].getxF()[f], fibres[i].getyF()[f], fibres[i].getzF()[f], i));
	    rVoxels[floor(fibres[i].getxF()[f])][floor(fibres[i].getyF()[f])][floor(fibres[i].getzF()[f])].addVertex(fibres[i].getxF()[f], fibres[i].getyF()[f], fibres[i].getzF()[f], i, cells.size() - 1);

            if(f > 0){
                cells[cells.size() - 1].addcon(cells.size() - 2);
                cells[cells.size() - 2].addcon(cells.size() - 1);
            }
        }

        for(int b = 1; b < fibres[i].getxB().size(); b+=1){

            cells.push_back(cell(fibres[i].getxB()[b], fibres[i].getyB()[b], fibres[i].getzB()[b], i));
	    rVoxels[floor(fibres[i].getxB()[b])][floor(fibres[i].getyB()[b])][floor(fibres[i].getzB()[b])].addVertex(fibres[i].getxB()[b], fibres[i].getyB()[b], fibres[i].getzB()[b], i, cells.size() - 1);
	
            if (b == 1){
                cells[cells.size()-1].addcon(firstCell);
                cells[firstCell].addcon(cells.size()-1);

            }
            if(b > 1){
                cells[cells.size() - 1].addcon(cells.size() - 2);
                cells[cells.size() - 2].addcon(cells.size() - 1);
            }
        }
    }
	

    // cout all measures regarding fibre map

    cout << "Coverage " << nOnes/(nOnes + nZeros) << endl;
    cout << "Number of fibres " << fibres.size() << endl;
    cout << "Mean nodes per fibre " << (float)totalLength/(float)fibres.size() << endl;



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

    aoba::SaveArrayAsNumpy("SAVE PATH FOR FLAT .NPY FIBRE MAP HERE", cellsData);
	
    return 0;
}




