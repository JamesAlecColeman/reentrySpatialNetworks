#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include "numpy.hpp"
#include <math.h>
#include "Cell.h"
#include <unistd.h>
#include <algorithm>

using namespace std;

 const double g = 3; // Coarse graining factor

 // Healthy sheep
 //const int Lx = (int)ceil(325/g);
 //const int Ly = (int)ceil(949/g);
 //const int Lz = (int)ceil(611/g);

 // Human
 //const int Lx = (int)ceil(303/g);
 //const int Ly = (int)ceil(275/g);
 //const int Lz = (int)ceil(244/g);

 // HF sheep
 const int Lx = (int)ceil(265/g);
 const int Ly = (int)ceil(1035/g);
 const int Lz = (int)ceil(623/g);

 double ch = 0.4400; // Characteristic distance (sigmoid parameter)
 string chch = "04400/";
 
 string directory = "SAVE PATH HERE";
 
 // CMP MODEL PARAMETERS
 int SA = 1000; // SA Period
 int T = 83; // Refractory Period
 const float delta = 1; // Fraction of cells susceptible to conduction block
 float epsilon = 0.0025; // Probability of conduction block
 int stopTime = 100000; // Number of time steps without re-entry after which to stop the model
 int nCircuits = 15; // Number of circuits to collect before terminating

 int nExcited = 0;
 int reentrant = 0;
 int saveIDint;
 string saveIDstring;
 int counter = 0;
 bool AF = 0;
 bool running = 1;

 vector<int> SAcells; // Indices of cells in SA node
 vector<cell> cells;

 // Coordinates of cells (for saving)
 vector<float> xCells;
 vector<float> yCells;
 vector<float> zCells;

 vector<int> frequency; // Number of re-entrant circuits at each cell

 vector<int> timeToReentry; // Stores number of timesteps at which re-entry occurs
 vector<int> circuitLength; // Stores lengths of re-entrant circuits

 vector<int> active; // Cells that are active
 vector<int> activeNew; // Cells that will be active in the next time step

 vector<int> cellsToExcite;

 int seedIndex;
 int minTime;

// Returns random number between 0 and 1
float rnd(){return static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * 1;}

// Allocates cells susceptible to conduction block
void doDefectiveCells(){
    for(int i = 0; i < cells.size(); i++){
        if(rnd() <= delta){cells[i].setdef(1);}
        else{cells[i].setdef(0);}
    }
}

void doTransverseConns(){
     double dist;
     double sigProb;
     double s = 7;

    // Remove connections between cells not on the same fibre
    for(int i = 0; i < cells.size(); i++){
        for(int p = 0; p < cells[i].getcons().size(); p++){
            if(cells[cells[i].getcons()[p]].getiFib() == cells[i].getiFib() && cells[cells[i].getcons()[p]].getjFib() == cells[i].getjFib() && cells[cells[i].getcons()[p]].getkFib() == cells[i].getkFib()){}
            else{cells[i].removeCon(p);p--;}
        }
    }
    // Allocates transverse connections (If two cells are proximal, only one of them has the other stored as a proximal cell)
    for(int i = 0; i < cells.size(); i++){
        for(int b = 0; b < cells[i].getProximalCells().size(); b++){
            dist = cells[i].getCellDistances()[b];
            sigProb = 1/(1 + exp(s * (dist - ch)));
            if(rnd() < sigProb){
                cells[i].addcon(cells[i].getProximalCells()[b]);
                cells[cells[i].getProximalCells()[b]].addcon(i);
            }
        }
    }
}

// Decides whether a cell connected to an active cell should be excited
void tryActivate(int index){

    if(cells[index].getstate() == 0){
        if(cells[index].getdef() == 1){
            if(rnd() >= epsilon){
                cells[index].setstate(T + 1);
                activeNew.push_back(index);
                cells[index].addTimeExcited(counter);
                nExcited++;
            }
        }
        else{
            cells[index].setstate(T + 1);
            activeNew.push_back(index);
            cells[index].addTimeExcited(counter);
            nExcited++;
        }
    }
}

void step(){
    cellsToExcite.clear();
    activeNew.clear();

    // Finds cells connected to active cells
    for(int m = 0; m < active.size(); m++){
        if(cells[active[m]].getstate() == T + 1){
            vector<int> cons = cells[active[m]].getcons();
            for(int c = 0; c < cons.size(); c++){
                cellsToExcite.push_back(cons[c]);
            }
        }
    }

    // Tries to excite cells connected to active cells
    sort(cellsToExcite.begin(), cellsToExcite.end());
    for (int m = 0; m < cellsToExcite.size(); m++){
        if (m==0){
            tryActivate(cellsToExcite[0]);
        }
        else if (cellsToExcite[m]!=cellsToExcite[m-1]){
            tryActivate(cellsToExcite[m]);
        }
    }

    // Changes the state of cells -1 and adds cells that are still excited to activeNew
    for(int m = 0; m < active.size(); m++){
        cells[active[m]].setstate(cells[active[m]].getstate() - 1);
        if(cells[active[m]].getstate() != 0){activeNew.push_back(active[m]);}
    }

	// Periodically excites SA cells
    if(counter % SA == 0){
        for(int v = 0; v < SAcells.size(); v++){
            if(cells[SAcells[v]].getstate() == 0){activeNew.push_back(SAcells[v]);}
            cells[SAcells[v]].setstate(T + 1);
            cells[SAcells[v]].addTimeExcited(counter);
            nExcited++;
        }
    }

	// Updates active cells
    active = activeNew;

	// At the end of every heartbeat, checks if re-entry has occurred
    if(counter % SA == (SA-1)){
       
        nExcited = 0;
        AF = 0;
        for(int i = 0; i < cells.size(); i++){
            if(cells[i].getTimesExcited().size()>1) {
                AF = 1;
            }
        }

		// If re-entry has occured, finds the first cell to be excited twice
        if(AF){
 	    nExcited = 0;
	    reentrant = reentrant + 1;
	    timeToReentry.push_back(counter);

	    
            minTime = counter+1;
            for(int i = 0; i < cells.size(); i++){
                if(cells[i].getTimesExcited().size()>1) {
                    if(cells[i].getTimesExcited()[1] < minTime){
                        minTime = cells[i].getTimesExcited()[1];
                        seedIndex = i;
                    }
                }
            }
            circuitLength.push_back(cells[seedIndex].getTimesExcited()[1] - cells[seedIndex].getTimesExcited()[0]);    
	    
            frequency[seedIndex]++;
            running = 0;
        }

        for(int i = 0; i < cells.size(); i++){
            cells[i].clearTimesExcited();
        }
    }

    // Frame counter
    counter = counter + 1;

}


int main (int argc, char **argv) {
	

    // Generates a seed
    srand (time(NULL) + getpid());

    saveIDint = (int)((rnd()) * 100000);
    stringstream temp;
	temp << saveIDint;
    saveIDstring = temp.str();

	// Stores the raw data containing all information about the cells
    vector<int> flatDataMeta;
    vector<double> flatData;

	// Loads raw data
    aoba::LoadArrayFromNumpy(directory+"YOUR FIBRE MAP FILE HERE.npy", flatDataMeta, flatData);

    int cellStart = 0;
    int nCells = 0;
    bool going = 1;
    int nCons;
    int nProximal;

	// Unpacks the raw data to reconstruct the cells
    while(going){

        cells.push_back(cell(flatData[cellStart], flatData[cellStart+1],flatData[cellStart+2],flatData[cellStart+3],flatData[cellStart+4],flatData[cellStart+5]));
        nCons = flatData[cellStart + 6];

        for(int p = cellStart + 7; p < cellStart + 7 + nCons; p++){
            cells[nCells].addcon(flatData[p]);
        }

        nProximal = flatData[cellStart + 7 + nCons];

        for(int k = cellStart + 7 + nCons + 1; k < cellStart + 7 + nCons + 1 + nProximal; k++){
            cells[nCells].addProximalCell(flatData[k + nProximal + 1], flatData[k]);
        }

        cellStart = cellStart + 7 + nCons + 1 + nProximal + nProximal + 1;

        if (cellStart > flatData.size()){
            going = 0;
        }

        nCells++;
    }

    cout << cells.size() << "<-----" << endl;

	// Sets the SA node cells
	
	// oldSheep SA
	/*for(int i = 0; i < cells.size(); i++){
		if(pow(cells[i].getx() - 120/g, 2) + pow(cells[i].gety() - 750/g, 2) + pow(cells[i].getz() - 466/g, 2) < 324/pow(g,2)){
			SAcells.push_back(i);
		}
	}*/
	
	
	// human SA
	/*(for(int i = 0; i < cells.size(); i++){
		if(pow(cells[i].getx() - 147/g, 2) + pow(cells[i].gety() - 213/g, 2) + pow(cells[i].getz() - 200/g, 2) < 324/pow(g,2)){
			SAcells.push_back(i);
		}
	}*/

	// newSheep SA
	for(int i = 0; i < cells.size(); i++){
		if(pow(cells[i].getx() - 120/g, 2) + pow(cells[i].gety() - 300/g, 2) + pow(cells[i].getz() - 150/g, 2) < 324/pow(g,2)){
			SAcells.push_back(i);
		}
	}




	// Stores the coordinates of the cells and makes a vector to store the frequency of re-entrant circuits at each cell
    for(int f = 0; f < cells.size(); f++){
        xCells.push_back(cells[f].getx());
        yCells.push_back(cells[f].gety());
        zCells.push_back(cells[f].getz());
        frequency.push_back(0);
    }

	// Saves the cell coordinates

    aoba::SaveArrayAsNumpy(directory+chch+"xs.npy", xCells);
    aoba::SaveArrayAsNumpy(directory+chch+"ys.npy", yCells);
    aoba::SaveArrayAsNumpy(directory+chch+"zs.npy", zCells);

	// Runs the model and updates the frequency and times to re-entry file
    while(reentrant < nCircuits){

        counter = 0;
        running = 1;

		// Allocates transverse connections and defective cells
        doTransverseConns();
        doDefectiveCells();

		// Sets states to 0, clears times excited and active cells
        for(int j = 0; j < cells.size(); j++){
            cells[j].setstate(0);
            cells[j].clearTimesExcited();
        }
        active.clear();

        while(running){
            step();
			// Stops the model if re-entry has not occurred after a chosen time
            if(counter > stopTime){
				running = 0;
				timeToReentry.push_back(-1);
			}
        }

		// Saves the frequencies and times to re-entry
        aoba::SaveArrayAsNumpy(directory+chch+"freqs"+saveIDstring+".npy", frequency);
	aoba::SaveArrayAsNumpy(directory+chch+"timestoreentry"+saveIDstring+".npy", timeToReentry);
	aoba::SaveArrayAsNumpy(directory+ chch + "circuitlengths"+saveIDstring+".npy", circuitLength);
    }

    return 0;
}
