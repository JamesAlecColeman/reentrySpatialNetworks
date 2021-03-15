#include "Cell.h"
#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

std::vector<int> cell::getcons(){return cell::cons;}
std::vector<double> cell::getCellDistances(){return cell::cellDistances;}
std::vector<int> cell::getProximalCells(){return cell::proximalCells;}
std::vector<int> cell::getTimesExcited(){return cell::timesExcited;}
void cell::removeCon(int ind){cell::cons.erase(cell::cons.begin() + ind);}


bool cell::getdef(){return cell::defective;}
int cell::getstate(){return cell::state;}

void cell::addTimeExcited(int time){cell::timesExcited.push_back(time);}

void cell::clearTimesExcited(){cell::timesExcited.clear();}
void cell::clearCons(){cell::cons.clear();}

int cell::getFib(){return cell::fib;}

double cell::getx(){return cell::x;}
double cell::gety(){return cell::y;}
double cell::getz(){return cell::z;}
void cell::setstate(int state){cell::state = state;}
void cell::setdef(bool defective){cell::defective=defective;}
void cell::addcon(float con){cell::cons.push_back(con);}


void cell::addProximalCell(int index, double distance){
    cell::proximalCells.push_back(index);
    cell::cellDistances.push_back(distance);
}

cell::cell(double xx, double yy, double zz, int i){
    cell::x=xx;
    cell::y=yy;
    cell::z=zz;
    cell::state=0;
    cell::defective=0;
    cell::fib=i;
}


