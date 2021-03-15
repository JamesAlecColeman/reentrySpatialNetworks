#include <stdlib.h>
#include <iostream>
#include <vector>

class cell {

private:
    int state;
    bool defective;
    std::vector<int> cons;

    double x;
    double y;
    double z;

    int fib;

    std::vector<double> cellDistances;
    std::vector<int> proximalCells;

    std::vector<int> timesExcited;


public:
    std::vector<int> getcons();
    std::vector<double> getCellDistances();
    std::vector<int> getProximalCells();
    std::vector<int> getTimesExcited();

    void removeCon(int ind);

    bool getdef();
    int getstate();

    void addTimeExcited(int time);

    void clearTimesExcited();
    void clearCons();

    int getFib();

    void addProximalCell(int index, double distance);

    double getx();
    double gety();
    double getz();
    void setstate(int state);
    void setdef(bool defective);
    void addcon(float con);

    cell(double xx, double yy, double zz, int i);
};
