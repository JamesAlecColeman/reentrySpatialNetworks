#include <stdlib.h>
#include <iostream>
#include <vector>

class voxel {

private:
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    std::vector<int> fibIndex;
    std::vector<int> cellIndex;

public:
    std::vector<double> getxs();
    std::vector<double> getys();
    std::vector<double> getzs();
    std::vector<int> getFibIndex();
    std::vector<int> getCellIndex();

    void addVertex(float xVertex, float yVertex, float zVertex, int fibreIndex, int cellNumber);

    voxel();
};
