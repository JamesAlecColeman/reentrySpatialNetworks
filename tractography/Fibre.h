#include <stdlib.h>
#include <iostream>
#include <vector>

class fibre {

private:

    std::vector<double> xF;
    std::vector<double> yF;
    std::vector<double> zF;

    std::vector<double> xB;
    std::vector<double> yB;
    std::vector<double> zB;

    double fibreLength;

public:

    void addPosF(double xxF, double yyF, double zzF);
    void addPosB(double xxB, double yyB, double zzB);

    std::vector<double> getxF();
    std::vector<double> getyF();
    std::vector<double> getzF();

    std::vector<double> getxB();
    std::vector<double> getyB();
    std::vector<double> getzB();

    void addLength(double length);

    double getFibreLength();

    void clearFibre();

    fibre();
};
