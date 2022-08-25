#include <iostream>

class YeeGrid
{

public:
    double nx;
    double ny;
    double nz;

    double originx;
    double originy;
    double originz;

    ThreeVec E;
    ThreeVec EPlusOne;

    ThreeVec BMinusHalf;
    ThreeVec BPlusHalf;
    ThreeVec BPlusThreeHalf;

    ThreeVec J;

    double rho;

    YeeGrid(int cellx, int celly, int cellz)
    {
        nx = cellx;
        ny = celly;
        nz = cellz;
    }

    void incrementRho(double rhoIncrement)
    {
        rho += rhoIncrement;
    }
};
