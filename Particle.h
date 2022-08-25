#include <iostream>
#include <math.h>
#include "threevector.h"

class Particle
{

public:
    // Index of which Yee Grid Particle is located in
    int partType; // 0 = neutral, 1 = ion, 2 = electron

    int yeeGridNx;
    int yeeGridNy;
    int yeeGridNz;

    // Index of which Yee Grid Particle is located in after push
    int yeeGridNxPlusOne;
    int yeeGridNyPlusOne;
    int yeeGridNzPlusOne;

    // E and B fields weighted to particle @ each time step
    ThreeVec E;
    ThreeVec B;

    // Particle Position with respect to main origin
    double x;
    double y;
    double z;

    double xPlusOne;
    double yPlusOne;
    double zPlusOne;

    double deltaX;
    double deltaY;
    double deltaZ;

    // Used for ZigZag method
    double addXs;
    double addYs;
    double addZs;

    // Used for ZigZag method, store location before adjusting for periodic BC (Can be out of domain)
    double preAdjustmentX;
    double preAdjustmentY;
    double preAdjustmentZ;

    // Particle position with respect to local Yee Grid

    // TODO: Replace local coords with ThreeVec
    double localX;
    double localY;
    double localZ;

    double localXPlusOne;
    double localYPlusOne;
    double localZPlusOne;

    // Particle Momentums
    ThreeVec pMinusHalf;
    ThreeVec pPlusHalf;

    // Particle Velocities
    ThreeVec vMinusHalf;
    ThreeVec vPlusHalf;

    double q;
    double m;
    int isMobile; // Weibel Instability: 0 = no, 1 = yes

    Particle(double xInit, double yInit, double zInit, int type, int mobility)
    {
        partType = type;

        x = xInit;
        y = yInit;
        z = zInit;

        // Initialization function for particle's local yee grid
        yeeGridNx = floor(x / dx);
        yeeGridNy = floor(y / dy);
        yeeGridNz = floor(z / dz);

        localX = fmod(x, dx);
        localY = fmod(y, dy);
        localZ = fmod(z, dz);

        if (mobility == 0)
        {
            isMobile = 0;
            yeeGridNxPlusOne = yeeGridNx;
            yeeGridNyPlusOne = yeeGridNy;
            yeeGridNzPlusOne = yeeGridNz;

            xPlusOne = x;
            yPlusOne = y;
            zPlusOne = z;

            localXPlusOne = localX;
            localYPlusOne = localY;
            localZPlusOne = localZ;
        }
        else
        {
            isMobile = 1;
        }

        // Initialization function for particle's local displacement from grid origin
    }

    void updateAfterPush(double xNew, double yNew, double zNew)
    {
        if (isMobile == 1)
        {
            xPlusOne = xNew;
            yPlusOne = yNew;
            zPlusOne = zNew;

            // Initialization function for particle's local yee grid
            yeeGridNxPlusOne = floor(xNew / dx);

            yeeGridNyPlusOne = floor(yNew / dy);
            yeeGridNzPlusOne = floor(zNew / dz);

            // Initialization function for particle's local displacement from grid origin
            localXPlusOne = fmod(xNew, dx);
            localYPlusOne = fmod(yNew, dy);
            localZPlusOne = fmod(zNew, dz);
        }
    }

    void updateAfterLoop()
    {
        // reset all variables
        if (isMobile == 1)
        {
            x = xPlusOne;
            xPlusOne = 0;
            y = yPlusOne;
            yPlusOne = 0;
            z = zPlusOne;
            zPlusOne = 0;

            yeeGridNx = yeeGridNxPlusOne;
            yeeGridNxPlusOne = 0;
            yeeGridNy = yeeGridNyPlusOne;
            yeeGridNyPlusOne = 0;
            yeeGridNz = yeeGridNzPlusOne;
            yeeGridNzPlusOne = 0;

            localX = localXPlusOne;
            localXPlusOne = 0;
            localY = localYPlusOne;
            localYPlusOne = 0;
            localZ = localZPlusOne;
            localZPlusOne = 0;

            deltaX = 0;
            deltaY = 0;
            deltaZ = 0;

            addXs = 0;
            addYs = 0;
            addZs = 0;

            preAdjustmentX = 0;
            preAdjustmentY = 0;
            preAdjustmentZ = 0;

            pMinusHalf.setX(pPlusHalf.getX());
            pPlusHalf.setX(0);
            pMinusHalf.setY(pPlusHalf.getY());
            pPlusHalf.setY(0);
            pMinusHalf.setZ(pPlusHalf.getZ());
            pPlusHalf.setZ(0);

            vMinusHalf.setX(vPlusHalf.getX());
            vPlusHalf.setX(0);
            vMinusHalf.setY(vPlusHalf.getY());
            vPlusHalf.setY(0);
            vMinusHalf.setZ(vPlusHalf.getZ());
            vPlusHalf.setZ(0);

            E.setX(0);
            E.setY(0);
            E.setZ(0);

            B.setX(0);
            B.setY(0);
            B.setZ(0);
        }
    }
};
