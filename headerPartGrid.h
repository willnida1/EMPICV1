#include <iostream>

#include "header.h"
#include "Particle.h"
#include "YeeGrid.h"
#include <vector>

YeeGrid ***yeegrid = NULL;
Particle *particles = NULL;

void initializePart()
{
    // initialize particle in memory

    particles = (Particle *)malloc(numMacroParticles * sizeof(Particle));

    double beamWidth = (xmax - xmin) * 1 / 5;
    double beamSeparation = ((xmax - xmin) - (2 * beamWidth)) / 16;

    // Initialize Particles
    double yPartE = ymin;
    double xPartE = ((xmax - xmin) - (beamWidth * 2 + beamSeparation)) / 2; // initialize half of beam separation from bottom of domain for periodic BC
    double zPartE = (zmax - zmin) / 2 - beamWidth / 2;

    double partSeparationY = dy / sqrt(numMacroPerCell);
    double partSeparationX = (beamWidth) / (sqrt(numMacroPerCell) * Nx); //(dy / sqrt(numMacroPerCell) / beamWidth / 2);
    double partSeparationZ;

    if (numDimensions == 3)
    {
        partSeparationY = dy / cbrt(numMacroPerCell);
        partSeparationX = 2 * (beamWidth) / (cbrt(numMacroPerCell) * Nx); //(dy / sqrt(numMacroPerCell) / beamWidth / 2);
        partSeparationZ = 2 * (beamWidth) / (cbrt(numMacroPerCell) * Nz);
    }
    // WEIBEL: Initialize Ions and electrons on top of one another, make Ions stationary

    for (int n = 0; n < numMacroElectrons; n++)
    {

        if (yPartE > ymax)
        {
            yPartE = ymin;
            xPartE += partSeparationX;
        }

        if (numDimensions == 3 && xPartE > ((xmax - xmin) + beamSeparation) / 2)
        {
            xPartE = beamSeparation / 2;
            zPartE += partSeparationZ;
        }

        int type = 2;     // electron value
        int isMobile = 1; // Mobile

        // Maxwellian distribution for thermal velocity of electron beams

        double vXThermal = sqrt(fabs(-2 * log((double)rand() / RAND_MAX) * vth_e * cos(2 * PI0 * ((double)rand() / RAND_MAX))));

        double vYThermal = 4 * sqrt(fabs(-2 * log((double)rand() / RAND_MAX) * vth_e * sin(2 * PI0 * ((double)rand() / RAND_MAX))));

        double gammaLorentzN = sqrt(1 + (pow((0.8 * cph + vXThermal), 2) / pow(cph, 2)));

        if (n < numMacroElectrons / 2)
        {

            // initialize bottom half of particles moving right
            if (numDimensions < 3)
                particles[n] = Particle(xPartE, yPartE, 0, type, isMobile);
            else if (numDimensions == 3)
                particles[n] = Particle(xPartE, yPartE, zPartE, type, isMobile); // Mobile Species
            particles[n].q = ech * -1 * numPartPerMacro;
            particles[n].m = mass_e * numPartPerMacro;
            // Lower Half, Rightward streaming Electron Beam
            particles[n].pMinusHalf.setX(particles[n].m * vXThermal);
            particles[n].pMinusHalf.setY(particles[n].m * gammaLorentzN * (0.8 * cph + vYThermal));
            if (numDimensions == 3)
                particles[n].pMinusHalf.setZ(particles[n].m * vYThermal);
        }
        else if (n == (int)(numMacroElectrons / 2) || n == (int)(numMacroElectrons / 2) + 1 || n == (int)(numMacroElectrons / 2) - 1)
        {
            yPartE = ymin;
            xPartE = ((xmax - xmin) - (beamWidth * 2 + beamSeparation)) / 2;
            zPartE = (zmax - zmin) / 2 - beamWidth / 2;
        }
        else
        {
            // initialize top half of particles moving left
            if (numDimensions < 3)
                particles[n] = Particle(xPartE, yPartE, 0, type, isMobile);
            if (numDimensions == 3)
                particles[n] = Particle(xPartE, yPartE, zPartE, type, isMobile); // Mobile Species
            particles[n].q = ech * -1 * numPartPerMacro;
            particles[n].m = mass_e * numPartPerMacro;
            // Upper Half, Leftward streaming Electron Beam
            particles[n].pMinusHalf.setX(particles[n].m * vXThermal);
            particles[n].pMinusHalf.setY(particles[n].m * -1 * gammaLorentzN * (0.8 * cph + vYThermal));
            if (numDimensions == 3)
                particles[n].pMinusHalf.setZ(particles[n].m * vYThermal);
        }

        yPartE += partSeparationY;
    }

    // Initialize Ions on top of electrons, uniform distribution

    double xPartI = xmin;
    double yPartI = ymin;
    double zPartI = zmin;

    for (int j = numMacroElectrons; j < numMacroParticles; j++)
    {
        if (xPartI >= xmax)
        {
            xPartI = xmin;
            yPartI += partSeparationY;
        }

        if (numDimensions == 3 && yPartI >= ymax)
        {
            yPartI = ymin;
            zPartI += partSeparationZ;
        }

        int type = 1;     // ion value
        int isMobile = 0; // non Mobile species
        particles[j] = Particle(xPartI, yPartI, zPartI, type, isMobile);
        particles[j].q = ech * numPartPerMacro;
        particles[j].m = mass_i * numPartPerMacro;

        xPartI += partSeparationX;
    }
}

void initializeGrid()
{
    // initialize array in memory
    yeegrid = (YeeGrid ***)malloc(Nx * sizeof(YeeGrid));
    for (int i = 0; i < Nx; i++)
    {
        yeegrid[i] = (YeeGrid **)malloc(Ny * sizeof(YeeGrid));
        for (int j = 0; j < Ny; j++)
        {
            yeegrid[i][j] = (YeeGrid *)malloc(Nz * sizeof(YeeGrid));
        }
    }

    // Initialize Yee Grids
    for (int ix = 0; ix < Nx; ix++)
    {
        for (int iy = 0; iy < Ny; iy++)
        { // 0 to 511, 512 N
            for (int iz = 0; iz < Nz; iz++)
            {

                // Cell Indicies
                yeegrid[ix][iy][iz].nx = ix;
                yeegrid[ix][iy][iz].ny = iy;
                yeegrid[ix][iy][iz].nz = iz;

                // Cell Origin Coordinates
                yeegrid[ix][iy][iz].originx = yeegrid[ix][iy][iz].nx * dx;
                yeegrid[ix][iy][iz].originy = yeegrid[ix][iy][iz].ny * dy;
                yeegrid[ix][iy][iz].originz = yeegrid[ix][iy][iz].nx * dz;

                // Initialize all values to zero to avoid null
                // Set test distribution fn for Ex
                double m = 8;
                double n = 6;
                double a = (m * PI0) / (double)(xmax - xmin);
                double b = (n * PI0) / (double)(ymax - ymin);

                double r = sqrt(pow((ix + 0.5) * dx, 2) + pow((iy + 0.5) * dy, 2));

                // yeegrid[ix][iy][iz].E.setX(exp(-25 * fabs((iy * dy) - (0.5 * dy))));

                // if (ix == 0 && iy == 0)
                //  std::cout << exp(-25 * fabs((iy * dy) - (0.5 * dy))) << std::endl;

                yeegrid[ix][iy][iz].E.setX(0);
                yeegrid[ix][iy][iz].E.setY(0);
                yeegrid[ix][iy][iz].E.setZ(0);

                yeegrid[ix][iy][iz].EPlusOne.setX(0);
                yeegrid[ix][iy][iz].EPlusOne.setY(0);
                yeegrid[ix][iy][iz].EPlusOne.setZ(0);

                yeegrid[ix][iy][iz].J.setX(0);
                yeegrid[ix][iy][iz].J.setY(0);
                yeegrid[ix][iy][iz].J.setZ(0);

                yeegrid[ix][iy][iz].BMinusHalf.setX(0);
                yeegrid[ix][iy][iz].BMinusHalf.setY(0);
                // yeegrid[ix][iy][iz].BMinusHalf.setZ(0);
                yeegrid[ix][iy][iz].BMinusHalf.setZ(cos((2 * PI0 * (0.5 + iy) * dy) / (ymax - ymin)));
                //   yeegrid[ix][iy][iz].BMinusHalf.setZ(exp(-25 * pow(r, 2)));
                //  yeegrid[ix][iy][iz].BMinusHalf.setZ(sin(((0.5 + ix) * dx * a)) * cos((0.5 + iy) * dy * b));
                //  yeegrid[ix][iy][iz].BMinusHalf.setZ(cos(((0.5 + ix) * dx * PI0 * m) / (xmax - xmin)) * cos(((0.5 + iy) * dy * PI0 * n) / (ymax - ymin)));

                yeegrid[ix][iy][iz].BPlusHalf.setX(0);
                yeegrid[ix][iy][iz].BPlusHalf.setY(0);
                yeegrid[ix][iy][iz].BPlusHalf.setZ(cos((2 * PI0 * (0.5 + iy) * dy) / (ymax - ymin)));

                // yeegrid[ix][iy][iz].BPlusHalf.setZ(0);
                //   yeegrid[ix][iy][iz].BPlusHalf.setZ(exp(-25 * pow(r, 2)));
                //   yeegrid[ix][iy][iz].BPlusHalf.setZ(sin(((0.5 + ix) * dx * a)) * cos((0.5 + iy) * dy * b));
                //   yeegrid[ix][iy][iz].BPlusHalf.setZ(cos(((0.5 + ix) * dx * PI0 * m) / (xmax - xmin)) * cos(((0.5 + iy) * dy * PI0 * n) / (ymax - ymin)));

                yeegrid[ix][iy][iz].BPlusThreeHalf.setX(0);
                yeegrid[ix][iy][iz].BPlusThreeHalf.setY(0);
                yeegrid[ix][iy][iz].BPlusThreeHalf.setZ(0);
            }
        }
    }
    /*
        std::vector<double> bValues;
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                bValues.push_back(yeegrid[i][j][0].BPlusHalf.getZ());
            }
        }

        double *bFinal = &bValues[0];

        heatPlotResults(dx, dy, Nx, Ny, bFinal, Nx, Ny);
        */
}