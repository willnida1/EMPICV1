/*Header to initialize numerical values*/

#ifndef HEADER_H
#define HEADER_H

#include <iostream>

#define PI0 3.14159265359
#define mass_e 9.10938356e-31
#define mass_p 1.6726219e-27
#define eps 8.85418782e-12
#define ech 1.60217662e-19
#define cph 3.0e8
#define boltz 1.38064852e-23
#define mu0 1.25663706212e-6 // H/m

#endif

double mass_i;

// Domain Information
double dt, tmax;
int t;
double dx, xmin, xmax;
double dy, ymin, ymax;
double dz, zmin, zmax;
int boundaryConditionX; // 1 = Periodic BC, 2 = Wall BCs
int boundaryConditionY; // 1 = Periodic BC, 2 = Wall BCs
int boundaryConditionZ; // 1 = Periodic BC, 2 = Wall BCs
double gridArea;
double gridVolume;
int Nx, Ny, Nz;

double numDimensions;

// Particles and Grid
//----------------------------
// moved to headerPartGrid.h

// Particle # Information
//----------------------------
int numMacroElectrons; // Number of Electron Macroparticles
int numMacroIons;      // Number of Ion Macroparticles
int numMacroParticles; // number of macroparticles
int numPartPerMacro;   // particles per macroparticle
int numMacroPerCell;

double densityIon;
double densityElectron;

double debyeLength;

// Particle Thermal Information
//----------------------------
double Te_init; // Thermal values for elec and ions eV
double Ti_init;
double T_background;

double vth_e; // Electron Thermal Vel
double vth_i; // Ion Thermal Vel

// Plasma Information
//----------------------------

double freqPlasma;

// misc
//----------------------------
int nrel, nmax, nout, nout_print, num_avg;
int restart;

void initialize()
{

    T_background = 300; // K

    vth_i = 0.01 * cph;
    vth_e = 0.01 * cph;

    // Temperatures
    Ti_init = mass_i * pow(vth_i, 2) / boltz;
    Te_init = mass_e * pow(vth_e, 2) / boltz;

    densityIon = densityElectron = 1e16; //  #/m^3
    debyeLength = sqrt((boltz * eps * Te_init) / (densityElectron * pow(ech, 2)));

    // [0 - start from init | 1 - restart from file]
    restart = 0;

    Nx = 64;
    Ny = 64;
    Nz = 1;

    numDimensions = 2;

    xmin = 0;
    xmax = 2 * PI0 * debyeLength;
    ymin = 0;
    ymax = 2 * PI0 * debyeLength;
    zmin = 0;
    zmax = 0.02;

    if (numDimensions == 3)
    {
        zmax = 2 * PI0 * debyeLength;
        Nz = 64;
    }

    dx = (xmax - xmin) / ((double)Nx);
    dy = (ymax - ymin) / ((double)Ny);
    dz = (zmax - zmin) / ((double)Nz);
    gridArea = dx * dy; // * dz;
    gridVolume = dx * dy * dz;

    boundaryConditionX = 1; // 1 = periodic, 2 = wall
    boundaryConditionY = 1;
    boundaryConditionZ = 1;

    numMacroPerCell = 25;
    numMacroParticles = 2 * numMacroPerCell * Nx * Ny * Nz;
    numMacroElectrons = numMacroParticles / 2;
    numMacroIons = numMacroParticles / 2;
    numPartPerMacro = densityElectron / numMacroElectrons;

    freqPlasma = sqrt((densityElectron * pow(ech, 2)) / (eps * mass_e));

    dt = 0.125 / freqPlasma / 500;
    tmax = dt * 100000;
    nmax = (int)(tmax / dt);
}

// functions for calculating distributed E field at any coordinates
double ExFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;
    double Ex = 0;
    return Ex;
}
double EyFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;
    double Ey = 0;
    return Ey;
}
double EzFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;

    double Ez = 0;
    return Ez;
}

double BxFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;
    double Bx = 0; // T
    return Bx;
}
double ByFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;
    double By = 0; // T
    return By;
}

double BzFieldDistribution(double x1, double y1, double z1)
{
    x1;
    y1;
    z1;
    double Bz = 0;
    return Bz;
}
