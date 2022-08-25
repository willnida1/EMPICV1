// Main File to run

#include <iostream>
#include "math.h"
#include "gnuplotHelper.h"
#include "picGrid.h"
#include <algorithm>
#include <vector>
//#include "include/engine.h"

using namespace std;

// Gnuplot results

// run initialize

int main()
{
    std::cout << std::fixed << std::setprecision(17);

    initialize();
    initializePart();

    initializeGrid();

    cout << "C * dt: " << cph * dt << endl;
    cout << "dx: " << dx << endl;
    cout << freqPlasma << endl;

    vector<double> xValues(1);
    vector<double> yValues(1);
    vector<double> zValues(1);

    vector<double> bValues;
    vector<double> tValues;

    // cout << particles[0].yeeGridNx << endl;'
    double errorTotal = 0;

    double timeReal = 0;
    double rhoTotal = 0;

    for (t = 0; t < nmax; t++)
    {

        double bTotal = 0;

        fieldGatherer();

        particlePusher();

        currentDeposer();

        fieldSolver();

        // double dEdx = (yeegrid[10][10][0].EPlusOne.getX() - yeegrid[9][10][0].EPlusOne.getX()) / dx;
        // double dEdy = (yeegrid[10][10][0].EPlusOne.getY() - yeegrid[10][9][0].EPlusOne.getY()) / dy;
        // cout << yeegrid[10][10][0].rho / eps << endl;

        // cout << yeegrid[10][20][0].BPlusThreeHalf.getZ() << endl;

        // yValues.push_back(yeegrid[10][20][0].BPlusThreeHalf.getZ());
        /*
                double error = 0;

                double m = 8;
                double n = 6;
                double a = (m * PI0) / (xmax - xmin);
                double b = (n * PI0) / (ymax - ymin);
                double w = cph * sqrt(pow(a, 2) + pow(b, 2));
                */
        for (int k = 0; k < Nz; k++)
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++)
                {
                    bTotal += fabs(yeegrid[i][j][0].BPlusThreeHalf.getZ());
                    // if (t == 200)
                    //   bValues.push_back(yeegrid[i][j][0].BPlusThreeHalf.getZ());

                    // error += fabs((sin(((0.5 + i) * dx * a)) * cos((0.5 + j) * dy * b) * cos(w * (t * dt) + (3 * dt / 2))) - yeegrid[i][j][0].BPlusThreeHalf.getZ());
                    // For TE Mode, Wall BC: error += fabs((cos(((0.5 + i) * dx * PI0 * m) / (xmax - xmin)) * cos(((0.5 + j) * dy * PI0 * n) / (ymax - ymin)) * cos(w * (t * dt) + (3 * dt / 2))) - yeegrid[i][j][0].BPlusThreeHalf.getZ());
                }
            }
        }

        // TE Mode Diagnostics
        // cout << error / Nx / Ny << endl;
        //   cout << "=====" << endl;
        // errorTotal += error / Nx / Ny;

        if (t % 20 == 0)
        {

            for (int ip = 0; ip < numMacroElectrons; ip++)
            {

                xValues.push_back(particles[ip].x);
                yValues.push_back(particles[ip].y);
                zValues.push_back(particles[ip].z);
            }

            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < Nx; i++)
                {
                    // bTotal += fabs(yeegrid[i][j][0].BPlusThreeHalf.getZ());
                    // bValues.push_back(yeegrid[i][j][0].BPlusThreeHalf.getZ());
                }
            }

            // xValues.push_back(dt * t * freqPlasma);
            // yValues.push_back(bTotal / (Nx * Ny));

            cout << "B Average:  " << bTotal / (Nx * Ny * Nz) << endl;
            cout << "Frequency Adjusted Time: " << dt * t * freqPlasma << endl;

            bValues.push_back(bTotal / (Nx * Ny * Nz));
            tValues.push_back(dt * t * freqPlasma);

            double *bFinal = &bValues[0];
            double *tFinal = &tValues[0];
            double *xFinal = &xValues[0];
            double *yFinal = &yValues[0];
            double *zFinal = &zValues[0];

            if (numDimensions < 3)
                plotResults(xFinal, yFinal, numMacroElectrons);

            else
                plotResults3D(xFinal, yFinal, zFinal, numMacroElectrons);

            linePlotResults(tFinal, bFinal, t);

            // heatPlotResults(dx, dy, Nx, Ny, bFinal, dt); //, error / Nx / Ny);
        }

        xValues.clear();
        yValues.clear();
        zValues.clear();
        pushUpdateObjectVariables();

        timeReal += dt;
    }
}
// ffmpeg -i particles%04d.png WeibelPartStable64PPC.mpeg
// ffmpeg -i field%04d.png WeibelFieldStable144PPC.mpeg
