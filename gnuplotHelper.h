

// Plotting functions for visualization of data
#include <chrono>
#include <thread>
#include <iostream>

int io = 0;
int jo = 0;

static void plotResults3D(double *xData, double *yData, double *zData, int dataSize)
{
    FILE *gnuplotPipe, *tempDataFile;
    char *tempDataFileName;
    double x;
    double y;
    double z;

    int i;
    tempDataFileName = "tempData";
    gnuplotPipe = popen("gnuplot", "w");

    if (gnuplotPipe)
    {

        tempDataFile = fopen(tempDataFileName, "w");

        for (i = 0; i <= dataSize; i++)
        {
            x = xData[i];
            y = yData[i];
            z = zData[i];
            fprintf(tempDataFile, "%lf %lf %lf\n", x, y, z);
        }

        fclose(tempDataFile);
        fprintf(gnuplotPipe, "set xrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set yrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set zrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set zrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set view 70, 60, 1, 1;\n", tempDataFileName);
        fprintf(gnuplotPipe, "set terminal png size 1600,1600; set output 'Particles/particles%04d.png';\n", io, tempDataFileName);
        io++;
        fprintf(gnuplotPipe, "splot \"%s\" \n", tempDataFileName);
        fflush(gnuplotPipe);
        // fclose(tempDataFile);
        // printf("press enter to continue...");
        // getchar();
        // remove(tempDataFileName);
        // fprintf(gnuplotPipe, "exit \n");
    }
    else
    {
        printf("gnuplot not found...");
    }
}

static void plotResults(double *xData, double *yData, int dataSize)
{
    FILE *gnuplotPipe, *tempDataFile;
    char *tempDataFileName;
    double x;
    double y;

    int i;
    tempDataFileName = "tempData";
    gnuplotPipe = popen("gnuplot", "w");

    if (gnuplotPipe)
    {

        tempDataFile = fopen(tempDataFileName, "w");

        for (i = 0; i <= dataSize; i++)
        {
            x = xData[i];
            y = yData[i];
            fprintf(tempDataFile, "%lf %lf\n", x, y);
        }

        fclose(tempDataFile);
        fprintf(gnuplotPipe, "set xrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set yrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set terminal png size 1600,1600; set output 'Particles/particles%04d.png';\n", io, tempDataFileName);
        io++;
        fprintf(gnuplotPipe, "plot \"%s\" \n", tempDataFileName);
        fflush(gnuplotPipe);
        // fclose(tempDataFile);
        // printf("press enter to continue...");
        // getchar();
        // remove(tempDataFileName);
        // fprintf(gnuplotPipe, "exit \n");
    }
    else
    {
        printf("gnuplot not found...");
    }
}

static void linePlotResults(double *xData, double *yData, int dataSize)
{
    FILE *gnuplotPipe, *tempDataFile;
    char *tempDataFileName;
    double x;
    double y;

    int i;
    tempDataFileName = "fieldLineData";
    gnuplotPipe = popen("gnuplot", "w");

    if (gnuplotPipe)
    {

        tempDataFile = fopen(tempDataFileName, "w");

        for (i = 0; i <= dataSize; i++)
        {
            x = xData[i];
            y = yData[i];
            fprintf(tempDataFile, "%lf %lf\n", x, y);
        }

        fclose(tempDataFile);
        fprintf(gnuplotPipe, "set xrange [0 : 2]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set yrange [0 : 12]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set terminal png size 1600,1600; set output 'Fields/bField%04d.png';\n", io, tempDataFileName);
        io++;
        fprintf(gnuplotPipe, "plot \"%s\" \n", tempDataFileName);
        fflush(gnuplotPipe);
        // fclose(tempDataFile);
        // printf("press enter to continue...");
        // getchar();
        // remove(tempDataFileName);
        // fprintf(gnuplotPipe, "exit \n");
    }
    else
    {
        printf("gnuplot not found...");
    }
}

static void heatPlotResults(double dx, double dy, int numCol, int numRow, double *val, double dt) //, double error)
{
    FILE *gnuplotPipe, *tempDataFile;
    char *tempDataFileName;
    double x, y;
    int i;
    tempDataFileName = "fieldData";
    gnuplotPipe = popen("gnuplot", "w");

    if (gnuplotPipe)
    {
        tempDataFile = fopen(tempDataFileName, "w");

        int k = 0;
        for (int j = 0; j < numRow; j++)
        {
            for (int l = 0; l < numCol; l++)
            {
                fprintf(tempDataFile, "%lf %lf %.30f\n", ((double)(l)*dx), ((double)(j)*dy), val[k]);
                k++;
            }
        }
        fclose(tempDataFile);

        fprintf(gnuplotPipe, "set xrange [0 : 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set yrange [0: 0.0034]\n", tempDataFileName);
        fprintf(gnuplotPipe, "set terminal png size 1600,1600; set output 'field%04d.png'; set title 'Weibel Instability -- Time: %04f nanoseconds' font 'Arial{,15}'; set xlabel 'X (m)'; set ylabel 'Y (m)'; set cblabel 'B (T)';\n", jo, jo * dt * 1e9 * 20, tempDataFileName);
        jo++;
        fprintf(gnuplotPipe, "plot \"%s\" with image \n", tempDataFileName);
        fflush(gnuplotPipe);
    }

    else
    {
        printf("gnuplot not found...");
    }
}
