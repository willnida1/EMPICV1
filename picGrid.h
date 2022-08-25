/*
                                                    ===========
                                                    picGrid.h
                                                    ===========

This file contains all functions pertaining to calculating and updating values of both the particles and
the Yee Grids.

- William Nida, 2022

*/

/*
---------------
Main TODO Items
_______________

1. Go through code to include things like exterior imposed E or B field --> most likely finished
2. Implement initialization function for when code first starts to prepare for boris push, decide values of E, B, etc.
    - Where to start? likely assign initial momentum to particles, push, and continue

*/

#include <iostream>
#include "math.h"
#include "headerPartGrid.h"
#include <algorithm>

void fieldGatherer()
{
    // Gathers field data and weights it to each macro particle for t = n, using En, and averaging B(n + 1/2) & B(n - 1/2)
    // TODO: 2D implementation for now, 3D later

    //============
    // Gather E
    //============

    // IMPORTANT: Could simplify:
    //  Check what quadrant particle is in, determine boxes to pull values from based on location

    //===============
    //      Ex
    //===============

    for (int ip = 0; ip < numMacroParticles; ip++)
    {

        if (particles[ip].isMobile == 1)
        {

            int PartGridXLeft = particles[ip].yeeGridNx;
            int PartGridYBottom = particles[ip].yeeGridNy;
            int PartGridZ = 0;
            if (numDimensions == 3)
            {
                PartGridZ = particles[ip].yeeGridNz;
            }

            int PartGridXRight = particles[ip].yeeGridNx;
            int PartGridYTop = particles[ip].yeeGridNy + 1;
            int PartGridZInc = 0;

            if (numDimensions == 3)
            {
                PartGridZInc = particles[ip].yeeGridNz + 1;
                if (PartGridZ == Nz - 1)
                    PartGridZInc = 0;
            }

            double adjustedX = particles[ip].localX;
            double adjustedY = particles[ip].localY;
            double adjustedZ = particles[ip].localZ;

            if (particles[ip].localX >= (dx / 2))
            {

                adjustedX -= (dx / 2);
                PartGridXLeft = particles[ip].yeeGridNx;
                PartGridXRight = particles[ip].yeeGridNx + 1;
                if (PartGridXLeft == Nx - 1)
                    PartGridXRight = 0;
            }

            else if (particles[ip].localX < (dx / 2))
            {

                adjustedX += (dx / 2);
                PartGridXLeft = particles[ip].yeeGridNx - 1;
                PartGridXRight = particles[ip].yeeGridNx;
                if (PartGridXRight == 0)
                    PartGridXLeft = Nx - 1;
            }

            if (PartGridYBottom >= Ny - 1)
                PartGridYTop = 0;

            double NetExParticle;

            if (numDimensions < 3)
            {
                NetExParticle = (adjustedX * adjustedY / gridArea) * yeegrid[PartGridXRight][PartGridYTop][PartGridZ].E.getX()                   /* bottom left box allocated to top right field */
                                + (dx - adjustedX * adjustedY / gridArea) * yeegrid[PartGridYTop][PartGridXLeft][PartGridZ].E.getX()             /* bottom right box allocated to top left field */
                                + (adjustedX * (dy - adjustedY) / gridArea) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZ].E.getX()       /* top left box to bottom right */
                                + ((dx - adjustedX) * (dy - adjustedY) / gridArea) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZ].E.getX() /* top right box to bottom left */
                                + ExFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);                                        // Add field Distribution
            }
            else if (numDimensions == 3)
            {

                // std::cout << PartGridZ << std::endl;
                // std::cout << PartGridZInc << std::endl;
                // std::cout << ip << std::endl;

                //   std::cout << "a" << std::endl;

                // std::cout << PartGridXLeft << std::endl;
                // std::cout << PartGridXRight << std::endl;
                //   std::cout << "b" << std::endl;

                // std::cout << PartGridYBottom << std::endl;
                // std::cout << PartGridYTop << std::endl;
                // std::cout << "c" << std::endl;

                NetExParticle = (adjustedX * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZInc].E.getX()                    /* bottom left box allocated to top right field */
                                + (dx - adjustedX * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZInc].E.getX()              /* bottom right box allocated to top left field */
                                + (adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZInc].E.getX()        /* top left box to bottom right */
                                + ((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZInc].E.getX(); /* top right box to bottom left */

                // Lower Level Additions

                NetExParticle += (adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZ].E.getX(); /* bottom left box allocated to top right field */

                NetExParticle += (dx - adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZ].E.getX(); /* bottom right box allocated to top left field */

                NetExParticle += (adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZ].E.getX(); /* top left box to bottom right */

                NetExParticle += ((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZ].E.getX() /* top right box to bottom left */

                                 + ExFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }

            particles[ip].E.setX(NetExParticle);
        }
    }

    //===============
    //      Ey
    //===============

    for (int ip = 0; ip < numMacroParticles; ip++)
    {

        if (particles[ip].isMobile == 1)
        {

            int PartGridXLeft = particles[ip].yeeGridNx;
            int PartGridYBottom = particles[ip].yeeGridNy;
            int PartGridZ = 0;
            if (numDimensions == 3)
            {
                PartGridZ = particles[ip].yeeGridNz;
            }

            int PartGridXRight = particles[ip].yeeGridNx + 1;
            int PartGridYTop = particles[ip].yeeGridNy;
            int PartGridZInc = 0;

            if (numDimensions == 3)
            {
                PartGridZInc = particles[ip].yeeGridNz + 1;
            }

            double adjustedX = particles[ip].localX;
            double adjustedY = particles[ip].localY;
            double adjustedZ = particles[ip].localZ;

            if (particles[ip].localY >= (dy / 2))
            {
                adjustedY -= (dy / 2);
                PartGridYTop = particles[ip].yeeGridNy + 1;
                PartGridYBottom = particles[ip].yeeGridNy;
                if (PartGridYBottom == Ny - 1)
                    PartGridYTop = 0;
            }

            else if (particles[ip].localY < (dy / 2))
            {

                adjustedY += (dy / 2);
                PartGridYTop = particles[ip].yeeGridNy;
                PartGridYBottom = particles[ip].yeeGridNy - 1;
                if (PartGridYTop == 0)
                    PartGridYBottom = Ny - 1;
            }
            if (PartGridXLeft == Nx - 1)
                PartGridXRight = 0;

            double NetEyParticle;
            if (numDimensions < 3)
            {
                NetEyParticle = (adjustedX * adjustedY / gridArea) * yeegrid[PartGridXRight][PartGridYTop][PartGridZ].E.getY()                   /* bottom left box allocated to top right field */
                                + ((dx - adjustedX) * adjustedY / gridArea) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZ].E.getY()           /* bottom right box allocated to top left field */
                                + (adjustedX * (dy - adjustedY) / gridArea) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZ].E.getY()       /* top left box to bottom right */
                                + ((dx - adjustedX) * (dy - adjustedY) / gridArea) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZ].E.getY() /* top right box to bottom left */
                                + EyFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
                // TODO: Implement double NetEzParticle and modify prior two
            }

            else if (numDimensions == 3)
            {
                NetEyParticle = (adjustedX * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZInc].E.getY()                   /* bottom left box allocated to top right field */
                                + ((dx - adjustedX) * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZInc].E.getY()           /* bottom right box allocated to top left field */
                                + (adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZInc].E.getY()       /* top left box to bottom right */
                                + ((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZInc].E.getY() /* top right box to bottom left */

                                // Lower Level Additions

                                + (adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZ].E.getY()                 /* bottom left box allocated to top right field */
                                + ((dx - adjustedX) * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZ].E.getY()           /* bottom right box allocated to top left field */
                                + (adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZ].E.getY()       /* top left box to bottom right */
                                + ((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZ].E.getY() /* top right box to bottom left */

                                + EyFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }

            particles[ip].E.setY(NetEyParticle);
        }
    }

    //==================
    // Ez Field
    //==================
    if (numDimensions == 3)
    {

        for (int ip = 0; ip < numMacroParticles; ip++)
        {

            if (particles[ip].isMobile == 1)
            {

                int PartGridXLeft = particles[ip].yeeGridNx;
                int PartGridXRight = particles[ip].yeeGridNx + 1;

                // Apply Periodic BC
                if (particles[ip].yeeGridNx == Nx - 1)
                {
                    PartGridXRight = 0;
                }

                int PartGridYBottom = particles[ip].yeeGridNy;
                int PartGridYTop = particles[ip].yeeGridNy + 1;
                if (particles[ip].yeeGridNy == Ny - 1)
                {
                    PartGridYTop = 0;
                }

                double adjustedX = particles[ip].localX;
                double adjustedY = particles[ip].localY;
                double adjustedZ = particles[ip].localZ;

                int PartGridZBottom;
                int PartGridZTop;

                if (particles[ip].localZ >= dz / 2)
                {
                    // If in upper half
                    PartGridZBottom = particles[ip].yeeGridNz;
                    PartGridZTop = particles[ip].yeeGridNz + 1;

                    // Periodic BC
                    if (PartGridZBottom == Nz - 1)
                    {
                        PartGridZTop = 0;
                    }

                    adjustedZ -= (dz / 2);
                }
                else if (particles[ip].localZ >= dz / 2)
                {
                    // If in lower half
                    PartGridZBottom = particles[ip].yeeGridNz - 1;
                    PartGridZTop = particles[ip].yeeGridNz;

                    // Periodic BC
                    if (PartGridZTop == 0)
                    {
                        PartGridZBottom = Nz - 1;
                    }

                    adjustedZ += (dz / 2);
                }

                double NetEzParticle = (adjustedX * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZTop].E.getZ()                   /* bottom left box allocated to top right field */
                                       + ((dx - adjustedX) * adjustedY * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZTop].E.getZ()           /* bottom right box allocated to top left field */
                                       + (adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZTop].E.getZ()       /* top left box to bottom right */
                                       + ((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZTop].E.getZ() /* top right box to bottom left */

                                       // Lower Level Additions

                                       + (adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYTop][PartGridZBottom].E.getZ()                 /* bottom left box allocated to top right field */
                                       + ((dx - adjustedX) * adjustedY * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYTop][PartGridZBottom].E.getZ()           /* bottom right box allocated to top left field */
                                       + (adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXRight][PartGridYBottom][PartGridZBottom].E.getZ()       /* top left box to bottom right */
                                       + ((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * yeegrid[PartGridXLeft][PartGridYBottom][PartGridZBottom].E.getZ() /* top right box to bottom left */

                                       + EyFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }
        }
    }

    //==================
    // Bz FIELD
    //==================

    for (int ip = 0; ip < numMacroParticles; ip++)
    {

        if (particles[ip].isMobile == 1)
        {

            int quadrant;

            int bottomIndex;
            int topIndex;
            int leftIndex;
            int rightIndex;

            int zBottomIndex = particles[ip].yeeGridNz;
            int zTopIndex = particles[ip].yeeGridNz + 1;
            // Periodic
            if (zBottomIndex == Nz - 1)
            {
                zTopIndex = 0;
            }

            double adjustedX;
            double adjustedY;
            double adjustedZ = particles[ip].localZ;

            int PartGridZ = 0;

            if (particles[ip].localX >= (dx / 2) && particles[ip].localY >= (dy / 2))
                quadrant = 1;
            else if (particles[ip].localX < (dx / 2) && particles[ip].localY >= (dy / 2))
                quadrant = 2;
            else if (particles[ip].localX < (dx / 2) && particles[ip].localY < (dy / 2))
                quadrant = 3;
            else if (particles[ip].localX >= (dx / 2) && particles[ip].localY < (dy / 2))
                quadrant = 4;

            if (quadrant == 1)
            {
                bottomIndex = particles[ip].yeeGridNy;
                topIndex = particles[ip].yeeGridNy + 1;
                if (bottomIndex == Ny - 1)
                    topIndex = 0;

                leftIndex = particles[ip].yeeGridNx;
                rightIndex = particles[ip].yeeGridNx + 1;
                if (leftIndex == Nx - 1)
                    rightIndex = 0;

                adjustedX = particles[ip].localX - (dx / 2);
                adjustedY = particles[ip].localY - (dy / 2);
            }
            else if (quadrant == 2)
            {
                bottomIndex = particles[ip].yeeGridNy;
                topIndex = particles[ip].yeeGridNy + 1;
                if (bottomIndex == Ny - 1)
                    topIndex = 0;

                leftIndex = particles[ip].yeeGridNx - 1;
                rightIndex = particles[ip].yeeGridNx;
                if (rightIndex == 0)
                    leftIndex = Nx - 1;

                adjustedX = particles[ip].localX + (dx / 2);
                adjustedY = particles[ip].localY - (dy / 2);
            }
            else if (quadrant == 3)
            {
                bottomIndex = particles[ip].yeeGridNy - 1;
                topIndex = particles[ip].yeeGridNy;
                if (topIndex == 0)
                    bottomIndex = Ny - 1;

                leftIndex = particles[ip].yeeGridNx - 1;
                rightIndex = particles[ip].yeeGridNx;
                if (rightIndex == 0)
                    leftIndex = Nx - 1;

                adjustedX = particles[ip].localX + (dx / 2);
                adjustedY = particles[ip].localY + (dy / 2);
            }
            else if (quadrant == 4)
            {
                bottomIndex = particles[ip].yeeGridNy - 1;
                topIndex = particles[ip].yeeGridNy;
                if (topIndex == 0)
                    bottomIndex = Ny - 1;

                leftIndex = particles[ip].yeeGridNx;
                rightIndex = particles[ip].yeeGridNx + 1;
                if (leftIndex == Nx - 1)
                    rightIndex = 0;

                adjustedX = particles[ip].localX - (dx / 2);
                adjustedY = particles[ip].localY + (dy / 2);
            }

            double NetBzParticle = ((adjustedX * adjustedY / gridArea) * (yeegrid[rightIndex][topIndex][PartGridZ].BMinusHalf.getZ() + yeegrid[rightIndex][topIndex][PartGridZ].BPlusHalf.getZ()) / 2)                     // bottom left to top right
                                   + (((dx - adjustedX) * adjustedY / gridArea) * (yeegrid[leftIndex][topIndex][PartGridZ].BMinusHalf.getZ() + yeegrid[leftIndex][topIndex][PartGridZ].BPlusHalf.getZ()) / 2)              // bottom right to top left
                                   + ((adjustedX * (dy - adjustedY) / gridArea) * (yeegrid[rightIndex][bottomIndex][PartGridZ].BMinusHalf.getZ() + yeegrid[rightIndex][bottomIndex][PartGridZ].BPlusHalf.getZ()) / 2)      // top left to bottom right
                                   + (((dx - adjustedX) * (dy - adjustedY) / gridArea) * (yeegrid[leftIndex][bottomIndex][PartGridZ].BMinusHalf.getZ() + yeegrid[leftIndex][bottomIndex][PartGridZ].BPlusHalf.getZ()) / 2) // top right to bottom left
                                   + BzFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);

            if (numDimensions == 3)
            {
                NetBzParticle = ((adjustedX * adjustedY * adjustedZ / gridVolume) * (yeegrid[rightIndex][topIndex][zTopIndex].BMinusHalf.getZ() + yeegrid[rightIndex][topIndex][zTopIndex].BPlusHalf.getZ()) / 2)                     // bottom left to top right
                                + (((dx - adjustedX) * adjustedY * adjustedZ / gridVolume) * (yeegrid[leftIndex][topIndex][zTopIndex].BMinusHalf.getZ() + yeegrid[leftIndex][topIndex][zTopIndex].BPlusHalf.getZ()) / 2)              // bottom right to top left
                                + ((adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[rightIndex][bottomIndex][zTopIndex].BMinusHalf.getZ() + yeegrid[rightIndex][bottomIndex][zTopIndex].BPlusHalf.getZ()) / 2)      // top left to bottom right
                                + (((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[leftIndex][bottomIndex][zTopIndex].BMinusHalf.getZ() + yeegrid[leftIndex][bottomIndex][zTopIndex].BPlusHalf.getZ()) / 2) // top right to bottom left

                                // Bottom addition

                                + ((adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[rightIndex][topIndex][zBottomIndex].BMinusHalf.getZ() + yeegrid[rightIndex][topIndex][zBottomIndex].BPlusHalf.getZ()) / 2)                   // bottom left to top right
                                + (((dx - adjustedX) * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[leftIndex][topIndex][zBottomIndex].BMinusHalf.getZ() + yeegrid[leftIndex][topIndex][zBottomIndex].BPlusHalf.getZ()) / 2)              // bottom right to top left
                                + ((adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[rightIndex][bottomIndex][zBottomIndex].BMinusHalf.getZ() + yeegrid[rightIndex][bottomIndex][zBottomIndex].BPlusHalf.getZ()) / 2)      // top left to bottom right
                                + (((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[leftIndex][bottomIndex][zBottomIndex].BMinusHalf.getZ() + yeegrid[leftIndex][bottomIndex][zBottomIndex].BPlusHalf.getZ()) / 2) // top right to bottom left

                                + BzFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }

            particles[ip].B.setX(0);
            particles[ip].B.setY(0);
            particles[ip].B.setZ(NetBzParticle);
        }
    }

    //==========
    // Bx Field
    //==========

    if (numDimensions == 3)
    {

        for (int ip = 0; ip < numMacroParticles; ip++)
        {

            if (particles[ip].isMobile == 1)
            {

                int quadrant;

                int yIndexBottom;
                int yIndexTop;
                int zIndexBottom;
                int zIndexTop;

                double adjustedX = particles[ip].localX;
                double adjustedY;
                double adjustedZ;

                int xIndexLeft = particles[ip].yeeGridNx;
                int xIndexRight = particles[ip].yeeGridNx + 1;
                if (xIndexLeft == Nx - 1)
                {
                    xIndexRight = 0;
                }

                // Imagine plane projected on left face of yee cube
                if (particles[ip].localY >= (dy / 2) && particles[ip].localZ >= (dz / 2))
                    quadrant = 1;
                else if (particles[ip].localY < (dy / 2) && particles[ip].localZ >= (dz / 2))
                    quadrant = 2;
                else if (particles[ip].localY < (dy / 2) && particles[ip].localZ < (dz / 2))
                    quadrant = 3;
                else if (particles[ip].localY >= (dy / 2) && particles[ip].localZ < (dz / 2))
                    quadrant = 4;

                if (quadrant == 1)
                {
                    yIndexBottom = particles[ip].yeeGridNy;
                    yIndexTop = particles[ip].yeeGridNy + 1;
                    if (yIndexBottom == Ny - 1)
                        yIndexTop = 0;

                    zIndexBottom = particles[ip].yeeGridNz;
                    zIndexTop = particles[ip].yeeGridNz + 1;
                    if (zIndexBottom == Nz - 1)
                        zIndexTop = 0;

                    adjustedY = particles[ip].localY - (dy / 2);
                    adjustedZ = particles[ip].localZ - (dz / 2);
                }
                else if (quadrant == 2)
                {
                    yIndexBottom = particles[ip].yeeGridNy - 1;
                    yIndexTop = particles[ip].yeeGridNy;
                    if (yIndexTop == 0)
                        yIndexBottom = Ny - 1;

                    zIndexBottom = particles[ip].yeeGridNz;
                    zIndexTop = particles[ip].yeeGridNz + 1;
                    if (zIndexBottom == Nz - 1)
                        zIndexTop = 0;

                    adjustedY = particles[ip].localY + (dy / 2);
                    adjustedZ = particles[ip].localZ - (dz / 2);
                }
                else if (quadrant == 3)
                {
                    yIndexBottom = particles[ip].yeeGridNy - 1;
                    yIndexTop = particles[ip].yeeGridNy;
                    if (yIndexTop == 0)
                        yIndexBottom = Ny - 1;

                    zIndexBottom = particles[ip].yeeGridNz - 1;
                    zIndexTop = particles[ip].yeeGridNz;
                    if (zIndexTop == 0)
                        zIndexBottom = Nz - 1;

                    adjustedY = particles[ip].localY + (dy / 2);
                    adjustedZ = particles[ip].localZ + (dz / 2);
                }
                else if (quadrant == 4)
                {
                    yIndexBottom = particles[ip].yeeGridNy;
                    yIndexTop = particles[ip].yeeGridNy + 1;
                    if (yIndexBottom == Ny - 1)
                        yIndexTop = 0;

                    zIndexBottom = particles[ip].yeeGridNz - 1;
                    zIndexTop = particles[ip].yeeGridNz;
                    if (zIndexTop == 0)
                        zIndexBottom = Nz - 1;

                    adjustedY = particles[ip].localY - (dy / 2);
                    adjustedZ = particles[ip].localZ + (dz / 2);
                }
                /*
                std::cout << xIndexRight << std::endl;
                std::cout << "===" << std::endl;
                std::cout << yIndexTop << std::endl;
                std::cout << "===" << std::endl;
                std::cout << zIndexTop << std::endl;
                std::cout << "===" << std::endl;
                */
                double NetBzParticle = ((adjustedX * adjustedY * adjustedZ / gridVolume) * (yeegrid[xIndexRight][yIndexTop][zIndexTop].BMinusHalf.getX() + yeegrid[xIndexRight][yIndexTop][zIndexTop].BPlusHalf.getX()) / 2)                     // bottom left to top right
                                       + (((dx - adjustedX) * adjustedY * adjustedZ / gridVolume) * (yeegrid[xIndexLeft][yIndexTop][zIndexTop].BMinusHalf.getX() + yeegrid[xIndexLeft][yIndexTop][zIndexTop].BPlusHalf.getX()) / 2)              // bottom right to top left
                                       + ((adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[xIndexRight][yIndexBottom][zIndexTop].BMinusHalf.getX() + yeegrid[xIndexRight][yIndexBottom][zIndexTop].BPlusHalf.getX()) / 2)      // top left to bottom right
                                       + (((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[xIndexLeft][yIndexBottom][zIndexTop].BMinusHalf.getX() + yeegrid[xIndexLeft][yIndexBottom][zIndexTop].BPlusHalf.getX()) / 2) // top right to bottom left

                                       // Bottom addition

                                       + ((adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexRight][yIndexTop][zIndexBottom].BMinusHalf.getX() + yeegrid[xIndexRight][yIndexTop][zIndexBottom].BPlusHalf.getX()) / 2)                   // bottom left to top right
                                       + (((dx - adjustedX) * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexLeft][yIndexTop][zIndexBottom].BMinusHalf.getX() + yeegrid[xIndexLeft][yIndexTop][zIndexBottom].BPlusHalf.getX()) / 2)              // bottom right to top left
                                       + ((adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexRight][yIndexBottom][zIndexBottom].BMinusHalf.getX() + yeegrid[xIndexRight][yIndexBottom][zIndexBottom].BPlusHalf.getX()) / 2)      // top left to bottom right
                                       + (((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexLeft][yIndexBottom][zIndexBottom].BMinusHalf.getX() + yeegrid[xIndexLeft][yIndexBottom][zIndexBottom].BPlusHalf.getX()) / 2) // top right to bottom left

                                       + BxFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }
        }
    }

    //==========
    // By Field
    //==========
    if (numDimensions == 3)
    {

        for (int ip = 0; ip < numMacroParticles; ip++)
        {

            if (particles[ip].isMobile == 1)
            {

                int quadrant;

                int xIndexLeft;
                int xIndexRight;
                int zIndexBottom;
                int zIndexTop;

                double adjustedX;
                double adjustedY = particles[ip].localY;
                double adjustedZ;

                int yIndexBottom = particles[ip].yeeGridNy;
                int yIndexTop = particles[ip].yeeGridNy + 1;
                if (yIndexBottom == Ny - 1)
                {
                    yIndexTop = 0;
                }

                // Imagine plane projected on back face of yee cube
                if (particles[ip].localX >= (dx / 2) && particles[ip].localZ >= (dz / 2))
                    quadrant = 1;
                else if (particles[ip].localX < (dx / 2) && particles[ip].localZ >= (dz / 2))
                    quadrant = 2;
                else if (particles[ip].localX < (dx / 2) && particles[ip].localZ < (dz / 2))
                    quadrant = 3;
                else if (particles[ip].localX >= (dx / 2) && particles[ip].localZ < (dz / 2))
                    quadrant = 4;

                if (quadrant == 1)
                {
                    xIndexLeft = particles[ip].yeeGridNx;
                    xIndexRight = particles[ip].yeeGridNx + 1;
                    if (xIndexLeft == Nx - 1)
                        xIndexRight = 0;

                    zIndexBottom = particles[ip].yeeGridNz;
                    zIndexTop = particles[ip].yeeGridNz + 1;
                    if (zIndexBottom == Nz - 1)
                        zIndexTop = 0;

                    adjustedX = particles[ip].localX - (dx / 2);
                    adjustedZ = particles[ip].localZ - (dz / 2);
                }
                else if (quadrant == 2)
                {
                    xIndexLeft = particles[ip].yeeGridNx - 1;
                    xIndexRight = particles[ip].yeeGridNx;
                    if (xIndexRight == 0)
                        xIndexLeft = Nx - 1;

                    zIndexBottom = particles[ip].yeeGridNz;
                    zIndexTop = particles[ip].yeeGridNz + 1;
                    if (zIndexBottom == Nz - 1)
                        zIndexTop = 0;

                    adjustedX = particles[ip].localX + (dx / 2);
                    adjustedZ = particles[ip].localZ - (dz / 2);
                }
                else if (quadrant == 3)
                {
                    xIndexLeft = particles[ip].yeeGridNx - 1;
                    xIndexRight = particles[ip].yeeGridNx;
                    if (xIndexRight == 0)
                        xIndexLeft = Nx - 1;

                    zIndexBottom = particles[ip].yeeGridNz - 1;
                    zIndexTop = particles[ip].yeeGridNz;
                    if (zIndexTop == 0)
                        zIndexBottom = Nz - 1;

                    adjustedX = particles[ip].localX + (dx / 2);
                    adjustedZ = particles[ip].localZ + (dz / 2);
                }
                else if (quadrant == 4)
                {
                    xIndexLeft = particles[ip].yeeGridNx;
                    xIndexRight = particles[ip].yeeGridNx + 1;
                    if (xIndexLeft == Nx - 1)
                        xIndexRight = 0;

                    zIndexBottom = particles[ip].yeeGridNz - 1;
                    zIndexTop = particles[ip].yeeGridNz;
                    if (zIndexTop == 0)
                        zIndexBottom = Nz - 1;

                    adjustedX = particles[ip].localX - (dx / 2);
                    adjustedZ = particles[ip].localZ + (dz / 2);
                }

                double NetBzParticle = ((adjustedX * adjustedY * adjustedZ / gridVolume) * (yeegrid[xIndexRight][yIndexTop][zIndexTop].BMinusHalf.getY() + yeegrid[xIndexRight][yIndexTop][zIndexTop].BPlusHalf.getY()) / 2)                     // bottom left to top right
                                       + (((dx - adjustedX) * adjustedY * adjustedZ / gridVolume) * (yeegrid[xIndexLeft][yIndexTop][zIndexTop].BMinusHalf.getY() + yeegrid[xIndexLeft][yIndexTop][zIndexTop].BPlusHalf.getY()) / 2)              // bottom right to top left
                                       + ((adjustedX * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[xIndexRight][yIndexBottom][zIndexTop].BMinusHalf.getY() + yeegrid[xIndexRight][yIndexBottom][zIndexTop].BPlusHalf.getY()) / 2)      // top left to bottom right
                                       + (((dx - adjustedX) * (dy - adjustedY) * adjustedZ / gridVolume) * (yeegrid[xIndexLeft][yIndexBottom][zIndexTop].BMinusHalf.getY() + yeegrid[xIndexLeft][yIndexBottom][zIndexTop].BPlusHalf.getY()) / 2) // top right to bottom left

                                       // Bottom addition

                                       + ((adjustedX * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexRight][yIndexTop][zIndexBottom].BMinusHalf.getY() + yeegrid[xIndexRight][yIndexTop][zIndexBottom].BPlusHalf.getY()) / 2)                   // bottom left to top right
                                       + (((dx - adjustedX) * adjustedY * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexLeft][yIndexTop][zIndexBottom].BMinusHalf.getY() + yeegrid[xIndexLeft][yIndexTop][zIndexBottom].BPlusHalf.getY()) / 2)              // bottom right to top left
                                       + ((adjustedX * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexRight][yIndexBottom][zIndexBottom].BMinusHalf.getY() + yeegrid[xIndexRight][yIndexBottom][zIndexBottom].BPlusHalf.getY()) / 2)      // top left to bottom right
                                       + (((dx - adjustedX) * (dy - adjustedY) * (dz - adjustedZ) / gridVolume) * (yeegrid[xIndexLeft][yIndexBottom][zIndexBottom].BMinusHalf.getY() + yeegrid[xIndexLeft][yIndexBottom][zIndexBottom].BPlusHalf.getY()) / 2) // top right to bottom left

                                       + ByFieldDistribution(particles[ip].x, particles[ip].y, particles[ip].z);
            }
        }
    }
}

void particlePusher()
{
    // Update particle momentum @ p(n+1/2) and particle position @ x(n+1)
    for (int ip = 0; ip < numMacroParticles; ip++)
    {
        if (particles[ip].isMobile == 1)
        {
            //==================
            // Begin Boris Push
            //==================

            ThreeVec pMinus = particles[ip].pMinusHalf + particles[ip].E.operator*(particles[ip].q * dt / 2);
            // Based off of (possibly incorrect) https://arxiv.org/pdf/1104.3163.pdf eqn for gamma(n)
            double gammaLorentzN = sqrt(1 + (pow(pMinus.mag(), 2) / pow(cph, 2)));

            ThreeVec T = particles[ip].B.operator*(particles[ip].q * dt / 2 / gammaLorentzN / particles[ip].m);
            ThreeVec S = T.operator*(2 / (1 + pow(T.mag(), 2)));

            ThreeVec pPrime = pMinus.operator+(pMinus.operator^(T));

            ThreeVec pPlus = pMinus.operator+(pPrime.operator^(S));

            particles[ip].pPlusHalf = pPlus.operator+(particles[ip].E.operator*(particles[ip].q * dt / 2));

            particles[ip].vPlusHalf = particles[ip].pPlusHalf / gammaLorentzN / particles[ip].m;

            double xNew = particles[ip].x + (particles[ip].vPlusHalf.getX() * dt);
            double yNew = particles[ip].y + (particles[ip].vPlusHalf.getY() * dt);
            double zNew = particles[ip].z + (particles[ip].vPlusHalf.getZ() * dt);

            // std::cout << particles[ip].vPlusHalf.getX() * dt << std::endl;

            particles[ip].deltaX = xNew - particles[ip].x;
            particles[ip].deltaY = yNew - particles[ip].y;
            particles[ip].deltaZ = zNew - particles[ip].z;

            particles[ip].preAdjustmentX = xNew;
            particles[ip].preAdjustmentY = yNew;
            particles[ip].preAdjustmentZ = zNew;

            particles[ip].addXs = xNew + particles[ip].x;
            particles[ip].addYs = yNew + particles[ip].y;
            particles[ip].addZs = zNew + particles[ip].z;
            // std::cout << xNew << std::endl;

            // Allocate periodic BC to particles leaving simulation domain
            if (boundaryConditionX == 1)
            {
                if (xNew >= xmax)
                {
                    xNew = fmod(xNew - xmin, (xmax - xmin)) + xmin;
                }
                if (xNew < 0)
                {
                    xNew = xmax - fabs(fmod(fabs(xNew - xmin), (xmax - xmin)));
                }
            }

            if (boundaryConditionY == 1)
            {
                if (yNew >= ymax)
                {
                    yNew = fmod(yNew - ymin, (ymax - ymin)) + ymin;
                }
                if (yNew < 0)
                {
                    yNew = ymax - fabs(fmod(fabs(yNew - ymin), (ymax - ymin)));
                }
            }

            if (numDimensions == 3)
            {
                if (boundaryConditionZ == 1)
                {
                    if (zNew >= zmax)
                    {
                        zNew = fmod(zNew - zmin, (zmax - zmin)) + zmin;
                    }
                    if (zNew < 0)
                    {
                        zNew = zmax - fabs(fmod(fabs(zNew - zmin), (zmax - zmin)));
                    }
                }
            }

            particles[ip].updateAfterPush(xNew, yNew, zNew);
        }
    }
}

void currentDeposer()
{
    for (int ip = 0; ip < numMacroParticles; ip++)
    {
        if (particles[ip].isMobile == 1)
        {

            // Useful Vars

            double deltaX = particles[ip].deltaX;
            double deltaY = particles[ip].deltaY;
            double deltaZ = particles[ip].deltaZ;

            // Implement Zig-Zag Scheme

            // If the particle did cross box boundary
            if (particles[ip].yeeGridNx != particles[ip].yeeGridNxPlusOne || particles[ip].yeeGridNy != particles[ip].yeeGridNyPlusOne || particles[ip].yeeGridNz != particles[ip].yeeGridNzPlusOne)
            {

                // Fix for Periodic BC
                // Calculate Relay Point for Zig-Zag Scheme
                //=====IMPORTANT=====
                // Plan:
                // Let position measurements span beyond boundary, but weight values to proper locations
                // nextGridIndex values represent the grid index of the particle at the end if the boundary conditions werent imposed
                // CHECK
                //===================

                int nextGridIndexX = particles[ip].yeeGridNxPlusOne;
                int nextGridIndexY = particles[ip].yeeGridNyPlusOne;
                int nextGridIndexZ = particles[ip].yeeGridNzPlusOne;

                // condition of going from max to min for x
                if (particles[ip].yeeGridNx == (Nx - 1) && particles[ip].yeeGridNxPlusOne == 0)
                {
                    nextGridIndexX = Nx;
                }

                // condition of going from min to max for x
                else if (particles[ip].yeeGridNx == 0 && particles[ip].yeeGridNxPlusOne == (Nx - 1))
                {
                    nextGridIndexX = -1;
                }

                // condition of going from max to min for y
                if (particles[ip].yeeGridNy == (Ny - 1) && particles[ip].yeeGridNyPlusOne == 0)
                {
                    nextGridIndexY = Ny;
                }

                // condition of going from min to max for y
                else if (particles[ip].yeeGridNx == 0 && particles[ip].yeeGridNxPlusOne == (Nx - 1))
                {
                    nextGridIndexX = -1;
                }

                // condition of going from max to min for z
                if (particles[ip].yeeGridNz == (Nz - 1) && particles[ip].yeeGridNzPlusOne == 0)
                {
                    nextGridIndexZ = Nz;
                }

                // condition of going from min to max for z
                else if (particles[ip].yeeGridNz == 0 && particles[ip].yeeGridNzPlusOne == (Nz - 1))
                {
                    nextGridIndexZ = -1;
                }

                double xr = fmin(fmin(particles[ip].yeeGridNx * dx, nextGridIndexX * dx) + dx, fmax(fmax(particles[ip].yeeGridNx * dx, nextGridIndexX * dx), particles[ip].addXs / 2));
                double yr = fmin(fmin(particles[ip].yeeGridNy * dy, nextGridIndexY * dy) + dy, fmax(fmax(particles[ip].yeeGridNy * dy, nextGridIndexY * dy), particles[ip].addYs / 2));
                double zr = fmin(fmin(particles[ip].yeeGridNz * dz, nextGridIndexZ * dz) + dz, fmax(fmax(particles[ip].yeeGridNz * dz, nextGridIndexZ * dz), particles[ip].addZs / 2));

                // Calculate Charge Flux
                double Fx1 = (xr - particles[ip].x) * particles[ip].q / dt;
                double Fx2 = (particles[ip].preAdjustmentX - xr) * particles[ip].q / dt;

                double Fy1 = (yr - particles[ip].y) * particles[ip].q / dt;
                double Fy2 = (particles[ip].preAdjustmentY - yr) * particles[ip].q / dt;

                double Fz1 = (zr - particles[ip].z) * particles[ip].q / dt;
                double Fz2 = (particles[ip].preAdjustmentZ - zr) * particles[ip].q / dt;

                // Calculate Shape Factor corresponding to linear weighting function
                double Wx1 = ((particles[ip].x + xr) / 2 / dx) - particles[ip].yeeGridNx;
                double Wx2 = ((particles[ip].preAdjustmentX + xr) / 2 / dx) - nextGridIndexX;

                double Wy1 = ((particles[ip].y + yr) / 2 / dy) - particles[ip].yeeGridNy;
                double Wy2 = ((particles[ip].preAdjustmentY + yr) / 2 / dy) - nextGridIndexY;

                double Wz1 = ((particles[ip].z + zr) / 2 / dz) - particles[ip].yeeGridNz;
                double Wz2 = ((particles[ip].preAdjustmentZ + zr) / 2 / dz) - nextGridIndexZ;

                // Assign Current Densities

                // First Box indicies
                int LeftIndex1 = particles[ip].yeeGridNx;
                int RightIndex1 = particles[ip].yeeGridNx + 1;
                int BottomIndex1 = particles[ip].yeeGridNy;
                int TopIndex1 = particles[ip].yeeGridNy + 1;

                // Assume is in same cell
                int LeftIndex2 = particles[ip].yeeGridNx;
                int RightIndex2 = particles[ip].yeeGridNx + 1;
                int BottomIndex2 = particles[ip].yeeGridNy;
                int TopIndex2 = particles[ip].yeeGridNy + 1;

                // If crosses boundary in X
                if (particles[ip].yeeGridNx != particles[ip].yeeGridNxPlusOne)
                {
                    LeftIndex2 = particles[ip].yeeGridNxPlusOne;
                    RightIndex2 = particles[ip].yeeGridNxPlusOne + 1;
                }

                // If crosses boundary in Y
                if (particles[ip].yeeGridNy != particles[ip].yeeGridNyPlusOne)
                {

                    BottomIndex2 = particles[ip].yeeGridNyPlusOne;
                    TopIndex2 = particles[ip].yeeGridNyPlusOne + 1;
                }

                //  adjust which yee grid boxes get values of current density based on imposed periodic BC
                if (boundaryConditionX == 1)
                {
                    // if going from max to min

                    if (particles[ip].yeeGridNx == (Nx - 1))
                    {
                        RightIndex1 = 0;
                    }

                    if (particles[ip].yeeGridNxPlusOne == (Nx - 1))
                    {
                        RightIndex2 = 0;
                    }

                    // if equivalent and at max, set to min
                    if (particles[ip].yeeGridNx == particles[ip].yeeGridNxPlusOne && particles[ip].yeeGridNx == (Nx - 1))
                    {
                        RightIndex1 = 0;
                        RightIndex2 = 0;
                    }
                }
                // adjust for y
                if (boundaryConditionY == 1)
                {

                    // if going from max to min

                    if (particles[ip].yeeGridNy == (Ny - 1))
                    {
                        TopIndex1 = 0;
                    }

                    if (particles[ip].yeeGridNyPlusOne == (Ny - 1))
                    {
                        TopIndex2 = 0;
                    }

                    // if equivalent and at max, set to min
                    if (particles[ip].yeeGridNy == particles[ip].yeeGridNyPlusOne && particles[ip].yeeGridNy == (Ny - 1))
                    {
                        TopIndex1 = 0;
                        TopIndex2 = 0;
                    }
                }

                // Update Current Density to each side of both boxes that particle resided in

                if (numDimensions < 3)
                {
                    // Bottom of First Box
                    yeegrid[LeftIndex1][BottomIndex1][0].J.setX(yeegrid[LeftIndex1][BottomIndex1][0].J.getX() + (Fx1 * (1 - Wy1) / dx / dy));
                    // Left of First Box
                    yeegrid[LeftIndex1][BottomIndex1][0].J.setY(yeegrid[LeftIndex1][BottomIndex1][0].J.getY() + (Fy1 * (1 - Wx1) / dx / dy));
                    // Top of First Box
                    yeegrid[LeftIndex1][TopIndex1][0].J.setX(yeegrid[LeftIndex1][TopIndex1][0].J.getX() + (Fx1 * Wy1 / dx / dy));
                    // Right of First Box
                    yeegrid[RightIndex1][BottomIndex1][0].J.setY(yeegrid[RightIndex1][BottomIndex1][0].J.getY() + (Fy1 * Wx1 / dx / dy));

                    // Bottom of Second Box
                    yeegrid[LeftIndex2][BottomIndex2][0].J.setX(yeegrid[LeftIndex2][BottomIndex2][0].J.getX() + (Fx2 * (1 - Wy2) / dx / dy));
                    // Left of Second Box
                    yeegrid[LeftIndex2][BottomIndex2][0].J.setY(yeegrid[LeftIndex2][BottomIndex2][0].J.getY() + (Fy2 * (1 - Wx2) / dx / dy));
                    // Top of Second Box
                    yeegrid[LeftIndex2][TopIndex2][0].J.setX(yeegrid[LeftIndex2][TopIndex2][0].J.getX() + (Fx2 * Wy2 / dx / dy));
                    // Right of Second Box
                    yeegrid[RightIndex2][BottomIndex2][0].J.setY(yeegrid[RightIndex2][BottomIndex2][0].J.getY() + (Fy2 * Wx2 / dx / dy));
                }

                else if (numDimensions == 3)
                {

                    int zBottom1 = particles[ip].yeeGridNz;
                    int zTop1 = particles[ip].yeeGridNz + 1;
                    int zBottom2 = particles[ip].yeeGridNz;
                    int zTop2 = particles[ip].yeeGridNz + 1;

                    // If crosses boundary in Z
                    if (particles[ip].yeeGridNz != particles[ip].yeeGridNzPlusOne)
                    {

                        zBottom2 = particles[ip].yeeGridNzPlusOne;
                        zTop2 = particles[ip].yeeGridNzPlusOne + 1;
                    }

                    if (boundaryConditionZ == 1)
                    {

                        // if going from max to min

                        if (particles[ip].yeeGridNz == (Nz - 1))
                        {
                            zTop1 = 0;
                        }

                        if (particles[ip].yeeGridNzPlusOne == (Nz - 1))
                        {
                            zTop2 = 0;
                        }

                        // if equivalent and at max, set to min
                        if (particles[ip].yeeGridNz == particles[ip].yeeGridNzPlusOne && particles[ip].yeeGridNz == (Nz - 1))
                        {
                            zTop1 = 0;
                            zTop2 = 0;
                        }

                        yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setX(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getX() + ((Fx1 * (1 - Wy1) * (1 - Wz1)) / dx / dy / dz));
                        yeegrid[LeftIndex1][TopIndex1][zBottom1].J.setX(yeegrid[LeftIndex1][TopIndex1][zBottom1].J.getX() + ((Fx1 * Wy1 * (1 - Wz1)) / dx / dy / dz));
                        yeegrid[LeftIndex1][BottomIndex1][zTop1].J.setX(yeegrid[LeftIndex1][BottomIndex1][zTop1].J.getX() + ((Fx1 * (1 - Wy1) * Wz1) / dx / dy / dz));
                        yeegrid[LeftIndex1][TopIndex1][zTop1].J.setX(yeegrid[LeftIndex1][TopIndex1][zTop1].J.getX() + ((Fx1 * Wy1 * Wz1) / dx / dy / dz));

                        yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setY(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getY() + ((Fy1 * (1 - Wx1) * (1 - Wz1)) / dx / dy / dz));
                        yeegrid[RightIndex1][BottomIndex1][zBottom1].J.setY(yeegrid[RightIndex1][BottomIndex1][zBottom1].J.getY() + ((Fy1 * Wx1 * (1 - Wz1)) / dx / dy / dz));
                        yeegrid[LeftIndex1][BottomIndex1][zTop1].J.setY(yeegrid[LeftIndex1][BottomIndex1][zTop1].J.getY() + ((Fy1 * (1 - Wx1) * Wz1) / dx / dy / dz));
                        yeegrid[RightIndex1][BottomIndex1][zTop1].J.setY(yeegrid[RightIndex1][BottomIndex1][zTop1].J.getY() + ((Fy1 * Wx1 * Wz1) / dx / dy / dz));

                        yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setZ(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getZ() + ((Fz1 * (1 - Wx1) * (1 - Wy1)) / dx / dy / dz));
                        yeegrid[RightIndex1][BottomIndex1][zBottom1].J.setZ(yeegrid[RightIndex1][BottomIndex1][zBottom1].J.getZ() + ((Fz1 * Wx1 * (1 - Wy1)) / dx / dy / dz));
                        yeegrid[LeftIndex1][TopIndex1][zBottom1].J.setZ(yeegrid[LeftIndex1][TopIndex1][zBottom1].J.getZ() + ((Fz1 * (1 - Wx1) * Wy1) / dx / dy / dz));
                        yeegrid[RightIndex1][TopIndex1][zBottom1].J.setZ(yeegrid[RightIndex1][TopIndex1][zBottom1].J.getZ() + ((Fz1 * Wx1 * Wy1) / dx / dy / dz));

                        // Second Grid

                        yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.setX(yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.getX() + ((Fx2 * (1 - Wy2) * (1 - Wz2)) / dx / dy / dz));
                        yeegrid[LeftIndex2][TopIndex2][zBottom2].J.setX(yeegrid[LeftIndex2][TopIndex2][zBottom2].J.getX() + ((Fx2 * Wy2 * (1 - Wz2)) / dx / dy / dz));
                        yeegrid[LeftIndex2][BottomIndex2][zTop2].J.setX(yeegrid[LeftIndex2][BottomIndex2][zTop2].J.getX() + ((Fx2 * (1 - Wy2) * Wz2) / dx / dy / dz));
                        yeegrid[LeftIndex2][TopIndex2][zTop2].J.setX(yeegrid[LeftIndex2][TopIndex2][zTop2].J.getX() + ((Fx2 * Wy2 * Wz2) / dx / dy / dz));

                        yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.setY(yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.getY() + ((Fy2 * (1 - Wx2) * (1 - Wz2)) / dx / dy / dz));
                        yeegrid[RightIndex2][BottomIndex2][zBottom2].J.setY(yeegrid[RightIndex2][BottomIndex2][zBottom2].J.getY() + ((Fy2 * Wx2 * (1 - Wz2)) / dx / dy / dz));
                        yeegrid[LeftIndex2][BottomIndex2][zTop2].J.setY(yeegrid[LeftIndex2][BottomIndex2][zTop2].J.getY() + ((Fy2 * (1 - Wx2) * Wz2) / dx / dy / dz));
                        yeegrid[RightIndex2][BottomIndex2][zTop2].J.setY(yeegrid[RightIndex2][BottomIndex2][zTop2].J.getY() + ((Fy2 * Wx2 * Wz2) / dx / dy / dz));

                        yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.setZ(yeegrid[LeftIndex2][BottomIndex2][zBottom2].J.getZ() + ((Fz2 * (1 - Wx2) * (1 - Wy2)) / dx / dy / dz));
                        yeegrid[RightIndex2][BottomIndex2][zBottom2].J.setZ(yeegrid[RightIndex2][BottomIndex2][zBottom2].J.getZ() + ((Fz2 * Wx2 * (1 - Wy2)) / dx / dy / dz));
                        yeegrid[LeftIndex2][TopIndex2][zBottom2].J.setZ(yeegrid[LeftIndex2][TopIndex2][zBottom2].J.getZ() + ((Fz2 * (1 - Wx2) * Wy2) / dx / dy / dz));
                        yeegrid[RightIndex2][TopIndex2][zBottom2].J.setZ(yeegrid[RightIndex2][TopIndex2][zBottom2].J.getZ() + ((Fz2 * Wx2 * Wy2) / dx / dy / dz));
                    }
                }
            }

            // Did not cross cell boundary
            else
            {

                // TODO: Make more understandable
                // Determine Yee Grid Indicies to weight various Js to
                // Bottom Jx, Left Jy will be weighted to current box, find J(x+1) to get right Jy, find J(y+1) to get top Jx

                int TopIndex1 = particles[ip].yeeGridNyPlusOne + 1;   // Y-Index for getting top Jx
                int RightIndex1 = particles[ip].yeeGridNxPlusOne + 1; // X-Index for getting right Jy
                int BottomIndex1 = particles[ip].yeeGridNyPlusOne;
                int LeftIndex1 = particles[ip].yeeGridNxPlusOne;

                int zBottom1 = particles[ip].yeeGridNzPlusOne;
                int zTop1 = particles[ip].yeeGridNzPlusOne + 1;

                // Calculate Charge Flux
                double Fx1 = deltaX * particles[ip].q / dt;
                double Fy1 = deltaY * particles[ip].q / dt;
                double Fz1 = deltaZ * particles[ip].q / dt;

                // Calculate Shape Factor corresponding to linear weighting function

                double Wx1 = (particles[ip].addXs / 2 / dx) - LeftIndex1;
                double Wy1 = (particles[ip].addYs / 2 / dy) - BottomIndex1;
                double Wz1 = (particles[ip].addZs / 2 / dz) - zBottom1;

                if (boundaryConditionX == 1)
                {
                    // Set Index of Grid's Top Jx to be at index 0, if at max x
                    if (particles[ip].yeeGridNxPlusOne == (Nx - 1))
                    {
                        RightIndex1 = 0;
                    }
                }
                if (boundaryConditionY == 1)
                {
                    // Set Index of Grid's right Jy to be at index 0, if at max y
                    if (particles[ip].yeeGridNyPlusOne == (Ny - 1))
                    {
                        TopIndex1 = 0;
                    }
                }

                if (boundaryConditionZ == 1)
                {
                    // Set Index of Grid's right Jy to be at index 0, if at max y
                    if (particles[ip].yeeGridNzPlusOne == (Nz - 1))
                    {
                        zTop1 = 0;
                    }
                }

                // stayed in box
                // Need to weight to bottom Jx, top Jx. Left Jy, right Jy
                // Deal w/periodic BC

                //!!!!!!!CHECK!!!!!!!
                // Probably Has errors

                if (numDimensions < 3)
                {
                    yeegrid[LeftIndex1][BottomIndex1][0].J.setX(yeegrid[LeftIndex1][BottomIndex1][0].J.getX() + (Fx1 * (1 - Wy1) / dx / dy));
                    yeegrid[LeftIndex1][BottomIndex1][0].J.setY(yeegrid[LeftIndex1][BottomIndex1][0].J.getY() + (Fy1 * (1 - Wx1) / dx / dy));
                    yeegrid[LeftIndex1][TopIndex1][0].J.setX(yeegrid[LeftIndex1][TopIndex1][0].J.getX() + (Fx1 * Wy1 / dx / dy));
                    yeegrid[RightIndex1][BottomIndex1][0].J.setY(yeegrid[RightIndex1][BottomIndex1][0].J.getY() + (Fy1 * Wx1 / dx / dy));
                }

                else if (numDimensions == 3)
                {
                    yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setX(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getX() + ((Fx1 * (1 - Wy1) * (1 - Wz1)) / dx / dy / dz));
                    yeegrid[LeftIndex1][TopIndex1][zBottom1].J.setX(yeegrid[LeftIndex1][TopIndex1][zBottom1].J.getX() + ((Fx1 * Wy1 * (1 - Wz1)) / dx / dy / dz));
                    yeegrid[LeftIndex1][BottomIndex1][zTop1].J.setX(yeegrid[LeftIndex1][BottomIndex1][zTop1].J.getX() + ((Fx1 * (1 - Wy1) * Wz1) / dx / dy / dz));
                    yeegrid[LeftIndex1][TopIndex1][zTop1].J.setX(yeegrid[LeftIndex1][TopIndex1][zTop1].J.getX() + ((Fx1 * Wy1 * Wz1) / dx / dy / dz));

                    yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setY(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getY() + ((Fy1 * (1 - Wx1) * (1 - Wz1)) / dx / dy / dz));
                    yeegrid[RightIndex1][BottomIndex1][zBottom1].J.setY(yeegrid[RightIndex1][BottomIndex1][zBottom1].J.getY() + ((Fy1 * Wx1 * (1 - Wz1)) / dx / dy / dz));
                    yeegrid[LeftIndex1][BottomIndex1][zTop1].J.setY(yeegrid[LeftIndex1][BottomIndex1][zTop1].J.getY() + ((Fy1 * (1 - Wx1) * Wz1) / dx / dy / dz));
                    yeegrid[RightIndex1][BottomIndex1][zTop1].J.setY(yeegrid[RightIndex1][BottomIndex1][zTop1].J.getY() + ((Fy1 * Wx1 * Wz1) / dx / dy / dz));

                    yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.setZ(yeegrid[LeftIndex1][BottomIndex1][zBottom1].J.getZ() + ((Fz1 * (1 - Wx1) * (1 - Wy1)) / dx / dy / dz));
                    yeegrid[RightIndex1][BottomIndex1][zBottom1].J.setZ(yeegrid[RightIndex1][BottomIndex1][zBottom1].J.getZ() + ((Fz1 * Wx1 * (1 - Wy1)) / dx / dy / dz));
                    yeegrid[LeftIndex1][TopIndex1][zBottom1].J.setZ(yeegrid[LeftIndex1][TopIndex1][zBottom1].J.getZ() + ((Fz1 * (1 - Wx1) * Wy1) / dx / dy / dz));
                    yeegrid[RightIndex1][TopIndex1][zBottom1].J.setZ(yeegrid[RightIndex1][TopIndex1][zBottom1].J.getZ() + ((Fz1 * Wx1 * Wy1) / dx / dy / dz));
                }
            }

            // Need to deal with particles leaving and experiencing periodic BC

            //============
            // Compute charge density (rho) at each grid origin
            //============

            int currentIndexX = particles[ip].yeeGridNxPlusOne;
            int currentIndexY = particles[ip].yeeGridNyPlusOne;
            int currentIndexZ = 0;

            double charge = particles[ip].q;

            double localXNew = particles[ip].localXPlusOne;
            double localYNew = particles[ip].localYPlusOne;

            int incrementedIndexX = particles[ip].yeeGridNxPlusOne + 1;
            int incrementedIndexY = particles[ip].yeeGridNyPlusOne + 1;

            if (boundaryConditionX == 1 && currentIndexX == Nx - 1)
            {
                incrementedIndexX = 0;
            }
            if (boundaryConditionY == 1 && currentIndexY == Ny - 1)
            {
                incrementedIndexY = 0;
            }

            // bottom left node receiving top right area
            yeegrid[currentIndexX][currentIndexY][currentIndexZ].incrementRho(charge * ((dx - localXNew) * (dy - localYNew) / gridArea)); // dx / dy;
            // bottom right node receiving top left area
            yeegrid[incrementedIndexX][currentIndexY][currentIndexZ].incrementRho(charge * (localXNew * (dy - localYNew) / gridArea)); // dx / dy;
            // top left node receiving bottom right area
            yeegrid[currentIndexX][incrementedIndexY][currentIndexZ].incrementRho(charge * ((dx - localXNew) * localYNew / gridArea)); // dx / dy;
            // top right node receiving bottom left area
            yeegrid[incrementedIndexX][incrementedIndexY][currentIndexZ].incrementRho(charge * (localXNew * localYNew / gridArea)); // dx / dy;
        }
    }
}

void fieldSolver()
{
    int z = 0;

    // Iterate through all yee grids and solve for Ex and Ey, accounting for boundary conditions
    if (numDimensions < 3)
    {
        for (int x = 0; x < Nx; x++)
        {
            for (int y = 0; y < Ny; y++)
            {

                // For BC
                double Ex = 0;
                double Ey = 0;

                // Periodic BC Y
                if (boundaryConditionY == 1)
                {

                    if (y == 0)
                    {
                        // compute Ex in case of BC imposed at minimum
                        Ex = (((yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x][Ny - 1][z].BPlusHalf.getZ()) * pow(cph, 2) / dy) - (yeegrid[x][y][z].J.getX() / eps)) * dt + yeegrid[x][y][z].E.getX();

                    } // Could just add change to base initial E values?
                    else
                    {
                        // compute Ex
                        Ex = (((yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x][y - 1][z].BPlusHalf.getZ()) * pow(cph, 2) / dy) - (yeegrid[x][y][z].J.getX() / eps)) * dt + yeegrid[x][y][z].E.getX();
                    }
                }

                // WallBC in Y
                if (boundaryConditionY == 2)
                {
                    if (y == 0) // || y == Ny - 1) //
                    {
                        // compute Ex in case of BC imposed at minimum
                        Ex = 0;

                    } // Could just add change to base initial E values?
                    else
                    {
                        // compute Ex
                        Ex = (((yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x][y - 1][z].BPlusHalf.getZ()) * (pow(cph, 2)) / dy) - (yeegrid[x][y][z].J.getX() / eps)) * dt + yeegrid[x][y][z].E.getX();
                    }
                }

                // Periodic BC X
                if (boundaryConditionX == 1)
                {
                    if (x == 0)
                    {
                        // compute Ey in case of BC imposed at minimum
                        Ey = ((-1 * (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[Nx - 1][y][z].BPlusHalf.getZ()) * (pow(cph, 2)) / dx) - (yeegrid[x][y][z].J.getY() / eps)) * dt + yeegrid[x][y][z].E.getY();
                    }
                    else
                    {
                        // compute Ey
                        Ey = ((-1 * (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x - 1][y][z].BPlusHalf.getZ()) * (pow(cph, 2)) / dx) - (yeegrid[x][y][z].J.getY() / eps)) * dt + yeegrid[x][y][z].E.getY();
                    }
                }

                // Wall BC X
                if (boundaryConditionX == 2)
                {
                    if (x == 0)
                    {
                        // compute Ey in case of BC imposed at minimum
                        Ey = 0;
                    }
                    else
                    {
                        // compute Ey
                        Ey = ((-1 * (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x - 1][y][z].BPlusHalf.getZ()) * (pow(cph, 2)) / dx) - (yeegrid[x][y][z].J.getY() / eps)) * dt + yeegrid[x][y][z].E.getY();
                        // if (x == 100 && y == 100)

                        //    std::cout << yeegrid[x - 1][y][z].BPlusHalf.getZ() << std::endl;
                    }
                }

                yeegrid[x][y][z].EPlusOne.setX(Ex);
                yeegrid[x][y][z].EPlusOne.setY(Ey);
                // if (x == xmax && y == ymax)
                // std::cout << yeegrid[x][y][z].EPlusOne.getX() << std::endl;
            }
        }

        // Iterate through all yee grids and solve for Bz, accounting for boundary conditions

        for (int x = 0; x < Nx; x++)
        {
            for (int y = 0; y < Ny; y++)
            {

                double Bz; // = yeegrid[x][y][z].BPlusThreeHalf.getZ(); // placeholder value

                // calculate next index to pull Ex or Ey from (lower Ex or Ey to calculate Bz will always be in the same grid index as Bz)
                int xNext = x + 1;
                int yNext = y + 1;

                // Placeholder value, in order to know if value will change at all, if doesnt change, impose normal conditions
                double ExBottom = 99.99;
                double ExTop = 99.99;
                double EyLeft = 99.99;
                double EyRight = 99.99;

                // Loop Back around for Periodic BC
                // INHERENT PERIODIC BC< MUST CHANGE

                if (boundaryConditionX == 1)
                {
                    if (x == Nx - 1)
                    {
                        EyRight = yeegrid[0][y][z].EPlusOne.getY();
                    }
                    if (EyRight == 99.99)
                    {
                        EyRight = yeegrid[xNext][y][z].EPlusOne.getY();
                    }
                    if (EyLeft == 99.99)
                    {
                        EyLeft = yeegrid[x][y][z].EPlusOne.getY();
                    }
                }
                else if (boundaryConditionX == 2)
                {
                    if (x == 0)
                    {
                        EyLeft = 0;
                    }
                    if (x == Nx - 1)
                    {
                        EyRight = 0;
                    }
                    if (EyLeft == 99.99)
                    {
                        EyLeft = yeegrid[x][y][z].EPlusOne.getY();
                    }
                    if (EyRight == 99.99)
                    {
                        EyRight = yeegrid[xNext][y][z].EPlusOne.getY();
                    }
                    /*
                    if (x == 1 && (y != 0 && y != Ny - 1))
                    {
                        Bz = ((((yeegrid[x][yNext][z].EPlusOne.getX() - yeegrid[x][y][z].EPlusOne.getX()) / dy) - ((yeegrid[xNext][y][z].EPlusOne.getY() - yeegrid[x][y][z].EPlusOne.getY()) / dx)) * dt) + yeegrid[x][y][z].BPlusHalf.getZ();
                        yeegrid[0][y][z].BPlusThreeHalf.setZ(Bz);
                    }

                    else if (x == Nx - 1 && (y != 0 && y != Ny - 1))
                    {
                        Bz = yeegrid[Nx - 2][y][z].BPlusThreeHalf.getZ();
                    }
                    */
                }

                if (boundaryConditionY == 1)
                {
                    if (y == Ny - 1)
                    {
                        ExTop = yeegrid[x][0][z].EPlusOne.getX();
                    }
                    if (ExTop == 99.99)
                    {
                        ExTop = yeegrid[x][yNext][z].EPlusOne.getX();
                    }
                    if (ExBottom == 99.99)
                    {
                        ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                    }
                }
                else if (boundaryConditionY == 2)
                {

                    if (y == 0)
                    {
                        ExBottom = 0;
                    }
                    if (y == Ny - 1)
                    {
                        ExTop = 0;
                    }
                    if (ExBottom == 99.99)
                    {
                        ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                    }
                    if (ExTop == 99.99)
                    {

                        ExTop = yeegrid[x][yNext][z].EPlusOne.getX();
                    }
                    /*
                    if (y == 1 && (x != 0 && x != Nx - 1))
                    {
                        Bz = ((((yeegrid[x][yNext][z].EPlusOne.getX() - yeegrid[x][y][z].EPlusOne.getX()) / dy) - ((yeegrid[xNext][y][z].EPlusOne.getY() - yeegrid[x][y][z].EPlusOne.getY()) / dx)) * dt) + yeegrid[x][y][z].BPlusHalf.getZ();
                        yeegrid[x][0][z].BPlusThreeHalf.setZ(Bz);
                    }
                    else if (y == Ny - 1 && (x != 0 && x != Nx - 1))
                    {
                        Bz = yeegrid[x][Ny - 2][z].BPlusThreeHalf.getZ();
                    }
                    */
                }

                // Find New field by calculating dBz,  ***commented out section ***subtracting previous bg field quantities and adding new ones
                Bz = ((((ExTop - ExBottom) / dy) - ((EyRight - EyLeft) / dx)) * dt) + yeegrid[x][y][z].BPlusHalf.getZ();
                yeegrid[x][y][z].BPlusThreeHalf.setZ(Bz);
            }
        }
    }

    else if (numDimensions == 3)
    {
        // For BC
        double Ex = 0;
        double Ey = 0;
        double Ez = 0;

        for (int x = 0; x < Nx; x++)
        {
            for (int y = 0; y < Ny; y++)
            {
                for (int z = 0; z < Nz; z++)
                {
                    int xIndex = x - 1;
                    int yIndex = y - 1;
                    int zIndex = z - 1;
                    double element1;
                    double element2;

                    if (boundaryConditionY == 1 && y == 0)
                    {
                        yIndex = Ny - 1;
                        element1 = (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x][yIndex][z].BPlusHalf.getZ()) / dy;
                    }
                    else if (boundaryConditionY == 2 && y == 0)
                    {
                        element1 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element1 = (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[x][yIndex][z].BPlusHalf.getZ()) / dy;
                    }

                    if (boundaryConditionZ == 1 && z == 0)
                    {
                        zIndex = Nz - 1;
                        element2 = (yeegrid[x][y][z].BPlusHalf.getY() - yeegrid[x][y][zIndex].BPlusHalf.getY()) / dz;
                    }
                    else if (boundaryConditionZ == 2 && z == 0)
                    {
                        element2 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element2 = (yeegrid[x][y][z].BPlusHalf.getY() - yeegrid[x][y][zIndex].BPlusHalf.getY()) / dz;
                    }

                    // compute Ex
                    Ex = (((element1 - element2) * pow(cph, 2)) - (yeegrid[x][y][z].J.getX() / eps)) * dt + yeegrid[x][y][z].E.getX();

                    if (boundaryConditionY == 2)
                    {
                        if (y == 0)
                        {
                            // compute Ex in case of BC imposed at minimum
                            Ex = 0;
                        }
                    }
                    if (boundaryConditionZ == 2)
                    {
                        if (z == 0)
                        {
                            // compute Ex in case of BC imposed at minimum
                            Ex = 0;
                        }
                    }

                    //===========
                    //    Ey
                    //===========

                    // Periodic BC X
                    xIndex = x - 1;
                    zIndex = z - 1;
                    element1 = 0;
                    element2 = 0;

                    if (boundaryConditionX == 1 && x == 0)
                    {
                        xIndex = Nx - 1;
                        element1 = (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[xIndex][y][z].BPlusHalf.getZ()) / dx;
                    }
                    else if (boundaryConditionX == 2 && x == 0)
                    {
                        element1 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element1 = (yeegrid[x][y][z].BPlusHalf.getZ() - yeegrid[xIndex][y][z].BPlusHalf.getZ()) / dx;
                    }

                    if (boundaryConditionZ == 1 && z == 0)
                    {
                        zIndex = Nz - 1;
                        element2 = (yeegrid[x][y][z].BPlusHalf.getX() - yeegrid[x][y][zIndex].BPlusHalf.getX()) / dz;
                    }
                    else if (boundaryConditionZ == 2 && z == 0)
                    {
                        element2 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element2 = (yeegrid[x][y][z].BPlusHalf.getX() - yeegrid[x][y][zIndex].BPlusHalf.getX()) / dz;
                    }

                    // compute Ey
                    Ey = (((element2 - element1) * pow(cph, 2)) - (yeegrid[x][y][z].J.getY() / eps)) * dt + yeegrid[x][y][z].E.getY();

                    if (boundaryConditionX == 2)
                    {
                        if (x == 0)
                        {
                            // compute Ey in case of BC imposed at minimum
                            Ey = 0;
                        }
                    }
                    if (boundaryConditionZ == 2)
                    {
                        if (z == 0)
                        {
                            // compute Ey in case of BC imposed at minimum
                            Ey = 0;
                        }
                    }

                    //==========
                    //    Ez
                    //==========
                    xIndex = x - 1;
                    yIndex = y - 1;
                    element1 = 0;
                    element2 = 0;

                    if (boundaryConditionX == 1 && x == 0)
                    {
                        xIndex = Nx - 1;
                        element1 = (yeegrid[x][y][z].BPlusHalf.getY() - yeegrid[xIndex][y][z].BPlusHalf.getY()) / dx;
                    }
                    else if (boundaryConditionX == 2 && x == 0)
                    {
                        element1 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element1 = (yeegrid[x][y][z].BPlusHalf.getY() - yeegrid[xIndex][y][z].BPlusHalf.getY()) / dx;
                    }

                    if (boundaryConditionY == 1 && y == 0)
                    {
                        yIndex = Ny - 1;
                        element2 = (yeegrid[x][y][z].BPlusHalf.getX() - yeegrid[x][yIndex][z].BPlusHalf.getX()) / dy;
                    }
                    else if (boundaryConditionZ == 2 && z == 0)
                    {
                        element2 = 0; // placeholder for calculations
                    }
                    else
                    {
                        element2 = (yeegrid[x][y][z].BPlusHalf.getX() - yeegrid[x][yIndex][z].BPlusHalf.getX()) / dy;
                    }

                    // compute Ez
                    Ez = (((element1 - element2) * pow(cph, 2)) - (yeegrid[x][y][z].J.getZ() / eps)) * dt + yeegrid[x][y][z].E.getZ();

                    if (boundaryConditionX == 2)
                    {
                        if (x == 0)
                        {
                            // compute Ez in case of BC imposed at minimum
                            Ez = 0;
                        }
                    }
                    if (boundaryConditionY == 2)
                    {
                        if (y == 0)
                        {
                            // compute Ez in case of BC imposed at minimum
                            Ez = 0;
                        }
                    }

                    yeegrid[x][y][z].EPlusOne.setX(Ex);
                    yeegrid[x][y][z].EPlusOne.setY(Ey);
                    yeegrid[x][y][z].EPlusOne.setY(Ez);
                    // if (x == xmax && y == ymax)
                    // std::cout << yeegrid[x][y][z].EPlusOne.getX() << std::endl;
                }
            }
        }

        // Iterate through all yee grids and solve for Bx, By, Bz, accounting for boundary conditions
        for (int x = 0; x < Nx; x++)
        {
            for (int y = 0; y < Ny; y++)
            {
                for (int z = 0; z < Nz; z++)
                {
                    //==========
                    // Bx
                    //==========
                    double Bx; // = yeegrid[x][y][z].BPlusThreeHalf.getZ(); // placeholder value

                    // calculate next index to pull Ex or Ey from (lower Ex or Ey to calculate Bz will always be in the same grid index as Bz)
                    int xNext = x + 1;
                    int yNext = y + 1;
                    int zNext = z + 1;

                    // Placeholder value, in order to know if value will change at all, if doesnt change, impose normal conditions
                    // IMPORTANT FOR NAMING CONVENTION: Observing left side of cell from center

                    double EyBottom = 99.99;
                    double EyTop = 99.99;
                    double EzLeft = 99.99;
                    double EzRight = 99.99;

                    // Loop Back around for Periodic BC
                    // INHERENT PERIODIC BC< MUST CHANGE

                    if (boundaryConditionY == 1)
                    {
                        if (y == Ny - 1)
                        {
                            EzRight = yeegrid[x][0][z].EPlusOne.getZ();
                        }
                        if (EzRight == 99.99)
                        {
                            EzRight = yeegrid[x][yNext][z].EPlusOne.getZ();
                        }
                        if (EzLeft == 99.99)
                        {
                            EzLeft = yeegrid[x][y][z].EPlusOne.getZ();
                        }
                    }
                    else if (boundaryConditionY == 2)
                    {
                        if (y == 0)
                        {
                            EzLeft = 0;
                        }
                        if (y == Ny - 1)
                        {
                            EzRight = 0;
                        }
                        if (EzLeft == 99.99)
                        {
                            EzLeft = yeegrid[x][y][z].EPlusOne.getZ();
                        }
                        if (EzRight == 99.99)
                        {
                            EzRight = yeegrid[x][yNext][z].EPlusOne.getZ();
                        }
                    }

                    if (boundaryConditionZ == 1)
                    {
                        if (z == Nz - 1)
                        {
                            EyTop = yeegrid[x][y][0].EPlusOne.getY();
                        }
                        if (EyTop == 99.99)
                        {
                            EyTop = yeegrid[x][y][zNext].EPlusOne.getY();
                        }
                        if (EyBottom == 99.99)
                        {
                            EyBottom = yeegrid[x][y][z].EPlusOne.getY();
                        }
                    }
                    else if (boundaryConditionZ == 2)
                    {

                        if (z == 0)
                        {
                            EyBottom = 0;
                        }
                        if (z == Nz - 1)
                        {
                            EyTop = 0;
                        }
                        if (EyBottom == 99.99)
                        {
                            EyBottom = yeegrid[x][y][z].EPlusOne.getY();
                        }
                        if (EyTop == 99.99)
                        {

                            EyTop = yeegrid[x][y][zNext].EPlusOne.getY();
                        }
                    }

                    Bx = ((((EyTop - EyBottom) / dz) - ((EzLeft - EzRight) / dy)) * dt) + yeegrid[x][y][z].BPlusHalf.getX();

                    //==========
                    // By
                    //==========
                    double By; // = yeegrid[x][y][z].BPlusThreeHalf.getZ(); // placeholder value

                    // calculate next index to pull Ex or Ey from (lower Ex or Ey to calculate Bz will always be in the same grid index as Bz)
                    xNext = x + 1;
                    zNext = z + 1;

                    // Placeholder value, in order to know if value will change at all, if doesnt change, impose normal conditions
                    // IMPORTANT FOR NAMING CONVENTION: Observing front side of cell from outside

                    double ExBottom = 99.99;
                    double ExTop = 99.99;
                    EzLeft = 99.99;
                    EzRight = 99.99;

                    // Loop Back around for Periodic BC
                    // INHERENT PERIODIC BC< MUST CHANGE

                    if (boundaryConditionX == 1)
                    {
                        if (x == Nx - 1)
                        {
                            EzRight = yeegrid[0][y][z].EPlusOne.getZ();
                        }
                        if (EzRight == 99.99)
                        {
                            EzRight = yeegrid[xNext][y][z].EPlusOne.getZ();
                        }
                        if (EzLeft == 99.99)
                        {
                            EzLeft = yeegrid[x][y][z].EPlusOne.getZ();
                        }
                    }
                    else if (boundaryConditionX == 2)
                    {
                        if (x == 0)
                        {
                            EzLeft = 0;
                        }
                        if (x == Nx - 1)
                        {
                            EzRight = 0;
                        }
                        if (EzLeft == 99.99)
                        {
                            EzLeft = yeegrid[x][y][z].EPlusOne.getZ();
                        }
                        if (EzRight == 99.99)
                        {
                            EzRight = yeegrid[xNext][y][z].EPlusOne.getZ();
                        }
                    }

                    if (boundaryConditionZ == 1)
                    {
                        if (z == Nz - 1)
                        {
                            ExTop = yeegrid[x][y][0].EPlusOne.getX();
                        }
                        if (ExTop == 99.99)
                        {
                            ExTop = yeegrid[x][y][zNext].EPlusOne.getX();
                        }
                        if (ExBottom == 99.99)
                        {
                            ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                        }
                    }
                    else if (boundaryConditionZ == 2)
                    {

                        if (z == 0)
                        {
                            ExBottom = 0;
                        }
                        if (z == Nz - 1)
                        {
                            ExTop = 0;
                        }
                        if (ExBottom == 99.99)
                        {
                            ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                        }
                        if (ExTop == 99.99)
                        {
                            ExTop = yeegrid[x][y][zNext].EPlusOne.getX();
                        }
                    }

                    By = ((((EzLeft - EzRight) / dx) - ((ExTop - ExBottom) / dz)) * dt) + yeegrid[x][y][z].BPlusHalf.getY();

                    //==========
                    // Bz
                    //==========

                    double Bz; // = yeegrid[x][y][z].BPlusThreeHalf.getZ(); // placeholder value

                    // calculate next index to pull Ex or Ey from (lower Ex or Ey to calculate Bz will always be in the same grid index as Bz)
                    xNext = x + 1;
                    yNext = y + 1;

                    // Placeholder value, in order to know if value will change at all, if doesnt change, impose normal conditions
                    ExBottom = 99.99;
                    ExTop = 99.99;
                    double EyLeft = 99.99;
                    double EyRight = 99.99;

                    // Loop Back around for Periodic BC
                    // INHERENT PERIODIC BC< MUST CHANGE

                    if (boundaryConditionX == 1)
                    {
                        if (x == Nx - 1)
                        {
                            EyRight = yeegrid[0][y][z].EPlusOne.getY();
                        }
                        if (EyRight == 99.99)
                        {
                            EyRight = yeegrid[xNext][y][z].EPlusOne.getY();
                        }
                        if (EyLeft == 99.99)
                        {
                            EyLeft = yeegrid[x][y][z].EPlusOne.getY();
                        }
                    }
                    else if (boundaryConditionX == 2)
                    {
                        if (x == 0)
                        {
                            EyLeft = 0;
                        }
                        if (x == Nx - 1)
                        {
                            EyRight = 0;
                        }
                        if (EyLeft == 99.99)
                        {
                            EyLeft = yeegrid[x][y][z].EPlusOne.getY();
                        }
                        if (EyRight == 99.99)
                        {
                            EyRight = yeegrid[xNext][y][z].EPlusOne.getY();
                        }
                    }

                    if (boundaryConditionY == 1)
                    {
                        if (y == Ny - 1)
                        {
                            ExTop = yeegrid[x][0][z].EPlusOne.getX();
                        }
                        if (ExTop == 99.99)
                        {
                            ExTop = yeegrid[x][yNext][z].EPlusOne.getX();
                        }
                        if (ExBottom == 99.99)
                        {
                            ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                        }
                    }
                    else if (boundaryConditionY == 2)
                    {

                        if (y == 0)
                        {
                            ExBottom = 0;
                        }
                        if (y == Ny - 1)
                        {
                            ExTop = 0;
                        }
                        if (ExBottom == 99.99)
                        {
                            ExBottom = yeegrid[x][y][z].EPlusOne.getX();
                        }
                        if (ExTop == 99.99)
                        {

                            ExTop = yeegrid[x][yNext][z].EPlusOne.getX();
                        }
                    }

                    Bz = ((((ExTop - ExBottom) / dy) - ((EyRight - EyLeft) / dx)) * dt) + yeegrid[x][y][z].BPlusHalf.getZ();

                    yeegrid[x][y][z].BPlusThreeHalf.setX(Bx);
                    yeegrid[x][y][z].BPlusThreeHalf.setY(By);
                    yeegrid[x][y][z].BPlusThreeHalf.setZ(Bz);
                }
            }
        }
    }
}

void pushUpdateObjectVariables()
{
    // Push backwards values of all Yee Grids, remove old data
    int z = 0;

    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {

            yeegrid[x][y][z].E.setX(yeegrid[x][y][z].EPlusOne.getX());
            yeegrid[x][y][z].E.setY(yeegrid[x][y][z].EPlusOne.getY());
            yeegrid[x][y][z].E.setZ(yeegrid[x][y][z].EPlusOne.getZ());

            yeegrid[x][y][z].EPlusOne.setX(0);
            yeegrid[x][y][z].EPlusOne.setY(0);
            yeegrid[x][y][z].EPlusOne.setZ(0);

            yeegrid[x][y][z].BMinusHalf.setX(yeegrid[x][y][z].BPlusHalf.getX());
            yeegrid[x][y][z].BPlusHalf.setX(yeegrid[x][y][z].BPlusThreeHalf.getX());
            yeegrid[x][y][z].BPlusThreeHalf.setX(0);

            yeegrid[x][y][z].BMinusHalf.setY(yeegrid[x][y][z].BPlusHalf.getY());
            yeegrid[x][y][z].BPlusHalf.setY(yeegrid[x][y][z].BPlusThreeHalf.getY());
            yeegrid[x][y][z].BPlusThreeHalf.setY(0);

            yeegrid[x][y][z].BMinusHalf.setZ(yeegrid[x][y][z].BPlusHalf.getZ());
            yeegrid[x][y][z].BPlusHalf.setZ(yeegrid[x][y][z].BPlusThreeHalf.getZ());
            yeegrid[x][y][z].BPlusThreeHalf.setZ(0);

            yeegrid[x][y][z].rho = 0;

            yeegrid[x][y][z].J.setX(0);
            yeegrid[x][y][z].J.setY(0);
            yeegrid[x][y][z].J.setZ(0);
        }
    }

    for (int ip = 0; ip < numMacroParticles; ip++)
    {
        // call particle fn reset
        particles[ip].updateAfterLoop();
    }
}

// utility functions
double min(double a, double b)
{
    if (a > b)
    {
        return b;
    }
    return a;
}
double max(double a, double b)
{
    if (a < b)
    {
        return b;
    }
    return a;
}