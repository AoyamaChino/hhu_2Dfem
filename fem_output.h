#pragma once
#ifndef _FEM_OUTPUT
#define _FEM_OUTPUT
#include "fem_basic.h"
#include "fem_global.h"

int fem_after()
{
    for (int i = 0; i < Node_num; i++)
    {
        int dofx = dof_global[i2(i, 0, 2)];
        int dofy = dof_global[i2(i, 1, 2)];
        // printf(" %d %d\n", dofx, dofy);

        if (dofx != -1)
            dispx[i] = dof_disp[dofx];
        if (dofy != -1)
            dispy[i] = dof_disp[dofy];

        // printf(" %le %le\n", dispx[i], dispy[i]);
    }

    return 0;
}

int fem_output()
{
    // in simple .VTK file format
    // visualize using paraview

    FILE *outputfile;
    char _buff[256];

    outputfile = fopen("./output/outfile.vtk", "w");

    fprintf(outputfile, "# vtk DataFile Version 1.0\n"
                        "Unstructured Grid Example\n"
                        "ASCII\n");

    fprintf(outputfile, "\n"
                        "DATASET UNSTRUCTURED_GRID\n");

    fprintf(outputfile, "POINTS\t%dfloat\n", Node_num);

    for (int i = 0; i < Node_num; i++)
    {
        fprintf(outputfile, "\t%le\t%le\t%le\n", coordx[i] + dispx[i], coordy[i] + dispy[i], 0.0); /* code */
    }

    fprintf(outputfile, "CELLS\t%d\t%d\n", Element_num, 5 * Element_num);
    for (int i = 0; i < Element_num; i++)
    {
        fprintf(outputfile, "\t%d\t%d\t%d\t%d\t%d\n",
                4,
                element_4info[n1(i)],
                element_4info[n2(i)],
                element_4info[n3(i)],
                element_4info[n4(i)]);
    }

    fprintf(outputfile, "CELL_TYPES\t%d\n", Element_num);
    for (int i = 0; i < Element_num; i++)
    {
        fprintf(outputfile, "\t%d\n", 9);
    }

    fprintf(outputfile, "POINT_DATA\t%d\n"
                        "SCALARS scalars float\n",
            Node_num);

    fprintf(outputfile, " LOOKUP_TABLE default\n ");
    for (int i = 0; i < Node_num; i++)
    {
        fprintf(outputfile, "\t%le\n", 1.0 * i);
    }

    fprintf(outputfile, "VECTORS displacement float\n");
    for (int i = 0; i < Node_num; i++)
    {
        fprintf(outputfile, "\t%le\t%le\t%le\n", dispx[i], dispy[i], 0.0);
    }

    fprintf(outputfile, "VECTORS stress float\n");
    for (int i = 0; i < Node_num; i++)
    {
        fprintf(outputfile, "\t%le\t%le\t%le\n", stressx[i], stressy[i], 0.0);
    }

    fprintf(outputfile, "VECTORS strain float\n");
    for (int i = 0; i < Node_num; i++)
    {
        fprintf(outputfile, "\t%le\t%le\t%le\n", strainx[i], strainy[i], 0.0);
    }

    fclose(outputfile);

    // testing
    // for (int i = 0; i < _MAXDOF; i++)
    //{
    //   printf("%le\n", dof_disp[i]);
    //}

    return 0;
}

#endif