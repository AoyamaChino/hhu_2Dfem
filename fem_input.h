#pragma once
#ifndef _FEM_INPUT
#define _FEM_INPUT
#include "fem_basic.h"
#include "fem_global.h"

FILE *notefile;
FILE *inputfile;
char _buff_inp[255];
fpos_t inppos = 0;

int filegotouseful(FILE *file)
{
    while (_buff_inp[0] != '*')
    {
        fgets(_buff_inp, 255, file);
    }
    // printf("%s", _buff_inp);
    return 0;
}

void skipline(int num, FILE *file)
{
    for (int i = 0; i < num; i++)
    {
        fgets(_buff_inp, 255, file);
    }
}

int fem_input()
{
    filegotouseful(inputfile);
    // printf("%s", _buff_inp);
    sscanf(_buff_inp, "*Node2d,%d", &Node_num);
    fprintf(notefile, "node_num=%d\n", Node_num);

    coordx = (double *)calloc(Node_num, sizeof(double));
    coordy = (double *)calloc(Node_num, sizeof(double));

    skipline(1, inputfile);

    for (int i = 0; i < Node_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);

        int id;
        double px, py;
        sscanf(_buff_inp, "%d %lf %lf", &id, &px, &py);

        if (i != id)
        {
            printf("READ NODE INFO ERROR\n");
            exit(0);
        }

        coordx[id] = px;
        coordy[id] = py;
        fprintf(notefile, "\t%d\t%lf\t%lf\n", id, px, py);
    }
    filegotouseful(inputfile);
    sscanf(_buff_inp, "*Element2d,%d", &Element_num);
    fprintf(notefile, "element_num=%d\n", Element_num);

    //#EleId:N1, N2, N3, N4, StepId, MatId; info:6 * int;
    element_4info = (int *)calloc(ELEMENT_DATA_NUM * Element_num, sizeof(int));

    skipline(1, inputfile);

    for (int i = 0; i < Element_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int tid, tn1, tn2, tn3, tn4, tstepid, tmatid;
        sscanf(_buff_inp, "%d %d %d %d %d %d %d", &tid, &tn1, &tn2, &tn3, &tn4,
               &tstepid, &tmatid);
        // printf("%d %d %d %d %d %d %d\n",
        //       id, n1, n2, n3, n4, stepid, matid);
        if (i != tid)
        {
            printf("ELEMENT INFO ERROR\n");
            exit(0);
        }
        element_4info[n1(i)]     = tn1;
        element_4info[n2(i)]     = tn2;
        element_4info[n3(i)]     = tn3;
        element_4info[n4(i)]     = tn4;
        element_4info[stepid(i)] = tstepid;
        element_4info[matid(i)]  = tmatid;
    }

    // longging
    for (int i = 0; i < Element_num * 6; i++)
    {
        fprintf(notefile, "\t%d", element_4info[i]);
        if ((i + 1) % 6 == 0)
            fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    sscanf(_buff_inp, "*Material,%d", &Material_num);
    fprintf(notefile, "Material_num=%d\n", Material_num);
    skipline(1, inputfile);

    filegotouseful(inputfile);
    skipline(1, inputfile);

    mat_seep = (double *)calloc(MAT_SEEP_NUM * Material_num, sizeof(double));

    for (int i = 0; i < Material_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int id;
        double d[MAT_SEEP_NUM];
        sscanf(_buff_inp, "%d,%lf,%lf,%lf,%lf,%lf", &id, &d[0], &d[1], &d[2], &d[3],
               &d[4]);

        // printf("%d,%le,%le,%le,%le,%le\n",
        //        id, d[0], d[1], d[2], d[3], d[4]);

        if (id != i)
        {
            printf("SEEP READ ERROR\n");
            exit(0);
        }

        // logging
        for (int j = 0; j < MAT_SEEP_NUM; j++)
        {
            mat_seep[i2(i, j, MAT_SEEP_NUM)] = d[j];
            fprintf(notefile, "%le ", mat_seep[i2(i, j, MAT_SEEP_NUM)]);
        }
        fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    skipline(1, inputfile);

    mat_common = (double *)calloc(MAT_COMMON_NUM * Material_num, sizeof(double));

    for (int i = 0; i < Material_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int id;
        double d[MAT_COMMON_NUM];

        sscanf(_buff_inp, "%d,%lf,%lf,%lf,%lf,%lf", &id, &d[0], &d[1], &d[2], &d[3],
               &d[4]);

        if (id != i)
        {
            printf("COMMON READ ERROR\n");
            exit(0);
        }

        // logging
        for (int j = 0; j < MAT_COMMON_NUM; j++)
        {
            mat_common[i2(i, j, MAT_COMMON_NUM)] = d[j];
            fprintf(notefile, "%le ", mat_common[i2(i, j, MAT_COMMON_NUM)]);
        }
        fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    skipline(1, inputfile);

    mat_slope_dry =
        (double *)calloc(MAT_SLOPE_NUM * Material_num, sizeof(double));

    for (int i = 0; i < Material_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        // printf("%s", _buff_inp);

        int id;
        double d[MAT_SLOPE_NUM];

        sscanf(_buff_inp, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &id, &d[0],
               &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9]);

        if (id != i)
        {
            printf("%d / %d:", id, i);
            printf("SLOPE DRY ERROR\n");
            exit(0);
        }

        // logging
        for (int j = 0; j < MAT_SLOPE_NUM; j++)
        {
            mat_slope_dry[i2(i, j, MAT_SLOPE_NUM)] = d[j];
            fprintf(notefile, "%lf ", mat_slope_dry[i2(i, j, MAT_SLOPE_NUM)]);
        }
        fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    skipline(1, inputfile);

    mat_slope_wet =
        (double *)calloc(MAT_SLOPE_NUM * Material_num, sizeof(double));

    for (int i = 0; i < Material_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int id;
        double d[MAT_SLOPE_NUM];

        // printf("%s", _buff_inp);
        sscanf(_buff_inp, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &id, &d[0],
               &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9]);

        if (id != i)
        {
            printf("SLOPE WET ERROR");
            exit(0);
        }

        // logging
        for (int j = 0; j < MAT_SLOPE_NUM; j++)
        {
            mat_slope_wet[i2(i, j, MAT_SLOPE_NUM)] = d[j];
            fprintf(notefile, "%lf ", mat_slope_wet[i2(i, j, MAT_SLOPE_NUM)]);
        }
        fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    skipline(1, inputfile);

    deform_info = (double *)calloc(DEF_INFO_NUM * Material_num, sizeof(double));

    for (int i = 0; i < Material_num; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int id;
        double d[DEF_INFO_NUM];

        sscanf(_buff_inp, "%d,%lf,%lf", &id, &d[0], &d[1]);

        if (id != i)
        {
            printf("DEFORM INFO ERROR");
            exit(0);
        }

        // logging
        for (int j = 0; j < DEF_INFO_NUM; j++)
        {
            deform_info[i2(i, j, DEF_INFO_NUM)] = d[j];
            fprintf(notefile, "%lf ", deform_info[i2(i, j, DEF_INFO_NUM)]);
        }
        fprintf(notefile, "\n");
    }

    filegotouseful(inputfile);
    skipline(1, inputfile);

    // BoundaryCondition

    filegotouseful(inputfile);
    fgetpos(inputfile, &inppos);
    skipline(2, inputfile);

    // printf("%s", _buff_inp);

    _DISPCONS_COUNT = 0;
    // get the amount of constraint nodes
    while (_buff_inp[0] == ' ')
    {
        _DISPCONS_COUNT++;
        // printf("%s", _buff_inp);
        fgets(_buff_inp, 255, inputfile);
    }
    // logging
    fprintf(notefile, "cons_count=%d\n", _DISPCONS_COUNT);
    fsetpos(inputfile, &inppos);

    fgets(_buff_inp, 255, inputfile);

    node_constraint_info = (int *)calloc(NOD_CON_NUM * _DISPCONS_COUNT, sizeof(int));

    for (int i = 0; i < _DISPCONS_COUNT; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        int d[6];
        sscanf(_buff_inp, "%d %d %d %d %d %d", &d[0], &d[1], &d[2], &d[3], &d[4],
               &d[5]);

        for (int j = 0; j < NOD_CON_NUM; j++)
        {
            node_constraint_info[i2(i, j, NOD_CON_NUM)] = d[j];

            // logging
            fprintf(notefile, "\t%d",
                    node_constraint_info[i2(i, j, NOD_CON_NUM)]);
        }

        fprintf(notefile, "\n", _buff_inp);
    }

    // printf("%s", _buff_inp);

    filegotouseful(inputfile);
    skipline(2, inputfile);

    load_count = 0;
    fgetpos(inputfile, &inppos);
    fgets(_buff_inp, 255, inputfile);
    while ((_buff_inp[0] == ' ') && (!feof(inputfile)))
    {
        load_count++;
        // printf("%s", _buff_inp);
        fgets(_buff_inp, 255, inputfile);
    }

    // logging
    fprintf(notefile, "load_count=%d\n", load_count);
    fsetpos(inputfile, &inppos);

    load_nodes = (int *)calloc(load_count, sizeof(int));
    load_info  = (double *)calloc(load_count * LOD_INFO_NUM, sizeof(double));

    for (int i = 0; i < load_count; i++)
    {
        fgets(_buff_inp, 255, inputfile);
        // printf("%s", _buff_inp);
        sscanf(_buff_inp, "%d %lf %lf",
               &load_nodes[i],
               &load_info[i2(i, 0, LOD_INFO_NUM)],
               &load_info[i2(i, 1, LOD_INFO_NUM)]);

        fprintf(notefile, "%d %lf %lf\n",
                load_nodes[i],
                load_info[i2(i, 0, LOD_INFO_NUM)],
                load_info[i2(i, 1, LOD_INFO_NUM)]);
    }

    return 0;
}

#endif