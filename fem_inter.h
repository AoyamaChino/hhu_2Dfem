#ifndef _FEM_INTER
#define _FEM_INTER

#include "fem_basic.h"
#include "fem_global.h"

int fem_inter()
{
    // memory alloc
    dispx   = (double *)calloc(Node_num, sizeof(double));
    dispy   = (double *)calloc(Node_num, sizeof(double));
    stressx = (double *)calloc(Node_num, sizeof(double));
    stressy = (double *)calloc(Node_num, sizeof(double));
    strainx = (double *)calloc(Node_num, sizeof(double));
    strainy = (double *)calloc(Node_num, sizeof(double));

    // displacement constrait
    ifxfixd = (int *)calloc(Node_num, sizeof(int));
    ifyfixd = (int *)calloc(Node_num, sizeof(int));

    for (int i = 0; i < _DISPCONS_COUNT; i++)
    {
        int nodeid      = node_constraint_info[i2(i, 0, NOD_CON_NUM)];
        ifxfixd[nodeid] = node_constraint_info[i2(i, 1, NOD_CON_NUM)];
        ifyfixd[nodeid] = node_constraint_info[i2(i, 2, NOD_CON_NUM)];
    }

    // dof counter
    _MAXDOF = 0;

    // test
    /*
    for (int i = 0; i < Node_num; i++)
    {
        if (ifxfixd[i] != 0)
        {
            printf("%d,%d,%d\n", i, ifxfixd[i], ifyfixd[i]);
        }
    }
    */
    return 0;
}

#endif