#pragma once
#ifndef _FEM_MEM
#define _FEM_MEM
#include "fem_basic.h"
#include "fem_global.h"

int memlocate();
int fem_memfree()
{
    free(coordx);
    free(coordy);

    free(stressx);
    free(stressy);

    free(strainx);
    free(strainy);
    free(element_4info);

    free(mat_seep);
    free(mat_common);
    free(mat_slope_dry);
    free(mat_slope_wet);

    free(deform_info);
    free(boundary_info);
    free(node_constraint_info);

    free(load_nodes);
    free(load_info);

    free(dispx);
    free(dispy);
    return 0;
}

#endif