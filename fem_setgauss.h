#pragma once
#ifndef _FEM_GAUSS_SET
#define _FEM_GAUSS_SET
#include "fem_basic.h"
#include "fem_global.h"

int fem_setgauss(int rank)
{
    if (rank == 2)
    {
        _GAUSS_POINTS = 2;

        _GAUSS_POS[0] = -0.577350269189626;
        _GAUSS_POS[1] = 0.577350269189626;

        _GAUSS_Height[0] = 1.0;
        _GAUSS_Height[1] = 1.0;
    }

    if (rank == 3)
    {
        _GAUSS_POINTS = 3;

        _GAUSS_POS[0] = 0.0000000000000000;
        _GAUSS_POS[1] = -0.7745966692414834;
        _GAUSS_POS[2] = 0.7745966692414834;

        _GAUSS_Height[0] = 0.8888888888888888;
        _GAUSS_Height[1] = 0.5555555555555556;
        _GAUSS_Height[2] = 0.5555555555555556;
    }

    if (rank == 4)
    {
        _GAUSS_POINTS = 4;

        _GAUSS_POS[0] = -0.3399810435848563;
        _GAUSS_POS[1] = 0.3399810435848563;
        _GAUSS_POS[2] = -0.8611363115940526;
        _GAUSS_POS[3] = 0.8611363115940526;

        _GAUSS_Height[0] = 0.6521451548625461;
        _GAUSS_Height[1] = 0.6521451548625461;
        _GAUSS_Height[2] = 0.3478548451374538;
        _GAUSS_Height[3] = 0.3478548451374538;
    }
    printf("Gauss int points = %d\n", _GAUSS_POINTS);
    return 0;
}

#endif