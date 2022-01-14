#pragma once
#ifndef _FEM_SOLVE_
#define _FEM_SOLVE_

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "fem_basic.h"
#include "fem_global.h"

double get_K_num(double *global_stiffness,
                 int row, int col,
                 int *diag_pos)
{
    int rrow = row;
    int rcol = col;

    if (col > row)
    {
        rrow = col;
        rcol = row;
    }

    int pos = diag_pos[rrow] - rrow + rcol - 1;

    if ((pos <= diag_pos[rrow - 1] - 1) && (rrow > 0))
        return 0.0;

    double value = global_stiffness[pos];
    return value;
}

int get_load(double *loadnodexy)
{
    for (int i = 0; i < _MAXDOF; i++)
    {
        for (int j = 0; j < load_count; j++)
        {
            // printf("load count %d\n", load_count);
            // printf("loadnodes:%d\n", load_nodes[j]);
            // printf("\t%e\t%e\n",
            //        load_info[i2(j, 0, 2)],
            //        load_info[i2(j, 1, 2)]);

            int lnode = load_nodes[j];
            int dofx  = dof_global[i2(lnode, 0, 2)];
            int dofy  = dof_global[i2(lnode, 1, 2)];

            // printf("\t%d\t%d\n", dofx, dofy);

            if (dofx != -1)
                loadnodexy[dofx] = load_info[i2(j, 0, 2)];

            if (dofy != -1)
                loadnodexy[dofy] = load_info[i2(j, 1, 2)];
        }
    }

    return 0;
}

//没啥用的东西
int solve_axb_CG()
{
    loadnodexy = (double *)calloc(_MAXDOF, sizeof(double));
    get_load(loadnodexy);

    dof_disp = (double *)calloc(_MAXDOF, sizeof(double));

    // B = loadnodexy
    // x = dof_disp
    // A = global_stiffness

    double *b = loadnodexy;
    double *x = (double *)calloc(_MAXDOF, sizeof(double));
    double *r = (double *)calloc(_MAXDOF, sizeof(double));

    /*
    for (int i = 0; i < _MAXDOF; i++)
    {
        printf("%d %f\n", i, b[i]);
    }
    */

    for (int i = 0; i < _MAXDOF; i++)
    {
        dof_disp[i] = x[i];
        // printf("%e\n", dof_disp[i]);
    }

    return 0;
}

int solve_eigen()
{
    loadnodexy = (double *)calloc(_MAXDOF, sizeof(double));
    get_load(loadnodexy);

    dof_disp = (double *)calloc(_MAXDOF, sizeof(double));
    // B = loadnodexy
    // x = dof_disp
    // A = global_stiffness

    Eigen::SparseMatrix<double, Eigen::RowMajor> A(_MAXDOF, _MAXDOF);
    Eigen::VectorXd b(_MAXDOF);
    Eigen::VectorXd x(_MAXDOF);

    int _sc = 0;
    for (int i = 0; i < _MAXDOF; i++)
    {
        for (int j = 0; j < _MAXDOF; j++)
        {
            double value =
                get_K_num(global_stiffness, i, j, diag_pos);
            if (abs(value) < 1e-10)
                continue;
            //
            A.insert(i, j) = value;
            //
            _sc++;
            // printf("\t%d\t%d\t%e\t%d\n", i, j, value, _sc);
        }
    }

    for (int i = 0; i < _MAXDOF; i++)
    {
        b(i) = loadnodexy[i];
        x(i) = 0;
        // printf("\t%d\t%e\n", i, b.coeffRef(i));
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                             Eigen::Lower>
        cg;
    cg.compute(A);
    x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;
    // update b, and solve again
    x = cg.solve(b);

    for (int i = 0; i < _MAXDOF; i++)
    {
        dof_disp[i] = x[i];
        // printf("%e\n", dof_disp[i]);
    }

    return 0;
}
#endif