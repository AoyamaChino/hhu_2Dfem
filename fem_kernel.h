#pragma once
#ifndef _FEM_KERNEL
#define _FEM_KERNEL

#include "fem_basic.h"
#include "fem_global.h"
#include "fem_setgauss.h"
#include "fem_solve.h"

int fem_get_pNip(double *nirs, double r, double s);
int fem_diagonal_pos();
int fem_global_stiffness();
int fem_element_stiffness(int element, double *element_stiff);

int fem_ass_sub(double *sub_stiff,
                double *element_stiff,
                int irow,
                int jcol);

int fem_cal_sub(double *sub_stiff,
                double *sub_dint);

int fem_get_elastic_d(double *ela_d);

int fem_sub_element_DINT(int element,
                         int irow, int jcol,
                         double *sub_dint);

int fem_findnixy(int element, double *jacobi, double *nixy, double *gauss_rs);
int fem_get_pNip(double *nirs, double r, double s);

int fem_dofserial()
{
    // serial number of Degree of Freedom
    // 2D problem
    // 整体自由度矩阵
    dof_global = (int *)calloc(2 * Node_num, sizeof(int));
    // let it size be 2 * Node_num

    for (int i = 0; i < 2 * Node_num; i++)
    {
        dof_global[i] = -1;
    }
    // 初始化自由度矩阵，-1代表是约束的

    _MAXDOF = 0;
    for (int i = 0; i < Node_num; i++)
    {
        if (ifxfixd[i] == 0)
        {
            dof_global[i2(i, 0, 2)] = _MAXDOF;
            _MAXDOF++;
        }
        if (ifyfixd[i] == 0)
        {
            dof_global[i2(i, 1, 2)] = _MAXDOF;
            _MAXDOF++;
        }
        if ((ifxfixd[i] + ifyfixd[i] > 0))
        {
            fprintf(notefile, "node %d constrait\n", i);
        }
    }
    // 总自由度数目
    printf("Total DOF %d\n", _MAXDOF);

    // logging
    for (int i = 0; i < Node_num; i++)
    {
        fprintf(notefile, "Node:%d\t|\t%d\t%d\n",
                i,
                dof_global[i2(i, 0, 2)],
                dof_global[i2(i, 1, 2)]);
    }

    // 找到单元内部的自由度序号
    // xy=*2, 对于每个单元是八个数据
    dof_element = (int *)calloc(2 * Element_num * ELEMENT_NODE_NUM, sizeof(int));

    for (int i = 0; i < Element_num; i++)
    {
        int tn[4];
        // 第i单元
        // tn是整体节点序号
        tn[0] = element_4info[n1(i)];
        tn[1] = element_4info[n2(i)];
        tn[2] = element_4info[n3(i)];
        tn[3] = element_4info[n4(i)];

        for (int ie = 0; ie < 4; ie++)
        {
            // 单元内第ie点
            // printf("\t%d", tn[ie]);
            // 获取整体自由度序号
            int dofx = dof_global[i2(tn[ie], 0, 2)];
            int dofy = dof_global[i2(tn[ie], 1, 2)];
            // printf("\t%d\t%d\n", dofx, dofy);
            //  8 * i + 2 * ie
            dof_element[i2(i, 2 * ie, 2 * ELEMENT_NODE_NUM)]     = dofx;
            dof_element[i2(i, 2 * ie + 1, 2 * ELEMENT_NODE_NUM)] = dofy;
        }
        // printf("\n");
    }

    // logging

    for (int i = 0; i < Element_num; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            fprintf(notefile, "\t%d", dof_element[8 * i + j]);
        }
        fprintf(notefile, "\t|\tElement:%d\n", i);
    }

    return 0;
}

int fem_diagonal_pos()
{
    // 找到矩阵中的对角线元素的位置
    // Principle element matrix | 2D data in memory
    //
    diag_pos = (int *)calloc(_MAXDOF, sizeof(int));
    //

    for (int ielement = 0; ielement < Element_num; ielement++)
    {
        int mindof = _MAXDOF;

        for (int idof = 0; idof < 2 * 4; idof++)
        {
            int tdof = dof_element[i2(ielement, idof, 2 * 4)];
            if ((tdof < mindof) && (tdof > -1)) //-1 means fixd dof
                mindof = tdof;
            // logging
            fprintf(notefile, "\t%d", tdof);
        }
        fprintf(notefile, "\t|\t%d\n", mindof); // find the min degree of freedom

        for (int idof = 0; idof < 2 * 4; idof++)
        {
            int tdof      = dof_element[i2(ielement, idof, 2 * 4)];
            int tdof_diff = -1;

            if (tdof > -1)
                tdof_diff = tdof - mindof + 1; //实际就是半带宽

            if (tdof_diff > diag_pos[tdof])
                diag_pos[tdof] = tdof_diff;

            fprintf(notefile, "\t%d", tdof_diff);
        }
        fprintf(notefile, "\n\n");
    }

    // find diag_pos
    for (int i = 0; i < _MAXDOF; i++)
    {
        if (i > 0)
            diag_pos[i] += diag_pos[i - 1];
        fprintf(notefile, "\t%d\t%d\n", i, diag_pos[i]);
        _MAX_DIAG_POS = diag_pos[i];
    }
    printf("Size of stiffness %d\n", _MAX_DIAG_POS);

    return 0;
}

int fem_global_stiffness()
{
    // 下三角存储的半带宽
    // 全局刚度矩阵
    // 非零数据数量 = _MAX_DIAG_POS
    // 矩阵形式 : _MAXDOF * _MAXODF
    global_stiffness = (double *)calloc(_MAX_DIAG_POS, sizeof(double));
    //
    //单元刚度矩阵计算，大小为8*8，可以分为4个2*2区域
    for (int i_element = 0; i_element < Element_num; i_element++)
    {
        // printf(" %d |", i_element);

        double *element_stiff;
        element_stiff = (double *)calloc(8 * 8, sizeof(double));
        fem_element_stiffness(i_element, element_stiff); //获得此单元的单元刚度矩阵

        for (int ele_row = 0; ele_row < 2 * 4; ele_row++)
        {
            int global_pos_row =
                dof_element[i2(i_element, ele_row, 2 * 4)];

            if (global_pos_row == -1)
                continue; //此点被约束

            for (int ele_col = 0; ele_col < 2 * 4; ele_col++)
            {
                int global_pos_col =
                    dof_element[i2(i_element, ele_col, 2 * 4)];

                if (global_pos_col == -1)
                    continue; //此点被约束

                if (global_pos_col > global_pos_row) //上三角区域
                    continue;

                double t_stiff =
                    element_stiff[i2(ele_row, ele_col, 2 * 4)];

                // printf("\t%f", t_stiff);

                // 行列的选择：MA[i] - i + j - 1; -1 是为了找到内存中的位置

                int pos = diag_pos[global_pos_row] - global_pos_row + global_pos_col - 1;
                global_stiffness[pos] += t_stiff;

                // printf("\t|%d,%d|%d", global_pos_row, global_pos_col, pos);
            }
            // printf("\n");
        }
        // printf("\n");
    }
    //
    // logging
    /*
    for (int i = 0; i < _MAXDOF; i++)
    {
        for (int j = 0; j < _MAXDOF; j++)
        {
            if (j > i)
            {
                fprintf(notefile, "\t%0.2e", 0.0);
                continue;
            }
            int pos = diag_pos[i] - i + j - 1;

            if ((pos < diag_pos[i - 1] - 1) && (i > 0))
            {
                fprintf(notefile, "\t%0.2e", 0.0);
                continue;
            }
            fprintf(notefile, "\t%0.2e", global_stiffness[pos]);
        }
        fprintf(notefile, "\n");
    }
    */
    // logging
    /*
    for (int i = 0; i < _MAXDOF; i++)
    {
        for (int j = 0; j < _MAXDOF; j++)
        {
            fprintf(notefile, "\t%0.2e", get_K_num(global_stiffness, i, j, diag_pos));
        }
        fprintf(notefile, "\n");
        // printf("%e\n", get_K_num(global_stiffness, i, i, diag_pos));
    }
    */
    return 0;
}

int fem_element_stiffness(int element, double *element_stiff)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // 分成四个2*2区域计算
            double *sub_stiff;
            sub_stiff = (double *)calloc(2 * 2, sizeof(double));

            double *sub_dint;
            sub_dint = (double *)calloc(2 * 2, sizeof(double));

            // 获取分块的 DXX 等，存入矩阵 sub_dint
            fem_sub_element_DINT(element, i, j, sub_dint);

            // 局部刚度计算
            fem_cal_sub(sub_stiff, sub_dint);

            // 组装单元刚度矩阵
            fem_ass_sub(sub_stiff, element_stiff, i, j);
        }
    }

    // testing
    /*
    printf("\telement:%d\n", element);

    for (int i = 0; i < 8 * 8; i++)
    {
        if (i % 8 == 0)
            printf("\n");
        printf("\t%f", element_stiff[i]);
    }
    */

    return 0;
}

int fem_ass_sub(double *sub_stiff,
                double *element_stiff,
                int irow,
                int jcol)
{
    // element_stiff = 8 * 8
    // testing

    /*
    for (int in = 0; in < 4; in++)
    {
        if (in % 2 == 0)
            printf("\n");
        printf("\t%f", sub_stiff[in]);
    }
    printf("\n");
    */

    int true_row = irow * 2;
    int true_col = jcol * 2;

    for (int subrow = 0; subrow < 2; subrow++)
    {
        for (int subcol = 0; subcol < 2; subcol++)
        {
            int trow = true_row + subrow;
            int tcol = true_col + subcol;

            // printf("\t%d\t%d\n", trow, tcol);

            element_stiff[i2(trow, tcol, 8)] =
                sub_stiff[i2(subrow, subcol, 2)];
        }
    }

    return 0;
}

int fem_cal_sub(double *sub_stiff,
                double *sub_dint)
{
    // sub_dint : 2*2
    // sub_stiff : 2*2
    double *ela_d;
    ela_d = (double *)calloc(3, sizeof(double)); //弹性量矩阵 D
    fem_get_elastic_d(ela_d);                    //获取这个矩阵中的三个有效数据

    sub_stiff[i2(0, 0, 2)] = ela_d[0] * sub_dint[i2(0, 0, 2)] +
                             ela_d[2] * sub_dint[i2(1, 1, 2)];

    sub_stiff[i2(0, 1, 2)] = ela_d[1] * sub_dint[i2(0, 1, 2)] +
                             ela_d[2] * sub_dint[i2(1, 0, 2)];

    sub_stiff[i2(1, 0, 2)] = ela_d[1] * sub_dint[i2(1, 0, 2)] +
                             ela_d[2] * sub_dint[i2(0, 1, 2)];

    sub_stiff[i2(1, 1, 2)] = ela_d[0] * sub_dint[i2(1, 1, 2)] +
                             ela_d[2] * sub_dint[i2(0, 0, 2)];

    /*
    for (int dx = 0; dx < 2; dx++)
    {
        for (int dy = 0; dy < 2; dy++)
        {
            sub_stiff[i2(dx, dy, 2)] = ela_d[dx ^ dy] * sub_dint[i2(dx, dy, 2)] +
                                       ela_d[2] * sub_dint[i2(1 - dx, 1 - dy, 2)];
        }
    }
    */

    return 0;
}

int fem_get_elastic_d(double *ela_d)
// ela_d: 3*3
{
    double E  = deform_info[i2(0, 0, 2)]; //# 弹模(KPa)
    double nu = deform_info[i2(0, 1, 2)]; //# 泊松比
    // printf("\t%f\t%f\n", E, nu);
    ela_d[0] = ((1 - nu) * E) /
               ((1 + nu) * (1 - 2 * nu));

    ela_d[1] = (nu * E) /
               ((1 + nu) * (1 - 2 * nu));

    ela_d[2] = (E) /
               (2 * (1 + nu));

    // printf("\t%f\t%f\t%f\n", ela_d[0], ela_d[1], ela_d[2]);

    return 0;
}

int fem_sub_element_DINT(int element,
                         int irow, int jcol,
                         double *sub_dint)
{
    // irow jcol是单元刚度矩阵里头的对应小块的编号
    // printf("\nele: %d sub: irow:%d jcol:%d\n", element, irow, jcol);
    // sub_dint 参见书P171页
    // 获取 sub_dint 矩阵

    // 这一个大循环就是获得了sub_dint矩阵

    for (int i = 0; i < _GAUSS_POINTS; i++)
    {
        // 对高斯积分点循环
        for (int j = 0; j < _GAUSS_POINTS; j++)
        {
            // 高斯积分点
            double gauss_rs[] = {
                _GAUSS_POS[i],
                _GAUSS_POS[j]};

            double gauss_h[] = {
                _GAUSS_Height[i],
                _GAUSS_Height[j]};

            /*
            printf("gauss_rs:\t%f\t%f\n", gauss_rs[0], gauss_rs[1]);
            printf("gauss_h:\t%f\t%f\t%f\n",
                   gauss_h[0], gauss_h[1],
                   gauss_h[0] * gauss_h[1]);
            */
            // 这一点处的jacobi矩阵
            double *tjacobi;
            tjacobi = (double *)calloc(2 * 2, sizeof(double));

            // 获取这一点的 Nixy 矩阵，以及 jacobi 矩阵;
            // Nixy形函数对整体坐标的导数
            double *nixy;
            nixy = (double *)calloc(4 * 2, sizeof(double));
            //
            fem_findnixy(element, tjacobi, nixy, gauss_rs);

            double det_jaco = m_det(tjacobi, 2);
            // printf("\t%f\n", det_jaco);

            if (det_jaco <= 0)
            {
                printf("JACOBI ERROR\n");
                exit(0);
            }

            //求得 sub_dint 矩阵
            for (int dx = 0; dx < 2; dx++)
            {
                for (int dy = 0; dy < 2; dy++)
                {
                    // printf("\t%d\t%d\t%d\t%d\n", irow, jcol, dx, dy);
                    /*
                    printf("\t%f\t%f\t%f\n",
                           nixy[i2(irow, dx, 2)],
                           nixy[i2(jcol, dy, 2)],
                           det_jaco);
                    */

                    sub_dint[i2(dx, dy, 2)] += nixy[i2(irow, dx, 2)] *
                                               nixy[i2(jcol, dy, 2)] *
                                               det_jaco *
                                               gauss_h[0] * gauss_h[1];
                }
            }

            /*
            for (int i = 0; i < 4; i++)
            {
                if (i % 2 == 0)
                    printf("\n");
                printf("\t%.2f", tjacobi[i]);
            }
            */

            // printf("\t%.2f", m_det(tjacobi, 2));
        }
    }
    return 0;
}

int fem_findnixy(int element, double *jacobi, double *nixy, double *gauss_rs)
{
    // 获取单元四个节点的坐标
    double x[] = {
        coordx[element_4info[n1(element)]],
        coordx[element_4info[n2(element)]],
        coordx[element_4info[n3(element)]],
        coordx[element_4info[n4(element)]]};

    double y[] = {
        coordy[element_4info[n1(element)]],
        coordy[element_4info[n2(element)]],
        coordy[element_4info[n3(element)]],
        coordy[element_4info[n4(element)]]};

    // gauss_rs是局部点坐标
    // printf("\t%.3f\t%.3f\n", x[2], y[2]);
    // printf("\t%.3f\t%.3f\n", gauss_rs[0], gauss_rs[1]);
    // nirs是对局部坐标的导数矩阵，带入r、s值，获取矩阵
    double *nirs;
    nirs = (double *)calloc(4 * 2, sizeof(double));
    // 此处是对应积分点的局部坐标导数
    fem_get_pNip(nirs, gauss_rs[0], gauss_rs[1]);

    /*
    for (int i = 0; i < 2 * 4; i++)
    {
        if (i % 2 == 0)
            printf("\n");
        printf("\t%0.7f", nirs[i]);
    }
    */

    //计算jacobi矩阵
    for (int i = 0; i < 4; i++)
    {
        jacobi[i2(0, 0, 2)] += nirs[i2(i, 0, 2)] * x[i];
        jacobi[i2(0, 1, 2)] += nirs[i2(i, 0, 2)] * y[i];
        jacobi[i2(1, 0, 2)] += nirs[i2(i, 1, 2)] * x[i];
        jacobi[i2(1, 1, 2)] += nirs[i2(i, 1, 2)] * y[i];
    }

    //计算jacobi矩阵的行列式的值
    double jacobi_det = m_det(jacobi, 2);

    // testing
    /*
    for (int i = 0; i < 4; i++)
    {
        if (i % 2 == 0)
            printf("\n");
        printf("\t%.2f", jacobi[i]);
    }

    printf("\t%.2f", jacobi_det);
    */

    // inverse matrix of the jacobi
    // 计算jacobi行列式的逆矩阵
    double *inv_jaco;
    inv_jaco = (double *)calloc(2 * 2, sizeof(double));

    inv_jaco[i2(0, 0, 2)] = jacobi[i2(1, 1, 2)] / jacobi_det;
    inv_jaco[i2(0, 1, 2)] = -jacobi[i2(0, 1, 2)] / jacobi_det;
    inv_jaco[i2(1, 0, 2)] = -jacobi[i2(1, 0, 2)] / jacobi_det;
    inv_jaco[i2(1, 1, 2)] = jacobi[i2(0, 0, 2)] / jacobi_det;

    // testing
    /*
    for (int i = 0; i < 4; i++)
    {
        if (i % 2 == 0)
            printf("\n");
        printf("\t%.2f", inv_jaco[i]);
    }
    printf("\n\n");
    */

    // nixy 是形函数Ni对整体坐标的导数矩阵
    // 行是形函数Ni，列是x,y

    // 相乘得到需要的形函数对整体坐标的导数
    for (int i = 0; i < 4; i++)
    {
        nixy[i2(i, 0, 2)] = inv_jaco[i2(0, 0, 2)] * nirs[i2(i, 0, 2)] +
                            inv_jaco[i2(0, 1, 2)] * nirs[i2(i, 1, 2)]; // i, x

        nixy[i2(i, 1, 2)] = inv_jaco[i2(1, 0, 2)] * nirs[i2(i, 0, 2)] +
                            inv_jaco[i2(1, 1, 2)] * nirs[i2(i, 1, 2)]; // i, y

        // printf("\t%f\t%f\n", nixy[i2(i, 0, 2)], nixy[i2(i, 1, 2)]);
    }

    return 0;
}

// r,s是局部坐标，实际就是高斯积分点
// 此函数获得局部导数值
int fem_get_pNip(double *nirs, double r, double s)
{
    // printf("\t%f\t%f\n", r, s);

    nirs[i2(0, 0, 2)] = -(1 - s) / 4.0; // 1, x
    nirs[i2(0, 1, 2)] = -(1 - r) / 4.0; // 1, y

    nirs[i2(1, 0, 2)] = (1 - s) / 4.0;
    nirs[i2(1, 1, 2)] = -(1 + r) / 4.0;

    nirs[i2(2, 0, 2)] = (1 + s) / 4.0;
    nirs[i2(2, 1, 2)] = (1 + r) / 4.0;

    nirs[i2(3, 0, 2)] = -(1 + s) / 4.0;
    nirs[i2(3, 1, 2)] = (1 - r) / 4.0;

    return 0;
}

#endif