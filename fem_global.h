#pragma once
#ifndef _FEM_GLOBAL
#define _FEM_GLOBAL
// space for global parameters and arrays
#include "fem_basic.h"

int
    _MAX_STEP,
    Node_num,
    Element_num,
    *node_constraint_info,
    load_count,
    _DISPCONS_COUNT,
    *element_4info, //#EleId:N1, N2, N3, N4, StepId, MatId; info; size of each box:6 * int;
    Material_num,
    *load_nodes,
    *ifxfixd, // if the node is fixed in direction x
    *ifyfixd,
    *dof_global,
    *dof_element,
    _MAXDOF = 0,
    *diag_pos,
    _MAX_DIAG_POS; // also the size of the stiffness matrix(in form of 1D memory)

double
    *loadnodexy,
    *coordx,
    *coordy,
    *dispx,
    *dispy,
    *stressx,
    *stressy,
    *strainx,
    *strainy,
    *global_stiffness, //整体刚度矩阵
    *mat_seep,         // id,Kx(m/s),Ky(m/s),Kz(m/s),孔隙率,给水度
    *mat_common,       // id,材料名,容重(KN/m3),浮容重(KN/m3),是新筑土层,采用非线性
    *mat_slope_dry,    //有效应力 c(KPa),有效应力 摩擦角(度),固结不排水 c(KPa),固结不排水 摩擦角(度),总应力 c(KPa),总应力 摩擦角(度),非线性Fi0(度),非线性dFi(度),施工期孔压系数,降落期孔压系数
    *mat_slope_wet,    // id,_有效应力 c(KPa),_有效应力 摩擦角(度),_固结不排水 c(KPa),_固结不排水 摩擦角(度),_总应力 c(KPa),_总应力 摩擦角(度),_非线性Fi0(度),_非线性dFi(度),_施工期孔压系数,_降落期孔压系数
    *deform_info,
    *boundary_info,
    *load_info,
    *dof_disp; //依据自由度编号的解空间

// gauss points
int _GAUSS_POINTS;
double _GAUSS_POS[50];
double _GAUSS_Height[50];

double m_det(double *m, int size)
{
    if (size != 2)
    {
        printf("matrix det error\n");
        exit(0);
    }
    double value =
        m[i2(0, 0, 2)] * m[i2(1, 1, 2)] - m[i2(0, 1, 2)] * m[i2(1, 0, 2)];

    return value;
}

#endif