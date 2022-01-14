#pragma once
#ifndef _FEM_BASIC
#define _FEM_BASIC

#define i2(i, j, n) (((i) * (n)) + (j))
#define n1(i) ((i)*6)
#define n2(i) (((i)*6) + 1)
#define n3(i) (((i)*6) + 2)
#define n4(i) (((i)*6) + 3)
#define stepid(i) (((i)*6) + 4)
#define matid(i) (((i)*6) + 5)

#define ELEMENT_DATA_NUM 6
#define ELEMENT_NODE_NUM 4 //矩形单元
#define MAT_SEEP_NUM 5
#define MAT_COMMON_NUM 5
#define MAT_SLOPE_NUM 10
#define DEF_INFO_NUM 2
#define NOD_CON_NUM 6
#define LOD_INFO_NUM 2
#define _PI 3.14159265358979

#include <cstdio>
#include <iostream>

#endif