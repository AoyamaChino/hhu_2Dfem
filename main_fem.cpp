#include "fem_basic.h"
#include "fem_global.h"
#include "fem_input.h"
#include "fem_inter.h"
#include "fem_kernel.h"
#include "fem_mem.h"
#include "fem_output.h"
#include "fem_setgauss.h"
#include "fem_solve.h"

/*
LIN FENG 2022
College of Mechanics and Materials,
Hohai University.
*/

int main(int argc, char const *argv[])
{
    /*
    double test[] = {1, 2,
                     2, 4};
    printf("%.2f\n", m_det(test, 2));
    */

    printf("\n"
           "Note that this program is only for 2D problem\n"
           "The element is fixed to 4 nodes element\n"
           "Only for consentrated force and displacement boundary\n"
           "\n");

    notefile = fopen("./note.txt", "w");
    //
    inputfile = fopen("./input/recbing.inp", "r");
    fem_input();       // input data
    fclose(inputfile); //
    //
    fem_inter();            // in file fem_inter
    fem_dofserial();        // in file fem_kernel
    fem_diagonal_pos();     // in file fem_kernel
    fem_setgauss(3);        // in file setgauss
    fem_global_stiffness(); // in file fem_kernel get global stiff mat
    //
    _MAX_STEP = 20000;
    // solve_axb_CG();
    solve_eigen();
    //
    fem_after();
    fem_output();
    // fem_memfree();
    //
    fclose(notefile);
    return 0;
}