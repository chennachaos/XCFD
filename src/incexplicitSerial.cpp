/*
! Program for 2D explicit finite element analysis of incompressible Navier-Stokes
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-May-2018
! Place : Swansea, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "ExplicitCFD.h"

using namespace std;



int main(int argc, char* argv[])
{
    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    string  infile = argv[1];


    ExplicitCFD  explcfd;

    explcfd.readInputData(infile);

    explcfd.readControlParameters();

    explcfd.prepareInputData();

    explcfd.solveExplicitStep();

    //explcfd.computeElementErrors(1);

    cout << " Program is successful \n " << endl;

    return 1;
}

