/*
! Module for computing basis functions and their derivatives
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 18-May-2018
! Place : Swansea, UK
!
!
*/

#include "BasisFunctionsBernstein.h"
#include "headersBasic.h"
#include "headersEigen.h"



int computeBasisFunctions2D(bool flag, int ELEM_TYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac)
{
    int  ii=0, jj=0, count=0, nlbf=0;

    if( ELEM_TYPE == 1 )
    {
      if(degree == 1)
      {
        count = 3;
        nlbf = count;
      }
      else if(degree == 2)
      {
        count = 6;
        nlbf = count;
      }
    }
    else if( ELEM_TYPE == 2 )
    {
      count = 9;
      nlbf = count;
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      exit(-2);
    }

    vector<double>  N1(count), N2(count), dN1(count), dN2(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf);
    double  xx, yy, detinv, B[2][2], Binv[2][2] ;

    if( ELEM_TYPE == 1 )
    {
      BernsteinBasisFunsTria(degree, param[0], param[1], &N[0], &dN_du1[0], &dN_du2[0]);
    }
    else if( ELEM_TYPE == 2 )
    {
      BernsteinBasisFunsQuad(degree, param[0], param[1], &N[0], &dN_du1[0], &dN_du2[0]);
    }
    else
    {
      //cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      //exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B[0][0] = B[1][0] = B[0][1] = B[1][1] = 0.0;
    for(ii=0; ii<nlbf; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];

      B[0][0] +=  (xx * dN_du1[ii]) ;
      B[1][0] +=  (xx * dN_du2[ii]) ;
      B[0][1] +=  (yy * dN_du1[ii]) ;
      B[1][1] +=  (yy * dN_du2[ii]) ;
    }

    Jac  = B[0][0]*B[1][1] - B[0][1]*B[1][0];

    //printf("Bmat  \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n\n", B[0][0], B[0][1], B[1][0], B[1][1], Jac);

    detinv = 1.0/Jac ;

    Binv[0][0] =  B[1][1] * detinv;
    Binv[0][1] = -B[0][1] * detinv;
    Binv[1][0] = -B[1][0] * detinv;
    Binv[1][1] =  B[0][0] * detinv;

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv[0][0] + dN_du2[ii] * Binv[0][1];
      dN_dy[ii] = dN_du1[ii] * Binv[1][0] + dN_du2[ii] * Binv[1][1];
    }

    return 0;
}



int computeBasisFunctions3D(bool flag, int ELEM_TYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac)
{
    int  ii=0, jj=0, count=0, nlbf=0;

    if( ELEM_TYPE == 1 )
    {
      if(degree == 1)
      {
        count = 4;
        nlbf = count;
      }
      else if(degree == 2)
      {
        count = 10;
        nlbf = count;
      }
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions3D... " << endl;
      exit(-2);
    }

    vector<double>  N1(count), N2(count), dN1(count), dN2(count), dN3(count);
    vector<double>  dN_du1(nlbf), dN_du2(nlbf), dN_du3(nlbf);
    double  xx, yy, zz, detinv;
    MatrixXd B(3,3), Binv(3,3) ;

    if( ELEM_TYPE == 1 )
    {
      BernsteinBasisFunsTetra(degree, param[0], param[1], param[2], &N[0], &dN_du1[0], &dN_du2[0], &dN_du3[0]);
    }
    else
    {
      //cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      //exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B.setZero();
    for(ii=0; ii<nlbf; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];
      zz = zNode[ii];

      B(0,0) +=  (xx * dN_du1[ii]) ;
      B(1,0) +=  (xx * dN_du2[ii]) ;
      B(2,0) +=  (xx * dN_du3[ii]) ;

      B(0,1) +=  (yy * dN_du1[ii]) ;
      B(1,1) +=  (yy * dN_du2[ii]) ;
      B(2,1) +=  (yy * dN_du3[ii]) ;

      B(0,2) +=  (zz * dN_du1[ii]) ;
      B(1,2) +=  (zz * dN_du2[ii]) ;
      B(2,2) +=  (zz * dN_du3[ii]) ;
    }

    //printMatrix(B);

    Jac  = B.determinant();
    Binv = B.inverse();

    //printMatrix(Binv);

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(ii=0; ii<nlbf; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1) + dN_du3[ii] * Binv(0,2);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1) + dN_du3[ii] * Binv(1,2);
      dN_dz[ii] = dN_du1[ii] * Binv(2,0) + dN_du2[ii] * Binv(2,1) + dN_du3[ii] * Binv(2,2);
    }

    return 0;
}



