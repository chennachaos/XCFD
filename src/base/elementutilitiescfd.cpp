/*
! Subroutines for calculating the stiffness
! and residual (or force) for the simple
! Lagrange family of elements for
!
! Incompressible Navier-Stokes using
! the stabilised formulation
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-May-2018
! Place : Swansea, UK
!
!
*/

#include "headersBasic.h"
#include "elementutilitiescfd.h"
#include "BasisFunctionsBernstein.h"


int  pointInverseTria6node(double* xNode, double* yNode, double* target, double* param)
{
    MatrixXd  J(2,2), Jinv(2,2);
    VectorXd  xi(2), f(2);

    double  fnorm;
    vector<double>  N(6), dN_du1(6), dN_du2(6);
    int ii, iter;

    xi.setZero();

    xi(0) = param[0];  xi(1) = param[1];
    for(iter=0; iter<10; iter++)
    {
      BernsteinBasisFunsTria(2, xi(0), xi(1), &N[0], &dN_du1[0], &dN_du2[0]);

      //J(0,0) = -2.0*xi3*xNode[0] + 2.0*xi1*xNode[1] + 2.0*(xi3-xi1)*xNode[3] + 2.0*xi2*xNode[4] - 2.0*xi2*xNode[5];
      //J(0,1) = -2.0*xi3*xNode[0] + 2.0*xi2*xNode[2] - 2.0*xi1*xNode[3] + 2.0*xi1*xNode[4] - 2.0*(xi3-xi2)*xNode[5];

      //J(1,0) = -2.0*xi3*yNode[0] + 2.0*xi1*yNode[1] + 2.0*(xi3-xi1)*yNode[3] + 2.0*xi2*yNode[4] - 2.0*xi2*yNode[5];
      //J(1,1) = -2.0*xi3*yNode[0] + 2.0*xi2*yNode[2] - 2.0*xi1*yNode[3] + 2.0*xi1*yNode[4] - 2.0*(xi3-xi2)*yNode[5];

      f.setZero();
      J.setZero();
      for(ii=0; ii<6; ii++)
      {
        f(0) += xNode[ii]*N[ii];
        f(1) += yNode[ii]*N[ii];

        J(0,0) +=  xNode[ii]*dN_du1[ii] ;
        J(0,1) +=  xNode[ii]*dN_du2[ii] ;
        J(1,0) +=  yNode[ii]*dN_du1[ii] ;
        J(1,1) +=  yNode[ii]*dN_du2[ii] ;
      }

      f(0) -= target[0];
      f(1) -= target[1];

      fnorm = f.norm();

      //cout << iter << '\t' << fnorm << '\t' << xi(0) << '\t' << xi(1) << endl;

      if( fnorm < 1.0e-8 )
      {
        break;
      }

      Jinv = J.inverse();

      xi -= Jinv*f;
    }

    //cout << " xcoord = " << target[0] << '\t' << xi(0) << endl;
    //cout << " ycoord = " << target[1] << '\t' << xi(1) << endl;

    param[0] = xi(0);
    param[1] = xi(1);

    return 0;
}




int  pointInverseTetra10node(double* xNode, double* yNode, double* zNode, double* target, double* param)
{
    MatrixXd  J(3,3), Jinv(3,3);
    VectorXd  xi(3), f(3);

    double  fnorm;
    vector<double>  N(10), dN_du1(10), dN_du2(10), dN_du3(10);
    int ii, iter;

    xi.setZero();

    xi(0) = param[0];  xi(1) = param[1];  xi(2) = param[2];
    for(iter=0; iter<10; iter++)
    {
      BernsteinBasisFunsTetra(2, xi(0), xi(1), xi(2), &N[0], &dN_du1[0], &dN_du2[0], &dN_du3[0]);

      f.setZero();
      J.setZero();
      for(ii=0; ii<10; ii++)
      {
        f(0) += xNode[ii]*N[ii];
        f(1) += yNode[ii]*N[ii];
        f(2) += zNode[ii]*N[ii];

        J(0,0) +=  xNode[ii]*dN_du1[ii] ;
        J(0,1) +=  xNode[ii]*dN_du2[ii] ;
        J(0,2) +=  xNode[ii]*dN_du3[ii] ;
        
        J(1,0) +=  yNode[ii]*dN_du1[ii] ;
        J(1,1) +=  yNode[ii]*dN_du2[ii] ;
        J(1,2) +=  yNode[ii]*dN_du3[ii] ;

        J(2,0) +=  zNode[ii]*dN_du1[ii] ;
        J(2,1) +=  zNode[ii]*dN_du2[ii] ;
        J(2,2) +=  zNode[ii]*dN_du3[ii] ;
      }

      f(0) -= target[0];
      f(1) -= target[1];
      f(2) -= target[2];

      fnorm = f.norm();

      //cout << iter << '\t' << fnorm << '\t' << xi(0) << '\t' << xi(1) << endl;

      if( fnorm < 1.0e-8 )
      {
        break;
      }

      Jinv = J.inverse();

      xi -= Jinv*f;
    }

    //cout << " xcoord = " << target[0] << '\t' << xi(0) << endl;
    //cout << " ycoord = " << target[1] << '\t' << xi(1) << endl;
    //cout << " zcoord = " << target[2] << '\t' << xi(2) << endl;

    param[0] = xi(0);
    param[1] = xi(1);
    param[2] = xi(2);

    return 0;
}
