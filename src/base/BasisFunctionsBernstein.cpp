
#include "BasisFunctionsBernstein.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>

using std::vector;

double  norm(vector<double>&  vec)
{
    double  val = 0.0;
    for(int ii=0;ii<vec.size();ii++)
        val += vec[ii]*vec[ii];

   return sqrt(val);
}

void  subtract2vecs(const vector<double>&  vec1, const vector<double>&  vec2, vector<double>&  vec)
{
    for(int ii=0;ii<vec1.size();ii++)
        vec[ii] = vec1[ii] - vec2[ii];
}


void Bernstein_BasisFuns1D(int p, double xi, double* N, double* dN_dxi)
{
  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi[0] = 0.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

          dN_dxi[0] = -0.5;
          dN_dxi[1] =  0.5;

      break;

      case 2:

          N[0] = 0.25*(1.0-xi)*(1.0-xi);
          N[1] = 0.50*(1.0-xi*xi);
          N[2] = 0.25*(1.0+xi)*(1.0+xi);

          dN_dxi[0] =  0.5*(xi-1.0);
          dN_dxi[1] = -xi;
          dN_dxi[2] =  0.5*(xi+1.0);

          //N[0] = (1.0-xi)*(1.0-xi);
          //N[1] =  2.0*xi*(1.0-xi);
          //N[2] =  xi*xi;

          //dN_dxi[0] = -2.0*(1.0-xi);
          //dN_dxi[1] =  2.0*(1.0-2.0*xi);
          //dN_dxi[2] =  2.0*xi;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void Bernstein_BasisFuns1D(int p, double xi, double* N)
{
  switch(p)
  {
      case 0:

          N[0] = 1.0;

      break;

      case 1:

          N[0] = 0.5*(1.0 - xi);
          N[1] = 0.5*(1.0 + xi);

      break;

      case 2:

          N[0] = (1.0-xi)*(1.0-xi);
          N[1] =  2.0*xi*(1.0-xi);
          N[2] =  xi*xi;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }

   return;
}



void BernsteinBasisFunsTria(int p, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
  double  xi3 = 1.0 - xi1 - xi2;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;

      break;

      case 1:

          N[0] = xi3;
          N[1] = xi1;
          N[2] = xi2;

          dN_dxi1[0] = -1.0;
          dN_dxi1[1] =  1.0;
          dN_dxi1[2] =  0.0;

          dN_dxi2[0] = -1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;

      break;

      case 2:

//
          N[0] = xi3*xi3;
          N[1] = xi1*xi1;
          N[2] = xi2*xi2;
          N[3] = 2.0*xi1*xi3;
          N[4] = 2.0*xi1*xi2;
          N[5] = 2.0*xi2*xi3;

          dN_dxi1[0] = -2.0*xi3;
          dN_dxi1[1] =  2.0*xi1;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  2.0*(xi3 - xi1);
          dN_dxi1[4] =  2.0*xi2;
          dN_dxi1[5] = -2.0*xi2;

          dN_dxi2[0] = -2.0*xi3;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  2.0*xi2;
          dN_dxi2[3] = -2.0*xi1;
          dN_dxi2[4] =  2.0*xi1;
          dN_dxi2[5] =  2.0*(xi3 - xi2);
//

/*
          N[0] = xi3*(2.0*xi3 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = 4.0*xi1*xi3;
          N[4] = 4.0*xi1*xi2;
          N[5] = 4.0*xi2*xi3;

          dN_dxi1[0] = -4.0*xi3 + 1.0;
          dN_dxi1[1] =  4.0*xi1 - 1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  4.0*(xi3 - xi1);
          dN_dxi1[4] =  4.0*xi2;
          dN_dxi1[5] = -4.0*xi2;

          dN_dxi2[0] = -4.0*xi3 + 1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  4.0*xi2 - 1.0;
          dN_dxi2[3] = -4.0*xi1;
          dN_dxi2[4] =  4.0*xi1;
          dN_dxi2[5] =  4.0*(xi3 - xi2);
*/
      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }
  return;
}



void BernsteinBasisFunsQuad(int p, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
    double  a1, a2, a3, b1, b2, b3;
    double  da1, da2, da3, db1, db2, db3;

    switch(p)
    {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;

      break;

      case 1:

          b1 = 0.5*(1.0 - xi2);
          b2 = 0.5*(1.0 + xi2);

          db1 = -0.5;
          db2 =  0.5;

          a1 = 0.5*(1.0 - xi1);
          a2 = 0.5*(1.0 + xi1);

          da1 = -0.5;
          da2 =  0.5;

          N[0] = b1*a1;
          N[1] = b1*a2;
          N[2] = b2*a2;
          N[3] = b2*a1;

          dN_dxi1[0] = b1*da1;
          dN_dxi1[1] = b1*da2;
          dN_dxi1[2] = b2*da2;
          dN_dxi1[3] = b2*da1;

          dN_dxi2[0] = db1*a1;
          dN_dxi2[1] = db1*a2;
          dN_dxi2[2] = db2*a2;
          dN_dxi2[3] = db2*a1;

      break;

      case 2:

          b1 = 0.25*(1.0-xi2)*(1.0-xi2);
          b2 = 0.50*(1.0-xi2*xi2);
          b3 = 0.25*(1.0+xi2)*(1.0+xi2);

          db1 =  0.5*(xi2-1.0);
          db2 = -xi2;
          db3 =  0.5*(xi2+1.0);

          a1 = 0.25*(1.0-xi1)*(1.0-xi1);
          a2 = 0.50*(1.0-xi1*xi1);
          a3 = 0.25*(1.0+xi1)*(1.0+xi1);

          da1 =  0.5*(xi1-1.0);
          da2 = -xi1;
          da3 =  0.5*(xi1+1.0);

          N[0] = b1*a1;
          N[1] = b1*a3;
          N[2] = b3*a3;
          N[3] = b3*a1;
          N[4] = b1*a2;
          N[5] = b2*a3;
          N[6] = b3*a2;
          N[7] = b2*a1;
          N[8] = b2*a2;

          dN_dxi1[0] = b1*da1;
          dN_dxi1[1] = b1*da3;
          dN_dxi1[2] = b3*da3;
          dN_dxi1[3] = b3*da1;
          dN_dxi1[4] = b1*da2;
          dN_dxi1[5] = b2*da3;
          dN_dxi1[6] = b3*da2;
          dN_dxi1[7] = b2*da1;
          dN_dxi1[8] = b2*da2;

          dN_dxi2[0] = db1*a1;
          dN_dxi2[1] = db1*a3;
          dN_dxi2[2] = db3*a3;
          dN_dxi2[3] = db3*a1;
          dN_dxi2[4] = db1*a2;
          dN_dxi2[5] = db2*a3;
          dN_dxi2[6] = db3*a2;
          dN_dxi2[7] = db2*a1;
          dN_dxi2[8] = db2*a2;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
    }

    return;
}






void BernsteinBasisFunsTetra(int p, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  xi4 = 1.0 - xi1 - xi2 - xi3;

  switch(p)
  {
      case 0:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;
          dN_dxi3[0] = 0.0;

      break;

      case 1:

          N[0] = xi4;
          N[1] = xi1;
          N[2] = xi2;
          N[3] = xi3;

          dN_dxi1[0] = -1.0;
          dN_dxi1[1] =  1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  0.0;

          dN_dxi2[0] = -1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;
          dN_dxi2[3] =  0.0;

          dN_dxi3[0] = -1.0;
          dN_dxi3[1] =  0.0;
          dN_dxi3[2] =  0.0;
          dN_dxi3[3] =  1.0;

      break;

      case 2: //xi4 = 1.0 - xi1 - xi2 - xi3;

          N[0] = xi4*xi4;
          N[1] = xi1*xi1;
          N[2] = xi2*xi2;
          N[3] = xi3*xi3;
          N[4] = 2.0*xi1*xi4;
          N[5] = 2.0*xi1*xi2;
          N[6] = 2.0*xi2*xi4;
          N[7] = 2.0*xi4*xi3;
          N[8] = 2.0*xi1*xi3;
          N[9] = 2.0*xi2*xi3;

          dN_dxi1[0] = -2.0*xi4;
          dN_dxi1[1] =  2.0*xi1;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  0.0;
          dN_dxi1[4] =  2.0*(xi4-xi1);
          dN_dxi1[5] =  2.0*xi2;
          dN_dxi1[6] = -2.0*xi2;
          dN_dxi1[7] = -2.0*xi3;
          dN_dxi1[8] =  2.0*xi3;
          dN_dxi1[9] =  0.0;

          dN_dxi2[0] = -2.0*xi4;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  2.0*xi2;
          dN_dxi2[3] =  0.0;
          dN_dxi2[4] = -2.0*xi1;
          dN_dxi2[5] =  2.0*xi1;
          dN_dxi2[6] =  2.0*(xi4-xi2);
          dN_dxi2[7] = -2.0*xi3;
          dN_dxi2[8] =  0.0;
          dN_dxi2[9] =  2.0*xi3;

          dN_dxi3[0] = -2.0*xi4;
          dN_dxi3[1] =  0.0;
          dN_dxi3[2] =  0.0;
          dN_dxi3[3] =  2.0*xi3;
          dN_dxi3[4] = -2.0*xi1;
          dN_dxi3[5] =  0.0;
          dN_dxi3[6] = -2.0*xi2;
          dN_dxi3[7] =  2.0*(xi4-xi3);
          dN_dxi3[8] =  2.0*xi1;
          dN_dxi3[9] =  2.0*xi2;

      break;

      default:

          printf("no basis functions defined for this degree = %5d \n", p);

      break;
  }
  return;
}




void BernsteinBasisFunsEdge2D(int p, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac)
{
   int  ii, jj, count, nlbf = p+1;

   double  du_dx, dx, dy ;
   vector<double>  dN1(nlbf);

   Bernstein_BasisFuns1D(p, param[0], N, &dN1[0]);

   if(p == 0)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = 0.0;
     for(ii=0; ii<nlbf; ii++)
     {
       dx +=  (xNode[ii] * dN1[ii]);
       dy +=  (yNode[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy);
   }

   du_dx = 1.0/Jac ;

  // Compute the normal
  normal[0] =  dy/Jac;
  normal[1] = -dx/Jac;

  return;
}



void BernsteinBasisFunsFace3D(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
  return;
}


void BernsteinBasisFunsFaceTria(int degree, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
   int  ii, jj, count, nlbf;
   if(degree == 1)
   {
     nlbf = 3;
   }
   else if(degree == 2)
   {
     nlbf = 6;
   }

   vector<double>  dN1(nlbf), dN2(nlbf);

   BernsteinBasisFunsTria(degree, param[0], param[1], N, &dN1[0], &dN2[0]);

   double  du[3], dv[3];
   du[0] = du[1] = du[2] = 0.0;
   dv[0] = dv[1] = dv[2] = 0.0;
   for(ii=0; ii<nlbf; ii++)
   {
     du[0] += (xNode[ii] * dN1[ii]);
     du[1] += (yNode[ii] * dN1[ii]);
     du[2] += (zNode[ii] * dN1[ii]);

     dv[0] += (xNode[ii] * dN2[ii]);
     dv[1] += (yNode[ii] * dN2[ii]);
     dv[2] += (zNode[ii] * dN2[ii]);
   }

   // Compute the normal
   normal[0] = du[1]*dv[2] - du[2]*dv[1];
   normal[1] = du[2]*dv[0] - du[0]*dv[2];
   normal[2] = du[0]*dv[1] - du[1]*dv[0];

   Jac = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

   // normalise -> unit normal
   normal[0] /= Jac;
   normal[1] /= Jac;
   normal[2] /= Jac;

   return;
}


void BernsteinBasisFunsFaceQuad(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
  return;
}



