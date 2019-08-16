
/* Node ordering for the quadratic triangular element
 4-----7-----3
 |     |     |
 |     |     |
 8-----9-----6
 |     |     |
 |     |     |
 1-----5-----2
*/

#include "BernsteinElem2DINSQuad9Node.h"
#include "BasisFunctionsBernstein.h"
#include "headersBasic.h"
#include "KimMoinFlow.h"
#include "elementutilitiescfd.h"

using namespace std;


BernsteinElem2DINSQuad9Node::BernsteinElem2DINSQuad9Node()
{
  ELEM_TYPE = 2;

  degree  = 2;
  npElem  = 9;

  nlbf    = 9;
  ndof    = 3;
  nsize   = nlbf*ndof;

  nGP = 9;
}


BernsteinElem2DINSQuad9Node::~BernsteinElem2DINSQuad9Node()
{

}


void BernsteinElem2DINSQuad9Node::prepareElemData(const vector<vector<double> >& node_coords)
{
    //printVector(nodeNums);

    // compute Volume and basis function derivatives

    double  dvol, Jac, param[2];

    int   ii, gp;

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    // Gauss point coordinates and weights
    vector<double>  gpts1(nGP), gpts2(nGP), gwts(nGP);

    gpts2[0] = -0.774596669241483;  gpts1[0] = -0.774596669241483;   gwts[0] = (5.0/9.0)*(5.0/9.0);
    gpts2[1] = -0.774596669241483;  gpts1[1] =  0.0;                 gwts[1] = (5.0/9.0)*(8.0/9.0);
    gpts2[2] = -0.774596669241483;  gpts1[2] =  0.774596669241483;   gwts[2] = (5.0/9.0)*(5.0/9.0);

    gpts2[3] =  0.0;                gpts1[3] = -0.774596669241483;   gwts[3] = (8.0/9.0)*(5.0/9.0);
    gpts2[4] =  0.0;                gpts1[4] =  0.0;                 gwts[4] = (8.0/9.0)*(8.0/9.0);
    gpts2[5] =  0.0;                gpts1[5] =  0.774596669241483;   gwts[5] = (8.0/9.0)*(5.0/9.0);

    gpts2[6] =  0.774596669241483;  gpts1[6] = -0.774596669241483;   gwts[6] = (5.0/9.0)*(5.0/9.0);
    gpts2[7] =  0.774596669241483;  gpts1[7] =  0.0;                 gwts[7] = (5.0/9.0)*(8.0/9.0);
    gpts2[8] =  0.774596669241483;  gpts1[8] =  0.774596669241483;   gwts[8] = (5.0/9.0)*(5.0/9.0);


    Nv.resize(nGP);
    dNvdx.resize(nGP);
    dNvdy.resize(nGP);

    Np.resize(nGP);
    dNpdx.resize(nGP);
    dNpdy.resize(nGP);

    for(int ii=0; ii<nGP; ii++)
    {
      Nv[ii].resize(nlbf);
      dNvdx[ii].resize(nlbf);
      dNvdy[ii].resize(nlbf);

      Np[ii].resize(nlbf);
      dNpdx[ii].resize(nlbf);
      dNpdy[ii].resize(nlbf);

      setZero(Nv[ii]);
      setZero(dNvdx[ii]);
      setZero(dNvdy[ii]);

      setZero(Np[ii]);
      setZero(dNpdx[ii]);
      setZero(dNpdy[ii]);
    }

    //printVector(nodeNums);

    elemVolGP.resize(nGP);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gpts1[gp];
      param[1] = gpts2[gp];

      computeBasisFunctions2D(false, ELEM_TYPE, degree, param, xNode, yNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], Jac);

      //cout << " Jac = " << Jac << endl;
      if(Jac < 0.0)
      {
        cout << " Negative Jacobian in 'BernsteinElem2DINSQuad9Node::prepareElemData' " << endl;
        cout << " Jac = " << Jac << endl;
        exit[1];
      }

      dvol = gwts[gp]*Jac;
      elemVolGP[gp] = dvol;

      elemVol += dvol;

      Np[gp][0] = 0.25*(1.0 - param[1])*(1.0 - param[0]);
      Np[gp][1] = 0.25*(1.0 - param[1])*(1.0 + param[0]);
      Np[gp][2] = 0.25*(1.0 + param[1])*(1.0 + param[0]);
      Np[gp][3] = 0.25*(1.0 + param[1])*(1.0 - param[0]);

    }//gp

    //charlen = sqrt(4.0*elemVol/PI);
    //charlen = 0.5*charlen;
    //charlen = charlen/3.0;

    charlen = (xNode[1]-xNode[0])*(xNode[1]-xNode[0])+(yNode[1]-yNode[0])*(yNode[1]-yNode[0]);
    charlen = min(charlen, (xNode[2]-xNode[1])*(xNode[2]-xNode[1])+(yNode[2]-yNode[1])*(yNode[2]-yNode[1]));
    charlen = min(charlen, (xNode[2]-xNode[3])*(xNode[2]-xNode[3])+(yNode[2]-yNode[3])*(yNode[2]-yNode[3]));
    charlen = min(charlen, (xNode[3]-xNode[0])*(xNode[3]-xNode[0])+(yNode[3]-yNode[0])*(yNode[3]-yNode[0]));

    charlen = 0.45*sqrt(charlen);

    //cout << " charlen = " << charlen << endl;

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return;
}


double  BernsteinElem2DINSQuad9Node::ResidualIncNavStokesAlgo1(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& veloVecPrev, const VectorXd& veloDotVec, const VectorXd& veloDotVecPrev, const VectorXd& presVec, const VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur)
{
    double  b1, b2, b3, b4, Da;
    double  velo[2], veloPrev[2], Du[2], dp[2], resi[2], rStab[2], gradTvel[2], pres;
    double  force[2], bforce[2], veloDot[2], divVelo;
    double  xx, yy, fact, fact2;
    double  grad[2][2], stress[2][2];

    int  gp, ii, jj, TI, TIp1;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    bforce[0] = elemData[2];
    bforce[1] = elemData[3];

    double  beta0  = elemData[5];
    double  betaSq = beta0*beta0;
    double  v_conv=-1.0, dtCric=1.0e10;
    double  beta=-1.0;
    double  v_diff = mu/rho/charlen;
    double  dtCric_visco = charlen*charlen*rho/mu/2.0;

    // time integration parameters
    //double  af   = timeData[1];
    //double  timefact = timeData[2];

    //double h2 = 4.0*elemVol/PI;
    //double tau = h2/(4.0*mu);

    //KimMoinFlow  analy(rho, mu, 0.0);

    //loop over Gauss points and compute element residual
    setZero(FlocalVelo);
    setZero(FlocalPres);
    force[0] = 0.0;    force[1] = 0.0;

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = 0.0;
          //veloPrev[0] = veloPrev[1] = 0.0;
          //veloDot[0] = veloDot[1] = 0.0;
          //Du[0]   = Du[1]   = 0.0;
          grad[0][0] = grad[0][1] = 0.0;
          grad[1][0] = grad[1][1] = 0.0;

          for(ii=0; ii<9; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            velo[0]    +=  b1*Nv[gp][ii];
            velo[1]    +=  b2*Nv[gp][ii];

            //b1 = 1.5*veloVec[TI]-0.5*veloVecPrev[TI];
            //b2 = 1.5*veloVec[TIp1]-0.5*veloVecPrev[TIp1];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[0][1] += b1*dNvdy[gp][ii];
            grad[1][0] += b2*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];

            //veloPrev[0] += veloVecPrev[TI]*Nv[gp][ii];
            //veloPrev[1] += veloVecPrev[TIp1]*Nv[gp][ii];

            //b1 = veloDotVec[TI];
            //b2 = veloDotVec[TIp1];

            //veloDot[0]    +=  b1*Nv[gp][ii];
            //veloDot[1]    +=  b2*Nv[gp][ii];
          }

          pres = 0.0;
          //dp[0] = dp[1] = 0.0;
          for(ii=0; ii<4; ii++)
          {
            //pres += (1.5*presVec[nodeNums[ii]]-0.5*presVecPrev[nodeNums[ii]])*Np[gp][ii];
            b1 = presVec[nodeNums[ii]];

            pres  += b1*Np[gp][ii];
            //dp[0] += b1*dNpdx[gp][ii];
            //dp[1] += b1*dNpdy[gp][ii];
          }

          //fact = 2.0*mu*(grad[0][0]+grad[1][1])/3.0;
          stress[0][0] = mu*grad[0][0] - pres;
          stress[0][1] = mu*grad[0][1];
          stress[1][0] = mu*grad[1][0];
          stress[1][1] = mu*grad[1][1] - pres;

          //force[0] = 0.0;
          //force[1] = 0.0;
          //force[0] = analy.computeForce(0, xx, yy, timeCur);
          //force[1] = analy.computeForce(1, xx, yy, timeCur);

          gradTvel[0] = velo[0]*grad[0][0] + velo[1]*grad[0][1];
          gradTvel[1] = velo[0]*grad[1][0] + velo[1]*grad[1][1];

          divVelo = grad[0][0]+grad[1][1];

          //gradTvel[0] = velo[0]*grad[0][0] + velo[1]*grad[0][1] + velo[0]*divVelo;
          //gradTvel[1] = velo[0]*grad[1][0] + velo[1]*grad[1][1] + velo[1]*divVelo;

          resi[0] = rho*(force[0] - gradTvel[0]) ;
          resi[1] = rho*(force[1] - gradTvel[1]) ;

          for(ii=0; ii<9; ii++)
          {
            TI   = ii*2;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*elemVolGP[gp];
            b2 = dNvdy[gp][ii]*elemVolGP[gp];
            b4 = Nv[gp][ii]*elemVolGP[gp];

            FlocalVelo[TI]   += (b4*resi[0] - b1*stress[0][0] - b2*stress[0][1]);
            FlocalVelo[TIp1] += (b4*resi[1] - b1*stress[1][0] - b2*stress[1][1]);
          }

          //v_conv = sqrt(veloPrev[0]*veloPrev[0]+veloPrev[1]*veloPrev[1]);
          //v_conv = sqrt(velo[0]*velo[0]+velo[1]*velo[1]);

          v_conv = max(v_conv, velo[0]*velo[0]+velo[1]*velo[1]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta = max(beta0, v_conv);
          //beta = beta0;
          //beta2 = beta*beta;

          fact = -elemVolGP[gp]*betaSq*divVelo;

          for(ii=0; ii<4; ii++)
          {
            FlocalPres[ii] += Np[gp][ii]*fact;
          }

          //dtCric = min(dtCric, charlen/(v_conv+sqrt(v_conv*v_conv+beta*beta)));
          //dtCric = min(dtCric, charlen/(v_conv+beta));
    }

    dtCric = charlen/(sqrt(v_conv)+sqrt(v_conv+betaSq));

    //return  min(dtCric, min(dtCric_visco, dtCric_pres));
    return  min(dtCric, dtCric_visco);
    //return  dtCric;
}



int  BernsteinElem2DINSQuad9Node::ResidualIncNavStokesAlgo2(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, VectorXd&  FlocalPres)
{
      double  velo[2], grad[2][2];
      double  xx, yy, fact, fact2, b1, b2, v_conv;

      int  gp, ii, jj, TI, TIp1;

      // material parameters
      double  rho = elemData[0];
      double  mu  = elemData[1];
      double  v_diff = mu/rho/charlen;

      double  beta2 = elemData[5]*elemData[5];
      double  beta=-1.0;
      double  beta0=elemData[5];

      // time integration parameters
      double  af   = timeData[1];
      double  timefact = timeData[2];

      //loop over Gauss points and compute element residual
      setZero(FlocalPres);

      for(gp=0; gp<nGP; gp++)
      {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = 0.0;
          grad[0][0] = grad[0][1] = grad[1][0] = grad[1][1] = 0.0;

          for(ii=0; ii<9; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];

            velo[0]    +=  b1*Nv[gp][ii];
            velo[1]    +=  b2*Nv[gp][ii];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];
          }

          v_conv = sqrt(velo[0]*velo[0]+velo[1]*velo[1]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta2 = beta*beta;

          fact = -elemVolGP[gp]*beta2*(grad[0][0]+grad[1][1]);

          for(ii=0; ii<4; ii++)
          {
            FlocalPres[ii] += Np[gp][ii]*fact;
          }
      }

    return 0;
}




double  BernsteinElem2DINSQuad9Node::calcCriticalTimeStep(const double* elemData, const double* timeData, const VectorXd&  veloVec)
{
    int ii, jj, gp, TI, TIp1;

    double  velo[2], b1, b2, v_conv, dtCric=1.0e10;
    double  beta=-1.0;
    double  beta0=elemData[5];
    double  rho=elemData[0];
    double  mu =elemData[1];
    double  v_diff = mu/rho/charlen;
    double  dtCric_visco = charlen*charlen*rho/mu/2.0;

    for(gp=0; gp<nGP; gp++)
    {
          velo[0] = velo[1] = 0.0;
          for(ii=0; ii<nlbf; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            velo[0]  +=  b1*Nv[gp][ii];
            velo[1]  +=  b2*Nv[gp][ii];
          }

          v_conv = sqrt(velo[0]*velo[0]+velo[1]*velo[1]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta = max(beta0, v_conv);
          beta = beta0;

          dtCric = min(dtCric, charlen/(v_conv+beta));
    }

    return  min(dtCric, dtCric_visco);
}


/*
int  BernsteinElem2DINSQuad9Node::StiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
{
    // For Laplace equation

    Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);
    //KimMoinFlow  analy(rho, mu, 1.0);

    int ii, jj, gp, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  b1, b2, b3, b4, b5, b6, b7, b8, val, force;
    double  dvol, pres, fact, fact1, fact2, bb1, bb2, param[2], dt, tCur;

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0; ii<npElem; ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    VectorXd  grad[2];

    double rho = elemData[0];
    double mu  = elemData[1];

    double  af = 1.0;//timeData[0];
    double  am = 1.0;//timeData[2];
    double  acceFact = rho*am*timeData[9];
    double  muTaf = mu*af;

    Klocal.setZero();
    Flocal.setZero();

    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP << endl;

    for(gp=0; gp<3; gp++)
    {
          dvol = elemVolGP[gp];

          //cout << " dvol = " << dvol << endl;
          //for(ii=0; ii<6; ii++)
            //cout << ii << '\t' << Nv[gp][ii] << '\t' << dNvdx[gp][ii] << '\t' << dNvdy[gp][ii] << endl;

          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          val = 0.0;
          grad.setZero();

          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            b1 = solnCur[nodeNums[ii]];

            val  +=  b1*Nv[gp][ii];

            grad[0] += b1*dNvdx[gp][ii];
            grad[1] += b1*dNvdy[gp][ii];
          }

          // this is pseudo-stress
          grad = mu*grad;

          force = 0.0;
          //force[0] = analy.computeForce(0, xx, yy, tCur);
          //cout << force[0] << '\t' << force[1] << endl;

          for(ii=0; ii<6; ii++)
          {
             b1 = dNvdx[gp][ii]*dvol;
             b2 = dNvdy[gp][ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;

             for(jj=0; jj<6; jj++)
             {
               Klocal(ii, jj) += ( b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj] );
             }

             Flocal[ii]  += (b4*force - b1*grad[0] - b2*grad[1] );
          }

    }//gp

    //printMatrix(Klocal); printf("\n\n\n");
    //printVector(Flocal); printf("\n\n\n");

    return 0;
}
*/


// Mass is assumed to be lumped.
// So, it is stored as a vector of diagonal vector
int BernsteinElem2DINSQuad9Node::MassMatrices(const vector<vector<double> >& node_coords, const double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2)
{
    int  ii;

    //cout << " VOLUME = " << dvol << endl;

    // Mass Matrix for the Velocity DOF
    // elemData[0] --> density
    double fact = elemData[0]*elemVol/9.0;

    for(ii=0; ii<9; ii++)
    {
      Mlocal1[ii] = fact;
    }

    // Mass Matrix for the Pressure DOF
    fact = elemVol/4.0;

    for(ii=0; ii<4; ii++)
    {
      Mlocal2[ii] = fact;
    }

    return 0;
}


//
int  BernsteinElem2DINSQuad9Node::StiffnessAndResidual(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
{
/*
    // Fully-implicit formulation

    //Stokes2DEx1 analy;
    Kovasznay  analy;
    analy.SetPressure(0.0);
    //KimMoinFlow  analy(rho, mu, 1.0);

    int ii, jj, gp, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;

    double  b1, b2, b3, b4, b5, b6, b7, b8, Da, Db;
    double  dvol, pres, fact, fact1, fact2, bb1, bb2, param[2], dt, tCur;

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0; ii<npElem; ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    double  res[3], vel[2], velDot[2], force[2], gradTvel[2];
    double  grad[2][2], gradN[2][2], stress[2][2];


    double rho = elemData[0];
    double mu  = elemData[1];

    double  af = 1.0;//timeData[0];
    double  am = 1.0;//timeData[2];
    double  acceFact = rho*am*timeData[9];
    double  muTaf = mu*af;

    setZero(Klocal);
    setZero(Flocal);

    //cout << " AAAAAAAAAA " << endl;
    //cout << nGP << endl;

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          vel[0] = vel[1] = 0.0;
          grad[0][0] = grad[0][1] = grad[1][0] = grad[1][1] = 0.0;
          pres = 0.0;
          for(ii=0; ii<npElem; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = solnCur[TI];
            b2 = solnCur[TIp1];
            b3 = solnCur[TIp2];

            vel[0]    +=  b1*Nv[gp][ii];
            vel[1]    +=  b2*Nv[gp][ii];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[0][1] += b1*dNvdy[gp][ii];
            grad[1][0] += b2*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];
          }

          for(ii=0; ii<4; ii++)
          {
            pres += solnCur[nodeNums[ii]*3+2]*Np[gp][ii];
          }

          velDot[0] = 0.0;
          velDot[1] = 0.0;

          // this is pseudo-stress
          stress[0][0] = mu*(grad[0][0]+grad[0][0]) - pres;
          stress[0][1] = mu*(grad[0][1]+grad[1][0]);
          stress[1][0] = stress[0][1];
          stress[1][1] = mu*(grad[1][1]+grad[1][1]) - pres;

          //cout << vel[0] << '\t' << vel[1] << endl;
          //cout << stress[0][0] << '\t' << stress[0][1] << '\t' << stress[1][0] << '\t' << stress[1][1] << endl;

          //force[0] = 0.0;
          //force[1] = 0.0;
          force[0] = analy.computeForce(0, xx, yy, tCur);
          force[1] = analy.computeForce(1, xx, yy, tCur);
          //cout << force[0] << '\t' << force[1] << endl;

          gradTvel = grad*vel;
          //gradTvel.setZero();

          res[0] = rho*(velDot[0] + gradTvel[0] - force[0]) ;
          res[1] = rho*(velDot[1] + gradTvel[1] - force[1]) ;
          res[2] = grad.trace();

          dvol = elemVolGP[gp];
          for(ii=0; ii<npElem; ii++)
          {
             TI   = 3*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;

             b1 = dNvdx[gp][ii]*dvol;
             b2 = dNvdy[gp][ii]*dvol;
             b4 = Nv[gp][ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b8 = Np[gp][ii]*dvol;

             for(jj=0; jj<npElem; jj++)
             {
               TJ   = 3*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*Nv[gp][jj];

               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj];

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               gradN = grad*(rho*Nv[gp][jj]);

               Db = rho*(vel[0]*dNvdx[gp][jj] + vel[1]*dNvdy[gp][jj]);

               gradN[0][0] += Db;
               gradN[1][1] += Db;

               Klocal(TI,   TJ)   += (b4*gradN[0][0]);
               Klocal(TI,   TJp1) += (b4*gradN[0][1]);
               Klocal(TIp1, TJ)   += (b4*gradN[1][0]);
               Klocal(TIp1, TJp1) += (b4*gradN[1][1]);

               // pressure term
               Klocal(TI,   TJp2) -= (b1*Np[gp][jj]);
               Klocal(TIp1, TJp2) -= (b2*Np[gp][jj]);

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dNvdx[gp][jj]);
               Klocal(TIp2, TJp1) -= (b8*dNvdy[gp][jj]);
               Klocal(TIp2, TJp2) += 0.0; //
             }

             Flocal(TI)   -= (b4*res[0] + b1*stress[0][0] + b2*stress[0][1] );
             Flocal(TIp1) -= (b4*res[1] + b1*stress[1][0] + b2*stress[1][1] );
             Flocal(TIp2) += (b8*res[2]);
          }

    }//gp

    //printVector(Flocal);
*/
    return 0;
}
//




int BernsteinElem2DINSQuad9Node::calcLoadVector(VectorXd& Flocal)
{
    return 0;
}



//
double  BernsteinElem2DINSQuad9Node::CalculateError(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, double timeCur, int index)
{
    double  b1, b2, b3, b4, fact, elemError;
    double  valNum[3], valExact[3], dp[2], gradNum[4], gradExact[4];

    int  nlbf=6, gp, ii, jj, TI, TIp1, TIp2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    //Stokes2DEx1 analy;
    //KimMoinFlow  analy(rho, mu, 0.0);
    Kovasznay  analy;
    //analy.SetPressure(0.0);

    elemError = 0.0;
    // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          valNum[0] = valNum[1] = valNum[2] = 0.0;
          gradNum[0] = gradNum[1] = gradNum[2] = gradNum[3] = 0.0;

          for(ii=0; ii<npElem; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];
            b3 = veloPrev[TIp2];

            valNum[0]  +=  b1*Nv[gp][ii];
            valNum[1]  +=  b2*Nv[gp][ii];
            valNum[2]  +=  b3*Np[gp][ii];

            gradNum[0] += b1*dNvdx[gp][ii];
            gradNum[2] += b1*dNvdy[gp][ii];
            gradNum[1] += b2*dNvdx[gp][ii];
            gradNum[3] += b2*dNvdy[gp][ii];
          }

          //valNum[2] = 0.0;
          //for(ii=0; ii<3; ii++)
            //valNum[2] += veloPrev[nodeNums[ii]*3+2]*Np[gp][ii];

          if(index < 3)
          {
            valExact[index] = analy.computeValue(index, xx, yy, timeCur);

            fact = valNum[index] - valExact[index];

            //if(index == 2)
              //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[2], valExact[2], fact);

            elemError += ( (fact*fact) * elemVolGP[gp] );
          }
          else
          {
            valExact[0] = analy.computeValue(0, xx, yy, timeCur);
            valExact[1] = analy.computeValue(1, xx, yy, timeCur);

            analy.computeDerivatives(xx, yy, gradExact);

            valNum[0] -= valExact[0];
            valNum[1] -= valExact[1];

            gradNum[0] -= gradExact[0];
            gradNum[1] -= gradExact[1];
            gradNum[2] -= gradExact[2];
            gradNum[3] -= gradExact[3];

            fact  = valNum[0]*valNum[0] + valNum[1]*valNum[1];
            fact += (gradNum[0]*gradNum[0]+gradNum[1]*gradNum[1]+gradNum[2]*gradNum[2]+gradNum[3]*gradNum[3]);

            elemError += ( fact * elemVolGP[gp] );
          }

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[index], valExact, fact);
    }

    return elemError;
}
//


/*
double  BernsteinElem2DINSQuad9Node::CalculateError(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, double timeCur, int index)
{
    double  b1, b2, b3, b4, fact, elemError, valExact[3];
    double  valNum[3], dp[2], gradNum[4], gradExact[4];

    int  gp, ii, jj, TI, TIp1, TIp2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
    }

    PoissonAnnulus analy;
    //PoissonEx3 analy;
    //Stokes2DEx1 analy;
    //KimMoinFlow  analy(rho, mu, 0.0);
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    elemError = 0.0;
    // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          valNum[0] = valNum[1] = valNum[2] = 0.0;
          gradNum[0] = gradNum[1] = gradNum[2] = gradNum[3] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            b1 = veloPrev[nodeNums[ii]];

            valNum[0]  +=  b1*Nv[gp][ii];

            gradNum[0] += b1*dNvdx[gp][ii];
            gradNum[1] += b1*dNvdy[gp][ii];
          }

          valExact[0] = analy.computeValue(index, xx, yy, timeCur);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[0], valExact[0], fact);

          valNum[0] -= valExact[0];

          fact = valNum[0]*valNum[0];

          if(index > 0)
          {
            analy.computeDerivatives(xx, yy, gradExact);

            //printf(" computed, exact \t %12.8f \t %12.8f \t %12.8f \t %12.8f  \t %12.8f \t %12.8f \n", valNum[0], valExact[0], gradNum[0], gradNum[1], gradExact[0], gradExact[1]);

            gradNum[0] -= gradExact[0];
            gradNum[1] -= gradExact[1];

            fact  += ( gradNum[0]*gradNum[0]+gradNum[1]*gradNum[1] );
          }

          elemError += ( fact * elemVolGP[gp] );
    }

    cout << " elemError = " << elemError << endl;

    return elemError;
}
*/



int BernsteinElem2DINSQuad9Node::toComputeInfSupCondition(const vector<vector<double> >& node_coords, const double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
/*
//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();

    // to compute the inf-sup constant

    int  sizep=9, ii, jj, kk, ll, gp, TI, TIp1, TJ, TJp1;
    double  bb1, bb2, bb3, bb4, nsize=18;

    double  fact, fact1, dvol;

    // resize local matrices and initialise them to zero
    if(Kuu.rows() != nsize)
    {
      Kuu.resize(nsize, nsize);
      Kup.resize(nsize, sizep);
      Kpp.resize(sizep, sizep);
    }
    Kuu.setZero();
    Kup.setZero();
    Kpp.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        dvol = elemVolGP[gp];

        // compute Kuu and Kup
        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dvol*dNvdx[gp][ii];
          bb2 = dvol*dNvdy[gp][ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          // compute Kuu
          for(jj=0; jj<nlbf; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            fact = bb1*dNvdx[gp][jj] + bb2*dNvdy[gp][jj] ;

            Kuu(TI,   TJ)   +=  fact ;
            Kuu(TIp1, TJp1) +=  fact ;
          }

          // compute Kup
          for(jj=0; jj<sizep; jj++)
          {
            Kup(TI,   jj) += bb1*Nv[gp][jj];
            Kup(TIp1, jj) += bb2*Nv[gp][jj];
          }
        }

        for(ii=0;ii<sizep;ii++)
        {
          bb3 = dvol*Nv[gp][ii];

          // compute Kup
          for(jj=0; jj<sizep; jj++)
          {
            Kpp(ii, jj) += bb3*Nv[gp][jj];
          }
        }
    }//gp
*/

    return 0;
}


