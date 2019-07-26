
/* Node ordering for the quadratic tetrahedral element
*/


#include "BernsteinElem3DINSTetra10Node.h"
#include "BasisFunctionsBernstein.h"
#include "headersBasic.h"
#include "KimMoinFlow.h"
#include "elementutilitiescfd.h"

using namespace std;


BernsteinElem3DINSTetra10Node::BernsteinElem3DINSTetra10Node()
{
  ELEM_TYPE = 1;

  degree  = 2;
  npElem  = 10;

  nlbf    = 10;
  ndof    = 3;
  nsize   = nlbf*ndof;

  nGP = 4;
}


BernsteinElem3DINSTetra10Node::~BernsteinElem3DINSTetra10Node()
{

}


void BernsteinElem3DINSTetra10Node::prepareElemData(vector<vector<double> >& node_coords)
{
    // compute Volume and basis function derivatives

    double  dvol, Jac, param[3];

    int   ii, gp;

    double xNode[npElem], yNode[npElem], zNode[npElem], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
      zNode[ii] = node_coords[nodeNums[ii]][2];
    }

    // Gauss point coordinates and weights
    vector<double>  gpts1(nGP), gpts2(nGP), gpts3(nGP), gwts(nGP);

    gpts1[0] = 0.1381966011250105;
    gpts1[1] = 0.5854101966249685;
    gpts1[2] = 0.1381966011250105;
    gpts1[3] = 0.1381966011250105;

    gpts2[0] = 0.1381966011250105;
    gpts2[1] = 0.1381966011250105;
    gpts2[2] = 0.5854101966249685;
    gpts2[3] = 0.1381966011250105;

    gpts3[0] = 0.1381966011250105;
    gpts3[1] = 0.1381966011250105;
    gpts3[2] = 0.1381966011250105;
    gpts3[3] = 0.5854101966249685;

    double  fact=1.0/24.0;
    gwts[0]  = fact;
    gwts[1]  = fact;
    gwts[2]  = fact;
    gwts[3]  = fact;


    Nv.resize(nGP);
    dNvdx.resize(nGP);
    dNvdy.resize(nGP);
    dNvdz.resize(nGP);

    Np.resize(nGP);
    dNpdx.resize(nGP);
    dNpdy.resize(nGP);
    dNpdz.resize(nGP);

    for(ii=0; ii<nGP; ii++)
    {
      Nv[ii].resize(nlbf);
      dNvdx[ii].resize(nlbf);
      dNvdy[ii].resize(nlbf);
      dNvdz[ii].resize(nlbf);

      Np[ii].resize(nlbf);
      dNpdx[ii].resize(nlbf);
      dNpdy[ii].resize(nlbf);
      dNpdz[ii].resize(nlbf);

      Nv[ii].setZero();
      dNvdx[ii].setZero();
      dNvdy[ii].setZero();
      dNvdz[ii].setZero();

      Np[ii].setZero();
      dNpdx[ii].setZero();
      dNpdy[ii].setZero();
      dNpdz[ii].setZero();
    }

    //printVector(nodeNums);

    elemVolGP.resize(nGP);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gpts1[gp];
      param[1] = gpts2[gp];
      param[2] = gpts3[gp];

      computeBasisFunctions3D(false, ELEM_TYPE, degree, param, xNode, yNode, zNode, &Nv[gp](0), &dNvdx[gp](0), &dNvdy[gp](0), &dNvdz[gp](0), Jac);

      //cout << " Jac = " << Jac << endl;
      if(Jac < 0.0)
      {
        cout << " Negative Jacobian in 'BernsteinEleme3DINSTetra10Node::prepareElemData' " << endl;
        cout << " Jac = " << Jac << endl;
        exit(1);
      }

      dvol = gwts[gp]*Jac;
      elemVolGP[gp] = dvol;

      elemVol += dvol;

      Np[gp][0] = 1.0-param[0]-param[1]-param[2];  Np[gp][1] = param[0];  Np[gp][2] = param[1];  Np[gp][3] = param[2];

      //computeBasisFunctions3D(false, ELEM_TYPE, 1, param, xNode, yNode, zNode, &Np[gp](0), &dNpdx[gp](0), &dNpdy[gp](0), &dNpdz[gp](0), Jac);
    }//gp

    charlen = pow(3.0*elemVol/PI/4.0, 1.0/3.0);
    charlen = 0.5*charlen;

    //charlen = (xNode[1]-xNode[0])*(xNode[1]-xNode[0])+(yNode[1]-yNode[0])*(yNode[1]-yNode[0]);
    //charlen = min(charlen, (xNode[2]-xNode[1])*(xNode[2]-xNode[1])+(yNode[2]-yNode[1])*(yNode[2]-yNode[1]));
    //charlen = min(charlen, (xNode[2]-xNode[0])*(xNode[2]-xNode[0])+(yNode[2]-yNode[0])*(yNode[2]-yNode[0]));

    //charlen = 0.5*sqrt(charlen);

    //cout << " charlen = " << charlen << endl;

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return;
}



double  BernsteinElem3DINSTetra10Node::calcCriticalTimeStep(double* elemData, double* timeData, VectorXd&  veloVec)
{
    int ii, jj, gp, TI, TIp1, TIp2;

    double  velo[3], b1, b2, b3, v_conv, dtCric=1.0e10;
    double  beta=-1.0;
    double  beta0=elemData[5];
    double  rho=elemData[0];
    double  mu =elemData[1];
    double  v_diff = mu/rho/charlen;
    double  dtCric_visco = charlen*charlen*rho/mu/3.0;

    for(gp=0; gp<nGP; gp++)
    {
          velo[0] = velo[1] = velo[2] = 0.0;
          for(ii=0; ii<10; ii++)
          {
            //xx += xNode[ii]*N(ii);
            //yy += yNode[ii]*N(ii);

            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];
            b3 = veloVec[TIp2];

            velo[0]  +=  b1*Nv[gp](ii);
            velo[1]  +=  b2*Nv[gp](ii);
            velo[2]  +=  b3*Nv[gp](ii);
          }

          v_conv = sqrt(velo[0]*velo[0] + velo[1]*velo[1] + velo[2]*velo[2]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta = max(beta0, v_conv);
          beta = beta0;

          dtCric = min(dtCric, charlen/(v_conv+beta));
    }

    return  min(dtCric, dtCric_visco);
}



double  BernsteinElem3DINSTetra10Node::ResidualIncNavStokesAlgo1(vector<vector<double> >& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& acceVec, VectorXd& acceVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur)
{
    double  b1, b2, b3, b4, Da;
    double  velo[3], Du[3], dp[3], resi[3], rStab[3], gradTvel[3], pres;
    double  force[3], acce[3];
    double  xx, yy, zz, fact, fact2;
    MatrixXd  grad(3,3), stress(3,3);

    int  gp, ii, jj, TI, TIp1, TIp2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];
    double  bforce[3] = {elemData[2], elemData[3], elemData[4]};
    double  beta0 = elemData[5];
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

    //loop over Gauss points and compute element residual
    Flocal1.setZero();
    Flocal2.setZero();

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = velo[2] = 0.0;
          acce[0] = acce[1] = acce[2] = 0.0;
          Du[0]   = Du[1]   = 0.0;
          grad.setZero();

          for(ii=0; ii<10; ii++)
          {
            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];
            b3 = veloVec[TIp2];

            velo[0]    +=  b1*Nv[gp](ii);
            velo[1]    +=  b2*Nv[gp](ii);
            velo[2]    +=  b3*Nv[gp](ii);

            grad(0,0) += b1*dNvdx[gp](ii);
            grad(0,1) += b1*dNvdy[gp](ii);
            grad(0,2) += b1*dNvdz[gp](ii);

            grad(1,0) += b2*dNvdx[gp](ii);
            grad(1,1) += b2*dNvdy[gp](ii);
            grad(1,2) += b2*dNvdz[gp](ii);

            grad(2,0) += b3*dNvdx[gp](ii);
            grad(2,1) += b3*dNvdy[gp](ii);
            grad(2,2) += b3*dNvdz[gp](ii);

            b1 = acceVec[TI];
            b2 = acceVec[TIp1];
            b3 = acceVec[TIp2];

            acce[0]    +=  b1*Nv[gp](ii);
            acce[1]    +=  b2*Nv[gp](ii);
            acce[2]    +=  b3*Nv[gp](ii);
          }

          pres = 0.0;
          dp[0] = dp[1] = dp[2] = 0.0;
          for(ii=0; ii<4; ii++)
          {
            b1 = presVec[nodeNums[ii]];

            pres  += b1*Np[gp](ii);
            dp[0] += b1*dNpdx[gp](ii);
            dp[1] += b1*dNpdy[gp](ii);
            dp[2] += b1*dNpdz[gp](ii);
          }

          //fact = 2.0*mu*(grad[0][0]+grad[1][1])/3.0;
          //stress[0][0] = mu*(grad[0][0]+grad[0][0]) - fact - pres;
          //stress[0][1] = mu*(grad[0][1]+grad[1][0]);
          //stress[1][0] = mu*(grad[1][0]+grad[0][1]);
          //stress[1][1] = mu*(grad[1][1]+grad[1][1]) - fact - pres;

          stress = mu*grad;

          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force[0] = 0.0;
          force[1] = 0.0;
          force[2] = 0.0;

          gradTvel[0] = grad(0,0)*velo[0] + grad(0,1)*velo[1] + grad(0,2)*velo[2];
          gradTvel[1] = grad(1,0)*velo[0] + grad(1,1)*velo[1] + grad(1,2)*velo[2];
          gradTvel[2] = grad(2,0)*velo[0] + grad(2,1)*velo[1] + grad(2,2)*velo[2];

          resi[0] = rho*(force[0] - gradTvel[0]) ;
          resi[1] = rho*(force[1] - gradTvel[1]) ;
          resi[2] = rho*(force[2] - gradTvel[2]) ;

          //rStab[0] = rho*(acce[0] + gradTvel[0]) - mu*Du[0] + dp[0] - force[0];
          //rStab[1] = rho*(acce[1] + gradTvel[1]) - mu*Du[1] + dp[1] - force[1];
          //rStab[2] = rho*(acce[2] + gradTvel[2]) - mu*Du[2] + dp[2] - force[2];

          //cout << " ddddddddd " << endl;

          for(ii=0; ii<10; ii++)
          {
            TI   = ii*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dNvdx[gp](ii)*elemVolGP[gp];
            b2 = dNvdy[gp](ii)*elemVolGP[gp];
            b3 = dNvdz[gp](ii)*elemVolGP[gp];
            b4 = Nv[gp](ii)*elemVolGP[gp];

            Flocal1(TI)   += (b4*resi[0] - b1*stress(0,0) - b2*stress(0,1) - b3*stress(0,2));
            Flocal1(TIp1) += (b4*resi[1] - b1*stress(1,0) - b2*stress(1,1) - b3*stress(1,2));
            Flocal1(TIp2) += (b4*resi[2] - b1*stress(2,0) - b2*stress(2,1) - b3*stress(2,2));
          }

          v_conv = max(v_conv, velo[0]*velo[0] + velo[1]*velo[1] + velo[2]*velo[2]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta2 = beta*beta;

          fact = -elemVolGP[gp]*betaSq*(grad(0,0) + grad(1,1) + grad(2,2));

          for(ii=0; ii<4; ii++)
          {
            Flocal2(ii) += Np[gp](ii)*fact;
          }

          //beta = max(beta0, max(v_conv, v_diff));
          //beta = max(beta0, v_conv);
          //beta = beta0;
    }

    dtCric = charlen/(sqrt(v_conv)+sqrt(v_conv+betaSq));

    return  min(dtCric, dtCric_visco);
}



int  BernsteinElem3DINSTetra10Node::ResidualIncNavStokesAlgo2(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& accePrev, VectorXd& presPrev, VectorXd&  Flocal2)
{
      double  velo[3];
      double  xx, yy, zz, fact, fact2, b1, b2, b3, v_conv;
      MatrixXd  grad(3,3);

      int  gp, ii, jj, TI, TIp1, TIp2;

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
      Flocal2.setZero();

      for(gp=0; gp<nGP; gp++)
      {
          // compute the gradient of displacement first
          xx = 0.0; yy = zz = 0.0;
          velo[0] = velo[1] = velo[2] = 0.0;
          grad(0,0) = grad(1,1) = grad(2,2) = 0.0;

          for(ii=0; ii<10; ii++)
          {
            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];
            b3 = veloPrev[TIp2];

            velo[0]    +=  b1*Nv[gp](ii);
            velo[1]    +=  b2*Nv[gp](ii);
            velo[2]    +=  b3*Nv[gp](ii);

            grad(0,0) += b1*dNvdx[gp](ii);
            grad(1,1) += b2*dNvdy[gp](ii);
            grad(2,2) += b2*dNvdy[gp](ii);
          }

          v_conv = sqrt(velo[0]*velo[0] + velo[1]*velo[1] + velo[2]*velo[2]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta2 = beta*beta;

          fact = -elemVolGP[gp]*beta2*(grad(0,0) + grad(1,1) + grad(2,2));

          for(ii=0; ii<4; ii++)
          {
            Flocal2(ii) += Np[gp](ii)*fact;
          }
      }

    return 0;
}




//
int  BernsteinElem3DINSTetra10Node::StiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
{
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

    VectorXd  res(3), vel(2), velDot(2), force(2), gradTvel(2);
    MatrixXd  grad(2,2), gradN(2,2), stress(2,2);


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

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          vel[0] = vel[1] = 0.0;
          grad.setZero();
          pres = 0.0;
          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp](ii);
            yy += yNode[ii]*Nv[gp](ii);

            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = solnCur[TI];
            b2 = solnCur[TIp1];
            b3 = solnCur[TIp2];

            vel[0]    +=  b1*Nv[gp](ii);
            vel[1]    +=  b2*Nv[gp](ii);

            grad(0,0) += b1*dNvdx[gp](ii);
            grad(0,1) += b1*dNvdy[gp](ii);
            grad(1,0) += b2*dNvdx[gp](ii);
            grad(1,1) += b2*dNvdy[gp](ii);
          }

          for(ii=0; ii<3; ii++)
          {
            pres += solnCur[nodeNums[ii]*3+2]*Np[gp](ii);
          }

          velDot(0) = 0.0;
          velDot(1) = 0.0;

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;

          //cout << vel(0) << '\t' << vel(1) << endl;
          //cout << stress(0,0) << '\t' << stress(0,1) << '\t' << stress(1,0) << '\t' << stress(1,1) << endl;

          //force[0] = 0.0;
          //force[1] = 0.0;
          force(0) = analy.computeForce(0, xx, yy, tCur);
          force(1) = analy.computeForce(1, xx, yy, tCur);
          //cout << force(0) << '\t' << force(1) << endl;

          gradTvel = grad*vel;
          //gradTvel.setZero();

          res(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
          res(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
          res(2) = grad.trace();

          dvol = elemVolGP[gp];
          for(ii=0; ii<6; ii++)
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

             for(jj=0; jj<6; jj++)
             {
               TJ   = 3*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               fact2 = rho*acceFact*Nv[gp](jj);

               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += b5*dNvdx[gp](jj)+b6*dNvdy[gp](jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // convection term

               gradN = grad*(rho*Nv[gp](jj));

               Db = rho*(vel(0)*dNvdx[gp](jj) + vel(1)*dNvdy[gp](jj));

               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Klocal(TI,   TJ)   += (b4*gradN(0,0));
               Klocal(TI,   TJp1) += (b4*gradN(0,1));
               Klocal(TIp1, TJ)   += (b4*gradN(1,0));
               Klocal(TIp1, TJp1) += (b4*gradN(1,1));

               // pressure term
               Klocal(TI,   TJp2) -= (b1*Np[gp](jj));
               Klocal(TIp1, TJp2) -= (b2*Np[gp](jj));

               // continuity equation
               Klocal(TIp2, TJ)   -= (b8*dNvdx[gp](jj));
               Klocal(TIp2, TJp1) -= (b8*dNvdy[gp](jj));
               Klocal(TIp2, TJp2) += 0.0; //
             }

             Flocal(TI)   -= (b4*res(0) + b1*stress(0,0) + b2*stress(0,1) );
             Flocal(TIp1) -= (b4*res(1) + b1*stress(1,0) + b2*stress(1,1) );
             Flocal(TIp2) += (b8*res(2));
          }

    }//gp

    //printVector(Flocal);

    return 0;
}
//



// Mass is assumed to be lumped.
// So, it is stored as a vector of diagonal vector
int BernsteinElem3DINSTetra10Node::MassMatrices(vector<vector<double> >& node_coords, double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2)
{
    int  ii;

    //cout << " VOLUME = " << dvol << endl;

    // Mass Matrix for the Velocity DOF
    // elemData[0] --> density
    double fact = elemData[0]*elemVol/10.0;

    for(ii=0; ii<10; ii++)
    {
      Mlocal1(ii) = fact;
    }

    // Mass Matrix for the Pressure DOF
    fact = elemVol/4.0;

    for(ii=0; ii<4; ii++)
    {
      Mlocal2(ii) = fact;
    }

    return 0;
}


int BernsteinElem3DINSTetra10Node::calcLoadVector(VectorXd& Flocal)
{
    return 0;
}



//
double  BernsteinElem3DINSTetra10Node::CalculateError(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& accePrev, VectorXd& presPrev, double timeCur, int index)
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
          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp](ii);
            yy += yNode[ii]*Nv[gp](ii);

            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];
            b3 = veloPrev[TIp2];

            valNum[0]  +=  b1*Nv[gp](ii);
            valNum[1]  +=  b2*Nv[gp](ii);
            valNum[2]  +=  b3*Np[gp](ii);

            gradNum[0] += b1*dNvdx[gp](ii);
            gradNum[2] += b1*dNvdy[gp](ii);
            gradNum[1] += b2*dNvdx[gp](ii);
            gradNum[3] += b2*dNvdy[gp](ii);
          }

          //valNum[2] = 0.0;
          //for(ii=0; ii<3; ii++)
            //valNum[2] += veloPrev[nodeNums[ii]*3+2]*Np[gp](ii);

          if(index < 3)
          {
            //valExact[index] = analy.computeValue(index, xx, yy, timeCur);
            if(index == 0)
              valExact[index] = 4.0*yy*(1.0-yy);
            else if(index == 1)
              valExact[index] = 0.0;
            else
              valExact[index] = 8.0*mu*(2.0-xx);


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
double  BernsteinElem3DINSTetra10Node::CalculateError(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& accePrev, VectorXd& presPrev, double timeCur, int index)
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
            xx += xNode[ii]*Nv[gp](ii);
            yy += yNode[ii]*Nv[gp](ii);

            b1 = veloPrev[nodeNums[ii]];

            valNum[0]  +=  b1*Nv[gp](ii);

            gradNum[0] += b1*dNvdx[gp](ii);
            gradNum[1] += b1*dNvdy[gp](ii);
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




int BernsteinElem3DINSTetra10Node::toComputeInfSupCondition(vector<vector<double> >& node_coords, double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
//
//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();

    // to compute the inf-sup constant

    int  sizep=3, ii, jj, kk, ll, gp, TI, TIp1, TJ, TJp1;
    double  bb1, bb2, bb3, bb4, nsize=12;

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
            Kup(TI,   jj) += bb1*Np[gp][jj];
            Kup(TIp1, jj) += bb2*Np[gp][jj];
          }
        }

        for(ii=0;ii<sizep;ii++)
        {
          bb3 = dvol*Np[gp][ii];

          // compute Kup
          for(jj=0; jj<sizep; jj++)
          {
            Kpp(ii, jj) += bb3*Np[gp][jj];
          }
        }
    }//gp

    return 0;
}


int  BernsteinElem3DINSTetra10Node::CalculateForces(int side, vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& presVec, VectorXd&  Flocal1)
{
    //cout << " CalculateForces " << endl;
    // faces for the 4-noded tetrahedron
    // face 1 - 2,4,3
    // face 2 - 3,4,1  (xi1=0)
    // face 3 - 1,4,2  (xi2=0)
    // face 4 - 1,2,3  (xi3=0)
    // 
    // faces for the 10-noded tetrahedron
    // face 1 - 2,4,3,9,10,6 ()
    // face 2 - 3,4,1,10,8,7 (xi1=0)
    // face 3 - 1,4,2,8,9,5  (xi2=0)
    // face 4 - 1,2,3,5,6,7  (xi3=0)

    int face_node_map[4][6]={{1,3,2,8,9,5},{2,3,0,9,7,6},{0,3,1,7,8,4},{0,1,2,4,5,6}};

    double  Jac, dvol, b1, b2, b3, b4, pres;
    double  grad[3][3], stress[3][3], trac[3], normal[3], tarpoint[3], Ne[6];
    double  xx, yy, zz, fact, fact2;
    double  xNode[10], yNode[10], zNode[10], xNodeFace[6], yNodeFace[6], zNodeFace[6], param[3];

    VectorXd  Nv(10), dNvdx(10), dNvdy(10), dNvdz(10), Np(4);
    Nv.setZero();    dNvdx.setZero();    dNvdy.setZero();    dNvdz.setZero();
    Np.setZero();

    int  gp, ii, jj, TI, TIp1, TIp2, degree=2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    for(ii=0; ii<10; ii++)
    {
      xNode[ii] = node_coords[nodeNums[ii]][0];
      yNode[ii] = node_coords[nodeNums[ii]][1];
      zNode[ii] = node_coords[nodeNums[ii]][2];
    }

    for(ii=0; ii<6; ii++)
    {
      xNodeFace[ii] = xNode[face_node_map[side][ii]];
      yNodeFace[ii] = yNode[face_node_map[side][ii]];
      zNodeFace[ii] = zNode[face_node_map[side][ii]];
    }


    vector<double>  gpts1(3), gpts2(3), gwts(3);
    gpts1[0] = 1.0/6.0;        gpts2[0] = 1.0/6.0;        gwts[0] = 1.0/6.0;
    gpts1[1] = 4.0/6.0;        gpts2[1] = 1.0/6.0;        gwts[1] = 1.0/6.0;
    gpts1[2] = 1.0/6.0;        gpts2[2] = 4.0/6.0;        gwts[2] = 1.0/6.0;

    //cout << " CalculateForces " << endl;

    //loop over Gauss points and compute element residual
    for(gp=0; gp<3; gp++)
    {
        // compute the normal and the Jacobian
        param[0] = gpts1[gp];
        param[1] = gpts2[gp];

        BernsteinBasisFunsFaceTria(degree, param, xNodeFace, yNodeFace, zNodeFace, Ne, normal, Jac);

        dvol = gwts[gp]*Jac;

        //compute the basis functions of calculating the stress

        tarpoint[0] = tarpoint[1] = tarpoint[2] = 0.0;
        for(ii=0; ii<6; ii++)
        {
          tarpoint[0] += xNodeFace[ii]*Ne[ii];
          tarpoint[1] += yNodeFace[ii]*Ne[ii];
          tarpoint[2] += zNodeFace[ii]*Ne[ii];
        }

        param[0] = 0.5;  param[1] = 0.5;  param[2] = 0.5;
        pointInverseTetra10node(xNode, yNode, zNode, tarpoint, param);

        computeBasisFunctions3D(false, ELEM_TYPE, degree, param, xNode, yNode, zNode, &Nv(0), &dNvdx(0), &dNvdy(0), &dNvdz(0), Jac);

        Np[0] = 1.0-param[0]-param[1]-param[2];  Np[1] = param[0];  Np[2] = param[1];  Np[3] = param[2];

        // calculate the velocity gradient

        grad[0][0] = grad[0][1] = grad[0][2] = 0.0;
        grad[1][0] = grad[1][1] = grad[1][2] = 0.0;
        grad[2][0] = grad[2][1] = grad[2][2] = 0.0;
        for(ii=0; ii<6; ii++)
        {
            TI   = nodeNums[ii]*3;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];
            b3 = veloVec[TIp2];

            grad[0][0] += b1*dNvdx[ii];
            grad[0][1] += b1*dNvdy[ii];
            grad[0][2] += b1*dNvdz[ii];

            grad[1][0] += b2*dNvdx[ii];
            grad[1][1] += b2*dNvdy[ii];
            grad[1][2] += b2*dNvdz[ii];

            grad[2][0] += b3*dNvdx[ii];
            grad[2][1] += b3*dNvdy[ii];
            grad[2][2] += b3*dNvdz[ii];
        }

        pres = 0.0;
        for(ii=0; ii<4; ii++)
        {
          pres += presVec[nodeNums[ii]]*Np[ii];
        }

        stress[0][0] = mu*grad[0][0] - pres;
        stress[0][1] = mu*grad[0][1] ;
        stress[0][2] = mu*grad[0][2] ;

        stress[1][0] = mu*grad[1][0] ;
        stress[1][1] = mu*grad[1][1] - pres;
        stress[1][2] = mu*grad[1][2] ;

        stress[2][0] = mu*grad[2][0] ;
        stress[2][1] = mu*grad[2][1] ;
        stress[2][2] = mu*grad[2][2] - pres;

        trac[0] = stress[0][0]*normal[0] + stress[0][1]*normal[1] + stress[0][2]*normal[2];
        trac[1] = stress[1][0]*normal[0] + stress[1][1]*normal[1] + stress[1][2]*normal[2];
        trac[2] = stress[2][0]*normal[0] + stress[2][1]*normal[1] + stress[2][2]*normal[2];

        for(ii=0; ii<6; ii++)
        {
          b4 = Ne[ii]*dvol;

          Flocal1(0) += (b4*trac[0]);
          Flocal1(1) += (b4*trac[1]);
          Flocal1(2) += (b4*trac[2]);
        }
    }

    //cout << " CalculateForces " << endl;

    return 0;
}


