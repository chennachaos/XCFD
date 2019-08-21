
#include "ExplicitCFD.h"
#include "elementutilitiescfd.h"
#include "ElementBase.h"
#include <chrono>



void ExplicitCFD::setSolver(int slv, int *parm, bool cIO)
{
    prepareMatrixPattern();

    calcMassMatrixForExplicitDynamics();

    return;
}




int ExplicitCFD::prepareMatrixPattern()
{
    cout <<  "\n     ExplicitCFD::prepareMatrixPattern()  .... STARTED ...\n" <<  endl;

    totalDOF = nNode*ndof;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    vector<vector<int> >  IDvelo;
    vector<vector<int> >  IEN, IEN2, ID, LM;
    vector<vector<bool> >  NodeTypeVelo;
    vector<bool>  NodeTypePres(nNode,  false);
    vector<int>   IDpres(nNode,  -1);


    NodeTypeVelo.resize(nNode);
    IDvelo.resize(nNode);

    int  ee, ii, jj, kk;

    for(ii=0;ii<nNode;++ii)
    {
      NodeTypeVelo[ii].resize(ndim);
      IDvelo[ii].resize(ndim);

      for(jj=0;jj<ndim;++jj)
      {
        NodeTypeVelo[ii][jj] = false;
        IDvelo[ii][jj] = -1;
      }
    }

    // fix the pressure at all the mid nodes
    for(ii=0; ii<nNode; ++ii)
    {
      if(midNodeData[ii][0])
        NodeTypePres[ii] = true;
    }

    if( (ndim == 2) && (npElem == 9) )
    {
      for(ee=0; ee<nElem; ++ee)
      {
        NodeTypePres[elemConn[ee][8]] = true;
      }
    }

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC_Velo; ++ii)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      NodeTypeVelo[DirichletBCsVelo[ii][0]][DirichletBCsVelo[ii][1]] = true;
    }

    for(ii=0; ii<nDBC_Pres; ++ii)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      NodeTypePres[DirichletBCsPres[ii][0]] = true;
    }


    for(ii=0;ii<nNode;++ii)
    {
      for(jj=0;jj<ndim;++jj)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(!NodeTypeVelo[ii][jj])
        {
          IDvelo[ii][jj] = totalDOF_Velo++;
          assyForSolnVelo.push_back(ii*ndim+jj);
        }
      }

      if(!NodeTypePres[ii])
      {
        IDpres[ii] = totalDOF_Pres++;
        assyForSolnPres.push_back(ii);
      }
    }

    cout << " Mesh statistics .....\n" << endl;
    cout << " nElem          = " << '\t' << nElem << endl;
    cout << " nNode          = " << '\t' << nNode  << endl;
    cout << " npElem         = " << '\t' << npElem << endl;
    cout << " ndof           = " << '\t' << ndof << endl;
    cout << " Velocity DOF   = " << '\t' << totalDOF_Velo << endl;
    cout << " Pressure DOF   = " << '\t' << totalDOF_Pres << endl;


    double  xx, yy, zz;
    int nn, dof;

    /*
    LM.resize(nElem);

    for(ee=0;ee<nElem;++ee)
    {
      npElem = IEN[ee].size();

      //printVector(IEN[ee]);

      ind = ndof*npElem;
      LM[ee].resize(ind);

      for(ii=0;ii<npElem;++ii)
      {
        ind = ndof*ii;

        kk = IEN[ee][ii];

        for(jj=0;jj<ndof;++jj)
        {
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << ID[kk][jj] << endl;
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << LM[ee][ind+jj] << '\t' << ID[IEN[ee][ii]][jj] << endl;
          LM[ee][ind+jj] = ID[kk][jj];
          //cout << " IIIIIIIII " << endl;
        }
      }
    }

    assyForSoln.resize(totalDOF);

    count = 0;
    for(ii=0;ii<nNode;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        //cout << ii << '\t' << jj << '\t' << ID[ii][jj] << endl;
        if( ID[ii][jj] != -1)
          assyForSoln[count++] = ii*ndof + jj;
      }
    }

    //for(ii=0; ii<totalDOF; ++ii)
      //cout << ii << '\t' << assyForSoln[ii] << endl;

    //cout << " totalDOF " << totalDOF << endl;
    for(ee=0; ee<nElem; ++ee)
    {
      elems[ee].forAssyVec = LM[ee];
    }
    */

    bool pp=false;
    //pp=true;
    if(pp)
    {
       printf("   IEN array \n\n");
       for(ii=0;ii<nElem;++ii)
       {
          for(jj=0;jj<npElem;++jj)
            cout << '\t' << IEN[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   ID array \n\n");
       for(ii=0;ii<nNode;++ii)
       {
          for(jj=0;jj<ndof;++jj)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   LM array \n\n");
       for(ii=0;ii<nElem;++ii)
       {
          for(jj=0;jj<LM[ii].size();++jj)
            cout << '\t' << LM[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assyForSoln array \n\n");
       for(ii=0;ii<totalDOF;++ii)
       {
          cout << assyForSolnVelo[ii] << endl;
       }
       printf("\n\n\n");
    }

    printf("\n element DOF values initialised \n\n");

    // matrix pattern needs to be prepared only for the implicit solver

    // remove data objects

    //IEN.clear();
    //LM.clear();
    //ID.clear();
    //forAssyMat.clear();

    printf("\n     ExplicitCFD::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}





void ExplicitCFD::applyBoundaryConditions(double timeFact)
{
    //cout <<  " applying boundary conditions .... " << endl;

    int ii, jj, n1, n2;
    double  fact1 = 1.0/(gamm1*dt), fact2 = (gamm1-1.0)/gamm1;

#pragma acc parallel loop   
    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        n1 = DirichletBCsVelo[ii][0];
        n2 = DirichletBCsVelo[ii][1];

        jj = n1*ndim+n2;

        velo[jj]    = veloApplied[jj] * timeFact;
        veloDot[jj] = fact1 * (velo[jj]-veloPrev[jj]) + fact2 * veloDotPrev[jj];
    }

#pragma acc parallel loop   
    for(ii=0; ii<nDBC_Pres; ++ii)
    {
        jj = DirichletBCsPres[ii][0];

        pres[jj]    = DirichletBCsPres[ii][2] * timeFact;
        presDot[jj] = fact1 * (pres[jj]-presPrev[jj]) + fact2 * presDotPrev[jj];
    }

    return;
}



void ExplicitCFD::addExternalForces(double loadFact)
{
    int  nn, dof, ii, ind;
    double specVal=0.0;

    VectorXd  vecTemp, Flocal;
    setZero(vecTemp);

    // specified nodal forces
    for(ii=0;ii<nodeForcesData.size();++ii)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      ind = nn*ndof+dof;

      vecTemp[ind] += specVal;
    }
    //printVector(vecTemp);

    return;
}




int ExplicitCFD::calcMassMatrixForExplicitDynamics()
{
    cout << " ExplicitCFD::calcMassMatrixForExplicitDynamics() STARTED " << endl;

    int  dd, ee, ii, jj;

    VectorXd  FlocalVelo(npElem*ndim), FlocalPres(npElem);

    // Compute global mass matrices
    //The Mass is assumed to be lumped so that the mass matrix is diagonal

    setZero(globalMassVelo);
    setZero(globalMassPres);

    for(ee=0; ee<nElem; ++ee)
    {
        // compute mass matrix and assemble it
        elems[ee].MassMatrices(node_coords, elemData, FlocalVelo, FlocalPres);

        //Assemble the element mass to the global mass
        for(ii=0; ii<npElemVelo; ++ii)
        {
            jj = elemConn[ee][ii]*ndim ;

            for(dd=0; dd<ndim; dd++)
              globalMassVelo[jj+dd]   +=  FlocalVelo[ii];
        }

        for(ii=0; ii<npElemPres; ++ii)
        {
            globalMassPres[elemConn[ee][ii]] += FlocalPres[ii];
        }
    }

    //printVector(globalMassVelo);
    //printVector(globalMassPres);

    cout << " ExplicitCFD::calcMassMatrixForExplicitDynamics() FINISHED " << endl;

    return 0;
}




int  ExplicitCFD::solveExplicitStepDTS()
{
/*
    calcMassMatrixForExplicitDynamics();

    cout << " Solving the ExplicitCFD " << endl;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    velo.setZero();      veloPrev.setZero();      veloCur.setZero();
    veloDot.setZero();   veloDotPrev.setZero();   veloDotCur.setZero();
    pres.setZero();      presPrev.setZero();      presCur.setZero();
    presDot.setZero();   presDotPrev.setZero();   presDotCur.setZero();

    double dtCrit=1.0e10, fact, fact1, fact2;
    int pseduostepsCount = 0;
    int pseduostepsMax   = 1000, n1, n2, ee, ii, jj, kk, ind;
    int  stepsCompleted=0;
    double dtau = 0.01;
    double  norm_velo = 1.0;
    double  norm_pres = 1.0;
    double  timeNow=0.0;
    double  timeFact=0.0;

    VectorXd  FlocalVelo(npElem*ndof), FlocalPres(npElem*ndof), TotalForce(3);

    setInitialConditions();

    //Time loop
    while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )
    {
        if(stepsCompleted < 2000)
        {
          timeFact = 0.5*(1-cos(PI*stepsCompleted/2000));
        }
        else
        {
          timeFact = 1.0;
        }

        timeFact = 1.0;
        applyBoundaryConditions(timeFact);
        //printVector(velo);

        //dtCrit=1.0e10;
        //for(ee=0; ee<nElem; ++ee)
          //dtCrit = min(dtCrit, elems[ee].calcCriticalTimeStep(elemData, timeData, veloPrev));
        //dt = dtCrit*0.9;
        //cout << " dtCrit = " << dtCrit << endl;

        //cout << " aaaaaaaaaaa " << endl;

        veloM   = velo;
        veloMp1 = velo;
        veloDotM   = veloDot;
        veloDotMp1 = veloDot;

        presM   = pres;
        presMp1 = pres;
        presDotM   = presDot;
        presDotMp1 = presDot;

        fact1 = 1.0/gamm1/dt;
        fact2 = (1.0-gamm1)/gamm1;

        rhsVecVelo2.setZero();
        ind = nNode_Velo*ndim;
        for(ii=0; ii<ind; ++ii)
        {
          rhsVecVelo2[ii] = globalMassVelo[ii]*(fact1*velo[ii] - fact2*veloDot[ii]);
        }

        rhsVecPres2.setZero();
        for(ii=0; ii<nNode_Pres; ++ii)
        {
          rhsVecPres2[ii] = globalMassPres[ii]*(fact1*pres[ii] - fact2*presDot[ii]);
        }

        pseduostepsCount = 0;
        pseduostepsMax   = 1000;
        dtau = 0.005;
        while(pseduostepsCount < pseduostepsMax)
        {
          //Add specified Dirichlet BC
          for(ii=0; ii<nDBC_Velo; ++ii)
          {
            n1 = DirichletBCsVelo[ii][0];
            n2 = DirichletBCsVelo[ii][1];

            jj = n1*ndim+n2;

            veloMp1[jj] = veloApplied[jj]*timeFact;

            veloDotMp1[jj] = (veloMp1[jj]-veloM[jj])/(gamm1*dt) - (1.0-gamm1)*veloDotM[jj]/gamm1;
          }

          //Loop over elements and compute the RHS

          rhsVecVelo.setZero();
          rhsVecPres.setZero();

          for(ee=0; ee<nElem; ++ee)
          {
            fact = timeNow - dt;
            //Compute the element force vector, including residual force
            dtCrit = elems[ee].ResidualIncNavStokesAlgo1(node_coords, elemData, timeData, veloM, veloPrev, veloDotM, veloDotPrev, presM, presPrev, FlocalVelo, FlocalPres, fact);

            //Assemble the element vector
            for(ii=0; ii<npElemVelo; ++ii)
            {
              jj = ndim*ii;
              kk = ndim*elemConn[ee][ii];

              rhsVecVelo[kk]   += FlocalVelo[jj];
              rhsVecVelo[kk+1] += FlocalVelo[jj+1];
            }

            for(ii=0; ii<npElemPres; ++ii)
            {
              rhsVecPres[elemConn[ee][ii]]   += FlocalPres[ii];
            }
          } //LoopElem

          // Add specified nodal force 
          //SolnData.addNodalForces(timeFact);

          //Add contribution from the Mass matrix terms

          for(ii=0; ii<totalDOF_Velo; ++ii)
          {
            jj = assyForSolnVelo[ii];

            veloDotMp1[jj] = (rhsVecVelo2[jj] + rhsVecVelo[jj] - globalMassVelo[jj]*veloM[jj]/gamm1/dt)/globalMassVelo[jj];

            veloMp1[jj] = veloM[jj] + dtau*(gamm1*veloDotMp1[jj]+(1.0-gamm1)*veloDotM[jj]);

            if( std::isnan(abs(veloMp1[jj])) )
            {
              cerr << " NAN encountered in velocity ... " << endl;
              cerr << " Program has been terminated ... " << endl;
              cout << jj << '\t' << globalMassVelo[jj] << '\t' << rhsVecVelo[jj] << endl;
              exit(-1);
            }
          }

          for(ii=0; ii<totalDOF_Pres; ++ii)
          {
            jj = assyForSolnPres[ii];

            presDotMp1[jj] = (rhsVecPres2[jj] + rhsVecPres[jj] - globalMassPres[jj]*presM[jj]/gamm2/dt)/globalMassPres[jj];

            presMp1[jj]    = presM[jj] + dtau*(gamm2*presDotMp1[jj]+(1.0-gamm2)*presDotM[jj]);

            if( std::isnan(abs(presMp1[jj])) )
            {
                cerr << " NAN encountered in pressure ... " << endl;
                cerr << " Program has been terminated ... " << endl;
                exit(-1);
            }
          }

          //cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;

          veloDiff = veloMp1 - veloM;
          presDiff = presMp1 - presM;

          norm_velo = veloDiff.norm();//veloMp1.norm();
          norm_pres = presDiff.norm();//presMp1.norm();

          if(norm_velo < 1.0e-6)
          {
            cout << " pseduostepsCompleted = " << pseduostepsCount << '\t' << norm_velo << endl;

            break;
          }

          veloM    = veloMp1;
          veloDotM    = veloDotMp1;

          presM    = presMp1;
          presDotM = presDotMp1;

          pseduostepsCount++;
        } // while pseduosteps

        velo    = veloMp1;
        veloDot    = veloDotMp1;

        pres    = presMp1;
        presDot = presDotMp1;

        //veloDiff = velo - veloPrev;
        //presDiff = pres - presPrev;

        //norm_velo = veloDiff.norm()/velo.norm();
        //norm_pres = presDiff.norm()/pres.norm();

        //if(norm_velo < 1.0e-6)
        //{
          //cout << " Solution fout_convdataged below the specified tolerance " << endl;
          //break;
        //}

        if(stepsCompleted%outputFreq == 0)
        {
          cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;

          cout << " velocity difference norm = " << '\t' << norm_velo << endl;

          postProcess();

          setZero(TotalForce);
          for(ii=0; ii<outputEdges.size(); ++ii)
          {
            elems[outputEdges[ii][0]]->CalculateForces(outputEdges[ii][1], node_coords, elemData, timeData, velo, pres, TotalForce);
          }

          fact1 = TotalForce[0];  fact2 = TotalForce[1];

          fout_convdata << timeNow << '\t' << stepsCompleted << '\t' << norm_velo << '\t' << norm_pres ;
          fout_convdata << '\t' << fact1 << '\t' << fact2 << endl;
        }

        // store the variables
        //SolnData.timeUpdate();
        veloPrev3  = veloPrev2;
        veloPrev2  = veloPrev;
        veloPrev   = velo;

        veloDotPrev  = veloDot;

        presPrev3  = presPrev2;
        presPrev2  = presPrev;
        presPrev   = pres;
        presDotPrev  = presDot;

        stepsCompleted = stepsCompleted + 1;
        timeNow = timeNow + dt;

      } //Time loop
*/

    return 0;
}


void ExplicitCFD::computeElementErrors(int ind)
{
    cout << " Computing errors \n " << endl;

    VectorXd  solnVTK(nNode_Velo*(ndim+1));
    int  n1, n2, dd;
    for(int ii=0; ii<nNode; ++ii)
    {
      n1 = ii*ndof;
      n2 = ii*ndim;

      for(dd=0; dd<ndim; dd++)
        solnVTK[n1+dd] = velo[n2+dd];

      solnVTK[n1+ndim] = pres[ii];
    }

    double totalError = 0.0, timeNow;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(int ee=0; ee<nElem; ++ee)
      {
        totalError += elems[ee].CalculateError(node_coords, elemData, timeData, solnVTK, veloDot, pres, timeNow, index);
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error in velocity   = %12.6E \n\n " , totalError);
    }

    return;
}



