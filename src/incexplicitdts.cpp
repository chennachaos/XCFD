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
#include "elementutilitiescfd.h"
#include "SolutionData.h"
#include "BernsteinElem2DINSTria6Node.h"

using namespace std;



int main(int argc, char* argv[])
{
    double tstart, tend;

    int ndim=2, ndof=2, npElem=6;
    double fact, xNode[6], yNode[6];

    vector<vector<int> > ElemDofArray;

    int forAssyVec(12), elemDofGlobal(12);

    int  nElem, nNode, nNode_Velo, nNode_Pres, nDBC, nFBC;
    int  ee, ii, jj, kk, ind, count, row, col;
    int  n1, n2, n3, n4, n5, n6, nsize;

    string  infileNodes, infileElems, infileDBCs;
    string  infileFBCs, charTemp, outFileName;

    int  npElemVelo=6;
    int  npElemPres=3;

    //Set file names
    //The file names are specified as inputs from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr <<  "Number of input files is not sufficient " << endl;
        cerr <<  "You must enter names of THREE files" << endl;
        cerr <<  "a.) Node file, b.) Element file, and c.) Dirichlet BC file" << endl;
        cerr << "Aborting..." << endl;
    }
    else
    {
       infileNodes = argv[1];
       infileElems = argv[2];
       infileDBCs  = argv[3];

        if(argc == 5)
           infileFBCs  = argv[3];
    }


    // Read nodal data files
    /////////////////////////////////////

    std::ifstream  infile_nodes(infileNodes);
    std::ifstream  infile_elems(infileElems);
    std::ifstream  infile_DBCs(infileDBCs);
    std::ifstream  infile_FBCs(infileFBCs);


    if(infile_nodes.fail())
    {
       cout << " Could not open the input nodes file " << endl;
       exit(1);
    }

    double  val[10];
    int  val2[10];

    std::string line;

    // read nodal coordinates
    ////////////////////////////////////////////

    cout << " reading nodes " << endl;

    nNode = 0;
    while (std::getline(infile_nodes, line))
      ++nNode;

    vector<vector<double> >  node_coords(nNode, vector<double>(3));

    infile_nodes.clear();
    infile_nodes.seekg(0, infile_nodes.beg);

      ii=0;
      while(infile_nodes >> val[0] >> val[1] >> val[2] )
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

        node_coords[ii][0] = val[1];
        node_coords[ii][1] = val[2];
        node_coords[ii][2] = 0.0;

        ii++;
      }


    // read elements
    ////////////////////////////////////////////
    cout << " reading elements " << endl;

    if(infile_elems.fail())
    {
       cout << " Could not open the input elements file " << endl;
       exit(1);
    }


    nElem = 0;
    while (std::getline(infile_elems, line))
      ++nElem;

    cout << " nElem   " << nElem << endl;

    infile_elems.clear();
    infile_elems.seekg(0, infile_elems.beg);

    vector<vector<int> >  elemNodeConn(nElem, vector<int>(npElem));

      ee=0;
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] >> val2[5] >> val2[6] >> val2[7] >> val2[8] >> val2[9] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=0; ii<6; ii++)
          elemNodeConn[ee][ii] = val2[4+ii]-1;

        ee++;
      }

    //
    // Read Dirichlet BC data
    //
    ////////////////////////////////////////////
    cout << " reading DBCs " << endl;

    if(infile_DBCs.fail())
    {
       cout << " Could not open the input elements file " << endl;
       exit(1);
    }


    nDBC = 0;
    while (std::getline(infile_DBCs, line))
      ++nDBC;

    cout << " nDBC         = " << '\t' << nDBC << endl;

    infile_DBCs.clear();
    infile_DBCs.seekg(0, infile_DBCs.beg);

    vector<vector<double> >  DirichletBCs(nDBC, vector<double>(3));

    ee=0;
    while(infile_DBCs >> val[0] >> val[1] >> val[2] )
    {
        DirichletBCs[ee][0] = val[0]-1;
        DirichletBCs[ee][1] = val[1]-1;
        DirichletBCs[ee][2] = val[2];

      ee++;
    }


    infile_nodes.close();
    infile_elems.close();
    infile_DBCs.close();

    cout << " Input files have been read successfully \n\n " << endl;


    vector<vector<int> >  IDvelo;
    vector<vector<bool> >  NodeTypeVelo;
    vector<int>  assyForSolnVelo, assyForSolnPres, IDpres;
    vector<bool>  NodeTypePres;

    NodeTypeVelo.resize(nNode);
    NodeTypePres.resize(nNode);
    IDvelo.resize(nNode);
    IDpres.resize(nNode);
    for(ii=0;ii<nNode;ii++)
    {
      NodeTypeVelo[ii].resize(ndof);
      IDvelo[ii].resize(ndof);

      for(jj=0;jj<ndof;jj++)
      {
        NodeTypeVelo[ii][jj] = false;
        IDvelo[ii][jj] = -1;
      }

      NodeTypePres[ii] = false;
      IDpres[ii] = -1;
    }


    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      if(DirichletBCs[ii][1] == 2)
        NodeTypePres[DirichletBCs[ii][0]] = true;
      else
        NodeTypeVelo[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    int totalDOFvelo = 0;
    int totalDOFpres = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(!NodeTypeVelo[ii][jj])
        {
          IDvelo[ii][jj] = totalDOFvelo++;
          assyForSolnVelo.push_back(ii*ndof+jj);
        }
      }

      if(!NodeTypePres[ii])
      {
        IDpres[ii] = totalDOFpres++;
        assyForSolnPres.push_back(ii);
      }
    }


      cout << " Mesh statistics .....\n" << endl;
      cout << " nElem          = " << '\t' << nElem << endl;
      cout << " nNode          = " << '\t' << nNode  << endl;
      cout << " npElem         = " << '\t' << npElem << endl;
      cout << " ndof           = " << '\t' << ndof << endl;
      cout << " Velocity DOF   = " << '\t' << totalDOFvelo << endl;
      cout << " Pressure DOF   = " << '\t' << totalDOFpres << endl;

      vector<vector<int> >  midNodeData;

      midNodeData.resize(nNode);

      for(ii=0; ii<nNode; ii++)
      {
        midNodeData[ii].resize(3);

        midNodeData[ii][0] = 0;  midNodeData[ii][1] = 0;  midNodeData[ii][2] = 0;
      }

      for(ee=0; ee<nElem; ee++)
      {
        ii = elemNodeConn[ee][3];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][0];
        midNodeData[ii][2] = elemNodeConn[ee][1];

        ii = elemNodeConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][1];
        midNodeData[ii][2] = elemNodeConn[ee][2];

        ii = elemNodeConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][2];
        midNodeData[ii][2] = elemNodeConn[ee][0];
      }

      VectorXd  pres, presCur, presPrev, presPrev2, presPrev3;
      VectorXd  presDot, presDotPrev, presDotCur, presDiff;
      VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3;
      VectorXd  veloDot, veloDotPrev, acceCur;

      VectorXd  veloM, veloMp1, acceM, acceMp1;
      VectorXd  presM, presMp1, presDotM, presDotMp1;


      VectorXd  solnApplied, solnVTK;
      VectorXd  globalMassVelo, globalMassPres;
      VectorXd  rhsVecVelo, rhsVecPres, rhsVecVelo2, rhsVecPres2;
      VectorXd  Flocal1(12), Flocal2(6);

      SolutionData  SolnData;

      SolnData.ndim       = ndim;
      SolnData.nNode_Velo = nNode_Velo;
      SolnData.nNode_Pres = nNode_Pres;

      nNode_Velo = nNode;
      nNode_Pres = nNode;

      ind = nNode_Velo*ndim;

      velo.resize(ind);
      velo.setZero();

      veloPrev  = velo;
      veloCur   = velo;
      veloPrev2 = velo;
      veloPrev3 = velo;

      veloDot     = velo;
      veloDotPrev = veloDot;
      acceCur  = veloDot;

      pres.resize(nNode_Pres);
      pres.setZero();

      presPrev  = pres;
      presPrev2 = pres;
      presPrev3 = pres;
      presCur   = pres;

      presDot     = velo;
      presDotPrev = pres;
      presDotCur  = pres;


      veloM    = velo;
      veloMp1  = velo;
      acceM    = velo;
      acceMp1  = velo;

      presM      = pres;
      presMp1    = pres;
      presDotM   = pres;
      presDotMp1 = pres;


      ind = nNode_Velo*2;
      globalMassVelo.resize(ind);
      rhsVecVelo.resize(ind);
      rhsVecVelo2.resize(ind);

      globalMassPres.resize(nNode_Pres);
      rhsVecPres.resize(nNode_Pres);
      rhsVecPres2.resize(nNode_Pres);

      ////////////////////////////////////////////////
      // Computations start from here
      ////////////////////////////////////////////////

    char fname[200];
    //sprintf(fname,"(),"-explicit.dat");

    ofstream fout("convergence-data.dat");

    if(fout.fail())
    {
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(14);


      double  elemData[50], timeData[50];
      double  norm_velo=1000.0, norm_pres=1000.0;

      //time integration parameters
      timeData[1] = 1.0;   timeData[2] = 0.0;

      //density
      elemData[0] = 1.0;
      //viscosity
      elemData[1] = 0.01;
      //Body force in X-, Y- and Z- direction
      elemData[2] = 0.0;   elemData[3] = 0.0; elemData[4] = 0.0;
      //beta
      elemData[5] = 2.0;

      double  rhoInf;
      double  am = 1.0;
      double  gamm = 0.5+am;

      int  stepsMax = 4000000;
      int  stepsCompleted=0;
      int  outputFreq = 10;
      int  fileCount=1;

      int  pseduostepsCount= 0;
      int  pseduostepsMax= 100;

      double  dt = 0.01;
      double  dtau = 0.1;
      double  timeFact=0.0;
      double  timeNow=0.0;
      double  timeFinal = 4000.0;
      double  num, denom;
      double  fact1, fact2;

      // create elements and prepare element data
      BernsteinElem2DINSTria6Node   **elems;

      elems = new BernsteinElem2DINSTria6Node* [nElem];

      for(ee=0;ee<nElem;ee++)
      {
        elems[ee] = new BernsteinElem2DINSTria6Node;

        elems[ee]->nodeNums = elemNodeConn[ee];

        //elems[ee]->SolnData = &(SolnData);

        elems[ee]->prepareElemData(node_coords);
      }

      cout << " elements are created and prepated " << endl;

      // Compute global mass matrices
      //The Mass is assumed to be lumped so that the mass matrix is diagonal

      globalMassVelo.setZero();
      globalMassPres.setZero();

      for(ee=0; ee<nElem; ee++)
      {
          // compute mass matrix and assemble it
          elems[ee]->MassMatrices(node_coords, elemData, Flocal1, Flocal2);

          //Assemble the element mass to the global mass
          for(ii=0; ii<npElemVelo; ii++)
          {
            jj = elemNodeConn[ee][ii]*2 ;

            globalMassVelo(jj)   +=  Flocal1(ii);
            globalMassVelo(jj+1) +=  Flocal1(ii);
          }

          for(ii=0; ii<npElemPres; ii++)
          {
            globalMassPres(elemNodeConn[ee][ii]) += Flocal2(ii);
          }
      }

      //printVector(globalMassVelo);
      //printVector(globalMassPres);

      cout << " Computing the solution \n" << endl;

      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////

      velo.setZero();      veloPrev.setZero();      veloCur.setZero();
      veloDot.setZero();      veloDotPrev.setZero();      acceCur.setZero();
      pres.setZero();      presPrev.setZero();      presCur.setZero();
      presDot.setZero();   presDotPrev.setZero();   presDotCur.setZero();

      //Time loop
      while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )
      {
        /*
        if(stepsCompleted > 1000)
        //if(norm_velo < 1.0e-4)
        {
          for(ii=0; ii<totalDOFvelo; ii++)
          {
            jj = assyForSolnVelo[ii];

            num = velo(jj) - veloPrev(jj);
            num = num*num;

            denom = (velo(jj) - veloPrev(jj)) - (veloPrev(jj) - veloPrev2(jj));

            if(abs(denom) > 1.0e-10)
            {
              //cout << " denom = " << denom << endl;
              velo(jj) = velo(jj) - num/denom;
            }
          }

          //
          for(ii=0; ii<totalDOFpres; ii++)
          {
            jj = assyForSolnPres[ii];

            num = pres(jj) - presPrev(jj);
            num = num*num;

            denom = (pres(jj) - presPrev(jj)) - (presPrev(jj) - presPrev2(jj));

            if(abs(denom) > 1.0e-10)
              pres(jj) = pres(jj) - num/denom;
          }
          //
        }
        */

        if(stepsCompleted < 1000)
        {
          timeFact = 0.5*(1-cos(PI*stepsCompleted/1000));
        }
        else
        {
          timeFact = 1.0;
        }

        //timeFact = 1.0;

        //SolnData.applyDirichletBCs(timeFact);

        //veloCur = af*velo + (1.0-af)*veloPrev
        //acceCur = am*veloDot + (1.0-am)*veloDotPrev
        //presCur = af*pres + (1.0-af)*presPrev

        veloM   = velo;
        veloMp1 = velo;
        acceM   = veloDot;
        acceMp1 = veloDot;

        presM   = pres;
        presMp1 = pres;
        presDotM   = presDot;
        presDotMp1 = presDot;

        fact1 = 1.0/gamm/dt;
        fact2 = (1.0-gamm)/gamm;

        rhsVecVelo2.setZero();
        ind = nNode_Velo*ndim;
        for(ii=0; ii<ind; ii++)
        {
          rhsVecVelo2(ii) = globalMassVelo(ii)*(fact1*velo(ii) - fact2*veloDot(ii));
        }

        rhsVecPres2.setZero();
        for(ii=0; ii<nNode_Pres; ii++)
        {
          rhsVecPres2(ii) = globalMassPres(ii)*(fact1*pres(ii) - fact2*presDot(ii));
        }

        pseduostepsCount = 0;
        pseduostepsMax   = 10;
        dtau = 0.01;
        while(pseduostepsCount < pseduostepsMax)
        {
          //Add specified Dirichlet BC
          for(ii=0; ii<nDBC; ii++)
          {
            n1 = DirichletBCs[ii][0];
            n2 = DirichletBCs[ii][1];

            if( n2 < 2)
            {
              jj = n1*2+n2;

              veloMp1(jj) = DirichletBCs[ii][2]*timeFact;

              //xx = node_coords[n1][0];
              //yy = node_coords[n1][1];
              //if(n2 == 0)
                //velo(jj) = -cos(xx)*sin(yy)*sin(2.0*timeNow) ;
              //else
                //velo(jj) =  sin(xx)*cos(yy)*sin(2.0*timeNow) ;

              acceMp1(jj) = (veloMp1(jj)-veloM(jj))/(gamm*dt) - (1.0-gamm)*acceM(jj)/gamm;
            }
          }

          //Loop over elements and compute the RHS

          rhsVecVelo.setZero();
          rhsVecPres.setZero();

          for(ee=0; ee<nElem; ee++)
          {
            fact = timeNow - dt;
            //Compute the element force vector, including residual force
            //elems[ee]->ResidualIncNavStokesTria6nodeAlgo1(node_coords, elemData, timeData, veloMp1, acceMp1, presMp1, Flocal1, Flocal2, fact);

            //Assemble the element vector
            for(ii=0; ii<npElemVelo; ii++)
            {
              jj = 2*ii;
              kk = 2*elemNodeConn[ee][ii];

              rhsVecVelo(kk)   += Flocal1(jj);
              rhsVecVelo(kk+1) += Flocal1(jj+1);
            }

            for(ii=0; ii<npElemPres; ii++)
            {
              rhsVecPres(elemNodeConn[ee][ii])   += Flocal2(ii);
            }
          } //LoopElem

          // Add specified nodal force 
          //SolnData.addNodalForces(timeFact);

          //Add contribution from the Mass matrix terms

          for(ii=0; ii<totalDOFvelo; ii++)
          {
            jj = assyForSolnVelo[ii];

            acceMp1(jj) = (rhsVecVelo2(jj) + rhsVecVelo(jj) - globalMassVelo(jj)*veloM(jj)/gamm/dt)/globalMassVelo(jj);

            veloMp1(jj) = veloM(jj) + dtau*(gamm*acceMp1(jj)+(1.0-gamm)*acceM(jj));

            if( isnan(abs(veloMp1(jj))) )
            {
              cerr << " NAN encountered in velocity ... " << endl;
              cerr << " Program has been terminated ... " << endl;
              //cout << jj << '\t' << globalMassVelo(jj) << '\t' << rhsVecVelo(jj) << endl;
              exit(-1);
            }
          }

          for(ii=0; ii<totalDOFpres; ii++)
          {
            jj = assyForSolnPres[ii];

            if( midNodeData[jj][0] == 0 )
            {
              presDotMp1(jj) = (rhsVecPres2(jj) + rhsVecPres(jj) - globalMassPres(jj)*presM(jj)/gamm/dt)/globalMassPres(jj);

              presMp1(jj)    = presM(jj) + dtau*(gamm*presDotMp1(jj)+(1.0-gamm)*presDotM(jj));

              if( isnan(abs(presMp1(jj))) )
              {
                cerr << " NAN encountered in pressure ... " << endl;
                cerr << " Program has been terminated ... " << endl;
                exit(-1);
              }
            }
          }

          //cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;

          veloDiff = veloMp1 - veloM;
          presDiff = presMp1 - presM;

          norm_velo = veloDiff.norm();//veloMp1.norm();
          norm_pres = presDiff.norm();//presMp1.norm();

          if(norm_velo < 1.0e-5)
          {
            cout << " pseduostep loop velocity difference norm = " << '\t' << norm_velo << endl;
            break;
          }

          veloM    = veloMp1;
          acceM    = acceMp1;

          presM    = presMp1;
          presDotM = presDotMp1;

          pseduostepsCount++;
        } // while pseduosteps

        velo    = veloMp1;
        veloDot    = acceMp1;

        pres    = presMp1;
        presDot = presDotMp1;


        if(stepsCompleted%outputFreq == 0)
        {
          // write the solution to the VTK file
          //writevtk(true, node_coords, elemNodeConn, midNodeData, velo, pres, fileCount);

          fileCount = fileCount+1;

          cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;

          veloDiff = velo - veloPrev;
          presDiff = pres - presPrev;

          norm_velo = veloDiff.norm()/velo.norm();
          norm_pres = presDiff.norm()/pres.norm();

          //cout << " velocity difference norm = " << '\t' << norm_velo << endl;

          fout << timeNow << '\t' << stepsCompleted << '\t' << norm_velo << '\t' << norm_pres << endl;
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

    /*
    cout << " Computing errors \n " << endl;
    double totalError = 0.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<3; index++)
    {
        for(ee=0; ee<nElem; ee++)
        {
            veloElem.setZero(); acceElem.setZero(); presElem.setZero();

            for(ii=0; ii<npElem; ii++)
            {
              n1 = elemNodeConn[ee][ii];

              xNode[ii] = node_coords[n1][0];
              yNode[ii] = node_coords[n1][1];

              jj = ii*ndof;
              kk = n1*ndof;

              veloElem(jj)   = velo(kk);
              veloElem(jj+1) = velo(kk+1);

              presElem(ii)   = pres(n1);
            }

            //Compute the element force vector, including residual force
            totalError += CalculateErrorTria6node(EQUAL_ORDER, xNode, yNode, elemData, timeData, veloElem, acceElem, presElem, timeNow, index);
      }

      totalError = sqrt(totalError);

      if(ii == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(ii == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
    }
    */

    //tend = MPI_Wtime()
    //write(*,*) "That took ", (tend-tstart), "seconds"

    if(elems != NULL)
    {
      for(ii=0;ii<nElem;ii++)
        delete elems[ii];

      delete [] elems;
      elems = NULL;
    }


    cout << " Program is successful \n " << endl;

    return 1;
}







