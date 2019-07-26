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
#include "BernsteinElem2DINSQuad9Node.h"
#include "BernsteinElem3DINSTetra10Node.h"


using namespace std;


int main(int argc, char* argv[])
{
    double tstart, tend;

    int ndim=3, ndof=3;
    double fact, xNode[10], yNode[10], zNode[10];

    vector<vector<int> > ElemDofArray;

    int  nElem, nNode, nNode_Velo, nNode_Pres, nDBC_Velo, nDBC_Pres, nFBC;
    int  dd, ee, ii, jj, kk, ind, count, row, col;
    int  n1, n2, n3, n4, n5, n6, nn, nsize, dof;

    string  infileNodes, infileElems, infileDBCs, infileOutput;
    string  infileFBCs, charTemp, outFileName;

    int  npElem=10;
    int  npElemVelo=10;
    int  npElemPres=4;

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

       //if(argc == 5)
         //infileFBCs  = argv[3];

       if(argc == 5)
         infileOutput = argv[4];
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

    if(ndim == 2)
    {
      ii=0;
      while(infile_nodes >> val[0] >> val[1] >> val[2] )
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

        node_coords[ii][0] = val[1];
        node_coords[ii][1] = val[2];
        node_coords[ii][2] = 0.0;

        ii++;
      }
    }
    else
    {
      ii=0;
      while(infile_nodes >> val[0] >> val[1] >> val[2] >> val[3] )
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

        node_coords[ii][0] = val[1];
        node_coords[ii][1] = val[2];
        node_coords[ii][2] = val[3];

        ii++;
      }
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
    if(npElem == 6)
    {
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] >> val2[5] >> val2[6] >> val2[7] >> val2[8] >> val2[9] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=0; ii<npElem; ii++)
          elemNodeConn[ee][ii] = val2[4+ii]-1;

        ee++;
      }
    }
    else if(npElem == 9)
    {
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] >> val2[5] >> val2[6] >> val2[7] >> val2[8] >> val2[9] >> val2[10] >> val2[11] >> val2[12] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=0; ii<npElem; ii++)
          elemNodeConn[ee][ii] = val2[4+ii]-1;

        ee++;
      }
    }
    else if(npElem == 10)
    {
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] >> val2[5] >> val2[6] >> val2[7] >> val2[8] >> val2[9] >> val2[10] >> val2[11] >> val2[12] >> val2[13] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=0; ii<npElem; ii++)
          elemNodeConn[ee][ii] = val2[4+ii]-1;

        ee++;
      }
    }
    else
    {
      cerr << " Invalid npElem " << npElem << endl;
      exit(-1);
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

    vector<vector<double> >  DirichletBCsVelo, DirichletBCsPres;
    vector<double>  vecDblTemp(3);

    nDBC_Velo = 0;
    nDBC_Pres = 0;
    while(infile_DBCs >> val[0] >> val[1] >> val[2] )
    {
      vecDblTemp[0] = val[0]-1;
      vecDblTemp[1] = val[1]-1;
      vecDblTemp[2] = val[2];

      if( (int) val[1] > ndim )
      {
        DirichletBCsPres.push_back(vecDblTemp);
        nDBC_Pres++;
      }
      else
      {
        DirichletBCsVelo.push_back(vecDblTemp);
        nDBC_Velo++;
      }
    }

    cout << " nDBC_Velo      = " << '\t' << nDBC_Velo << endl;
    cout << " nDBC_Pres      = " << '\t' << nDBC_Pres << endl;

    infile_nodes.close();
    infile_elems.close();
    infile_DBCs.close();

    //
    // Read Output data
    //
    ////////////////////////////////////////////

    vector<int>  OutputData;

    if( !infileOutput.empty() )
    {
      std::ifstream  infile_output(infileOutput);

      cout << " reading output data " << endl;

      if(infile_output.fail())
      {
        cout << " Could not open the input elements file " << endl;
        exit(1);
      }

      ii = 0;
      while (std::getline(infile_output, line))
        ++ii;

      infile_output.clear();
      infile_output.seekg(0, infile_output.beg);

      OutputData.resize(ii);

      ee=0;
      while(infile_output >> val[0] )
      {
        OutputData[ee] = val[0]-1;
        ee++;
      }

      infile_output.close();
    }

    cout << " Input files have been read successfully \n\n " << endl;

    vector<vector<int> >  midNodeData;

    midNodeData.resize(nNode);

    for(ii=0; ii<nNode; ii++)
    {
      midNodeData[ii].resize(3);

      midNodeData[ii][0] = 0;  midNodeData[ii][1] = 0;  midNodeData[ii][2] = 0;
    }

    vector<int>  pressure_nodes;
    if(npElem == 6)
    {
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

        pressure_nodes.push_back(elemNodeConn[ee][0]);
        pressure_nodes.push_back(elemNodeConn[ee][1]);
        pressure_nodes.push_back(elemNodeConn[ee][2]);
      }
    }
    else if(npElem == 9)
    {
      for(ee=0; ee<nElem; ee++)
      {
        ii = elemNodeConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][0];
        midNodeData[ii][2] = elemNodeConn[ee][1];

        ii = elemNodeConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][1];
        midNodeData[ii][2] = elemNodeConn[ee][2];

        ii = elemNodeConn[ee][6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][2];
        midNodeData[ii][2] = elemNodeConn[ee][3];

        ii = elemNodeConn[ee][7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][3];
        midNodeData[ii][2] = elemNodeConn[ee][0];

        pressure_nodes.push_back(elemNodeConn[ee][0]);
        pressure_nodes.push_back(elemNodeConn[ee][1]);
        pressure_nodes.push_back(elemNodeConn[ee][2]);
        pressure_nodes.push_back(elemNodeConn[ee][3]);
      }
    }
    else if(npElem == 10)
    {
      for(ee=0; ee<nElem; ee++)
      {
        ii = elemNodeConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][0];
        midNodeData[ii][2] = elemNodeConn[ee][1];

        ii = elemNodeConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][1];
        midNodeData[ii][2] = elemNodeConn[ee][2];

        ii = elemNodeConn[ee][6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][0];
        midNodeData[ii][2] = elemNodeConn[ee][2];

        ii = elemNodeConn[ee][7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][3];
        midNodeData[ii][2] = elemNodeConn[ee][0];

        ii = elemNodeConn[ee][8];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][1];
        midNodeData[ii][2] = elemNodeConn[ee][3];

        ii = elemNodeConn[ee][9];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemNodeConn[ee][2];
        midNodeData[ii][2] = elemNodeConn[ee][3];

        pressure_nodes.push_back(elemNodeConn[ee][0]);
        pressure_nodes.push_back(elemNodeConn[ee][1]);
        pressure_nodes.push_back(elemNodeConn[ee][2]);
        pressure_nodes.push_back(elemNodeConn[ee][3]);
      }
    }

    findUnique(pressure_nodes);
    nNode_Pres = pressure_nodes.size();

    vector<int>  pressure_nodes_map(nNode,-1);

    for(ii=0; ii<nNode_Pres; ii++)
    {
      pressure_nodes_map[pressure_nodes[ii]] = ii;
    }






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

    // fix the pressure at all the mid nodes
    for(ii=0; ii<nNode; ii++)
    {
      if(midNodeData[ii][0])
        NodeTypePres[ii] = true;
    }

    if( (ndim == 2) && (npElem == 9) )
    {
      for(ee=0; ee<nElem; ee++)
      {
        NodeTypePres[elemNodeConn[ee][8]] = true;
      }
    }

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC_Velo; ii++)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      NodeTypeVelo[DirichletBCsVelo[ii][0]][DirichletBCsVelo[ii][1]] = true;
    }

    for(ii=0; ii<nDBC_Pres; ii++)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      NodeTypePres[DirichletBCsPres[ii][0]] = true;
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





      VectorXd  pres, presCur, presPrev, presPrev2, presPrev3;
      VectorXd  presDot, presDotPrev, presDotCur, presDiff;
      VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3, veloTemp;
      VectorXd  veloDot, veloDotPrev, acceCur;

      VectorXd  veloApplied, solnVTK(nNode*(ndim+1));
      VectorXd  globalMassVelo, globalMassPres;
      VectorXd  rhsVecVelo, rhsVecPres, Flocal1(npElem*ndim), Flocal2(npElem);

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
      veloApplied = velo;
      veloTemp = velo;

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

      ind = nNode_Velo*ndim;
      globalMassVelo.resize(ind);
      rhsVecVelo.resize(ind);

      globalMassPres.resize(nNode_Pres);
      rhsVecPres.resize(nNode_Pres);

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


      double  xx, yy, zz;
//
      // loop over the nodes and adjust nodal coordinates
      for(nn=0; nn<nNode; nn++)
      {
        if( midNodeData[nn][0] )
        {
          n1 = midNodeData[nn][1];
          n2 = midNodeData[nn][2];

          xx = 0.25*node_coords[n1][0] + 0.25*node_coords[n2][0];
          yy = 0.25*node_coords[n1][1] + 0.25*node_coords[n2][1];

          node_coords[nn][0] = 2.0*(node_coords[nn][0] - xx);
          node_coords[nn][1] = 2.0*(node_coords[nn][1] - yy);

          if(ndim == 3)
          {
            zz = 0.25*node_coords[n1][2] + 0.25*node_coords[n2][2];

            node_coords[nn][2] = 2.0*(node_coords[nn][2] - zz);
          }
        }
      }

     veloApplied.setZero();
     for(ii=0; ii<nDBC_Velo; ii++)
     {
          n1 = DirichletBCsVelo[ii][0];
          n2 = DirichletBCsVelo[ii][1];

          jj = n1*ndim+n2;

          veloApplied(jj) = DirichletBCsVelo[ii][2];
     }

     for(ii=0; ii<nDBC_Velo; ii++)
     {
        nn  = DirichletBCsVelo[ii][0];
        dof = DirichletBCsVelo[ii][1];

        if( midNodeData[nn][0] )
        {
          fact = 0.25*veloApplied(midNodeData[nn][1]*ndof+dof) + 0.25*veloApplied(midNodeData[nn][2]*ndof+dof);

          veloApplied[nn*ndof+dof] = 2.0*(veloApplied(nn*ndof+dof) - fact);
        }
     }
     //printVector(solnApplied);

      double  elemData[50], timeData[50];
      double  norm_velo=1000.0, norm_pres=1000.0;

      //time integration parameters
      timeData[1] = 1.0;   timeData[2] = 0.0;

      //density
      elemData[0] = 1.0;
      //viscosity
      elemData[1] = 0.01;
      //elemData[1] = 1.0/73.0;
      //elemData[1] = 1.0/229.0;
      //Body force in X-, Y- and Z- direction
      elemData[2] = 0.0;   elemData[3] = 0.0; elemData[4] = 0.0;
      //beta
      elemData[5] = 2.0;

      double  rhoInf;
      double  am = 1.0;
      double  gamm1 = 0.5+am;
      double  gamm2 = 0.5+am;
      //double  gamm2 = 1.0;

      int  stepsMax = 4000000;
      int  stepsCompleted=0;
      int  outputFreq = 10;
      int  fileCount=1;

      double  dt = 0.01;
      double  timeFact=0.0;
      double  timeNow=0.0;
      double  timeFinal = 500.0;
      double  num, denom;
      double  fact1, fact2;

      // create elements and prepare element data
      //BernsteinElem2DINSTria6Node   **elems;
      //BernsteinElem2DINSQuad9Node   **elems;
      BernsteinElem3DINSTetra10Node   **elems;

      //elems = new BernsteinElem2DINSTria6Node* [nElem];
      //elems = new BernsteinElem2DINSQuad9Node* [nElem];
      elems = new BernsteinElem3DINSTetra10Node* [nElem];

      for(ee=0;ee<nElem;ee++)
      {
        //elems[ee] = new BernsteinElem2DINSTria6Node;
        //elems[ee] = new BernsteinElem2DINSQuad9Node;
        elems[ee] = new BernsteinElem3DINSTetra10Node;

        elems[ee]->nodeNums = elemNodeConn[ee];

        elems[ee]->SolnData = &(SolnData);

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

          //printVector(Flocal1);
          //printVector(Flocal2);

          //Assemble the element mass to the global mass
          for(ii=0; ii<npElemVelo; ii++)
          {
            jj = elemNodeConn[ee][ii]*ndim ;

            for(dd=0; dd<ndim; dd++)
              globalMassVelo(jj+dd)   +=  Flocal1(ii);
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

      //for(ii=0; ii<nNode_Velo; ii++)
      //{
        //yy = node_coords[ii][1];
        //veloPrev(ii*2) = 2.0*yy*(3.0-yy)/3.0;
        //veloPrev(ii*2) = 1.0;
      //}
      //velo = veloPrev;

      cout << " Computing the critical time step \n" << endl;

      double dtCrit=1.0e10;

      //solnVTK.setZero();
      //writevtk(ndim, node_coords, elemNodeConn, midNodeData, solnVTK, fileCount);

      //Time loop
      while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )
      {
        dtCrit=1.0e10;
        for(ee=0; ee<nElem; ee++)
          dtCrit = min(dtCrit, elems[ee]->calcCriticalTimeStep(elemData, timeData, veloPrev));

        dt = dtCrit*0.9;
        cout << " dtCrit = " << dtCrit << endl;

        //
        if(stepsCompleted < 2000)
        {
          timeFact = 0.5*(1-cos(PI*stepsCompleted/2000));
        }
        else
        {
          timeFact = 1.0;
        }

        //timeFact = 1.0;

        //SolnData.applyDirichletBCs(timeFact);

        //Add specified Dirichlet BC
        for(ii=0; ii<nDBC_Velo; ii++)
        {
          n1 = DirichletBCsVelo[ii][0];
          n2 = DirichletBCsVelo[ii][1];

          jj = n1*ndim+n2;

          velo(jj) = veloApplied(jj)*timeFact;

          //xx = node_coords[n1][0];
          //yy = node_coords[n1][1];
          //if(n2 == 0)
            //velo(jj) = -cos(xx)*sin(yy)*sin(2.0*timeNow) ;
          //else
            //velo(jj) =  sin(xx)*cos(yy)*sin(2.0*timeNow) ;

          veloDot(jj) = (velo(jj)-veloPrev(jj))/(gamm1*dt) - (1.0-gamm1)*veloDotPrev(jj)/gamm1;
        }

        //printVector(velo);

        for(ii=0; ii<nDBC_Pres; ii++)
        {
          jj = DirichletBCsPres[ii][0];

          //cout << ii << '\t' << jj << '\t' << DirichletBCsPres[ii][2] << endl;

          pres(n1) = DirichletBCsPres[ii][2]*timeFact;

          presDot(jj) = (pres(jj)-presPrev(jj))/(gamm2*dt) - (1.0-gamm2)*presDotPrev(jj)/gamm2;
        }

        //cout << " aaaaaaaaaaa " << endl;
        //Loop over elements and compute the RHS

        rhsVecVelo.setZero();
        rhsVecPres.setZero();

        for(ee=0; ee<nElem; ee++)
        {
            fact = timeNow - dt;
            //Compute the element force vector, including residual force
            elems[ee]->ResidualIncNavStokesAlgo1(node_coords, elemData, timeData, velo, veloPrev, veloDot, veloDotPrev, pres, presPrev, Flocal1, Flocal2, fact);

            //printVector(Flocal1);
            //printVector(Flocal2);

            //Assemble the element vector
            for(ii=0; ii<npElemVelo; ii++)
            {
              jj = ndim*ii;
              kk = ndim*elemNodeConn[ee][ii];

              for(dd=0; dd<ndim; dd++)
                rhsVecVelo(kk+dd) += Flocal1(jj+dd);
            }

            for(ii=0; ii<npElemPres; ii++)
            {
              rhsVecPres(elemNodeConn[ee][ii])   += Flocal2(ii);
            }
        } //LoopElem

        //cout << " cccccccccccc " << endl;

        // Add specified nodal force 
        //SolnData.addNodalForces(timeFact);

        //Add contribution from the Mass matrix terms

        for(ii=0; ii<totalDOFvelo; ii++)
        {
          jj = assyForSolnVelo[ii];

          veloDot(jj) = rhsVecVelo(jj)/(am*globalMassVelo(jj)) - (1.0-am)*veloDotPrev(jj)/am;

          //veloDot(jj) = rhsVecVelo(jj)/globalMassVelo(jj);

          velo(jj) = veloPrev(jj) + dt*(gamm1*veloDot(jj)+(1.0-gamm1)*veloDotPrev(jj));

          if( isnan(abs(velo(jj))) )
          {
            cerr << " NAN encountered in velocity ... " << endl;
            cerr << " Program has been terminated ... " << endl;
            cout << jj << '\t' << globalMassVelo(jj) << '\t' << rhsVecVelo(jj) << endl;
            exit(-1);
          }
        }

        //cout << " bbbbbbbbbbb " << endl;

        //veloTemp = 1.5*velo - 0.5*veloPrev;
        //veloPrev3 = 0.5*(velo+veloPrev)+dt*veloDot;
        //veloTemp = 0.5*(velo+veloPrev);
        //veloTemp = (23.0/12.0)*velo - (4.0/3.0)*veloPrev + (5.0/12.0)*veloPrev2;

        /*
        gamm2 = 1.0;
        rhsVecPres.setZero();
        for(ee=0; ee<nElem; ee++)
        {
            fact = timeNow - dt;
            //Compute the element force vector, including residual force
            elems[ee]->ResidualIncNavStokesTria6nodeAlgo2(node_coords, elemData, timeData, veloTemp, veloDot, pres, Flocal2);

            for(ii=0; ii<npElemPres; ii++)
            {
              rhsVecPres(elemNodeConn[ee][ii])   += Flocal2(ii);
            }
        } //LoopElem
        */

        //cout << " aaaaaaaaaaa " << endl;

        for(ii=0; ii<totalDOFpres; ii++)
        {
            jj = assyForSolnPres[ii];

            if( midNodeData[jj][0] == 0 )
            {
              presDot(jj) = rhsVecPres(jj)/(am*globalMassPres(jj)) - (1.0-am)*presDotPrev(jj)/am;

              pres(jj)    = presPrev(jj) + dt*(gamm2*presDot(jj)+(1.0-gamm2)*presDotPrev(jj));

              if( isnan(abs(pres(jj))) )
              {
                cerr << " NAN encountered in pressure ... " << endl;
                cerr << " Program has been terminated ... " << endl;
                exit(-1);
              }
            }
        }

        //cout << " cccccccccccccc " << endl;

        if(stepsCompleted%outputFreq == 0)
        {
          solnVTK.resize(nNode_Velo*(ndim+1));
          solnVTK.setZero();
          for(ii=0; ii<nNode; ii++)
          {
             n1 = ii*(ndim+1);
             n2 = ii*ndim;

             for(dd=0; dd<ndim; dd++)
               solnVTK(n1+dd) = velo(n2+dd);

             solnVTK(n1+ndim) = pres(ii);
          }


          //cout << " ddddddddddddddd " << endl;

          // write the solution to the VTK file
          writevtk(ndim, node_coords, elemNodeConn, midNodeData, solnVTK, fileCount);

          fileCount = fileCount+1;

          cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;

          veloDiff = velo - veloPrev;
          presDiff = pres - presPrev;

          norm_velo = 1.0;
          norm_pres = 1.0;
          if(stepsCompleted > 0)
          {
            //norm_velo = veloDiff.norm();
            //norm_pres = presDiff.norm();

            norm_velo = veloDiff.norm()/velo.norm();
            norm_pres = presDiff.norm()/pres.norm();
          }
          //cout << " velocity difference norm = " << '\t' << norm_velo << endl;

          fact1 = 0.0;
          fact2 = 0.0;
          for(ii=0; ii<OutputData.size(); ii++)
          {
            jj = OutputData[ii]*2;

            fact1 += rhsVecVelo(jj);
            fact2 += rhsVecVelo(jj+1);
          }

          fout << timeNow << '\t' << stepsCompleted << '\t' << norm_velo << '\t' << norm_pres ;
          fout << '\t' << fact1 << '\t' << fact2 << endl;
        }

        if(norm_velo < 1.0e-6)
        {
          cout << " Solution converged below the specified tolerance " << endl;
          break;
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

    solnVTK.resize(nNode_Velo*(ndim+1));
    solnVTK.setZero();
    for(ii=0; ii<nNode; ii++)
    {
      n1 = ii*3;
      n2 = ii*2;

      solnVTK(n1)   = velo(n2);
      solnVTK(n1+1) = velo(n2+1);
      solnVTK(n1+2) = pres(ii);
    }

    double totalError = 0.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<3; index++)
    {
      totalError = 0.0;
      for(ee=0; ee<nElem; ee++)
      {
        //Compute the element force vector, including residual force
        totalError += elems[ee]->CalculateError(node_coords, elemData, timeData, solnVTK, veloDot, pres, timeNow, index);
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




/*
        //
        if(stepsCompleted > 1000)
        //if(norm_velo < 1.0e-4)
        {
          //
          for(ii=0; ii<totalDOFvelo; ii++)
          {
            jj = assyForSolnVelo[ii];

            num = velo(jj) - veloPrev(jj);
            num = num*num;

            denom = (velo(jj) - veloPrev(jj)) - (veloPrev(jj) - veloPrev2(jj));

            if(abs(denom) > 1.0e-4)
            {
              //cout << " denom = " << denom << endl;
              velo(jj) = velo(jj) - num/denom;
            }
          }
          //
          //
          for(ii=0; ii<totalDOFpres; ii++)
          {
            jj = assyForSolnPres[ii];

            if( midNodeData[jj][0] == 0 )
            {
              num = pres(jj) - presPrev(jj);
              num = num*num;

              denom = (pres(jj) - presPrev(jj)) - (presPrev(jj) - presPrev2(jj));
 
              if(abs(denom) > 1.0e-4)
                pres(jj) = pres(jj) - num/denom;
            }
          }
          //
        }
        //
*/



        /*
        if(norm_velo < 1.0e-3)
        {
          veloTemp = velo;
          for(int aa=0; aa<2; aa++)
          {
            cout << " Aitken Accelerator " << endl;
          //
          for(ii=0; ii<totalDOFvelo; ii++)
          {
            jj = assyForSolnVelo[ii];

            num = veloTemp(jj) - veloPrev(jj);
            num = num*num;

            denom = (veloTemp(jj) - veloPrev(jj)) - (veloPrev(jj) - veloPrev2(jj));

            if(abs(denom) > 1.0e-10)
            {
              //cout << " denom = " << denom << endl;
              veloTemp(jj) = veloTemp(jj) - num/denom;
            }
          }

          //printVector(veloPrev);
          //printVector(velo);
          //for(ii=0; ii<nNode_Velo; ii++)
            //cout << ii << '\t' << velo(ii*2) << '\t' << velo(ii*2+1) << '\t' << pres(ii) << endl;
          }
        }
        else
        {
        */