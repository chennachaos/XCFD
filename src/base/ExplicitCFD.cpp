
#include <algorithm>
#include <chrono>
#include "ExplicitCFD.h"
#include "KimMoinFlow.h"
#include "elementutilitiescfd.h"
#include "BernsteinElem2DINSTria6Node.h"
#include "BernsteinElem2DINSQuad9Node.h"
#include "BernsteinElem3DINSTetra10Node.h"
#include "omp.h"

using namespace std;


ExplicitCFD::ExplicitCFD()
{
    ndof = 0; nElem = 0; nNode = 0; npElem = 0; fileCount = 0;
    totalDOF_Velo = 0; totalDOF_Pres = 0; totalDOF = 0;

    //char convfilename[200];
    //sprintf(convfilename,"%s%s%06d%s", "convergence-data-", infilename.c_str(),".dat");

    fout_convdata.open("convergence-data.dat", ios::out | ios::trunc );

    if(fout_convdata.fail())
    {
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout_convdata.setf(ios::fixed);
    fout_convdata.setf(ios::showpoint);
    fout_convdata.precision(14);

}


ExplicitCFD::~ExplicitCFD()
{
}




void  ExplicitCFD::readInputData(string&  fname)
{
    cout << " Reading input data \n\n " << endl;

    infilename = fname;

    std::ifstream  infile( string(infilename+".dat") );

    if(infile.fail())
    {
       cout << " Could not open the input nodes file " << endl;
       exit(1);
    }

    string  line, stringVal, stringVec[10];
    int  ii, arrayInt[100];
    double  tempDbl;

    // read the dimension
    infile >> stringVal >> ndim;
    ndof = ndim+1;

    // read number of nodes per element
    infile >> stringVal >> npElem;
    npElemVelo = npElem;
    if(ndim == 2)
    {
      if(npElem == 6)
        npElemPres = 3;
      else
        npElemPres = 4;
    }
    {
      if(npElem == 10)
        npElemPres = 4;
    }

    // read number of nodes
    infile >> stringVal >> nNode;
    nNode_Velo = nNode;
    nNode_Pres = nNode;

    // read number of elements
    infile >> stringVal >> nElem;

    // read number of Dirichlet BCs
    infile >> stringVal >> nDBC;

    // read number of Force BCs
    infile >> stringVal >> nFBC;

    // read number of Output nodes
    infile >> stringVal >> nOutputFaceLoads;

    cout << " ndim              =  " << ndim << endl;
    cout << " nNode             =  " << nNode << endl;
    cout << " nElem             =  " << nElem << endl;
    cout << " npElem            =  " << npElem << endl;
    cout << " nDBC              =  " << nDBC << endl;
    cout << " nFBC              =  " << nFBC << endl;
    cout << " nOutputFaceLoads  =  " << nOutputFaceLoads << endl;

    // read nodal coordinates
    ////////////////////////////////////////////

    cout << " reading nodes " << endl;

    node_coords.resize(nNode);
    for(ii=0; ii<nNode; ++ii)
      node_coords[ii].resize(ndim);

    infile >> stringVal ;
    cout << " reading " << stringVal << endl;
    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] ;
      }
    }
    else
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] >> node_coords[ii][2];
      }
    }

    // read elements
    ////////////////////////////////////////////
    infile >> stringVal ;
    cout << " reading " << stringVal << '\t' << npElem << endl;

    elemConn.resize(nElem);

    if(npElem == 6)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] ;

        //printf("%6d \t %6d \t %6d \t %6d \n", arrayInt[4], arrayInt[5], arrayInt[6], arrayInt[7]);

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 9)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] >> arrayInt[10] >> arrayInt[11] >> arrayInt[12] ;

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 10)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        //infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] >> arrayInt[10] >> arrayInt[11] >> arrayInt[12] >> arrayInt[13] ;
        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] >> arrayInt[10] ;

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[1+ii]-1;

        //printVector(elemConn[ee]);
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

    vector<double>  vecDblTemp(3);

    infile >> stringVec[0] >> stringVec[1] >> stringVec[2] ;
    cout << " reading " << stringVec[0] << stringVec[1] << stringVec[2] << endl;

    nDBC_Velo = 0;
    nDBC_Pres = 0;
    for(ii=0; ii<nDBC; ++ii)
    {
      infile >> arrayInt[0] >> arrayInt[1] >> tempDbl ;

      vecDblTemp[0] = arrayInt[0]-1;
      vecDblTemp[1] = arrayInt[1]-1;
      vecDblTemp[2] = tempDbl;

      if( arrayInt[1] > ndim )
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

    //
    // Read Output data
    //
    ////////////////////////////////////////////

    infile >> stringVal ;
    cout << " reading " << stringVal << endl;

    outputEdges.resize(nOutputFaceLoads);

    for(ii=0; ii<nOutputFaceLoads; ++ii)
    {
      outputEdges[ii].resize(2);

      infile >> arrayInt[0] >> arrayInt[1];

      outputEdges[ii][0] = arrayInt[0]-1;
      outputEdges[ii][1] = arrayInt[1]-1;
    }

    infile.close();

    cout << " Input files have been read successfully \n\n " << endl;

    return;
}



void ExplicitCFD::prepareInputData()
{
    printf("\n     ExplicitCFD::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    assert(ndim > 0 && ndim < 4);

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    // create elements and prepare element data
    //elems = new ElementBase* [nElem];
    //elems = new BernsteinElem2DINSTria6Node* [nElem];
    elems.resize(nElem);

    for(ee=0;ee<nElem;++ee)
    {
      //if(npElem == 6)
      //elems = new BernsteinElem2DINSTria6Node;
      //else if(npElem == 9)
        //elems[ee] = new BernsteinElem2DINSQuad9Node;
      //else if(npElem == 10)
        //elems[ee] = new BernsteinElem3DINSTetra10Node;

      elems[ee].nodeNums = elemConn[ee];

      elems[ee].prepareElemData(node_coords);
    }

    cout << " elements are created and prepated " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // find mid nodes and modify the nodal coordinates and Dirichlet BCs
    //
    ///////////////////////////////////////////////////////////////////

    // create the arrays for proce midnodes
    // index 0 --- 1 - midnode, 0 - not midnode
    // index 1 --- connecting node 1
    // index 2 --- connecting node 2

    midNodeData.resize(nNode);
    for(ii=0; ii<nNode; ++ii)
    {
      midNodeData[ii].resize(3);

      // set the default value to 0 ('not a mid-node')
      midNodeData[ii][0] = 0;  midNodeData[ii][1] = 0;  midNodeData[ii][2] = 0;
    }

    if(npElem == 6)
    {
      for(ee=0; ee<nElem; ++ee)
      {
        ii = elemConn[ee][3];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][0];
        midNodeData[ii][2] = elemConn[ee][1];

        ii = elemConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][1];
        midNodeData[ii][2] = elemConn[ee][2];

        ii = elemConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][2];
        midNodeData[ii][2] = elemConn[ee][0];

        pressure_nodes.push_back(elemConn[ee][0]);
        pressure_nodes.push_back(elemConn[ee][1]);
        pressure_nodes.push_back(elemConn[ee][2]);
      }
    }
    else if(npElem == 9)
    {
      for(ee=0; ee<nElem; ++ee)
      {
        ii = elemConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][0];
        midNodeData[ii][2] = elemConn[ee][1];

        ii = elemConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][1];
        midNodeData[ii][2] = elemConn[ee][2];

        ii = elemConn[ee][6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][2];
        midNodeData[ii][2] = elemConn[ee][3];

        ii = elemConn[ee][7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][3];
        midNodeData[ii][2] = elemConn[ee][0];

        pressure_nodes.push_back(elemConn[ee][0]);
        pressure_nodes.push_back(elemConn[ee][1]);
        pressure_nodes.push_back(elemConn[ee][2]);
        pressure_nodes.push_back(elemConn[ee][3]);
      }
    }
    else if(npElem == 10)
    {
      for(ee=0; ee<nElem; ++ee)
      {
        ii = elemConn[ee][4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][0];
        midNodeData[ii][2] = elemConn[ee][1];

        ii = elemConn[ee][5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][1];
        midNodeData[ii][2] = elemConn[ee][2];

        ii = elemConn[ee][6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][0];
        midNodeData[ii][2] = elemConn[ee][2];

        ii = elemConn[ee][7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][3];
        midNodeData[ii][2] = elemConn[ee][0];

        ii = elemConn[ee][8];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][1];
        midNodeData[ii][2] = elemConn[ee][3];

        ii = elemConn[ee][9];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = elemConn[ee][2];
        midNodeData[ii][2] = elemConn[ee][3];

        pressure_nodes.push_back(elemConn[ee][0]);
        pressure_nodes.push_back(elemConn[ee][1]);
        pressure_nodes.push_back(elemConn[ee][2]);
        pressure_nodes.push_back(elemConn[ee][3]);
      }
    }

    findUnique(pressure_nodes);
    nNode_Pres = pressure_nodes.size();

    pressure_nodes_map.resize(nNode,-1);
    for(ii=0; ii<nNode_Pres; ++ii)
    {
      pressure_nodes_map[pressure_nodes[ii]] = ii;
    }

    cout << " nNode_Pres = " << nNode_Pres << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    nNode_Velo = nNode;
    nNode_Pres = nNode;

    ind = nNode_Velo*ndim;
    globalMassVelo.resize(ind);
    rhsVecVelo.resize(ind);
    setZero(rhsVecVelo);
    rhsVecVelo2 = rhsVecVelo;


    globalMassPres.resize(nNode_Pres);
    rhsVecPres.resize(nNode_Pres);
    setZero(rhsVecPres);
    rhsVecPres2 = rhsVecPres;


    ind = nNode_Velo*ndim;

    velo.resize(ind);
    setZero(velo);

    veloPrev  = velo;
    //veloCur   = velo;
    veloPrev2 = velo;
    //veloPrev3 = velo;
    veloApplied = velo;
    veloDiff    = velo;

    veloDot     = velo;
    veloDotPrev = veloDot;
    //veloDotCur  = veloDot;

    pres.resize(nNode_Pres);
    setZero(pres);

    presPrev  = pres;
    presPrev2 = pres;
    //presPrev3 = pres;
    //presCur   = pres;
    presApplied = pres;
    presDiff    = pres;

    presDot     = velo;
    presDotPrev = pres;
    //presDotCur  = pres;

    soln.resize(nNode*ndof);
    setZero(soln);

    solnInit = soln;

    double  xx, yy, zz, fact;

    cout << " ========== " << endl;

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

    cout << " aaaaaaaaaaaaaaa " << endl;

    Kovasznay analy;

    setZero(veloApplied);
    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        n1 = DirichletBCsVelo[ii][0];
        n2 = DirichletBCsVelo[ii][1];

        jj = n1*ndim+n2;

        //xx = node_coords[n1][0] ;
        //yy = node_coords[n1][1] ;

        //DirichletBCsVelo[ii][2] = analy.computeValue(n2, xx, yy);

        veloApplied[jj] = DirichletBCsVelo[ii][2];
    }

    cout << " ppppppppp " << endl;

    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        nn  = DirichletBCsVelo[ii][0];
        dof = DirichletBCsVelo[ii][1];

        if( midNodeData[nn][0] )
        {
          fact = 0.25*veloApplied[midNodeData[nn][1]*ndim+dof] + 0.25*veloApplied[midNodeData[nn][2]*ndim+dof];

          veloApplied[nn*ndim+dof] = 2.0*(veloApplied[nn*ndim+dof] - fact);
        }
    }
    //printVector(solnApplied);

    setZero(presApplied);
    for(ii=0; ii<nDBC_Pres; ++ii)
    {
        n1 = DirichletBCsPres[ii][0];

        jj = n1;

        //xx = node_coords[n1][0] ;
        //yy = node_coords[n1][1] ;

        //DirichletBCsPres[ii][2] = analy.computeValue(ndim, xx, yy);

        presApplied[jj] = DirichletBCsPres[ii][2];
    }

    cout << " aaaaaaaaaaaaaaa " << endl;

    prepareMatrixPattern();

    printf("     ExplicitCFD::prepareInputData()  .... FINISHED ...\n\n");

    return;
}


void ExplicitCFD::readControlParameters()
{
    cout << " ExplicitCFD::readControlParameters " << endl;

    ifstream  infile("control-parameters.dat");

    if(infile.fail())
    {
       cout << " Could not open 'control-parameters.dat' file " << endl;
       exit(-1);
    }

    //time integration parameters
    timeData[1] = 1.0;   timeData[2] = 0.0;

    string  stringVal;

    //density
    infile >> stringVal >> elemData[0];

    //viscosity
    infile >> stringVal >> elemData[1];

    //Body force in X-, Y- and Z- direction
    infile >> stringVal >> elemData[2] >> elemData[3] >> elemData[4];

    //beta
    infile >> stringVal >> elemData[5];

    //rhoInf = 0.0;
    //am = 1.0/(1.0+rhoInf);
    //gamm1 = 0.5+am;
    //gamm2 = 0.5+am;
    infile >> stringVal >> gamm1;
    infile >> stringVal >> gamm2;

    // CFL number
    infile >> stringVal >> CFL;

    // timestep
    infile >> stringVal >> dt;

    // final time
    infile >> stringVal >> timeFinal;

    // Maximum number of time steps
    infile >> stringVal >> stepsMax;

    // Output file frequency
    infile >> stringVal >> outputFreq;

    // convergence tolerance
    infile >> stringVal >> conv_tol;

    infile.close();

    cout << " Control parameters are successfully read " << endl;

    return;
}




void ExplicitCFD::printInfo()
{
/*
  printf("\n Background Fluid Grid information \n");
    printf("\n ----------------------------------------------------------------------------- \n");
    printf("\t   Problem   Dimension        =  %5d\n\n", DIM);
    printf("\t   Polynomial Degree          =  %5d\t%5d\t%5d\n\n", degree[0], degree[1], degree[2]);
    printf("\t   Number of Elements         =  %5d\t%5d\t%5d\n\n", nelem[0], nelem[1], nelem[2]);
    printf("\t   Origin of the grid         =  %12.6f\t%12.6f\t%12.6f\n\n", origin[0], origin[1], origin[2]);
    printf("\t   Grid Length                =  %12.6f\t%12.6f\t%12.6f\n\n", gridLEN[0], gridLEN[1], gridLEN[2]);
    printf("\t   MAXIMUM REFINEMENT LEVEL   =  %5d\n\n", MAX_LEVEL);
    printf("\t   DOF at each Control Point  =  %5d\n\n", ndof);
    printf("\t   Total number of DOF        =  %5d\n\n", totalDOF);
    printf("\n ----------------------------------------------------------------------------- \n");
*/
    return;
}


void ExplicitCFD::assignBoundaryConditions()
{

}




void ExplicitCFD::setInitialConditions()
{
    //double  xx=0.0, yy=0.0, zz=0.0, fact;
    //double  specVal;

    // adjust the velocity values for the midnoes
    //VectorXd  velTemp;

    /*
    for(nn=0; nn<nNode; nn++)
    {
      if( midNodeData[nn][0] ) // if midnode
      {
        for(dd=0; dd<ndof; dd++)
        {
          xx = 0.25*velTemp(midNodeData[nn][1]*ndof+dd) + 0.25*velTemp(midNodeData[nn][2]*ndof+dd);

          //SolnData.var1Dot[nn*ndof+dd] = 2.0*(velTemp(nn*ndof+dd) - xx);
        }
        //cout << velTemp(nn*ndof) << '\t' << SolnData.var1Dot(nn*ndof) << endl;
      }
    }
    */

    for(int ii=0; ii<nNode_Velo; ++ii)
    {
        //xx = node_coords[ii][0];
        //yy = node_coords[ii][1];
        //zz = node_coords[ii][2];

        //veloPrev[ii*2] = 2.0*yy*(3.0-yy)/3.0;
        veloPrev[ii*ndim] = 1.0;

        //veloPrev[ii*ndim] = 16.0*0.45*yy*zz*(0.41-yy)*(0.41-zz)/0.41/0.41/0.41/0.41;
    }
    velo = veloPrev;

    return;
}






void ExplicitCFD::setTimeParam()
{
  //SolnData.setTimeParam();

  return;
}



void ExplicitCFD::updateIterStep()
{
  //SolnData.updateIterStep();

  return;
}




void ExplicitCFD::writeNodalData()
{
  return;
}



int  ExplicitCFD::solveExplicitStep()
{
    // OpenMP parallelisation using a single rhs vector and critical construct

    calcMassMatrixForExplicitDynamics();

    cout << " Solving the ExplicitCFD " << endl;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    setZero(velo);      setZero(veloPrev);      setZero(veloCur);
    setZero(veloDot);   setZero(veloDotPrev);   setZero(veloDotCur);
    setZero(pres);      setZero(presPrev);      setZero(presCur);
    setZero(presDot);   setZero(presDotPrev);   setZero(presDotCur);

    int  stepsCompleted=0;
    int  nsize_velo = rhsVecVelo.size();
    int  nsize_pres = rhsVecPres.size();
    int  aa, bb, dd, ee, ii, jj, kk, count, row, col, ind;

    double  norm_velo=1.0, norm_pres=1.0;
    double  dtCrit=1.0e10, fact, fact1, fact2;
    double  dtgamma11, dtgamma12, norm_velo_diff, norm_pres_diff;
    double  timeNow=0.0, timeFact=0.0, dt=0.0;
    double  Re=2.0/elemData[1];

    cout << " Re = " << Re << endl;

    //VectorXd  Flocal1(npElem*ndof), Flocal2(npElem*ndof);
    VectorXd  TotalForce(3);
    //double  FlocalVelo[npElem*ndof], FlocalPres[npElem*ndof];

    int nthreads, thread_cur;

    //#pragma omp parallel
    //nthreads = omp_get_num_threads();
    //cout << " nthreads = " << nthreads << endl;

    double time1 = omp_get_wtime();
    //auto time1 = chrono::steady_clock::now();

    //setInitialConditions();
    //timeFact = 1.0;
    //Time loop
    while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )
    {
        if(stepsCompleted < 5000)
        {
          timeFact = 0.5*(1-cos(PI*stepsCompleted/5000));
          //timeFact = stepsCompleted/5000.0;
        }
        else
        {
          timeFact = 1.0;
        }
        //timeFact = 1.0;

          //#pragma omp sections
          //{
            //#pragma omp section
            for(ii=0; ii<nsize_velo; ii++)
            {
              rhsVecVelo[ii]=0.0;
            }

            //#pragma omp section
            for(ii=0; ii<nsize_pres; ii++)
            {
              rhsVecPres[ii]=0.0;
            }
          //}

        //#pragma omp parallel private(ee, ii, jj, kk, dd, fact) default(shared)
        //{
          //Loop over elements and compute the RHS and time step
            double  FlocalVelo[npElem*ndof], FlocalPres[npElem*ndof];
          dtCrit=1.0e10;
          //#pragma acc parallel loop private(fact)  reduction(min : dtCrit) //schedule(dynamic,100) //shared(ndim, nElem, rhsVecVelo, rhsVecPres, elemConn)
          for(int ee=0; ee<nElem; ee++)
          {
            fact = timeNow - dt;
            //Compute the element force vector, including residual force and time step

            dtCrit = min(dtCrit, elems[ee].ResidualIncNavStokesAlgo1(node_coords, elemData, timeData, velo, veloPrev, veloDot, veloDotPrev, pres, presPrev, FlocalVelo, FlocalPres, fact) );

            //printVector(FlocalVelo);
            //printVector(Flocal2);
            int ii, jj, kk, dd;
            //Assemble the element vector
            //#pragma omp critical
            //{
              for(ii=0; ii<npElemVelo; ii++)
              {
                jj = ndim*ii;
                kk = ndim*elemConn[ee][ii];

                for(dd=0; dd<ndim; dd++)
                {
                  //#pragma acc atomic
                  rhsVecVelo.data()[kk+dd] += FlocalVelo[jj+dd];
                }
              }
              for(ii=0; ii<npElemPres; ii++)
              {
                  //#pragma acc atomic
                  rhsVecPres.data()[elemConn[ee][ii]]   += FlocalPres[ii];
              }
            //}
          } //LoopElem
        //}

        //cout << dtCrit << endl;

          dt = dtCrit*CFL/gamm1;

          // Add specified nodal force 
          //addExternalForces(timeFact);

          // compute the solution at t_{n+1}
          dtgamma11 = dt*gamm1;
          dtgamma12 = dt*(1.0-gamm1);

          //#pragma omp for
          for(ii=0; ii<totalDOF_Velo; ii++)
          {
            jj = assyForSolnVelo[ii];

            veloDot[jj] = rhsVecVelo[jj]/globalMassVelo[jj];
            velo[jj]    = veloPrev[jj] + dtgamma11*veloDot[jj] + dtgamma12*veloDotPrev[jj];
          }

          //#pragma omp for
          for(ii=0; ii<totalDOF_Pres; ii++)
          {
            jj = assyForSolnPres[ii];

            presDot[jj] = rhsVecPres[jj]/globalMassPres[jj];
            pres[jj]    = presPrev[jj] + dtgamma11*presDot[jj] + dtgamma12*presDotPrev[jj];
          }

          //#pragma omp single
          applyBoundaryConditions(timeFact);
          //printVector(velo);

          // compute the norms and store the variables
          norm_velo = 0.0;
          norm_velo_diff = 0.0;
          //#pragma omp for reduction(+:norm_velo,norm_velo_diff)
          for(ii=0; ii<nsize_velo; ii++)
          {
            norm_velo += velo[ii]*velo[ii];

            fact1 = velo[ii]-veloPrev[ii];
            norm_velo_diff += fact1*fact1;

            // store the variables
            //veloPrev3  = veloPrev2;
            //veloPrev2  = veloPrev;
            veloPrev[ii]  = velo[ii];
            veloDotPrev[ii]  = veloDot[ii];
          }

          norm_pres = 0.0;
          norm_pres_diff = 0.0;
          //#pragma omp for reduction(+:norm_pres, norm_pres_diff)
          for(ii=0; ii<nsize_pres; ii++)
          {
            norm_pres += pres[ii]*pres[ii];

            fact1 = pres[ii]-presPrev[ii];
            norm_pres_diff += fact1*fact1;

            //presPrev3  = presPrev2;
            //presPrev2  = presPrev;
            presPrev[ii]     = pres[ii];
            presDotPrev[ii]  = presDot[ii];
          }

        if( std::isnan(norm_velo) || std::isnan(norm_pres) )
        {
          cerr << " NAN encountered in the solution ... " << endl;
          cerr << " Program has been terminated ... " << endl;
          exit(-1);
        }

        norm_velo = sqrt(norm_velo_diff/norm_velo);
        norm_pres = sqrt(norm_pres_diff/norm_pres);
        if(stepsCompleted == 0)
        {
          norm_velo = 1.0;
          norm_pres = 1.0;
        }

        if( (stepsCompleted%outputFreq == 0) || (norm_velo < conv_tol) )
        {
            cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;
            cout << " velocity difference norm = " << '\t' << norm_velo << endl;

            //postProcess();

            /*
            setZero(TotalForce);
            for(ii=0; ii<outputEdges.size(); ii++)
            {
              elems[outputEdges[ii][0]]->CalculateForces(outputEdges[ii][1], node_coords, elemData, timeData, velo, pres, TotalForce);
            }

            fout_convdata << timeNow << '\t' << stepsCompleted << '\t' << norm_velo << '\t' << norm_pres ;
            fout_convdata << '\t' << TotalForce[0] << '\t' << TotalForce[1] << '\t' << TotalForce[2] << endl;
            */
        }

        if(norm_velo < conv_tol)
        {
          cout << " Solution convdataged below the specified tolerance " << endl;
          break;
        }

        stepsCompleted = stepsCompleted + 1;

        timeNow = timeNow + dt;
    } //Time loop

    postProcess();

    //auto tend = chrono::steady_clock::now();
    double tend = omp_get_wtime();

    //auto duration = chrono::duration_cast<chrono::milliseconds>(tend-time1).count();
    double  duration = tend - time1;

    cout << " \n \n Total time taken = " << duration << " seconds \n\n" << endl;

    return 0;
}
//




