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
#include <Eigen/SuperLUSupport>


using namespace std;



int main(int argc, char* argv[])
{
    double tstart, tend;

    int ndim=2, ndof=3, npElem=6;
    double fact, xNode[50], yNode[50];

    vector<vector<int> > ElemDofArray;

    int  nElem, nNode, nDBC, nFBC;
    int  ee, ii, jj, kk, ind, count, row, col;
    int  n1, n2, n3, n4, n5, n6, nsize;
    int  nn, dof;

    string  infileNodes, infileElems, infileDBCs, infileOutput;
    string  infileFBCs, charTemp, outFileName;

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

    double  val[50];
    int  val2[50];

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

    vector<vector<double> >  DirichletBCs;
    vector<double>  vecDblTemp(3);

    nDBC = 0;
    while(infile_DBCs >> val[0] >> val[1] >> val[2] )
    {
      vecDblTemp[0] = val[0]-1;
      vecDblTemp[1] = val[1]-1;
      vecDblTemp[2] = val[2];

      //if( val[1] <= ndof )
      //{
        DirichletBCs.push_back(vecDblTemp);
        nDBC++;
      //}
    }

    cout << " nDBC  = " << '\t' << nDBC << endl;

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
      }
    }

    vector<vector<int> >  ID;
    vector<vector<bool> >  NodeType;
    vector<int>  assyForSoln;

    NodeType.resize(nNode);

    ID.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      NodeType[ii].resize(ndof);
      ID[ii].resize(ndof);

      for(jj=0;jj<ndof;jj++)
      {
        NodeType[ii][jj] = false;
        ID[ii][jj] = -1;
      }
    }

    // fix the pressure at all the mid nodes
    for(ii=0; ii<nNode; ii++)
    {
      if(midNodeData[ii][0])
        NodeType[ii][2] = true;
    }

    if(npElem == 9)
    {
      for(ee=0; ee<nElem; ee++)
      {
        NodeType[elemNodeConn[ee][8]][2] = true;
      }
    }

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC; ii++)
    {
      //cout << ii << '\t' << DirichletBCs[ii][0] << '\t' << DirichletBCs[ii][1] << endl;
      NodeType[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    //for(ii=0; ii<totalDOF; ii++)
      //cout << ii << '\t' << assyForSoln[ii] << endl;

    int totalDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << NodeType[ii][jj] << endl;
        if(!NodeType[ii][jj])
        {
          ID[ii][jj] = totalDOF++;
          assyForSoln.push_back(ii*ndof+jj);
        }
      }
    }

      cout << " Mesh statistics .....\n" << endl;
      cout << " nElem          = " << '\t' << nElem << endl;
      cout << " nNode          = " << '\t' << nNode  << endl;
      cout << " npElem         = " << '\t' << npElem << endl;
      cout << " ndof           = " << '\t' << ndof << endl;
      cout << " Total DOF      = " << '\t' << totalDOF << endl;


      vector<vector<int> >  LM, forAssyMat;

      LM.resize(nElem);

      for(ee=0;ee<nElem;ee++)
      {
        npElem = elemNodeConn[ee].size();

        //printVector(IEN[ee]);

        ind = ndof*npElem;
        LM[ee].resize(ind);

        for(ii=0;ii<npElem;ii++)
        {
          ind = ndof*ii;

          kk = elemNodeConn[ee][ii];

          for(jj=0;jj<ndof;jj++)
          {
            LM[ee][ind+jj] = ID[kk][jj];
          }
        }
      }

      printf("\n element DOF values initialised \n\n");
      printf("\n Preparing matrix pattern \n\n");

      // compute matrix pattern needs to be prepared only for the implicit solver

      forAssyMat.clear();
      forAssyMat.resize(totalDOF);

      int *tt, r, c;

    for(ee=0;ee<nElem;ee++)
    {
      tt = &(LM[ee][0]);
      nsize = LM[ee].size();

      for(ii=0;ii<nsize;ii++)
      {
        r = tt[ii];

        if(r != -1)
        {
          for(jj=0;jj<nsize;jj++)
          {
            if(tt[jj] != -1)
            {
              //printf("ii.... %5d \t %5d \t %5d \t %5d \n",ii, jj, r, tt[jj]);
              forAssyMat[r].push_back(tt[jj]);
            }
          }
        }
      }
    }

    printf("\n Preparing matrix pattern DONE \n\n");

    VectorXi  nnzVec(totalDOF);

    int nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(forAssyMat[ii]);

      nnzVec[ii] = forAssyMat[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;


    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   Number of non-zeros = %5d \n\n", nnz);
       printf("   dof to dof connectivity ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<forAssyMat[ii].size();jj++)
            cout << '\t' << forAssyMat[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");
    }

      /////////////////////////////////////////
      //
      // Eigen based solver
      //
      /////////////////////////////////////////

      cout << " Eigen based solver " << totalDOF << endl;

      //solver->rhsVec.resize(nRow);

      SparseMatrixXd  mtx(totalDOF, totalDOF);

      mtx.reserve(nnz);
      mtx.reserve(nnzVec);

      for(ii=0;ii<totalDOF;ii++)
      {
        for(jj=0;jj<forAssyMat[ii].size();jj++)
        {
          mtx.coeffRef(ii, forAssyMat[ii][jj]) = 0.0;
        }
      }

      mtx.makeCompressed();


////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////

      VectorXd  pres, presCur, presPrev, presPrev2, presPrev3;
      VectorXd  presDot, presDotPrev, presDotCur, presDiff;
      VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3;
      VectorXd  veloDot, veloDotPrev, acceCur;

      VectorXd  solnApplied, solnVTK, soln, solnCur, solnTemp, rhsVec;

      ind = nNode*ndim;

      velo.resize(ind);
      velo.setZero();

      veloPrev  = velo;
      veloCur   = velo;
      veloPrev2 = velo;
      veloPrev3 = velo;

      veloDot     = velo;
      veloDotPrev = veloDot;
      acceCur  = veloDot;

      pres.resize(nNode);
      pres.setZero();

      presPrev  = pres;
      presPrev2 = pres;
      presPrev3 = pres;
      presCur   = pres;

      presDot     = velo;
      presDotPrev = pres;
      presDotCur  = pres;


      soln.resize(nNode*ndof);
      soln.setZero();

      solnCur = soln;
      solnApplied = soln;

      solnTemp.resize(totalDOF);
      solnTemp.setZero();

      rhsVec = solnTemp;

      ////////////////////////////////////////////////
      // Computations start from here
      ////////////////////////////////////////////////

      char fname[200];
      //sprintf(fname,"(),"-explicit.dat");

      ofstream fout("convergence-data2.dat");

      if(fout.fail())
      {
        cout << " Could not open the Output file" << endl;
        exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(14);


      double  xx, yy;
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
        }
      }
//

      double  elemData[50], timeData[50];
      double  norm_rhs=1000.0, norm_pres=1000.0;

      //time integration parameters
      timeData[1] = 1.0;   timeData[2] = 0.0;

      //density
      elemData[0] = 1.0;
      //viscosity
      //elemData[1] = 1.0;
      //elemData[1] = 0.025;
      elemData[1] = 1.0/229.0;
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

      double  dt = 0.01;
      double  timeFact=0.0;
      double  timeNow=0.0;
      double  timeFinal = 4000.0;
      double  num, denom;
      double  fact1, fact2;

      // create elements and prepare element data
      BernsteinElem2DINSTria6Node   **elems;
      //BernsteinElem2DINSQuad9Node   **elems;

      elems = new BernsteinElem2DINSTria6Node* [nElem];
      //elems = new BernsteinElem2DINSQuad9Node* [nElem];

      for(ee=0;ee<nElem;ee++)
      {
        elems[ee] = new BernsteinElem2DINSTria6Node;
        //elems[ee] = new BernsteinElem2DINSQuad9Node;

        elems[ee]->nodeNums = elemNodeConn[ee];

        //elems[ee]->SolnData = &(SolnData);

        elems[ee]->prepareElemData(node_coords);

        elems[ee]->forAssyVec = LM[ee];
      }

      cout << " elements are created and prepated " << endl;
      cout << " Computing the solution \n" << endl;

      //SimplicialLDLT<SparseMatrix<double> > solver;
      SuperLU<SparseMatrixXd > solver;

      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////

      velo.setZero();      veloPrev.setZero();      veloCur.setZero();
      veloDot.setZero();      veloDotPrev.setZero();      acceCur.setZero();
      pres.setZero();      presPrev.setZero();      presCur.setZero();
      presDot.setZero();   presDotPrev.setZero();   presDotCur.setZero();

      //for(ii=0; ii<nNode_Velo; ii++)
      //{
        //veloPrev(ii*2) = 1.0;
      //}

      //velo = veloPrev;

        //SolnData.applyDirichletBCs(timeFact);

        solnApplied.setZero();

        //Add specified Dirichlet BC
        for(ii=0; ii<nDBC; ii++)
        {
          nn  = DirichletBCs[ii][0];
          dof = DirichletBCs[ii][1];

          jj = nn*ndof+dof;

          solnApplied(jj) = DirichletBCs[ii][2];
        }

        //printVector(soln);

        for(ii=0; ii<nDBC; ii++)
        {
          nn   = (int) (DirichletBCs[ii][0]);
          dof  = (int) (DirichletBCs[ii][1]);

          if( midNodeData[nn][0] )
          {
            fact = 0.25*solnApplied(midNodeData[nn][1]*ndof+dof) + 0.25*solnApplied(midNodeData[nn][2]*ndof+dof);

            solnApplied[nn*ndof+dof] = 2.0*(solnApplied(nn*ndof+dof) - fact);
          }
        }
        //printVector(solnApplied);

      //Time loop
      //while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )

      vector<int>  vecTempInt;

      ind = npElem*ndof;

      VectorXd  Flocal(ind);
      MatrixXd  Klocal(ind, ind);

      int iter=0, loadIncr;

      soln.setZero();

      for(loadIncr=0; loadIncr<10; loadIncr++)
      {
        //timeFact = 0.5*(1-cos(PI*stepsCompleted/2000));
        timeFact = 0.1*(loadIncr+1);

        for(ii=0; ii<nDBC; ii++)
        {
          nn   = (int) (DirichletBCs[ii][0]);
          dof  = (int) (DirichletBCs[ii][1]);

          jj = nn*ndof + dof;

          soln[jj] = timeFact*solnApplied(jj);
        }

        //cout << " aaaaaaaaaaa " << endl;

      for(iter=0; iter<20; iter++)
      {
        //Loop over elements and compute the RHS

        mtx *= 0.0;
        rhsVec.setZero();

        for(ee=0; ee<nElem; ee++)
        {
            fact = timeNow - dt;

            //Compute the element force vector, including residual force
            elems[ee]->StiffnessAndResidual(node_coords, elemData, timeData, soln, Klocal, Flocal, fact);

            vecTempInt = elems[ee]->forAssyVec;
            ind = vecTempInt.size();

            //Assemble the element vector
            for(ii=0; ii<ind; ii++)
            {
              r = vecTempInt[ii];
              if( r != -1)
              {
                rhsVec(r)  += Flocal(ii);

                for(jj=0; jj<ind; jj++)
                {
                  c = vecTempInt[jj];
                  if( c != -1)
                  {
                    mtx.coeffRef(r, c)  += Klocal(ii, jj);
                  }
                }
              }
            }
        } //LoopElem

        // Add specified nodal force 
        //SolnData.addNodalForces(timeFact);

        norm_rhs = rhsVec.norm();

        cout << " RHS norm = " << norm_rhs << endl;

        if(norm_rhs < 1.0e-8)
        {
          cout << " Solution converged below the specified tolerance " << endl;
          break;
        }
        else
        {
          solver.compute(mtx);

          solnTemp = solver.solve(rhsVec);

          for(ii=0; ii<totalDOF; ii++)
          {
            soln(assyForSoln[ii]) += solnTemp(ii);
          }

          //printVector(soln);
        }
      } // iteration loop

          // write the solution to the VTK file

          velo.resize(nNode*3);
          velo.setZero();
          pres.setZero();

          /*
          for(ii=0; ii<nNode; ii++)
          {
            n1 = ii*2;
            n2 = ii*3;

            velo(n1)   = soln(n2);
            velo(n1+1) = soln(n2+1);

            pres(ii)   = soln(n2+2);

            //pres(ii)   = soln(ii);
          }

          for(ii=0; ii<nNode; ii++)
          {
            if(midNodeData[ii][0])
            {
              jj = ii*2;
              n1 = midNodeData[ii][1]*2;
              n2 = midNodeData[ii][2]*2;

              velo(jj)    = 0.5*velo(jj)   + 0.25*velo(n1)   + 0.25*velo(n2);
              velo(jj+1)  = 0.5*velo(jj+1) + 0.25*velo(n1+1) + 0.25*velo(n2+1);

              pres(ii)    = 0.5*(pres(midNodeData[ii][1]) + pres(midNodeData[ii][2]) );
            }
          }
          */

          //writevtk(ndim, node_coords, elemNodeConn, midNodeData, soln, fileCount);

          fileCount = fileCount+1;

          //fact1 = 0.0;
          //fact2 = 0.0;
          //for(ii=0; ii<OutputData.size(); ii++)
          //{
            //jj = OutputData[ii]*2;

            //fact1 += rhsVec(jj);
            //fact2 += rhsVec(jj+1);
          //}

          //fout << timeNow << '\t' << stepsCompleted << '\t' << norm_rhs << '\t' << norm_pres ;
          //fout << '\t' << fact1 << '\t' << fact2 << endl;


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

      } //Time loop

    //
    cout << " Computing errors \n " << endl;
    double totalError = 0.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(ee=0; ee<nElem; ee++)
      {
        //Compute the element force vector, including residual force
        totalError += elems[ee]->CalculateError(node_coords, elemData, timeData, soln, veloDot, pres, timeNow, index);
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
    //

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







