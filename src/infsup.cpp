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

    int ndim=2, ndof=2, npElem=6;
    double fact, xNode[50], yNode[50];

    vector<vector<int> > ElemDofArray;

    int  nElem, nNode, nDBC, nFBC;
    int  ee, ii, jj, kk, ind, count, row, col;
    int  n1, n2, n3, n4, n5, n6, nsize;
    int  nn, dof;
    int  rr, cc;


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

    int nNode_Pres=0;
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

    findUnique(pressure_nodes);
    nNode_Pres = pressure_nodes.size();

    vector<int>  pressure_nodes_map(nNode,-1);

    for(ii=0; ii<nNode_Pres; ii++)
    {
      pressure_nodes_map[pressure_nodes[ii]] = ii;
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


      vector<vector<int> >  LM;

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

      /////////////////////////////////////////
      //
      // Eigen based solver
      //
      /////////////////////////////////////////

      cout << " Eigen based solver " << totalDOF << endl;

////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////


      double  xx, yy;

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


      double  elemData[50];

      //density
      elemData[0] = 1.0;
      //viscosity
      //elemData[1] = 1.0;
      elemData[1] = 0.025;
      //Body force in X-, Y- and Z- direction
      elemData[2] = 0.0;   elemData[3] = 0.0; elemData[4] = 0.0;
      //beta
      elemData[5] = 2.0;

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

      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////

     cout << " nNode_Pres = " << nNode_Pres << endl;

      vector<int>  vecTempInt1, vecTempInt2;

    // to compute inf-sup number

    //nNode_Pres = nNode;

    ind = npElem*ndim;
    MatrixXd  Kuu(ind,ind), Kup(ind,9), Kpp(9,9);

    MatrixXd  Kglobal, Mglobal, eigen_vectors;
    VectorXd  eigen_values, eigvec, Flocal, vecTemp;
    MatrixXd  Vmat, Qmat, Bmat;

    Vmat.resize(totalDOF, totalDOF);
    Bmat.resize(totalDOF, nNode_Pres);
    Qmat.resize(totalDOF, totalDOF);

    Kglobal.resize(totalDOF, totalDOF);
    Mglobal.resize(nNode_Pres, nNode_Pres);

    Kglobal.setZero();
    Mglobal.setZero();
    Vmat.setZero();
    Bmat.setZero();
    Qmat.setZero();

    /////////////////////////////////////////
    // compute and assemble matrices
    /////////////////////////////////////////

    for(ee=0; ee<nElem; ee++)  // loop over all the elements
    {
        //cout << "       elem... : " << (ee+1) << endl;

        elems[ee]->toComputeInfSupCondition(node_coords, elemData, Kuu, Kup, Kpp);

        vecTempInt1=elems[ee]->forAssyVec;
        vecTempInt2=elems[ee]->nodeNums;

        n1=vecTempInt1.size();
        n2=3;

        //printVector(vecTempInt1);
        //printVector(vecTempInt2);

        // Kuu
        for(ii=0; ii<n1; ii++)
        {
          rr = vecTempInt1[ii];
          if( rr != -1 )
          {
            for(jj=0; jj<n1; jj++)
            {
              cc = vecTempInt1[jj];
              if( cc != -1 )
              {
                 Kglobal(rr, cc) += Kuu(ii,jj);
              }
            }
          }
        }

        //cout << " aaaaaaaaaaa " << endl;

        // Kup
        for(ii=0; ii<n1; ii++)
        {
          rr = vecTempInt1[ii];
          if( rr != -1 )
          {
            for(jj=0; jj<n2; jj++)
            {
              cc = pressure_nodes_map[vecTempInt2[jj]];
              //cc = vecTempInt2[jj];

              Bmat(rr,cc) += Kup(ii,jj);
            }
          }
        }

        //cout << " aaaaaaaaaaa " << endl;

        // Kpp
        for(ii=0; ii<n2; ii++)
        {
          rr = pressure_nodes_map[vecTempInt2[ii]];
          //rr = vecTempInt2[ii];

          for(jj=0; jj<n2; jj++)
          {
            cc = pressure_nodes_map[vecTempInt2[jj]];
            //cc = vecTempInt2[jj];

            Mglobal(rr, cc) += Kpp(ii,jj);
          }
        }
    }

    cout << "  Solving eigenvalue problem ... " << endl;

    Kglobal = Kglobal.inverse();
    //printMatrix(Vmat);
    //printf("\n\n");
    //printMatrix(globalK);
    //printf("\n\n");
    //printMatrix(globalM);
    //printf("\n\n");
    Kglobal = (Bmat.transpose()*Kglobal)*Bmat;

    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(Kglobal, Mglobal, EigenvaluesOnly);
    eigen_values = es.eigenvalues();
    //eigen_vectors = es.eigenvectors();

    //EigenSolver<MatrixXd>  es(Kglobal, EigenvaluesOnly);
    //cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
    //eigen_values = es.eigenvalues().col(0);

    cout << "  Eigen analysis successfully completed ... " << endl;

    cout << " The first " << min(20, (int) eigen_values.rows()) << " eigenvalues are ... " << endl;

    for(int ii=0;ii<min(20, (int) eigen_values.rows());ii++)
      printf("\t %5d \t %12.10f \t %12.10f \n", (ii+1), eigen_values(ii), sqrt(abs(eigen_values(ii))));



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







