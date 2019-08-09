
#include "ExplicitCFD.h"
#include <fstream>

#ifndef VTK_LEGACY_FORMAT
#include <vtkSmartPointer.h>

#include <vtkPoints.h>
#include <vtkVertex.h>

#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>

#include <vtkQuadraticTriangle.h>
#include <vtkBiQuadraticQuad.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticHexahedron.h>

#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

using namespace std;



void  ExplicitCFD::postProcess()
{
    //cout << " ExplicitCFD::postProcess " << endl;

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////

#ifdef VTK_LEGACY_FORMAT

  postProcessVTKlegacy();

#else

    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    //vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    //vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    //vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    //vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK   =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};

    vecVTK->SetNumberOfTuples(nNode);
    vecVTK->SetNumberOfComponents(3);
    scaVTK->SetNumberOfTuples(nNode);

    vecVTK->SetName("velocity");
    scaVTK->SetName("pressure");

    vtkIdType pt[10];

    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(node_coords[ii][0], node_coords[ii][1], 0.0);

        if( midNodeData[ii][0] )
        {
          jj = ii*ndim;
          n1 = midNodeData[ii][1]*ndim;
          n2 = midNodeData[ii][2]*ndim;

          vec[0] = 0.5*velo[jj]   + 0.25*velo[n1]   + 0.25*velo[n2];
          vec[1] = 0.5*velo[jj+1] + 0.25*velo[n1+1] + 0.25*velo[n2+1];

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          val = 0.5*(pres[midNodeData[ii][1]] + pres[midNodeData[ii][2]]);

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
        else
        {
          kk = ii*ndim;

          vec[0] = velo[kk];
          vec[1] = velo[kk+1];
          val    = pres[ii];

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
      }

      for(ee=0; ee<nElem; ee++)
      {
        npElem = elemConn[ee].size();

        if(npElem == 6)
        {
          for(ii=0; ii<npElem; ii++)
            tria6VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if(npElem == 9)
        {
          vec[0]  = 0.25*velo[elemConn[ee][8]*ndim];
          vec[0] += 0.0625*(velo[elemConn[ee][0]*ndim]+velo[elemConn[ee][1]*ndim]+velo[elemConn[ee][2]*ndim]+velo[elemConn[ee][3]*ndim]);
          vec[0] += 0.1250*(velo[elemConn[ee][4]*ndim]+velo[elemConn[ee][5]*ndim]+velo[elemConn[ee][6]*ndim]+velo[elemConn[ee][7]*ndim]);

          vec[1]  = 0.25*velo[elemConn[ee][8]*ndim+1];
          vec[1] += 0.0625*(velo[elemConn[ee][0]*ndim+1]+velo[elemConn[ee][1]*ndim+1]+velo[elemConn[ee][2]*ndim+1]+velo[elemConn[ee][3]*ndim+1]);
          vec[1] += 0.1250*(velo[elemConn[ee][4]*ndim+1]+velo[elemConn[ee][5]*ndim+1]+velo[elemConn[ee][6]*ndim+1]+velo[elemConn[ee][7]*ndim+1]);

          vecVTK->InsertTuple(elemConn[ee][8], vec);

          val = 0.25*(pres[elemConn[ee][0]]+pres[elemConn[ee][1]]+pres[elemConn[ee][2]]+pres[elemConn[ee][3]]);

          scaVTK->SetTuple1(elemConn[ee][8], val);

          for(ii=0; ii<npElem; ii++)
            quad9VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<nNode; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(node_coords[ii][0], node_coords[ii][1], node_coords[ii][2]);

        if( midNodeData[ii][0] )
        {
          jj = ii*ndim;
          n1 = midNodeData[ii][1]*ndim;
          n2 = midNodeData[ii][2]*ndim;

          vec[0] = 0.5*velo[jj]   + 0.25*velo[n1]   + 0.25*velo[n2];
          vec[1] = 0.5*velo[jj+1] + 0.25*velo[n1+1] + 0.25*velo[n2+1];
          vec[2] = 0.5*velo[jj+2] + 0.25*velo[n1+2] + 0.25*velo[n2+2];

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          val   = 0.5*(pres[midNodeData[ii][1]] + pres[midNodeData[ii][1]]) ;

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
        else
        {
          kk = ii*ndim;

          vec[0] = velo[kk];
          vec[1] = velo[kk+1];
          vec[2] = velo[kk+2];
          val    = pres[ii];

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
      }

      for(ee=0; ee<nElem; ee++)
      {
        npElem = elemConn[ee].size();

        if(npElem == 10)
        {
          for(ii=0; ii<npElem; ii++)
            tetra10VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);

    //Write the file.
    char VTKfilename[200];
    sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtu");
    writerUGridVTK->SetFileName(VTKfilename);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    fileCount = fileCount+1;
#endif

    return;
}






void  ExplicitCFD::postProcessVTKlegacy()
{
    //  using the legacy VTK format

    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};

    char VTKfilename[200];
    sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtk");

    ofstream fout(VTKfilename);

    if(fout.fail())
    {
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(6);


    //Directives
    fout <<  "# vtk DataFile Version 4.0" << endl;
    fout <<  "IncNS example" << endl;

    // ASCII or Binary
    fout <<  "ASCII" << endl;

    // Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
    fout <<  "DATASET UNSTRUCTURED_GRID" << endl;


    // Coordinates of the points (nodes)
    fout <<  "POINTS " << '\t' << nNode << '\t' << " float" << endl;

    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ii++)
      {
        fout <<  node_coords[ii][0] << '\t' << node_coords[ii][1] << '\t' <<  0.0 <<  endl;
      }
    }
    else
    {
      for(ii=0; ii<nNode; ii++)
      {
        fout <<  node_coords[ii][0] << '\t' << node_coords[ii][1] << '\t' <<  node_coords[ii][2] <<  endl;
      }
    }

    // Element<->Nodes connectivity
    // In VTK terminology, Cell<->Points
    // <number of nodes> <n1> <n2> <n3> ...
    // Starting index is 0
    int  ind = nElem*(npElem+1);
    fout <<  "CELLS " << '\t' << nElem << '\t' << ind << endl;


    if(ndim == 2)
    {
      if(npElemVelo == 6)
      {
        for(ee=0; ee<nElem; ee++)
        {
          fout << npElemVelo      << '\t' << elemConn[ee][0] << '\t' << elemConn[ee][1] << '\t';
          fout << elemConn[ee][2] << '\t' << elemConn[ee][3] << '\t' << elemConn[ee][4] << '\t' ;
          fout << elemConn[ee][5] << endl;
        }
      }
      if(npElemVelo == 9)
      {
        for(ee=0; ee<nElem; ee++)
        {
          fout << npElemVelo      << '\t' << elemConn[ee][0] << '\t' << elemConn[ee][1] << '\t';
          fout << elemConn[ee][2] << '\t' << elemConn[ee][3] << '\t' << elemConn[ee][4] << '\t' ;
          fout << elemConn[ee][5] << '\t' << elemConn[ee][6] << '\t' << elemConn[ee][7] << '\t' ;
          fout << elemConn[ee][8] << endl;
        }
      }
    }
    else
    {
      if(npElemVelo == 10)
      {
        for(ee=0; ee<nElem; ee++)
        {
          fout << npElemVelo      << '\t' << elemConn[ee][0] << '\t' << elemConn[ee][1] << '\t';
          fout << elemConn[ee][2] << '\t' << elemConn[ee][3] << '\t' << elemConn[ee][4] << '\t' ;
          fout << elemConn[ee][5] << '\t' << elemConn[ee][6] << '\t' << elemConn[ee][7] << '\t' ;
          fout << elemConn[ee][8] << '\t' << elemConn[ee][9] << endl;
        }
      }
    }

    // Cell type, as per VTK
    fout << "CELL_TYPES" << '\t' << nElem << endl;
    if(ndim == 2)
    {
      if(npElemVelo == 6)
      {
        for (ee=0; ee<nElem; ee++)
          fout <<  22 <<  endl;
      }
      if(npElemVelo == 9)
      {
        for (ee=0; ee<nElem; ee++)
          fout <<  28 <<  endl;
      }
    }
    else
    {
      if(npElemVelo == 10)
      {
        for (ee=0; ee<nElem; ee++)
          fout <<  24 <<  endl;
      }
    }


    // Point data
    fout << "POINT_DATA" << '\t' << nNode << endl;
    fout <<  "VECTORS velocity float" <<  endl;

    for(ii=0; ii<nNode; ii++)
    {
      if( midNodeData[ii][0] )
      {
        jj = ii*ndim;
        n1 = midNodeData[ii][1]*ndim;
        n2 = midNodeData[ii][2]*ndim;

        vec[0] = 0.5*velo[jj]   + 0.25*velo[n1]   + 0.25*velo[n2];
        vec[1] = 0.5*velo[jj+1] + 0.25*velo[n1+1] + 0.25*velo[n2+1];
        if(ndim ==  3)
          vec[2] = 0.5*velo[jj+2] + 0.25*velo[n1+2] + 0.25*velo[n2+2];
      }
      else
      {
        kk = ii*ndim;

        vec[0] = velo[kk];
        vec[1] = velo[kk+1];
          if(ndim ==  3)
            vec[2] = velo[kk+2];
      }

      fout << vec[0] <<  '\t' << vec[1] << '\t' << vec[2] << endl;
    }

    if(npElemVelo == 9)
    {
      for(ee=0; ee<nElem; ee++)
      {
        vec[0]  = 0.25*velo[elemConn[ee][8]*ndim];
        vec[0] += 0.0625*(velo[elemConn[ee][0]*ndim]+velo[elemConn[ee][1]*ndim]+velo[elemConn[ee][2]*ndim]+velo[elemConn[ee][3]*ndim]);
        vec[0] += 0.1250*(velo[elemConn[ee][4]*ndim]+velo[elemConn[ee][5]*ndim]+velo[elemConn[ee][6]*ndim]+velo[elemConn[ee][7]*ndim]);

        vec[1]  = 0.25*velo[elemConn[ee][8]*ndim+1];
        vec[1] += 0.0625*(velo[elemConn[ee][0]*ndim+1]+velo[elemConn[ee][1]*ndim+1]+velo[elemConn[ee][2]*ndim+1]+velo[elemConn[ee][3]*ndim+1]);
        vec[1] += 0.1250*(velo[elemConn[ee][4]*ndim+1]+velo[elemConn[ee][5]*ndim+1]+velo[elemConn[ee][6]*ndim+1]+velo[elemConn[ee][7]*ndim+1]);

        fout << vec[0] <<  '\t' << vec[1] << '\t' << vec[2] << endl;
      }
    }

    fout << "SCALARS pressure float 1" << endl;
    fout << "LOOKUP_TABLE default" << endl;

    for(ii=0; ii<nNode; ii++)
    {
      if( midNodeData[ii][0] )
      {
        // pressure at the mid-nodes is computed as the average
        // only for the purpose of plotting
        val = 0.5*(pres[midNodeData[ii][1]] + pres[midNodeData[ii][2]]) ;
      }
      else
      {
        val    = pres[ii];
      }
      fout << val << endl;
    }

    // pressure at the face node for the 9-noded quadrilateral element
    if(npElemVelo == 9)
    {
      for(ee=0; ee<nElem; ee++)
      {
        val = 0.25*(pres[elemConn[ee][0]]+pres[elemConn[ee][1]]+pres[elemConn[ee][2]]+pres[elemConn[ee][3]]);

        fout << val << endl;
      }
    }

    fout.close();

    fileCount = fileCount+1;

    return;
}




