
#include "ExplicitCFD.h"
#include "headersVTK.h"


void  ExplicitCFD::postProcess()
{
    //cout << " ExplicitCFD::postProcess " << endl;

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


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

          vec[0] = 0.5*velo(jj)   + 0.25*velo(n1)   + 0.25*velo(n2);
          vec[1] = 0.5*velo(jj+1) + 0.25*velo(n1+1) + 0.25*velo(n2+1);

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          val = 0.5*pres(midNodeData[ii][1]) + 0.5*pres(midNodeData[ii][2]) ;

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
          vec[0]  = 0.25*velo(elemConn[ee][8]*ndim);
          vec[0] += 0.0625*(velo(elemConn[ee][0]*ndim)+velo(elemConn[ee][1]*ndim)+velo(elemConn[ee][2]*ndim)+velo(elemConn[ee][3]*ndim));
          vec[0] += 0.1250*(velo(elemConn[ee][4]*ndim)+velo(elemConn[ee][5]*ndim)+velo(elemConn[ee][6]*ndim)+velo(elemConn[ee][7]*ndim));

          vec[1]  = 0.25*velo(elemConn[ee][8]*ndim+1);
          vec[1] += 0.0625*(velo(elemConn[ee][0]*ndim+1)+velo(elemConn[ee][1]*ndim+1)+velo(elemConn[ee][2]*ndim+1)+velo(elemConn[ee][3]*ndim+1));
          vec[1] += 0.1250*(velo(elemConn[ee][4]*ndim+1)+velo(elemConn[ee][5]*ndim+1)+velo(elemConn[ee][6]*ndim+1)+velo(elemConn[ee][7]*ndim+1));

          vecVTK->InsertTuple(elemConn[ee][8], vec);

          val = 0.25*(pres(elemConn[ee][0])+pres(elemConn[ee][1])+pres(elemConn[ee][2])+pres(elemConn[ee][3]));

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

          vec[0] = 0.5*velo(jj)   + 0.25*velo(n1)   + 0.25*velo(n2);
          vec[1] = 0.5*velo(jj+1) + 0.25*velo(n1+1) + 0.25*velo(n2+1);
          vec[2] = 0.5*velo(jj+2) + 0.25*velo(n1+2) + 0.25*velo(n2+2);

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          val   = 0.5*pres(midNodeData[ii][1]) + 0.5*pres(midNodeData[ii][1]) ;

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
    sprintf(VTKfilename,"%s%d%s", "vtk-output-",fileCount,".vtu");
    writerUGridVTK->SetFileName(VTKfilename);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();


    fileCount = fileCount+1;

    return;
}








