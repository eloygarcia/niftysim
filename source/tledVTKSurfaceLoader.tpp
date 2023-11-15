// =========================================================================
// File:       tledVTKSurfaceLoader.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifdef _Visualisation_
template <class TSurface>
void tledVTKSurfaceLoader<TSurface>::ReadFile() {
  vtkSmartPointer<vtkPolyData> sp_vtkMesh;

  {
    vtkSmartPointer<vtkPolyDataReader> sp_reader;

    sp_reader = vtkSmartPointer<vtkPolyDataReader>::New();
    sp_reader->SetFileName(this->GetFilename().c_str());

    sp_reader->Update();
    if (sp_reader->GetErrorCode()) {
      tledLogErrorStream(tledHelper::FatalError() << "VTK error:\n" << vtkErrorCode::GetStringFromErrorCode(sp_reader->GetErrorCode()));
    }

    if (TSurface::Facet::NumberOfVertices == 3) {
      vtkSmartPointer<vtkTriangleFilter> sp_triangulator;

      sp_triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
      tledVTK6CompatSetInput(sp_triangulator, sp_reader->GetOutput());
      sp_triangulator->Update();

      if (sp_triangulator->GetErrorCode()) {
	tledLogErrorStream(tledHelper::FatalError() << "VTK error:\n" << vtkErrorCode::GetStringFromErrorCode(sp_triangulator->GetErrorCode()));
      }

      sp_vtkMesh = sp_triangulator->GetOutput();
    } else {
      sp_vtkMesh = sp_reader->GetOutput();
    }
  }

  {
    float *p_nodeCd;
    int pInd, cInd;

    this->GetOutput().SetNumberOfNodes(sp_vtkMesh->GetNumberOfPoints());
    for (pInd = 0, p_nodeCd = this->GetOutput().GetAllNodeCoordinates(); pInd < sp_vtkMesh->GetNumberOfPoints(); pInd++) {
      for (cInd = 0; cInd < 3; cInd++, p_nodeCd++) *p_nodeCd = sp_vtkMesh->GetPoint(pInd)[cInd];
    }
  } /* if do node merging else ... */

  {
    const int numCellNodes = Surface::Facet::NumberOfVertices;

    int numFacets = 0;

    for (int eInd = 0; eInd < sp_vtkMesh->GetNumberOfCells(); eInd++) numFacets += (sp_vtkMesh->GetCell(eInd)->GetNumberOfPoints() == numCellNodes);
    this->GetOutput().SetNumberOfFacets(numFacets);
    for (int eInd = 0, fInd = 0; eInd < sp_vtkMesh->GetNumberOfCells(); eInd++) {
      if (sp_vtkMesh->GetCell(eInd)->GetNumberOfPoints() == numCellNodes) {
	int vtkCNodeInd;
	  
	for (vtkCNodeInd = 0; vtkCNodeInd < numCellNodes; vtkCNodeInd++) this->GetOutput().GetFacet(fInd).NodeIndices[vtkCNodeInd] = sp_vtkMesh->GetCell(eInd)->GetPointId(vtkCNodeInd);
	fInd++;
      }	  	
    }
 
    if (this->GetOutput().GetAllFacets().size() == 0) {
      tledLogErrorStream(tledHelper::FatalError() << "No " << (numCellNodes == 3? "triangle" : numCellNodes == 4? "quad" : "unsupported") << " elements found in file " << this->GetFilename());
    }    
  }  
}
#else 
template <class TSurface>
void tledVTKSurfaceLoader<TSurface>::ReadFile() {
  tledFatalFeatureNotEnabledError;
}
#endif
