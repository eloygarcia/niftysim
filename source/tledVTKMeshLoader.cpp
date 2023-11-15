// =========================================================================
// File:       tledVTKMeshLoader.cpp
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
#include "tledVTKMeshLoader.h"

#include <iostream>
#include <algorithm>
#include <cassert>
#ifdef _Visualisation_
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCell.h>
#include <vtkSmartPointer.h>
#include <vtkErrorCode.h>
#endif

using namespace std;

tledVTKMeshLoader::tledVTKMeshLoader(void) {
  SetNumberOfDOFs(3);
}

#ifdef _Visualisation_
static inline float _CastDouble2Float(const double &x) {
  return float(x);
}

void tledVTKMeshLoader::ReadFile() {
  vtkSmartPointer<vtkUnstructuredGrid> sp_vtkMesh;

  {
    vtkSmartPointer<vtkUnstructuredGridReader> sp_reader;

    sp_reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    sp_reader->SetFileName(this->GetFilename().c_str());
    sp_reader->Update();

    if (sp_reader->GetErrorCode()) {
      tledLogErrorStream(tledHelper::FatalError() << "VTK error:\n" << vtkErrorCode::GetStringFromErrorCode(sp_reader->GetErrorCode()));
    }

    sp_vtkMesh = sp_reader->GetOutput();
  }

  this->GetOutput().SetNumberOfNodes(sp_vtkMesh->GetNumberOfPoints(), 3);
  for (int pInd = 0; pInd < sp_vtkMesh->GetNumberOfPoints(); pInd++) {    
    std::transform(sp_vtkMesh->GetPoint(pInd), sp_vtkMesh->GetPoint(pInd) + 3, this->GetOutput().GetAllNodeCds() + 3*pInd, _CastDouble2Float);
  }

  {
    std::vector<int> elements;      

    elements.reserve(sp_vtkMesh->GetNumberOfCells());
    for (int eInd = 0; eInd < sp_vtkMesh->GetNumberOfCells(); eInd++) {
      if (sp_vtkMesh->GetCell(eInd)->GetNumberOfPoints() == this->GetNumberOfElementNodes()) {
	for (int vtkCNodeInd = 0; vtkCNodeInd < this->GetNumberOfElementNodes(); vtkCNodeInd++) elements.push_back(sp_vtkMesh->GetCell(eInd)->GetPointId(vtkCNodeInd));
      }	  	
    }
 
    if (elements.size() == 0) {
      tledLogErrorStream(tledHelper::FatalError() << "No " << this->GetMeshType() << " elements found in file " << this->GetFilename());
    }    

    this->GetOutput().SetNumberOfElements(elements.size()/this->GetNumberOfElementNodes(), this->GetMeshType());
    std::copy(elements.begin(), elements.end(), this->GetOutput().GetAllElNodeInds());
  }  
}
#else
void tledVTKMeshLoader::ReadFile() {
  tledFatalFeatureNotEnabledError;
} 
#endif /* __TLED_VTK_MESH_LOADING */
