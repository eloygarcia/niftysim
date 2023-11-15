// =========================================================================
// File:       tledVTKMeshSource.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledVTKMeshSource.h"
#ifdef _Visualisation_

template <>
const float* tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetPointCoordinates(const int ptInd) const {
  return this->GetInput().GetAllNodeCds() + 3*ptInd;
}

template <>
int tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetNumberOfPoints() const {
  return this->GetInput().GetNumNodes();
}

template <>
const int* tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetElementBegin(const int elementIndex) const {
  return this->GetInput().GetAllElNodeInds() + this->GetInput().GetNodesPerEl()*elementIndex;
}

template <>
const int* tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetElementEnd(const int elementIndex) const {
  return this->GetElementBegin(elementIndex) + this->GetInput().GetNodesPerEl();
}

template <>
int tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetNumberOfElements() const {
  return this->GetInput().GetNumEls();
}

template <>
int tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::GetVTKElementType(void) const {
  return (std::string(this->GetInput().GetElType()) == "T4" || std::string(this->GetInput().GetElType()) == "T4ANP")? VTK_TETRA : VTK_HEXAHEDRON == 3? VTK_TETRA : VTK_HEXAHEDRON;
}

template <>
void tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid>::Update();
#endif
