// =========================================================================
// File:       tledVTKSurfaceSource.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledVTKSurfaceSource.h"

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetVTKElementType(void) const {
  return this->GetInput().GetNumberOfFacetVertices() == 3? VTK_TRIANGLE : VTK_QUAD;
}

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetNumberOfPoints(void) const {
  return this->GetInput().GetNumberOfNodes();
}

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetNumberOfElements(void) const {
  return this->GetInput().GetNumberOfFacets();
}

template <>
const float* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetPointCoordinates(const int ptInd) const {
  return this->GetInput().GetNodeCoordinates(ptInd);
}

template <>
const int* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetElementBegin(const int elementIndex) const {
  return this->GetInput().GetFacetNodeIndices(elementIndex);
}

template <>
const int* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetElementEnd(const int elementIndex) const {
  return this->GetInput().GetFacetNodeIndices(elementIndex) + this->GetInput().GetNumberOfFacetVertices();
}
