// =========================================================================
// File:       tledVTKUnstructuredContactSurfaceSource.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledVTKUnstructuredContactSurfaceSource.h"

template <>
int tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetVTKElementType(void) const {
  return GetInput().GetNumberOfFacetVertices() == 3? VTK_TRIANGLE : VTK_QUAD;
}

template <>
int tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetNumberOfPoints(void) const {
  return GetInput().GetNumberOfNodes();
}

template <>
int tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetNumberOfElements(void) const {
  return GetInput().GetNumberOfFacets();
}

template <>
const float* tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetPointCoordinates(const int ptInd) const {
  return GetInput().GetNodeCoordinates(ptInd);
}

template <>
const int* tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetElementBegin(const int elementIndex) const {
  return GetInput().GetFacetNodeIndices(elementIndex);
}

template <>
const int* tledGenericVTKMeshSource<tledContactSurface, vtkPolyData>::GetElementEnd(const int elementIndex) const {
  return GetInput().GetFacetNodeIndices(elementIndex) + GetInput().GetNumberOfFacetVertices();
}
