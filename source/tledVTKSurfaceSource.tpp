// =========================================================================
// File:       tledVTKSurfaceSource.tpp
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

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetVTKElementType(void) const;

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetNumberOfPoints(void) const;

template <>
int tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetNumberOfElements(void) const;

template <>
const float* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetPointCoordinates(const int ptInd) const;

template <>
const int* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetElementBegin(const int elementIndex) const;

template <>
const int* tledGenericVTKMeshSource<tledSurface, vtkPolyData>::GetElementEnd(const int elementIndex) const;
