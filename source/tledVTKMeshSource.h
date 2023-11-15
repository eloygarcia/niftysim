// =========================================================================
// File:       tledVTKMeshSource.h
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
#ifndef tledVTKMeshSource_H
#define tledVTKMeshSource_H
#ifdef _Visualisation_

#include "tledMesh.h"
#include "tledGenericVTKMeshSource.h"
#include "tledHelper.h"

typedef tledGenericVTKMeshSource<tledMesh, vtkUnstructuredGrid> tledVTKMeshSource;

#else

class tledVTKMeshSource {};

#endif
#endif
