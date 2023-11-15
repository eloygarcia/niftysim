// =========================================================================
// File:       tledVTKMeshExporter.h
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
#ifndef tledVTKMeshExporter_H
#define tledVTKMeshExporter_H

#include "tledGenericVTKMeshExporter.h"

#ifdef _Visualisation_
#include <vtkUnstructuredGrid.h>

/**
 * \brief tledMesh-Exporter for writing of VTK unstructured grids
 * \ingroup export
 *
 *
 * Allows for making of more sophisticated animations in Paraview.
 */
typedef tledGenericVTKMeshExporter<tledMesh, vtkUnstructuredGrid> tledVTKMeshExporter;
#else
typedef tledGenericVTKMeshExporter<tledMesh> tledVTKMeshExporter;
#endif
#endif
