// =========================================================================
// File:       tledVTKSurfaceExporter.h
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
#ifndef tledVTKSurfaceExporter_H
#define tledVTKSurfaceExporter_H

#include "tledGenericVTKMeshExporter.h"
#include "tledSurface.h"

#ifdef _Visualisation_
#include <vtkPolyData.h>

/**
 * \brief VTK F/S exporter for tledSurface-type meshes
 * \ingroup export
 *
 *
 * Can be used for writing unstructured contact surfaces and membranes to VTK files.
 */
typedef tledGenericVTKMeshExporter<tledSurface, vtkPolyData> tledVTKSurfaceExporter;
#else
typedef tledGenericVTKMeshExporter<tledSurface> tledVTKSurfaceExporter;
#endif
#endif
