// =========================================================================
// File:       tledVTKUnstructuredContactSurfaceSource.h
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
#ifndef tledVTKUnstructuredContactSurfaceSource_H
#define tledVTKUnstructuredContactSurfaceSource_H

#ifdef _Visualisation_

#include "tledVTKSurfaceSource.h"

/** 
 * \brief Converts NiftySim unstructured contact surfaces to VTK PolyData objects.
 * \ingroup surface
 * \ingroup contact
 */
typedef tledVTKSurfaceSource tledVTKUnstructuredContactSurfaceSource;

#else 

class tledVTKUnstructuredContactSurfaceSource {};

#endif

#endif
