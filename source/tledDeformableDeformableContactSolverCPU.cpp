// =========================================================================
// File:       tledDeformableDeformableContactSolverCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledDeformableContactSurfaceCPU.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledDeformableDeformableContactSolverCPU.h"
#include "tledHelper.h"

tledDeformableDeformableContactSolverCPU* tledDeformableDeformableContactSolverCPU::CreateContactSolver(tledUnstructuredContactManager &r_manager) {
  if (r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() != 3) {
    tledFatalNotYetImplementedError;
  }

  if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_manager.GetDeformableSurface<tledSurface>()) == NULL) {
    return new tledDeformableDeformableContactSolverImplCPU<tledDeformableContactSurfaceT3CPU>(r_manager);
  } else {
    return new tledDeformableDeformableContactSolverImplCPU<tledDeformableMembraneContactSurfaceT3CPU>(r_manager);
  }
}
