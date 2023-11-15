// =========================================================================
// File:       tledDeformableDeformableContactSolverGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableDeformableContactSolverGPU_CU
#define tledDeformableDeformableContactSolverGPU_CU

#include "tledDeformableDeformableContactSolverGPU.h"
#include "tledDeformableDeformableContactSolverImplGPU.h"
#include "tledDeformableContactSurfaceGPU.h"

tledDeformableDeformableContactSolver* tledDeformableDeformableContactSolverGPU::CreateContactSolver(tledUnstructuredContactManager &r_manager) {
  tledDeformableDeformableContactSolver *p_solver = NULL;

  if (r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() == 3) {
    p_solver = new tledDeformableDeformableContactSolverImplGPU<tledDeformableContactSurfaceT3GPU>(r_manager);
  } else {
    tledFatalError("Setup unsupported.");
  }

  return p_solver;
}

#endif
