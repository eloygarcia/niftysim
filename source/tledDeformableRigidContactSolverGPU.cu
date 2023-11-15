// =========================================================================
// File:       tledDeformableRigidContactSolverGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidContactSolverGPU_CU
#define tledDeformableRigidContactSolverGPU_CU

#include "tledDeformableRigidContactSolverGPU.h"
#include "tledDeformableRigidContactSolverImplGPU.h"

tledDeformableRigidContactSolverGPU* tledDeformableRigidContactSolverGPU::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  tledDeformableRigidContactSolverGPU *p_solver = NULL;

  if (r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() == 3 && r_manager.GetRigidSurface<tledRigidContactSurface>(surfaceIndex).GetNumberOfFacetVertices() == 3) {
    p_solver = new tledDeformableRigidContactSolverImplGPU<tledDeformableContactSurfaceT3GPU, tledRigidContactSurfaceT3GPU>(r_manager, surfaceIndex);
  } else {
    tledFatalError("Setup unsupported.");
  }

  return p_solver;
}


#endif
