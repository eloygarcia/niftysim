// =========================================================================
// File:       tledDeformableMovingRigidContactSolverGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableMovingRigidContactSolverGPU_CU
#define tledDeformableMovingRigidContactSolverGPU_CU

#include "tledDeformableMovingRigidContactSolverGPU.h"
#include "tledDeformableMovingRigidContactSolverImplGPU.h"

tledDeformableRigidContactSolverGPU* tledDeformableMovingRigidContactSolverGPU::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  tledDeformableRigidContactSolverGPU *p_solver = NULL;

  if (r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() == 3 && r_manager.GetRigidSurface<tledRigidContactSurface>(surfaceIndex).GetNumberOfFacetVertices() == 3) {
    p_solver = new tledDeformableMovingRigidContactSolverImplGPU<tledDeformableContactSurfaceT3GPU, tledMovingRigidContactSurfaceT3GPU>(r_manager, surfaceIndex);
  } else {
    tledFatalError("Setup unsupported.");
  }

  return p_solver;  
}

#endif
