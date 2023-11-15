// =========================================================================
// File:       tledDeformableMovingRigidContactSolverCPU.cpp
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
#include "tledDeformableMovingRigidContactSolverCPU.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledHelper.h"

tledDeformableRigidContactSolverCPU* tledDeformableMovingRigidContactSolverCPU::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  const tledMovingRigidContactSurfaceCPU &rigSurf = r_manager.GetRigidSurface<tledMovingRigidContactSurfaceCPU>(surfaceIndex);

  if (rigSurf.GetNumberOfFacetVertices() != 3 || r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() != 3) {
    tledFatalNotYetImplementedError;
  }

  if (dynamic_cast<tledDeformableMembraneContactSurfaceCPU*>(&r_manager.GetDeformableSurface<tledSurface>()) == NULL) {
    return new tledDeformableMovingRigidContactSolverImplCPU<tledDeformableContactSurfaceT3CPU, tledMovingRigidContactSurfaceT3CPU>(r_manager, surfaceIndex);
  } else {
     return new tledDeformableMovingRigidContactSolverImplCPU<tledDeformableMembraneContactSurfaceT3CPU, tledMovingRigidContactSurfaceT3CPU>(r_manager, surfaceIndex);
  }
}
