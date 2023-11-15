// =========================================================================
// File:       tledDeformableRigidContactSolverCPU.cpp
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

#include "tledVectorArithmetic.h"
#include "tledHelper.h"
#include "tledDeformableRigidContactSolverCPU.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledDeformableMovingRigidContactSolverCPU.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"

tledDeformableRigidContactSolverCPU* tledDeformableRigidContactSolverCPU::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  if (r_manager.GetRigidSurface<tledRigidContactSurface>(surfaceIndex).GetNumberOfFacetVertices() != 3 || r_manager.GetDeformableSurface<tledSurface>().GetNumberOfFacetVertices() != 3) {
    tledFatalNotYetImplementedError;
  }

  if (dynamic_cast<tledDeformableMembraneContactSurface*>(&r_manager.GetDeformableSurface<tledSurface>()) == NULL) {
    return new tledDeformableRigidContactSolverImplCPU<tledDeformableContactSurfaceT3CPU, tledRigidContactSurfaceT3CPU>(r_manager, surfaceIndex);
  } else {
    return new tledDeformableRigidContactSolverImplCPU<tledDeformableMembraneContactSurfaceT3CPU, tledRigidContactSurfaceT3CPU>(r_manager, surfaceIndex);
  }  
}
