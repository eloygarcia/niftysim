// =========================================================================
// File:       tledDeformableRigidContactSolver.cpp
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

#include "tledDeformableRigidContactSolver.h"
#include "tledDeformableRigidContactSolverCPU.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledDeformableMovingRigidContactSolverCPU.h"
#include "tledHelper.h"
#ifdef GPU_GP_CONTACT
#include "tledDeformableRigidContactSolverGPU.h"
#include "tledDeformableMovingRigidContactSolverGPU.h"
#endif

#include <iostream>
#include <cstdlib>

tledDeformableRigidContactSolver* tledDeformableRigidContactSolver::CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex) {
  if (r_manager.GetRigidSurface<tledRigidContactSurfaceCPU>(surfaceIndex).IsMoving()) {
    if (!r_manager.UseGPU()) {
      return tledDeformableMovingRigidContactSolverCPU::CreateContactSolver(r_manager, surfaceIndex);
    } else {
#ifdef GPU_GP_CONTACT
      return tledDeformableMovingRigidContactSolverGPU::CreateContactSolver(r_manager, surfaceIndex);
#else
      tledFatalFeatureNotEnabledError;
#endif      
    } 
  } else {
    if (r_manager.UseGPU()) {
#ifdef GPU_GP_CONTACT
      return tledDeformableRigidContactSolverGPU::CreateContactSolver(r_manager, surfaceIndex);
#else
      tledFatalFeatureNotEnabledError;
#endif
    } else {
      return tledDeformableRigidContactSolverCPU::CreateContactSolver(r_manager, surfaceIndex);
    }  
  }
  tledFatalError("Setup not supported.");

  return NULL;
}
