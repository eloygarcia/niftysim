// =========================================================================
// File:       tledDeformableMovingRigidContactSolverGPU.h
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
#ifndef tledDeformableMovingRigidContactSolverGPU_H
#define tledDeformableMovingRigidContactSolverGPU_H

#include "tledDeformableRigidContactSolverGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledDeformableMovingRigidBVHTraverserGPU.h"

/**
 * \name Interface for the modelling of contact between deformable and moving rigid bodies (GPU version).
 * \ingroup contact
 */
class tledDeformableMovingRigidContactSolverGPU : public tledDeformableRigidContactSolverGPU {
  /**
   * \name Construction and Destruction
   * @{
   */
public:
  static tledDeformableRigidContactSolverGPU* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);

  tledDeformableMovingRigidContactSolverGPU(void) {}
  virtual ~tledDeformableMovingRigidContactSolverGPU(void) {}
  /** @} */
};

#endif
