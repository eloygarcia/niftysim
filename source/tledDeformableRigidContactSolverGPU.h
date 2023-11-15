// =========================================================================
// File:       tledDeformableRigidContactSolverGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidContactSolverGPU_H
#define tledDeformableRigidContactSolverGPU_H

#ifdef _GPU_
#include "tledDeformableRigidContactSolver.h"
#include "tledUnstructuredContactManager.h"
#include "tledContactSolverGPU.h"

/**
 * \brief Deformable-rigid contact handler.
 * \ingroup contact
 */
class tledDeformableRigidContactSolverGPU : public tledDeformableRigidContactSolver, public tledContactSolverGPU {
  /**
   * \name Instantiation
   * @{
   */
public:
  static tledDeformableRigidContactSolverGPU* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDeformableRigidContactSolverGPU(void) {}
  
public:  
  virtual ~tledDeformableRigidContactSolverGPU(void) {}
  /** @} */
};

#endif
#endif
