// =========================================================================
// File:       tledDeformableDeformableContactSolverGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableDeformableContactSolverGPU_H
#define tledDeformableDeformableContactSolverGPU_H

#include "tledDeformableDeformableContactSolver.h"
#include "tledUnstructuredContactManager.h"
#include "tledContactSolverGPU.h"

/**
 * \brief Deformable-deformable contact handler.
 * \ingroup contact
 */
class tledDeformableDeformableContactSolverGPU : public tledDeformableDeformableContactSolver, public tledContactSolverGPU {
  /**
   * \name Instantiation
   * @{
   */
public:
  static tledDeformableDeformableContactSolver* CreateContactSolver(tledUnstructuredContactManager &r_manager);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDeformableDeformableContactSolverGPU(void) {}
  
public:
  virtual ~tledDeformableDeformableContactSolverGPU(void) {}
  /** @} */
};

#endif
