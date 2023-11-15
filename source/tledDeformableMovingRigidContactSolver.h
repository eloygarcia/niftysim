// =========================================================================
// File:       tledDeformableMovingRigidContactSolver.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableMovingRigidContactSolver_H
#define tledDeformableMovingRigidContactSolver_H

#include "tledMovingRigidContactSurface.h"
#include "tledDeformableRigidContactSolver.h"

class tledDeformableMovingRigidContactSolver : public tledDeformableRigidContactSolver {  
  /**
   * \name Construction and Destruction
   * @{
   */
public:
  static tledDeformableRigidContactSolver* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);

  tledDeformableMovingRigidContactSolver(void) {}
  virtual ~tledDeformableMovingRigidContactSolver(void) {}
  /** @} */
};

#endif
