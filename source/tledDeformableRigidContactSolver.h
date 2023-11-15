// =========================================================================
// File:       tledDeformableRigidContactSolver.h
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

#ifndef tledDeformableRigidContactSolver_H
#define tledDeformableRigidContactSolver_H

#include "tledModel.h"
#include "tledContactSolver.h"
#include "tledContactSolverCPU.h"
#include "tledUnstructuredContactManager.h"
#include "tledDeformableRigidBVHTraverserCPU.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

class tledDeformableRigidContactSolver : public tledContactSolver {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSolver Superclass;
  /** @} */

   /**
   * \name Associated Rigid Surface
   * @{
   */
private:
  int m_SurfaceIndex;

public:
  void SetRigidSurfaceIndex(const int sInd) { m_SurfaceIndex = sInd; } 
  int GetRigidSurfaceIndex(void) const { return m_SurfaceIndex; } 
  /** @} */  

  /**
   * \name Construction and Destruction
   * @{
   */
public:
  static tledDeformableRigidContactSolver* CreateContactSolver(tledUnstructuredContactManager &r_manager, const int surfaceIndex);

  tledDeformableRigidContactSolver(void) {}
  virtual ~tledDeformableRigidContactSolver(void) {}
  /** @} */
};

#endif
