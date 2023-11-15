// =========================================================================
// File:       tledDeformableDeformableContactSolver.h
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

#ifndef tledDeformableDeformableContactSolver_H
#define tledDeformableDeformableContactSolver_H

#include "tledContactSolver.h"
#include "tledVectorArithmetic.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>

/**
 * \brief Solver for deformable-deformable contacts (including self-collisions)
 * \ingroup contact
 */
class tledDeformableDeformableContactSolver : public tledContactSolver {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSolver Superclass;
  /** @} */

  /**
   * \name Contact Search Parameters
   * @{
   */
private:
  bool m_DoSelfCollision;

public:
  /** 
   * \brief Detection and handling of self collision, i.e. collisions of topologically connected geometry.
   *
   * Turned on by default.
   */
  bool DoSelfCollision(void) const { return m_DoSelfCollision; }

  void ToggleDoSelfCollision(void) { m_DoSelfCollision = !m_DoSelfCollision; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDeformableDeformableContactSolver* CreateContactSolver(tledUnstructuredContactManager &r_manager);

  tledDeformableDeformableContactSolver(void) : m_DoSelfCollision(true) {}
  virtual ~tledDeformableDeformableContactSolver(void) {}
  /** @} */  
}; 

#endif
