// =========================================================================
// File:       tledTimeStepperCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTimeStepperCPU_H
#define tledTimeStepperCPU_H

#include "tledTimeStepper.h"

/**
 * \brief CPU Time ODE solver
 * \ingroup solver
 */
class tledTimeStepperCPU : public tledTimeStepper {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledTimeStepper Superclass;
  /** @} */

  /**
   * \name Solutions
   * @{
   */
public:
  virtual void SetCurrentDisplacements(const float U[]);
  /** @} */

  /**
   * \name Displacement Evolution and Related Computation
   * @{
   */
public:
  /** 
   * \brief Computes the next time step displacements. 
   *
   * Solution is not committed until FinishTimeStep is called.
   */
  virtual void EvolveDisplacements(const float effectiveF[]) = 0;

  /** Signals the end of a time step */
  virtual void FinishTimeStep(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledTimeStepperCPU(const int numNodes);
  virtual ~tledTimeStepperCPU(void) {}
  /** @} */
};

#endif
