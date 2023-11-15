// =========================================================================
// File:       tledCentralDifferenceTimeStepper.h
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
#ifndef tledCentralDifferenceTimeStepper_H
#define tledCentralDifferenceTimeStepper_H

#include "tledTimeStepper.h"

#include <cassert>

/**
 * \brief Explicit Central Difference ODE solver (CPU, GPU base class). 
 * \ingroup solver
 */
template <class TTimeStepperBase>
class tledCentralDifferenceTimeStepper : public TTimeStepperBase {
  /**
   * \name Types
   * @{
   */
public:
  typedef TTimeStepperBase Superclass;
  /** @} */

  /**
   * \name Solutions
   */
private:
  float *mp_PreviousDisplacements;

protected:
  /** Right rotation of Previous, Current, Next buffers, p_inBuffer != NULL is not supported! */
  virtual float* RotateDisplacementBuffers(float *p_inBuffer);

public:
  float* GetPreviousDisplacements(void) { return mp_PreviousDisplacements; }
  const float* GetPreviousDisplacements(void) const { return mp_PreviousDisplacements; }

  virtual float* GetCurrentDeltaDisplacements(float *p_du) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledCentralDifferenceTimeStepper(const int numNodes);
  virtual ~tledCentralDifferenceTimeStepper(void);
  /** @} */
};

#include "tledCentralDifferenceTimeStepper.tpp"
#endif
