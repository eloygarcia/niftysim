// =========================================================================
// File:       tledNewmarkTimeStepper.h
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
#ifndef tledNewmarkTimeStepper_H
#define tledNewmarkTimeStepper_H

#include "tledTimeStepper.h"

#include <cassert>

/**
 * \brief Explicit Newmark ODE solver (CPU, GPU base class). 
 * \ingroup solver
 */
template <class TTimeStepperBase>
class tledNewmarkTimeStepper : public TTimeStepperBase {
  /**
   * \name Types
   * @{
   */
public:
  typedef TTimeStepperBase Superclass;
  /** @} */

  /**
   * \name Computation
   * @{
   */
private:
  const float mc_Dt;

public:
  float GetTimeStep(void) const { return mc_Dt; }
  virtual float* GetCurrentDeltaDisplacements(float *p_du) const;
  /** @} */

  /**
   * \name Solutions
   * @{
   */
private:
  float *mp_CurrentVelocities, *mp_CurrentAccelerations;

protected:
  float* ReplaceCurrentVelocityBuffer(float *p_newBuffer);
  float* ReplaceCurrentAccelerationBuffer(float *p_newBuffer);

  float* GetCurrentVelocities(void) { return mp_CurrentVelocities; }
  float* GetCurrentAccelerations(void) { return mp_CurrentAccelerations; }

public:
  const float* GetCurrentVelocities(void) const { return mp_CurrentVelocities; }
  const float* GetCurrentAccelerations(void) const { return mp_CurrentAccelerations; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledNewmarkTimeStepper(const int numNodes, const float dt);
  virtual ~tledNewmarkTimeStepper(void);
  /** @} */
};

#include "tledNewmarkTimeStepper.tpp"
#endif
