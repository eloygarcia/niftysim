// =========================================================================
// File:       tledNewmarkTimeStepperCPU.h
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
#ifndef tledNewmarkTimeStepperCPU_H
#define tledNewmarkTimeStepperCPU_H

#include "tledTimeStepperCPU.h"
#include "tledNewmarkTimeStepper.h"

/**
 * \brief CPU Newmark time integration
 * \ingroup solver
 */
class tledNewmarkTimeStepperCPU : public tledNewmarkTimeStepper<tledTimeStepperCPU> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledNewmarkTimeStepper<tledTimeStepperCPU> Superclass;
  /** @} */

  /**
   * \name Displacement Evolution
   * @{
   */
private:
  const float *mpc_M;  
  const float mc_Alpha;
  float *mp_PreviousAccelerations, *mp_PreviousVelocities;

public:
  const float* GetPreviousVelocities(void) const { return mp_PreviousVelocities; }
  const float* GetPreviousAccelerations(void) const { return mp_PreviousAccelerations; }

  virtual void EvolveDisplacements(const float effectiveF[]);
  virtual void FinishTimeStep(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Requires M to be persistent */
  tledNewmarkTimeStepperCPU(const int numNodes, const float dt, const float alpha, const float M[]);
  virtual ~tledNewmarkTimeStepperCPU(void);
  /** @} */
};

#endif
