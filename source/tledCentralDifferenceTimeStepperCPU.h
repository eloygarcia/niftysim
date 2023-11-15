// =========================================================================
// File:       tledCentralDifferenceTimeStepperCPU.h
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
#ifndef tledCentralDifferenceTimeStepperCPU_H
#define tledCentralDifferenceTimeStepperCPU_H

#include "tledTimeStepperCPU.h"
#include "tledCentralDifferenceTimeStepper.h"

/**
 * \brief CPU central difference time integration
 * \ingroup solver
 */
class tledCentralDifferenceTimeStepperCPU : public tledCentralDifferenceTimeStepper<tledTimeStepperCPU> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledCentralDifferenceTimeStepper<tledTimeStepperCPU> Superclass;
  /** @} */

  /**
   * \name Displacement Evolution
   * @{
   */
private:
  float *mp_A, *mp_B, *mp_C;  

public:
  virtual void EvolveDisplacements(const float effectiveF[]);
  virtual void FinishTimeStep(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledCentralDifferenceTimeStepperCPU(const int numNodes, const float dt, const float alpha, const float M[]);
  virtual ~tledCentralDifferenceTimeStepperCPU(void);
  /** @} */
};

#endif
