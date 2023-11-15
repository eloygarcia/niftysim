// =========================================================================
// File:       tledParallelSolverCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallelSolverCPU_H
#define tledParallelSolverCPU_H

#include "tledSolverCPU.h"
#include "tledParallel.h"

#ifdef _USE_BOOST_

/**
 * \brief Main TLED solver with parallel CPU execution
 * \ingroup solver
 */
class tledParallelSolverCPU : public tledSolverCPU, public tledParallel {  
  /**
   * \name Displacement Evolution
   */
private:
  float *mp_ThreadResultBuffer;

protected:
  virtual void ComputeNewForces();
  virtual void ComputeNewForcesANP();
  /** @} */

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
public:
  virtual void Init(tledModel* p_model);

  tledParallelSolverCPU(void) : tledParallel(tledParallel::GetNumberOfHardwareThreads()), mp_ThreadResultBuffer(NULL) {}
  virtual ~tledParallelSolverCPU(void);
  /** @} */
};

#endif
#endif
