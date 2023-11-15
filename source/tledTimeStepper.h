// =========================================================================
// File:       tledTimeStepper.h
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
#ifndef tledTimeStepper_H
#define tledTimeStepper_H

#include "tledHelper.h"

/**
 * \brief Base class for time ODE integrators.
 * \ingroup solver
 */
class tledTimeStepper {
  /**
   * \name Solution Access
   * @{
   */
private:
  const int mc_NumNodes;
  float *mp_CurrentDisplacements, *mp_NextDisplacements;  

protected:
  /** 
   * If inBuffer == NULL, Current and Next buffers are simply swapped and NULL is returned.<br />
   * Otherwise, Next -> Current, inBuffer -> Next, pointer to previous Current is returned.
   */
  virtual float* RotateDisplacementBuffers(float *p_inBuffer);

public:
  int GetNumberOfNodes(void) const { return mc_NumNodes; }
  
  float* GetCurrentDisplacements(void) { return mp_CurrentDisplacements; }
  const float* GetCurrentDisplacements(void) const { return mp_CurrentDisplacements; }

  float* GetNextDisplacements(void) { return mp_NextDisplacements; }
  const float* GetNextDisplacements(void) const { return mp_NextDisplacements; }

  virtual void SetCurrentDisplacements(const float U[]) = 0;
  /** @} */
  
  /**
   * \name Displacement Evolution and Related Computation
   * @{
   */
public:
  /** Signals the end of a time step */
  virtual void FinishTimeStep(void) = 0;

  /** Computes the nodal displacement increments for output purposes */
  virtual float* GetCurrentDeltaDisplacements(float *p_du) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledTimeStepper(const int numNodes);
  virtual ~tledTimeStepper(void);
  /** @} */
};

#endif
