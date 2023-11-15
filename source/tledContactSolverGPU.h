// =========================================================================
// File:       tledContactSolverGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef tledContactSolverGPU_H
#define tledContactSolverGPU_H

#ifdef _GPU_
#include "tledCUDAHelpers.h"
#include "tledContactSolver.h"

/**
 * \brief Interface for GPU contact solvers 
 * \ingroup contact
 */
class tledContactSolverGPU {
  /**
   * \name Contact Response
   * @{
   */
public:
  /** Nodal contact force data structure, holding the target index and force */ 
  struct ContactResponse {
    int NodeIndex;
    float3 Force;
  };

public:
  /**
   * \brief Computes the contact responses based on current, predicted next iteration displacements.
   */
  virtual bool ComputeContactResponses(float4 *dp_R, const float4 *dpc_uNexts, const float4 *dpc_uCurrs) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactSolverGPU(void) {}
  virtual ~tledContactSolverGPU(void) {}
  /** @} */
};

#else

#error "CUDA file included despite CUDA-support not being enabled at compile time."

#endif /* _GPU_ */
#endif
