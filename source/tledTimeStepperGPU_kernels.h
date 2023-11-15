// =========================================================================
// File:       tledTimeStepperGPU_kernels.h
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
#ifndef tledTimeStepperGPU_kernels_H
#define tledTimeStepperGPU_kernels_H

/**
 * \brief API for GPU on-device time integration operations
 */
namespace tledTimeStepper_kernels {
  /** Evolves the displacement for one specific node */
  template <class TSolutionVariables>
  __device__ float3 ComputeNodalDisplacement(TSolutionVariables *p_sols, const int nodeInd, const float3 &nodeEffectiveF);
}
#endif
