// =========================================================================
// File:       tledSolverGPU_kernels.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSolverGPU_kernels_H
#define tledSolverGPU_kernels_H
#ifdef _GPU_

/**
 * \brief GPU solver functions
 * \ingroup solver
 */
namespace tledSolverGPU_kernels {
  __device__ float4 GetCurrentDisplacement(const int nodeIndex);
}
#endif
#endif
