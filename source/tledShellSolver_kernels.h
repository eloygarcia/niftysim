// =========================================================================
// File:       tledShellSolver_kernels.h
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
#ifndef tledShellSolver_kernels_H
#define tledShellSolver_kernels_H

/**
 * \brief Namespace of GPU shell/membrane functions
 * \ingroup shell
 */
namespace tledShellSolver_kernels {
  /**
   * Entry point for GPU shell/membrane internal force computation
   * \ingroup shell
   */
  template <class TElement, class TMaterial>
  __global__ void ComputeNewForces(float4 *gp_F, typename TElement::GPUElement *gp_elements, const typename TMaterial::GPUMaterial *gpc_materialParameters, const int numElements);

  __global__ void ComputeNodalForces(float4 *gp_f, const int numberOfNodes);
}

#endif
