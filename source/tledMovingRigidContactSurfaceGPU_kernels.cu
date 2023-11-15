// =========================================================================
// File:       tledMovingRigidContactSurfaceGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMovingRigidContactSurfaceGPU_kernels_CU
#define tledMovingRigidContactSurfaceGPU_kernels_CU

namespace tledMovingRigidContactSurfaceGPU_kernels {
  __global__ void TranslateNodesKernel(float3 *p_currCds, float3 *p_oldCds, const float3 *x0s, const int numNodes, const float3 tCurr, const float3 tOld) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    if (tid < numNodes) {
      p_currCds[tid] = x0s[tid] + tCurr;
      p_oldCds[tid] = x0s[tid] + tOld;
    }
  } 
}

#endif
