// =========================================================================
// File:       tledDeformableContactSurfaceGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableContactSurfaceGPU_kernels_CU
#define tledDeformableContactSurfaceGPU_kernels_CU

#include "tledContactSurfaceGPU_kernels.h"
#include "tledDeformableContactSurfaceGPU.h"

namespace tledDeformableContactSurfaceGPU_kernels {
  __global__ void UpdateNodePositionsKernel(float3 *p_surfaceNodes, const float3 *surfaceNodes0, const float4 *uNexts, const int *surfaceToVolumeNodeIndexMap, const int numNodes) {
    const int tInd = threadIdx.x + blockIdx.x*blockDim.x;

    if (tInd < numNodes) {
      const int sInd = tInd;
      const int gInd = surfaceToVolumeNodeIndexMap[sInd];

      p_surfaceNodes[sInd] = surfaceNodes0[sInd] + uNexts[gInd];
    }
  }
}

#endif
