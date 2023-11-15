// =========================================================================
// File:       tledBottomUpBVHUpdaterGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBottomUpBVHUpdaterGPU_kernels_CU
#define tledBottomUpBVHUpdaterGPU_kernels_CU
#include "tledBV_kernels.h"
#include "tledCUDAHelpers.h"

#include <thrust/functional.h>

namespace tledBottomUpBVHUpdaterGPU_kernels {
  template <class TBoundingVolume, class TMesh>
  __global__ void RefitLeafLevelKernel(TBoundingVolume *p_bvs, const int *leafInds, const int2 *indexBounds, const TMesh *pc_mesh, const float margin) {
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    const int2 currBounds = indexBounds[0];

    if (currBounds.x + tid < currBounds.y) {
      tledCudaAssert(p_bvs[leafInds[currBounds.x+tid]].PrimitiveIndex >= 0);
      tledBV_kernels::RefitLeaf(p_bvs[leafInds[currBounds.x+tid]], *pc_mesh, margin);
    }
  }

  template <class TBoundingVolume>
  __global__ void RefitInteriorLevelKernel(TBoundingVolume *p_bvs, const int *levels, const int2 *levelIndexBounds, const int levelInd) {
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    const int2 currLevel = levelIndexBounds[levelInd];
    
    if (tid + currLevel.x < currLevel.y) {
      const int bvInd = levels[currLevel.x+tid];

      tledCudaAssert(p_bvs[bvInd].PrimitiveIndex < 0);
      tledBV_kernels::RefitInterior(p_bvs, bvInd);      
    }        
  }
}

#endif
