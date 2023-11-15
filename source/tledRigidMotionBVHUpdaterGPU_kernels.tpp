// =========================================================================
// File:       tledRigidMotionBVHUpdaterGPU_kernels.cu
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
#ifndef tledRigidMotionBVHUpdaterGPU_kernels_CU
#define tledRigidMotionBVHUpdaterGPU_kernels_CU

#include "tledBV_kernels.h"

namespace tledRigidMotionBVHUpdaterGPU_kernels {
  template <class TBV>
  __global__ void TranslateBVHKernel(TBV *p_bvs, const int numBVs, const float3 t) {
    const int tid = threadIdx.x + blockDim.x*blockIdx.x;

    if (tid < numBVs) {
      tledBV_kernels::Translate(p_bvs[tid], t);
    }
  }
}

#endif
