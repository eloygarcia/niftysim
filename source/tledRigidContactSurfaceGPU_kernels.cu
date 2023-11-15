// =========================================================================
// File:       tledRigidContactSurfaceGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledRigidContactSurfaceGPU_kernels_CU
#define tledRigidContactSurfaceGPU_kernels_CU

namespace tledRigidContactSurfaceGPU_kernels { 
  template <class TSurface>
  __device__ const float4* GetFacetProjectionOperator(const TSurface &mesh, const int facetIndex) {
    return mesh.FacetProjectionOperators + 3*facetIndex;
  }
}
#endif
