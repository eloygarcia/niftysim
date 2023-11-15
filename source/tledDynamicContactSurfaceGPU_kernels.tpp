// =========================================================================
// File:       tledDynamicContactSurfaceGPU_kernels.tpp
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
#ifndef tledDynamicContactSurfaceGPU_kernels_TPP
#define tledDynamicContactSurfaceGPU_kernels_TPP

#include "tledCUDA_operators.cu"

namespace tledDynamicContactSurfaceGPU_kernels {
  template <class TSurface>
  __device__ float3 ComputeOldNormal(const TSurface &mesh, const int fInd) {
    float3 n;

    n = cross(mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[1]] - mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[0]], mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[2]] - mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[0]]);    

    return n/norm(n);
  }
}

#endif
