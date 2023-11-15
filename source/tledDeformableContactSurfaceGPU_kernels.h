// =========================================================================
// File:       tledDeformableContactSurfaceGPU_kernels.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableContactSurfaceGPU_kernels_H
#define tledDeformableContactSurfaceGPU_kernels_H

namespace tledDeformableContactSurfaceGPU_kernels {
  template <class TGPUSurface>
  __device__ bool IsAdjacent(const TGPUSurface &mesh, const int f0Ind, const int f1Ind);

  template <class TGPUSurface>
  __device__ float3 ComputeOldNormal(const TGPUSurface &mesh, const int fInd);

  template <class TGPUSurface>
  __device__ float3 GetNodeNormal(const int nInd, const TGPUSurface &mesh);
}
#endif
