// =========================================================================
// File:       tledContactSurfaceGPU_kernels.h
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
#ifndef tledContactSurfaceGPU_kernels_H
#define tledContactSurfaceGPU_kernels_H
namespace tledContactSurfaceGPU_kernels {
  template <class TSurface>
  __device__ void ComputeShapeValues(float *p_out, const float3 &xi);

  template <class TSurface>
  __device__ void ComputeProjectionOperator(float4 *p_projOp, const TSurface &mesh, const int fInd);
  __device__ bool ProjectOntoFacet(float3 &r_xi, const float3 &x, const float4 projOp[]);
  
  template <class TSurface>
  __device__ const float3& GetNodeNormal(const TSurface &mesh, const int nodeIndex);
};
#endif
