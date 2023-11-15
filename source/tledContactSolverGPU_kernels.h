// =========================================================================
// File:       tledContactSolverGPU_kernels.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSolverGPU_kernels_H
#define tledContactSolverGPU_kernels_H
namespace tledContactSolverGPU_kernels {
  struct ContactResponseUnion {
    __device__ tledContactSolverGPU::ContactResponse operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const;
  };

  struct ContactResponseOrdering {
    __device__ bool operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const;
  };
  
  struct ContactResponseSameNodePredicate {
    __device__ bool operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const;
  };

  template <class TProjection, class TXFormOperator>
  __global__ void ProjectionToResponseTransformKernel(tledContactSolverGPU::ContactResponse *p_responses, const TProjection *projections, const int numProjections, const TXFormOperator xFormOp, const int numResponsesPerProjection);

  __global__ void ConvertResponsesToExternalForce(float4 *p_R, const tledContactSolverGPU::ContactResponse *responses, const int numResponses, const int *dpc_surface2VolumeMap);
  __global__ void ComputeNodeVelocityKernel(float3 *dp_v, const float4 *uNexts, const float4 *uCurrs, const int *surface2VolumeMap, const int numNodes);
};
#endif
