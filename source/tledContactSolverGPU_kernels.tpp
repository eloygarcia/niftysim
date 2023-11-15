// =========================================================================
// File:       tledContactSolverGPU_kernels.tpp
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
#ifndef tledContactSolverGPU_kernels_TPP
#define tledContactSolverGPU_kernels_TPP

#include "tledContactSolverGPU_kernels.h"
#include "tledContactSolverGPU.h"
#include "tledCUDA_operators.cu"

namespace tledContactSolverGPU_kernels {
  __device__ tledContactSolverGPU::ContactResponse ContactResponseUnion::operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const {
    tledContactSolverGPU::ContactResponse r01;
    float nf, pf;
    
    r01.Force = r0.Force;
    r01.NodeIndex = r0.NodeIndex;
    nf = norm(r1.Force);    

    tledCudaPrintf(2, "Consolidating forces on %d: %f, %f, %f - %f, %f, %f\n", r01.NodeIndex, r0.Force.x, r0.Force.y, r0.Force.z, r1.Force.x, r1.Force.y, r1.Force.z);
    if (nf != 0 && (pf = dot(r0.Force, r1.Force)/nf) < nf) {
      r01.Force += (1 - pf/nf)*r1.Force;
    }
    tledCudaPrintf(2, "Consolidated forces on %d: %f, %f, %f\n", r01.NodeIndex, r01.Force.x, r01.Force.y, r01.Force.z);

    return r01;
  }

  __device__ bool ContactResponseOrdering::operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const {
    return r0.NodeIndex < r1.NodeIndex;
  }  
  
  __device__ bool ContactResponseSameNodePredicate::operator()(const tledContactSolverGPU::ContactResponse &r0, const tledContactSolverGPU::ContactResponse &r1) const {
    return r0.NodeIndex == r1.NodeIndex;
  }  

  /* Only difference, client must take care of division by dt (at least implicitly) */
  __global__ void ComputeNodeVelocityKernel(float3 *dp_v, const float4 *uNexts, const float4 *uCurrs, const int *surface2VolumeMap, const int numNodes) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    if (tid < numNodes) {
      const int nInd = tid;
      const int gnInd = surface2VolumeMap[nInd];

      float4 v4 = uNexts[gnInd] - uCurrs[gnInd];
      
      dp_v[nInd] = make_float3(v4.x, v4.y, v4.z);
    }
  } 

  __global__ void ConvertResponsesToExternalForce(float4 *p_R, const tledContactSolverGPU::ContactResponse *responses, const int numResponses, const int *dpc_surface2VolumeMap) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    if (tid < numResponses) {
      const int nInd = responses[tid].NodeIndex;
      const float3 respF = responses[tid].Force;

      float3 avgF = respF;
      float fProj;
      int rangeStart, rangeEnd;

      tledCudaAssert(responses[tid].NodeIndex >= 0);
      tledCudaAssert(tid == 0 || responses[tid-1].NodeIndex <= responses[tid].NodeIndex);
      for (rangeEnd = tid + 1; rangeEnd < numResponses && responses[rangeEnd].NodeIndex == nInd; rangeEnd++) {
	avgF += responses[rangeEnd].Force;
      }      

      for (rangeStart = tid - 1; rangeStart >= 0 && responses[rangeStart].NodeIndex == nInd; rangeStart--) {
	avgF += responses[rangeStart].Force;
      }      
      rangeStart += 1;
      tledCudaAssert(rangeStart >= 0 && responses[rangeStart].NodeIndex == nInd);

      if (norm(avgF) > 0.f) {
	avgF = avgF/norm(avgF);
	fProj = dot(avgF, respF);	
	for (; rangeStart < rangeEnd; rangeStart++) {
	  const float cFProj = dot(responses[rangeStart].Force, avgF);

	  if (cFProj > fProj) break;
	  else if (cFProj == fProj && rangeStart < tid) break;
	}
	
	if (rangeStart == rangeEnd) {
	  tledCudaPrintf(0, "Applying response to %d/%d: %f %f %f\n", nInd, dpc_surface2VolumeMap[nInd], (fProj*avgF).x, (fProj*avgF).y, (fProj*avgF).z);
	  p_R[dpc_surface2VolumeMap[nInd]] -= fProj*avgF;	  
	}
      }
    }
  }

  template <class TProjection, class TXFormOperator>
  __global__ void ProjectionToResponseTransformKernel(tledContactSolverGPU::ContactResponse *p_responses, const TProjection *projections, const int numProjections, const TXFormOperator xFormOp, const int numResponsesPerProjection) {
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    const int projInd = tid/numResponsesPerProjection;

    if (projInd < numProjections) {      
      tledCudaPrintf(3, "processing projection %d - %d\n", projInd, tid%numResponsesPerProjection);
      p_responses[tid] = xFormOp(projections[projInd], tid%numResponsesPerProjection);
    }
  }
}
#endif
