// =========================================================================
// File:       tledCentralDifferenceTimeStepperGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledCentralDifferenceTimeStepperGPU_kernels_CU
#define tledCentralDifferenceTimeStepperGPU_kernels_CU

#include "tledCentralDifferenceTimeStepperGPU.h"
#include "tledSolverGPU_kernels.h"
#include "tledTimeStepperGPU_kernels.h"
#include "tledCUDAHelpers.h"

#include "tledCentralDifferenceTimeStepperGPU_kernels.cu"

namespace tledCentralDifferenceTimeStepper_kernels {
  texture<float4, 1, cudaReadModeElementType> tx_CD;
}

namespace tledTimeStepper_kernels {
  template <>
  __device__ float3 ComputeNodalDisplacement<tledCentralDifferenceTimeStepperGPU::SolutionVariables>(tledCentralDifferenceTimeStepperGPU::SolutionVariables *p_sols, const int nodeInd, const float3 &nodeEffectiveF) {
    const float4 CD = tex1Dfetch(tledCentralDifferenceTimeStepper_kernels::tx_CD, nodeInd);
    const float4 ucurr = tledSolverGPU_kernels::GetCurrentDisplacement(nodeInd);

    float3 unext;
    
    unext = CD.x*nodeEffectiveF; 
    unext = unext + CD.y*make_float3(ucurr.x, ucurr.y, ucurr.z);
    unext = unext + CD.z*make_float3(p_sols->PreviousDisplacements[nodeInd].x, p_sols->PreviousDisplacements[nodeInd].y, p_sols->PreviousDisplacements[nodeInd].z);    

    return unext;
  }
}

#endif
