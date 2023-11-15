// =========================================================================
// File:       tledNewmarkTimeStepperGPU_kernels.cu
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
#ifndef tledNewmarkTimeStepperGPU_kernels_CU
#define tledNewmarkTimeStepperGPU_kernels_CU

#include "tledTimeStepperGPU_kernels.h"
#include "tledSolverGPU_kernels.h"
#include "tledCUDA_operators.h"
#include "tledCUDAHelpers.h"

namespace tledNewmarkTimeStepper_kernels {
  texture<float, 1, cudaReadModeElementType> tx_M;

  __device__ __constant__ float c_Dt;
  __device__ __constant__ float c_Alpha;
}

namespace tledTimeStepper_kernels {
  template <>
  __device__ float3 ComputeNodalDisplacement<tledNewmarkTimeStepperGPU::SolutionVariables>(tledNewmarkTimeStepperGPU::SolutionVariables *p_sols, const int nodeInd, const float3 &nodeEffectiveF) {
    using namespace tledNewmarkTimeStepper_kernels;

    const float m = tex1Dfetch(tx_M, nodeInd);
    const float4 ucurr = tledSolverGPU_kernels::GetCurrentDisplacement(nodeInd);    
    const float4 aprev = p_sols->PreviousAccelerations[nodeInd];
    const float4 vprev = p_sols->PreviousVelocities[nodeInd];

    float4 &r_acurr = p_sols->CurrentAccelerations[nodeInd];
    float4 &r_vcurr = p_sols->CurrentVelocities[nodeInd];
    float4 unext;
    
    //tledCudaPrintf(0, "a = %f dt = %f apx = %f apy = %f apz = %f\n", tledNewmarkTimeStepper_kernels::c_Alpha, tledNewmarkTimeStepper_kernels::c_Dt, aprev.x, aprev.y, aprev.z);
    //tledCudaPrintf(0, "m = %f vpx = %f vpy = %f vpz = %f\n", m, vprev.x, vprev.y, vprev.z);
    r_acurr = (make_float4(nodeEffectiveF.x/m, nodeEffectiveF.y/m, nodeEffectiveF.z/m, 0) - c_Alpha*(vprev + c_Dt/2*aprev))/(1 + c_Alpha*c_Dt/2);
    r_vcurr = vprev + c_Dt/2*(r_acurr + aprev);
    unext = ucurr + c_Dt*(r_vcurr + c_Dt/2*r_acurr);

    return make_float3(unext.x, unext.y, unext.z);
  }
}

#endif
