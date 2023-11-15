// =========================================================================
// File:       tledElementMembraneNonLinear_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledElementMembraneNonLinear_kernels_CU
#define tledElementMembraneNonLinear_kernels_CU

#include "tledShellSolver_kernels.h"
#include "tledSolverGPU_kernels.h"

namespace tledElementMembraneNonLinear_kernels {
  __device__ void _ComputeCG(float *p_dst, const float3 J[]) {
    for (int r = 0; r < 2; r++) for (int c = 0; c < 2; c++) p_dst[r*2+c] = dot(J[r], J[c]);
  }

  __device__ void _CopyC0(float *p_dst, const float2 C0[]) {
    for (int i = 0; i < 2; i++) {
      p_dst[2*i] = C0[i].x;
      p_dst[2*i+1] = C0[i].y;
    }
  }

  __device__ void _PSKToF(float4 &r_f, const float stress[], const float2 &dN, const float3 J[], const float t, const float A) {
    float pskDN[2];
    float3 tmpF;

    for (int r = 0; r < 2; r++) pskDN[r] = dot(make_float2(stress[2*r], stress[2*r+1]), dN);
    tmpF = (pskDN[0]*J[0] + pskDN[1]*J[1])*t*A;
    r_f = make_float4(tmpF.x, tmpF.y, tmpF.z, 0.0f);
  }
}

namespace tledElementMembrane_kernels {  
  template <>
  __device__ void ComputeStrain<tledElementMembraneNonLinear<3> >(float *p_strainDst, tledElementMembraneNonLinear<3>::GPUElement *gp_elements) {
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const float4 uV0 = tledSolverGPU_kernels::GetCurrentDisplacement(gp_elements[elInd].ElementNodeIndices.x);

    tledElementMembraneNonLinear<3>::GPUElement &r_element = gp_elements[elInd];

    r_element.J[0] = r_element.InitialJacobian[0] + (tledSolverGPU_kernels::GetCurrentDisplacement(r_element.ElementNodeIndices.y) - uV0);
    r_element.J[1] = r_element.InitialJacobian[1] + (tledSolverGPU_kernels::GetCurrentDisplacement(r_element.ElementNodeIndices.z) - uV0);
    tledElementMembraneNonLinear_kernels::_ComputeCG(p_strainDst + 2*2, r_element.J);
    tledElementMembraneNonLinear_kernels::_CopyC0(p_strainDst, r_element.InitialCauchyGreenTensor);   
  }

  template <>
  __device__ void ComputeForces<tledElementMembraneNonLinear<3> >(float4 *p_fDst, tledElementMembraneNonLinear<3>::GPUElement *gp_elements, const float stress[]) {
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const tledElementMembraneNonLinear<3>::GPUElement &element = gp_elements[elInd];

    tledElementMembraneNonLinear_kernels::_PSKToF(p_fDst[3*elInd], stress, make_float2(-1, -1), element.J, element.Area, element.Thickness);
    tledElementMembraneNonLinear_kernels::_PSKToF(p_fDst[3*elInd+1], stress, make_float2(1, 0), element.J, element.Area, element.Thickness);
    tledElementMembraneNonLinear_kernels::_PSKToF(p_fDst[3*elInd+2], stress, make_float2(0, 1), element.J, element.Area, element.Thickness);
  }
}
#endif
