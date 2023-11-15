// =========================================================================
// File:       tledElementMembraneSimpleLinear_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledElementMembraneSimpleLinear_kernels_CU
#define tledElementMembraneSimpleLinear_kernels_CU

#include "tledElementMembrane_kernels.h"
#include "tledElementMembraneSimpleLinear.h"
#include "tledSolverGPU_kernels.h"
#include "tledCUDA_operators.h"

namespace tledElementMembraneSimpleLinear_kernels {
  template <const int t_numFacetVtcs>
  __device__ void _ComputeForces(float4 *p_fDst, const float stress[], const float2 phiInvT[], const float3 elBasis[], const float t, const float A) {
    for (int vInd = 0; vInd < t_numFacetVtcs; vInd++) {
      float2 planeF;
      float3 spaceF;

      switch (vInd) {
      case 0:
	planeF.x = -(phiInvT[0].x + phiInvT[0].y)*stress[0] - (phiInvT[1].x + phiInvT[1].y)*stress[2];
	planeF.y = -(phiInvT[1].x + phiInvT[1].y)*stress[1] - (phiInvT[0].x + phiInvT[0].y)*stress[2];
	break;

      case 1:
	planeF.x = phiInvT[0].x*stress[0] + phiInvT[1].x*stress[2];
	planeF.y = phiInvT[1].x*stress[1] + phiInvT[0].x*stress[2];
	break;

      case 2:
	planeF.x = phiInvT[0].y*stress[0] + phiInvT[1].y*stress[2];
	planeF.y = phiInvT[1].y*stress[1] + phiInvT[0].y*stress[2];
	break;
      }

      spaceF = planeF.x*elBasis[0] + planeF.y*elBasis[1];
      p_fDst[vInd] = make_float4(spaceF.x, spaceF.y, spaceF.z, 0);
    }
  }
}

namespace tledElementMembrane_kernels {
  template <>
  __device__ void ComputeStrain<tledElementMembraneSimpleLinear<3> >(float *p_strainDst, tledElementMembraneSimpleLinear<3>::GPUElement *gp_elements) {
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const tledElementMembraneSimpleLinear<3>::GPUElement &element = gp_elements[elInd];
    const float4 uV0 = tledSolverGPU_kernels::GetCurrentDisplacement(element.ElementNodeIndices.x);

    float2 planeUs[2];
    
    for (int vInd = 0; vInd < 2; vInd++) {
      float4 localU;

      localU = tledSolverGPU_kernels::GetCurrentDisplacement(vInd == 0? element.ElementNodeIndices.y : element.ElementNodeIndices.z) - uV0;
      planeUs[vInd].x = dot(element.ElementBasis[0], localU);
      planeUs[vInd].y = dot(element.ElementBasis[1], localU);
    }

    p_strainDst[0] = element.PhiInvT[0].x*planeUs[0].x + element.PhiInvT[0].y*planeUs[1].x;
    p_strainDst[1] = element.PhiInvT[1].x*planeUs[0].y + element.PhiInvT[1].y*planeUs[1].y;
    p_strainDst[2] = element.PhiInvT[1].x*planeUs[0].x + element.PhiInvT[0].x*planeUs[0].y;
    p_strainDst[2] += element.PhiInvT[1].y*planeUs[1].x + element.PhiInvT[0].y*planeUs[1].y;
  }

  template <>
  __device__ void ComputeForces<tledElementMembraneSimpleLinear<3> >(float4 *p_fDst, tledElementMembraneSimpleLinear<3>::GPUElement *gp_elements, const float stress[]) {
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const tledElementMembraneSimpleLinear<3>::GPUElement &element = gp_elements[elInd];

    tledElementMembraneSimpleLinear_kernels::_ComputeForces<3>(p_fDst + 3*elInd, stress, element.PhiInvT, element.ElementBasis, element.Thickness, element.Area);
  }

  template <>
  __device__ void ComputeStrain<tledElementMembraneSimpleLinear<4> >(float *p_strainDst, tledElementMembraneSimpleLinear<4>::GPUElement *gp_elements) {
  }

  template <>
  __device__ void ComputeForces<tledElementMembraneSimpleLinear<4> >(float4 *p_fDst, tledElementMembraneSimpleLinear<4>::GPUElement *gp_elements, const float stress[]) {    
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;
    const tledElementMembraneSimpleLinear<4>::GPUElement &element = gp_elements[elInd];

    tledElementMembraneSimpleLinear_kernels::_ComputeForces<4>(p_fDst + 4*elInd, stress, element.PhiInvT, element.ElementBasis, element.Thickness, element.Area);
  }
}
#endif
