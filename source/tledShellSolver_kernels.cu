// =========================================================================
// File:       tledShellSolver_kernels.cu
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
#ifndef tledShellSolver_kernels_CU
#define tledShellSolver_kernels_CU

#include "tledShellSolver_kernels.h"
#include "tledElementMembrane_kernels.h"
#include "tledShellMaterial_kernels.h"
#include "tledCUDAHelpers.h"

namespace tledShellSolver_kernels {
  texture<int2, 1, cudaReadModeElementType> tx_NodeElementVertexLookupBaseIndex;
  texture<int, 1, cudaReadModeElementType> tx_NodeIndexLookupTable;
  texture<int, 1, cudaReadModeElementType> tx_NodeElementVertexLookupTable;
  texture<float4, 1, cudaReadModeElementType> tx_ElementVertexForces;  

  template <class TElement, class TMaterial>
  __global__ void ComputeNewForces(float4 *gp_F, typename TElement::GPUElement *gp_elements, const typename TMaterial::GPUMaterial *gpc_materialParameters, const int numElements) {
    const int elInd = blockIdx.x*blockDim.x + threadIdx.x;

    tledCudaAssert(numElements > 0);
    if (elInd < numElements) {
      float strain[TMaterial::NumberOfStrainComponents];
      float stress[TMaterial::NumberOfStressComponents];

      tledElementMembrane_kernels::ComputeStrain<TElement>(strain, gp_elements);
      tledShellMaterial_kernels::ComputeStress<TMaterial>(stress, strain, gpc_materialParameters);
      tledElementMembrane_kernels::ComputeForces<TElement>(gp_F, gp_elements, stress);
    }
  }

  __global__ void ComputeNodalForces(float4 *gp_f, const int numberOfNodes) {
    const int nInd = blockIdx.x*blockDim.x + threadIdx.x;

    if (nInd < numberOfNodes) {
      float4 f;
      
      {
	const int2 indBnds = tex1Dfetch(tx_NodeElementVertexLookupBaseIndex, nInd);
	
	f = make_float4(0, 0, 0, 0);
	for (int ltbInd = indBnds.x; ltbInd < indBnds.y; ltbInd++) {
	  const int vtxInd = tex1Dfetch(tx_NodeElementVertexLookupTable, ltbInd);

	  f += tex1Dfetch(tx_ElementVertexForces, vtxInd);
	}
      }

      tledCudaAssert(f.w == 0);

      {
	const int globalNInd = tex1Dfetch(tx_NodeIndexLookupTable, nInd);

	tledCudaAssert(globalNInd >= 0 && globalNInd < c_NumNodes);
	gp_f[globalNInd] -= f;
      }
    }
  }
}

#endif
