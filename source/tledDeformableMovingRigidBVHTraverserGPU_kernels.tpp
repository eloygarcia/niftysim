// =========================================================================
// File:       tledDeformableMovingRigidBVHTraverserGPU_kernels.tpp
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
#ifndef tledDeformableMovingRigidBVHTraverserGPU_kernels_CU
#define tledDeformableMovingRigidBVHTraverserGPU_kernels_CU

#include "tledCUDAHelpers.h"
#include "tledDynamicContactSurfaceGPU_kernels.h"
#include "tledContactSurfaceGPU_kernels.h"
#include "tledBVHTraverserGPU_kernels.h"

namespace tledDeformableMovingRigidBVHTraverserGPU_kernels {
  template <const int t_blockSize, class TMasterSurface, class TSlaveSurface, class TNodeFacetResult>
  __global__ void NodeFacetNarrowPhaseStage1(TNodeFacetResult *p_nodeList, int *p_numNodes, const TMasterSurface *pc_masterMesh, const TSlaveSurface *pc_slaveMesh, const int2 *narrowPhasePairs, const int numPairs) {
    typedef tledBVHTraverserGPU_kernels::NodeFacetProjectionIndexOrdering<TNodeFacetResult, TMasterSurface::NumberOfFacetVertices> __OutputOrdering;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const TMasterSurface &masterMesh = *pc_masterMesh;
    const TSlaveSurface &slaveMesh = *pc_slaveMesh;

    __shared__ TNodeFacetResult results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];

    tledCudaAssert(t_blockSize == blockDim.x);
    numResults[threadIdx.x] = 0;
    if (tInd < numPairs) {
      const int masterFInd = narrowPhasePairs[tInd].y;
      const int slaveNInd = narrowPhasePairs[tInd].x;
      const float3 oldV0 = masterMesh.OldNodeCoordinates[masterMesh.Facets[masterFInd].NodeIndices[0]];
      const float3 oldN = tledDynamicContactSurfaceGPU_kernels::ComputeOldNormal(masterMesh, masterFInd);
      
      float4 projOp[3];      
      
      tledContactSurfaceGPU_kernels::ComputeProjectionOperator(projOp, masterMesh, masterFInd);

      if (dot(slaveMesh.OldNodeCoordinates[slaveNInd] - oldV0, oldN) >= -tledUnstructuredContactManager_kernels::GetCloseDistance() 
	  && tledContactSurfaceGPU_kernels::ProjectOntoFacet(results[threadIdx.x].Xi, slaveMesh.NodeCoordinates[slaveNInd], projOp) 
	  && results[threadIdx.x].Xi.z < tledUnstructuredContactManager_kernels::GetCloseDistance()*1.25f) {
	results[threadIdx.x].ContactNodeIndices[0] = slaveNInd;
	for (int v = 0; v < TMasterSurface::NumberOfFacetVertices; v++) results[threadIdx.x].ContactNodeIndices[1+v] = masterMesh.Facets[masterFInd].NodeIndices[v];
	numResults[threadIdx.x] += 1;
      }
    }

    tledCUDAHelpers::MakeSortedUnique<t_blockSize, 1, TNodeFacetResult, __OutputOrdering>(results, numResults, __OutputOrdering(), true);
    tledCUDAHelpers::ParallelFlush<TNodeFacetResult>(p_nodeList, p_numNodes, results, numResults[0]);
  }

  template <const int t_blockSize, class TEdgeEdgeFinalProjection, class TEdgeEdgeIntermediateProjection, class TMasterSurface, class TSlaveSurface>
  __global__ void EdgeEdgeNarrowPhaseStage2(TEdgeEdgeFinalProjection *p_results, int *p_numResults, const TEdgeEdgeIntermediateProjection *intermediates, const int numIntermediates, const TMasterSurface *pc_masterMesh, const TSlaveSurface *pc_slaveMesh, const float safetyMargin) {
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionOrdering<TEdgeEdgeFinalProjection> __OutputOrdering;
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionSameSlavePredicate<TEdgeEdgeFinalProjection> __OutputEquality;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const float3 *masterNormals = pc_masterMesh->NodeNormals;
    const float3 *slaveNormals = pc_slaveMesh->NodeNormals;
    const float3 *masterX0 = pc_masterMesh->OldNodeCoordinates;
    const float3 *slaveX0 = pc_slaveMesh->OldNodeCoordinates;

    __shared__ TEdgeEdgeFinalProjection results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];    

    numResults[threadIdx.x] = 0;
    if (tInd < numIntermediates) {
      const TEdgeEdgeIntermediateProjection inItem = intermediates[tInd];

      float3 A0, B0, C0, D0;

      results[threadIdx.x].MasterEdge = inItem.MasterEdge;
      results[threadIdx.x].SlaveEdge = inItem.SlaveEdge;
      results[threadIdx.x].Xi.y = inItem.Xi.x;
      results[threadIdx.x].Xi.z = inItem.Xi.y;
      
      A0 = slaveX0[inItem.SlaveEdge.x];
      B0 = slaveX0[inItem.SlaveEdge.y];
      C0 = masterX0[inItem.MasterEdge.x];
      D0 = masterX0[inItem.MasterEdge.y];

      tledCudaPrintf(3, "%d/%d: Performing narrow-phase stage 2 on edges (%d, %d) - (%d, %d)\n", blockIdx.x, threadIdx.x, inItem.SlaveEdge.x, inItem.SlaveEdge.y, inItem.MasterEdge.x, inItem.MasterEdge.y);
      if (tledBVHTraverserGPU_kernels::ComputeEdgeEdgePenetrationDepth<true, true>(results[threadIdx.x].Xi.x, results[threadIdx.x].Normal, 
										   A0, B0, inItem.A1, inItem.B1, slaveNormals[inItem.SlaveEdge.x], slaveNormals[inItem.SlaveEdge.y], results[threadIdx.x].Xi.y,
										   C0, D0, inItem.C1, inItem.D1, masterNormals[inItem.MasterEdge.x], masterNormals[inItem.MasterEdge.y], results[threadIdx.x].Xi.z,
										   safetyMargin)) {
	numResults[threadIdx.x] = 1;
      }
    }

    tledCUDAHelpers::MergeSort<t_blockSize, 1, TEdgeEdgeFinalProjection, __OutputOrdering>(results, numResults, __OutputOrdering());
    tledCUDAHelpers::Unique<t_blockSize>(results, numResults[0], __OutputEquality());
    tledCUDAHelpers::ParallelFlush(p_results, p_numResults, results, numResults[0]);
  }
}

#endif
