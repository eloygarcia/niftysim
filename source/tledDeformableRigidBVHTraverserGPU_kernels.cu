// =========================================================================
// File:       tledDeformableRigidBVHTraverserGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableRigidBVHTraverserGPU_kernels_CU
#define tledDeformableRigidBVHTraverserGPU_kernels_CU

#include "tledContactSurfaceGPU_kernels.h"
#include "tledDeformableContactSurfaceGPU_kernels.h"
#include "tledRigidContactSurfaceGPU_kernels.h"

namespace tledDeformableRigidBVHTraverserGPU_kernels {
  template <const int t_blockSize, class TRigidMesh, class TDeformableMesh, class TResult>
  __global__ void NodeFacetNarrowPhaseStage1RigidMaster(TResult *p_results, int *p_numResults, const TRigidMesh *pc_mMesh, const TDeformableMesh *pc_sMesh, const int2 *items, const int numItems) {
    typedef tledBVHTraverserGPU_kernels::NodeFacetProjectionIndexOrdering<TResult, TRigidMesh::NumberOfFacetVertices> __OutputOrdering;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const TRigidMesh &mMesh = *pc_mMesh;
    const TDeformableMesh &sMesh = *pc_sMesh;

    __shared__ TResult results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];

    tledCudaAssert(t_blockSize == blockDim.x);
    numResults[threadIdx.x] = 0;
    if (tInd < numItems) {
      const int masterFInd = items[tInd].y;
      const int slaveNInd = items[tInd].x;
      const float3 V0 = mMesh.NodeCoordinates[mMesh.Facets[masterFInd].NodeIndices[0]];
      const float3 n = mMesh.FacetNormals[masterFInd];

      tledCudaAssert(masterFInd < mMesh.NumberOfFacets && masterFInd >= 0);
      tledCudaAssert(slaveNInd >= 0 && slaveNInd < sMesh.NumberOfNodes);
      if (dot(sMesh.OldNodeCoordinates[slaveNInd] - V0, n) >= -tledUnstructuredContactManager_kernels::GetCloseDistance() 
	  && tledContactSurfaceGPU_kernels::ProjectOntoFacet(results[threadIdx.x].Xi, sMesh.NodeCoordinates[slaveNInd], tledRigidContactSurfaceGPU_kernels::GetFacetProjectionOperator<TRigidMesh>(mMesh, masterFInd)) 
	  && results[threadIdx.x].Xi.z < tledUnstructuredContactManager_kernels::GetCloseDistance()*1.25f) {
	results[threadIdx.x].ContactNodeIndices[0] = slaveNInd;
	for (int v = 0; v < TRigidMesh::NumberOfFacetVertices; v++) results[threadIdx.x].ContactNodeIndices[1+v] = mMesh.Facets[masterFInd].NodeIndices[v];
	numResults[threadIdx.x] += 1;
      }
    }

    tledCUDAHelpers::MakeSortedUnique<t_blockSize, 1, TResult, __OutputOrdering>(results, numResults, __OutputOrdering(), true);
    tledCUDAHelpers::ParallelFlush<TResult>(p_results, p_numResults, results, numResults[0]);
  }

  template <const int t_blockSize, class TDeformableMesh, class TRigidMesh, class TResult>
  __global__ void NodeFacetNarrowPhaseStage1DeformableMaster(TResult *p_results, int *p_numResults, const TDeformableMesh *pc_mMesh, const TRigidMesh *pc_sMesh, const int2 *items, const int numItems) {
    typedef tledBVHTraverserGPU_kernels::NodeFacetProjectionIndexOrdering<TResult, TDeformableMesh::NumberOfFacetVertices> __OutputOrdering;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const TRigidMesh &sMesh = *pc_sMesh;
    const TDeformableMesh &mMesh = *pc_mMesh;

    __shared__ TResult results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];

    tledCudaAssert(t_blockSize == blockDim.x);
    numResults[threadIdx.x] = 0;
    if (tInd < numItems) {
      const int masterFInd = items[tInd].y;
      const int slaveNInd = items[tInd].x;
      const float3 oldV0 = mMesh.OldNodeCoordinates[mMesh.Facets[masterFInd].NodeIndices[0]];
      const float3 oldN = tledDeformableContactSurfaceGPU_kernels::ComputeOldNormal(mMesh, masterFInd);
      
      float4 projOp[3];      
      
      tledContactSurfaceGPU_kernels::ComputeProjectionOperator<TDeformableMesh>(projOp, mMesh, masterFInd);
      if (dot(sMesh.NodeCoordinates[slaveNInd] - oldV0, oldN) >= -tledUnstructuredContactManager_kernels::GetCloseDistance() 
	  && tledContactSurfaceGPU_kernels::ProjectOntoFacet(results[threadIdx.x].Xi, sMesh.NodeCoordinates[slaveNInd], projOp) 
	  && results[threadIdx.x].Xi.z < tledUnstructuredContactManager_kernels::GetCloseDistance()*1.25f) {
	results[threadIdx.x].ContactNodeIndices[0] = slaveNInd;
	for (int v = 0; v < TDeformableMesh::NumberOfFacetVertices; v++) results[threadIdx.x].ContactNodeIndices[1+v] = mMesh.Facets[masterFInd].NodeIndices[v];
	numResults[threadIdx.x] += 1;
      }
    }

    tledCUDAHelpers::MakeSortedUnique<t_blockSize, 1, TResult, __OutputOrdering>(results, numResults, __OutputOrdering(), true);
    tledCUDAHelpers::ParallelFlush<TResult>(p_results, p_numResults, results, numResults[0]);
  }

  template <const int t_blockSize, class TEdgeEdgeFinalProjection, class TEdgeEdgeIntermediateProjection, class TMasterMesh, class TSlaveMesh>
  __global__ void EdgeEdgeNarrowPhaseStage2RigidMaster(TEdgeEdgeFinalProjection *p_results, int *p_numResults, const TEdgeEdgeIntermediateProjection *intermediates, const int numIntermediates, const TMasterMesh *pc_mMesh, const TSlaveMesh *pc_sMesh, const float safetyMargin) {
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionOrdering<TEdgeEdgeFinalProjection> __OutputOrdering;
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionSameSlavePredicate<TEdgeEdgeFinalProjection> __OutputEquality;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const float3 *sNormals = pc_sMesh->NodeNormals;
    const float3 *mNormals = pc_mMesh->NodeNormals;
    const float3 *oldNodeCoordinates = pc_sMesh->OldNodeCoordinates;

    __shared__ TEdgeEdgeFinalProjection results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];    

    numResults[threadIdx.x] = 0;
    if (tInd < numIntermediates) {
      const TEdgeEdgeIntermediateProjection inItem = intermediates[tInd];

      float3 A0, B0;

      results[threadIdx.x].MasterEdge = inItem.MasterEdge;
      results[threadIdx.x].SlaveEdge = inItem.SlaveEdge;
      results[threadIdx.x].Xi.y = inItem.Xi.x;
      results[threadIdx.x].Xi.z = inItem.Xi.y;
      
      A0 = oldNodeCoordinates[inItem.SlaveEdge.x];
      B0 = oldNodeCoordinates[inItem.SlaveEdge.y];

      if (tledBVHTraverserGPU_kernels::ComputeEdgeEdgePenetrationDepth<true, false>(results[threadIdx.x].Xi.x, results[threadIdx.x].Normal, 
										    A0, B0, inItem.A1, inItem.B1, sNormals[inItem.SlaveEdge.x], sNormals[inItem.SlaveEdge.y], results[threadIdx.x].Xi.y,
										    inItem.C1, inItem.D1, inItem.C1, inItem.D1, mNormals[inItem.MasterEdge.x], mNormals[inItem.MasterEdge.y], results[threadIdx.x].Xi.z,
										    safetyMargin)) {
	numResults[threadIdx.x] = 1;
      }
    }

    tledCUDAHelpers::MergeSort<t_blockSize, 1, TEdgeEdgeFinalProjection, __OutputOrdering>(results, numResults, __OutputOrdering());
    tledCUDAHelpers::Unique<t_blockSize>(results, numResults[0], __OutputEquality());
    tledCUDAHelpers::ParallelFlush(p_results, p_numResults, results, numResults[0]);
  }

  template <const int t_blockSize, class TEdgeEdgeFinalProjection, class TEdgeEdgeIntermediateProjection, class TMasterMesh, class TSlaveMesh>
  __global__ void EdgeEdgeNarrowPhaseStage2DeformableMaster(TEdgeEdgeFinalProjection *p_results, int *p_numResults, const TEdgeEdgeIntermediateProjection *intermediates, const int numIntermediates, const TMasterMesh *pc_mMesh, const TSlaveMesh *pc_sMesh, const float safetyMargin) {
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionOrdering<TEdgeEdgeFinalProjection> __OutputOrdering;
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionSameSlavePredicate<TEdgeEdgeFinalProjection> __OutputEquality;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const float3 *sNormals = pc_sMesh->NodeNormals;
    const float3 *mNormals = pc_mMesh->NodeNormals;
    const float3 *oldNodeCoordinates = pc_mMesh->OldNodeCoordinates;

    __shared__ TEdgeEdgeFinalProjection results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];    

    numResults[threadIdx.x] = 0;
    if (tInd < numIntermediates) {
      const TEdgeEdgeIntermediateProjection inItem = intermediates[tInd];

      float3 C0, D0;

      results[threadIdx.x].MasterEdge = inItem.MasterEdge;
      results[threadIdx.x].SlaveEdge = inItem.SlaveEdge;
      results[threadIdx.x].Xi.y = inItem.Xi.x;
      results[threadIdx.x].Xi.z = inItem.Xi.y;
      
      C0 = oldNodeCoordinates[inItem.MasterEdge.x];
      D0 = oldNodeCoordinates[inItem.MasterEdge.y];

      if (tledBVHTraverserGPU_kernels::ComputeEdgeEdgePenetrationDepth<false, true>(results[threadIdx.x].Xi.x, results[threadIdx.x].Normal, 
										    inItem.A1, inItem.B1, inItem.A1, inItem.B1, sNormals[inItem.SlaveEdge.x], sNormals[inItem.SlaveEdge.y], results[threadIdx.x].Xi.y,
										    C0, D0, inItem.C1, inItem.D1, mNormals[inItem.MasterEdge.x], mNormals[inItem.MasterEdge.y], results[threadIdx.x].Xi.z,
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
