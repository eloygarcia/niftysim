// =========================================================================
// File:       tledBVHTraverserGPU_kernels.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserGPU_kernels_TPP
#define tledBVHTraverserGPU_kernels_TPP

#include "tledBVHTraverserGPU_kernels.h"
#include "tledContactSurfaceGPU_kernels.h"
#include "tledUnstructuredContactManager_kernels.h"
#include "tledBV_kernels.h"
#include "tledCUDAHelpers.h"

#include "tledCUDA_operators.cu"

namespace tledBVHTraverserGPU_kernels {
  template <class TProjection, class TSurface>
  __device__ bool ProjectOntoSurfaceC0(TProjection &r_proj, const float3 &x, const TSurface &mesh, const float3 normals[]) {
    const float tauRelChange = 1e-5f;
    const float facetMargin = 1e-4f;
    const int maxProjIts = 10;

    float3 a, b, r, xi;
    
    tledCudaAssert(TSurface::NumberOfFacetVertices == 3);
    xi.x = r_proj.ShapeValues[1], xi.y = r_proj.ShapeValues[2], xi.z = r_proj.GapValue;         

    {
      float3 fVtcs[3];

      for (int v = 0; v < 3; v++) fVtcs[v] = mesh.NodeCoordinates[r_proj.ContactNodeIndices[1+v]];
      a = fVtcs[1] - fVtcs[0];
      b = fVtcs[2] - fVtcs[0];
      r = x - fVtcs[0];
    }

    for (int it = 0; it < maxProjIts; it++) {
      float3 n, tmpXi, shapeVals;

      shapeVals.x = (1 - xi.x - xi.y), shapeVals.y = xi.x, shapeVals.z = xi.y;      
      n = normals[0]*shapeVals.x + normals[1]*shapeVals.y + normals[2]*shapeVals.z;
      n /= norm(n);
      tledCUDAMaths::SolveEquation3x3(tmpXi.x, tmpXi.y, tmpXi.z, a, b, n, r);
      if (tmpXi.x < -facetMargin || tmpXi.y < -facetMargin || tmpXi.x > 1 + facetMargin || tmpXi.y > 1 + facetMargin || tmpXi.x + tmpXi.y > 1 + facetMargin) return false;
      else {
	float sum;

	if (tmpXi.x < 0) tmpXi.x = 0;
	if (tmpXi.y < 0) tmpXi.y = 0;
	if ((sum = tmpXi.x + tmpXi.y) > 1) {
	  tmpXi.x /= sum;
	  tmpXi.y /= sum;
	}

	if (norm(tmpXi - xi) <= tauRelChange*norm(xi)) {
	  r_proj.Normal = n;
	  r_proj.ShapeValues[0] = shapeVals.x, r_proj.ShapeValues[1] = shapeVals.y, r_proj.ShapeValues[2] = shapeVals.z;
	  r_proj.GapValue = xi.z;

	  return true;
	} else xi = tmpXi;
      }
    }

    r_proj.ShapeValues[0] = 1 - xi.x - xi.y;
    r_proj.ShapeValues[1] = xi.x;
    r_proj.ShapeValues[2] = xi.y;

    r_proj.GapValue = xi.z;

    r_proj.Normal = normals[0]*r_proj.ShapeValues[0];
    for (int i = 1; i < 3; i++) r_proj.Normal += normals[i]*r_proj.ShapeValues[i];
    r_proj.Normal /= norm(r_proj.Normal);

    return true;
  }

  template <const int t_blockSize, class TMasterBV, class TSlaveBV>
  __global__ void FilterNonIntersectingBroadPhaseBVPairsKernel(int2 *p_outList, int *p_outCtr, const int2 *inList, const int *pc_inCtr, const TMasterBV *masterBVs, const TSlaveBV *slaveBVs) {
    const int numIn = *pc_inCtr;
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ bool doKeep[t_blockSize];
    __shared__ int baseInd, numItems;
    int offOut;

    doKeep[threadIdx.x] = false;
    if (tid < numIn) {
      doKeep[threadIdx.x] = tledBV_kernels::DoIntersect(masterBVs[inList[tid].x], slaveBVs[inList[tid].y]);
    }
    
    if (__syncthreads_or(doKeep[threadIdx.x])) {
      offOut = tledCUDAHelpers::ComputeOffsetFromFlag<t_blockSize>(doKeep);
      if (threadIdx.x == t_blockSize - 1) {
    	tledCudaAssert(offOut < t_blockSize);
    	numItems = doKeep[threadIdx.x] + offOut;
    	baseInd = atomicAdd(p_outCtr, numItems);
      }
      __syncthreads();
      tledCudaBlockPrintf(1, "block %d BV-pair filtering: %d output BV pairs.\n", blockIdx.x, numItems);
    
      if (doKeep[threadIdx.x]) {
    	p_outList[baseInd+offOut] = inList[tid];
    	tledCudaAssert(tledBV_kernels::DoIntersect(masterBVs[p_outList[baseInd+offOut].x], slaveBVs[p_outList[baseInd+offOut].y]));
      }
    } else {
      tledCudaAssert(!doKeep[threadIdx.x]);
      tledCudaBlockPrintf(1, "block %d BV-pair filtering: no intersecting BV-pairs found.\n", blockIdx.x);
    }
  }

  template <const int t_blockSize, class TMasterBV>
  __global__ void DescendInMasterBVHKernel(int2 *p_outList, int *p_outCtr, const int2 *inList, const int *pc_inCtr, const TMasterBV *masterBVs) {
    const int numIn = *pc_inCtr;
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ int2 pairs[t_blockSize*TMasterBV::NumberOfChildBVs];
    __shared__ unsigned short numItems[t_blockSize];
    
    numItems[threadIdx.x] = 0;
    if (tid < numIn) {
      const int2 &inPair = inList[tid];
      const TMasterBV &masterParent = masterBVs[inPair.x];

      tledCudaAssert(inPair.x >= 0 && inPair.y >= 0);
      if (tledBV_kernels::IsLeaf<TMasterBV>(masterParent)) {
    	tledCudaPrintf(3, "%d/%d: leaf reached %d\n", threadIdx.x, blockIdx.x, inPair.x);
    	pairs[threadIdx.x*TMasterBV::NumberOfChildBVs] = inPair;
    	numItems[threadIdx.x] = 1;
      } else {
    	tledCudaPrintf(3, "%d/%d: descending in %d\n", threadIdx.x, blockIdx.x, inPair.x);
    	for (int c = 0; c < TMasterBV::NumberOfChildBVs; c++) {
    	  if (TMasterBV::NumberOfChildBVs == 2 || masterParent.ChildIndices[c] >= 0) {
    	    pairs[threadIdx.x*TMasterBV::NumberOfChildBVs+numItems[threadIdx.x]].x = masterParent.ChildIndices[c];
    	    pairs[threadIdx.x*TMasterBV::NumberOfChildBVs+numItems[threadIdx.x]].y = inList[tid].y;

    	    tledCudaPrintf(3, "%d/%d: new master BV %d\n", threadIdx.x, blockIdx.x, pairs[threadIdx.x*TMasterBV::NumberOfChildBVs+numItems[threadIdx.x]].x);
    	    tledCudaAssert(pairs[threadIdx.x*TMasterBV::NumberOfChildBVs+numItems[threadIdx.x]].x > inPair.x);

    	    numItems[threadIdx.x] += 1;
    	  } 
    	}
      }
      tledCudaAssert(numItems[threadIdx.x] <= TMasterBV::NumberOfChildBVs);
    }
    __syncthreads();
    
    tledCUDAHelpers::CollectResults<t_blockSize>(pairs, numItems, TMasterBV::NumberOfChildBVs);
    tledCudaBlockPrintf(1, "block %d: %d output next-level BV-pairs.\n", blockIdx.x, numItems[0]);
    tledCUDAHelpers::ParallelFlush<int2>(p_outList, p_outCtr, pairs, numItems[0]);
  }

  template <class TMasterBV, class TSlaveBV>
  __global__ void ConvertBVToPrimitivePairsKernel(int2 *p_list, const int *pc_inCtr, const TMasterBV *masterBVs, const TSlaveBV *slaveBVs) {
    const int numIn = *pc_inCtr;
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    
    if (tid < numIn) {
      tledCudaAssert(tledBV_kernels::IsLeaf<TSlaveBV>(slaveBVs[p_list[tid].y]));
      tledCudaAssert(tledBV_kernels::IsLeaf<TMasterBV>(masterBVs[p_list[tid].x]));

      p_list[tid].x = masterBVs[p_list[tid].x].PrimitiveIndex;
      p_list[tid].y = slaveBVs[p_list[tid].y].PrimitiveIndex;
    }
  }

  template <const int t_numFacetVertices, const int t_blockSize, const bool t_doMaster>
  __global__ void AddNarrowPhaseNodeFacetTests(int2 *p_nodeFacetPairs, int *p_numNodeFacetPairs, const int *meshFacetVertexIndices, const int2 *narrowPhasePairs, const int *dpc_numPairs) {
    const int numPairs = *dpc_numPairs;
    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ int2 nodeFacetPairs[t_numFacetVertices*t_blockSize];      
    __shared__ unsigned short numNodeFacetPairs[t_blockSize];      

    tledCudaAssert(blockDim.x == t_blockSize);

    numNodeFacetPairs[threadIdx.x] = 0;
    if (tInd < numPairs) {
      const int masterFInd = t_doMaster? narrowPhasePairs[tInd].x : narrowPhasePairs[tInd].y;
      const int slaveFInd = t_doMaster? narrowPhasePairs[tInd].y : narrowPhasePairs[tInd].x;
      
      numNodeFacetPairs[threadIdx.x] = t_numFacetVertices;
      for (int v = 0; v < t_numFacetVertices; v++) {
	nodeFacetPairs[t_numFacetVertices*threadIdx.x+v].x = meshFacetVertexIndices[slaveFInd*t_numFacetVertices+v];
	nodeFacetPairs[t_numFacetVertices*threadIdx.x+v].y = masterFInd;
      }	
    }
      
    tledCUDAHelpers::MakeSortedUnique<t_blockSize, t_numFacetVertices, int2, BroadPhaseResultOrdering>(nodeFacetPairs, numNodeFacetPairs, BroadPhaseResultOrdering(), false);
    tledCUDAHelpers::ParallelFlush(p_nodeFacetPairs, p_numNodeFacetPairs, nodeFacetPairs, numNodeFacetPairs[0]);
  }

  template <const int t_numFacetEdges, const int t_blockSize, const bool t_doMaster>
  __global__ void AddNarrowPhaseEdgeEdgeTests(int2 *p_edgeEdgePairs, int *p_numEdgeEdgePairs, const int *slaveFacetEdgeIndices, const int *masterFacetEdgeIndices, const int2 *narrowPhasePairs, const int *dpc_numPairs) {
    const int numPairs = *dpc_numPairs;
    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const int itemInd = tInd/t_numFacetEdges;
    const int edgeInd = tInd%t_numFacetEdges;

    __shared__ int2 edgeEdgePairs[t_numFacetEdges*t_blockSize];      
    __shared__ unsigned short numEdgeEdgePairs[t_blockSize];      
      
    numEdgeEdgePairs[threadIdx.x] = 0;
    if (itemInd < numPairs) {
      const int masterFInd = t_doMaster? narrowPhasePairs[itemInd].x : narrowPhasePairs[itemInd].y;
      const int slaveFInd = t_doMaster? narrowPhasePairs[itemInd].y : narrowPhasePairs[itemInd].x;
      const int slaveEdgeInd = slaveFacetEdgeIndices[slaveFInd*t_numFacetEdges+edgeInd];

      numEdgeEdgePairs[threadIdx.x] = t_numFacetEdges;
      for (int me = 0; me < t_numFacetEdges; me++) {
	edgeEdgePairs[t_numFacetEdges*threadIdx.x+me].x = slaveEdgeInd;
	edgeEdgePairs[t_numFacetEdges*threadIdx.x+me].y = masterFacetEdgeIndices[masterFInd*t_numFacetEdges+me];
      }	
    }
      
    tledCUDAHelpers::MakeSortedUnique<t_blockSize, t_numFacetEdges, int2, BroadPhaseResultOrdering>(edgeEdgePairs, numEdgeEdgePairs, BroadPhaseResultOrdering(), false);
    tledCUDAHelpers::ParallelFlush(p_edgeEdgePairs, p_numEdgeEdgePairs, edgeEdgePairs, numEdgeEdgePairs[0]);
  }

  template <const int t_slaveNumFacetVtcs, const int t_masterNumFacetVtcs>
  __global__ void PrepareNarrowPhaseKernel(int2 *p_nodeFacetPairs, int *p_numNodeFacet, int2 *p_edgePairs, int *p_numEdgeEdge, const int2 *narrowPhasePrimitivePairs, const int *pc_numPairs, const int *slaveFacetVertices, const int *masterFacetEdges, const int *slaveFacetEdges, const bool doMaster) {
    const int numPairs = *pc_numPairs;
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;

    if (tid < numPairs) {
      const int masterF = doMaster? narrowPhasePrimitivePairs[tid].x : narrowPhasePrimitivePairs[tid].y;
      const int slaveF = doMaster? narrowPhasePrimitivePairs[tid].y : narrowPhasePrimitivePairs[tid].x;
      const int nfBase = atomicAdd(p_numNodeFacet, 3);
      const int eeBase = atomicAdd(p_numEdgeEdge, 9);

      for (int v = 0; v < 3; v++) {
	p_nodeFacetPairs[nfBase+v].x = slaveFacetVertices[slaveF*t_slaveNumFacetVtcs+v];
	p_nodeFacetPairs[nfBase+v].y = masterF;
      }

      for (int es = 0; es < 3; es++) for (int em = 0; em < 3; em++) {
	  p_edgePairs[eeBase+3*es+em].x = slaveFacetEdges[slaveF*t_slaveNumFacetVtcs+es];
	  p_edgePairs[eeBase+3*es+em].y = masterFacetEdges[masterF*t_masterNumFacetVtcs+em];
	}
    }
  }

  inline __device__ bool ComputeEdgeEdgeClosestPointParameters(float &r_r, float &r_q, const float3 &A, const float3 &B, const float3 &C, const float3 &D) {
    const float safetyMargin = 1e-2f;

    float3 u, v, w;
    float a, b, c, d, e, denom;

    u = B - A;
    v = D - C;
    w = A - C;

    a = dot(u, u);
    b = dot(u, v);
    c = dot(v, v);
    d = dot(u, w);
    e = dot(v, w);
    tledCudaAssert(a != a || a > 0.f);
    tledCudaAssert(c != c || c > 0.f);

    denom = a*c - b*b;
    tledCudaAssert(denom != denom || denom >= -1e-4f*(a + c));
    tledCudaAssert((a + c) != (a + c) || (a > 0.f && c > 0.f));
    if (denom <= 1e-7f*(a + c)) {
      r_r = 0.5f;
      r_q = (d + a/2)/b;
      if (r_q < -safetyMargin) {
	r_q = 0.0f;
	r_r = -d/a;
      } else if (r_q > 1 + safetyMargin) {
	r_q = 1.0f;
	r_r = (c - e)/b;
      }
    } else {
      r_r = (b*e - c*d)/denom;
      r_q = (a*e - b*d)/denom;
    }

    if (r_r >= -safetyMargin && r_r <= 1.0f + safetyMargin && r_q >= -safetyMargin && r_q <= 1.0f + safetyMargin) {
      const float thr = (a + c)/1024;

      float3 ab, cd, d;

      tledCudaAssert((a != a || a >= 0) && (c != c || c >= 0));
      ab = interpolate(A, B, r_r);
      cd = interpolate(C, D, r_q);
      d = ab - cd;

      tledCudaAssert(fabsf(dot(u, d)) < thr*128 && fabsf(dot(v, d)) < thr*128);
      r_r = fmaxf(0.0f, fminf(r_r, 1.0f));
      r_q = fmaxf(0.0f, fminf(r_q, 1.0f));
      if (!(fabsf(dot(u, d)) < thr && fabsf(dot(v, d)) < thr)) return false;
      else return true;
    } else return false;
  }

  template <const bool t_isMovingSlave, const bool t_isMovingMaster>
  __device__ bool ComputeEdgeEdgePenetrationDepth(float &r_g, float3 &r_n, 
						  const float3 &slave0T0, const float3 &slave1T0, const float3& slave0T1, const float3 &slave1T1, const float3 &slaveN0, const float3 &slaveN1, const float r, 
						  const float3 &master0T0, const float3 &master1T0, const float3 &master0T1, const float3 &master1T1, const float3 &masterN0, const float3 &masterN1, const float q,
						  const float safetyMargin) {
    float3 slaveX0, masterX0, d0;

    slaveX0 = interpolate(slave0T0, slave1T0, r);
    masterX0 = interpolate(master0T0, master1T0, q);
    r_n = interpolate(masterN0, masterN1, q);
    r_n /= norm(r_n);
    d0 = slaveX0 - masterX0;
    if (dot(r_n, d0) >= -safetyMargin) {
      float3 *p_slaveX1, *p_masterX1, d1, slaveX1, masterX1;
      
      if (t_isMovingSlave) {
	p_slaveX1 = &slaveX1;
	slaveX1 = interpolate(slave0T1, slave1T1, r);
      } else {
	p_slaveX1 = &slaveX0;
      }
      
      if (t_isMovingMaster) {
	masterX1 = interpolate(master0T1, master1T1, q);
	p_masterX1 = &masterX1;    
      } else {
	p_masterX1 = &masterX0;
      }
      
      d1 = *p_slaveX1 - *p_masterX1;
      r_g = dot(r_n, d1);
      if (r_g < safetyMargin) {
	float3 slaveN;
	
	slaveN = interpolate(slaveN0, slaveN1, r);
	
	return dot(r_n, slaveN) < 0 && fabsf(r_g) > 0.5f*norm(d1);
      }
    } 

    return false;    
  }

  template <const int t_blockSize, class TMasterSurface, class TSlaveSurface, class TEdgeEdgeResult>
  __global__ void EdgeEdgeNarrowPhaseStage1(TEdgeEdgeResult *p_edgeList, int *p_numEdges, const TMasterSurface *pc_mMesh, const TSlaveSurface *pc_sMesh, const int2 *narrowPhasePairs, const int numPairs) {
    typedef tledBVHTraverserGPU_kernels::EdgeEdgeProjectionIndexOrdering<TEdgeEdgeResult> __OutputOrdering;

    const int tInd = blockIdx.x*blockDim.x + threadIdx.x;
    const TMasterSurface &mMesh = *pc_mMesh;
    const TSlaveSurface &sMesh = *pc_sMesh;

    __shared__ TEdgeEdgeResult results[t_blockSize];
    __shared__ unsigned short numResults[t_blockSize];

    numResults[threadIdx.x] = 0;
    tledCudaAssert(t_blockSize == blockDim.x);
    if (tInd < numPairs) {
      const int slaveInd = narrowPhasePairs[tInd].x;
      const int masterInd = narrowPhasePairs[tInd].y;

      int2 mEdge, sEdge;
      
      tledCudaAssert(slaveInd < sMesh.NumberOfEdges && slaveInd >= 0);
      tledCudaAssert(masterInd < mMesh.NumberOfEdges && masterInd >= 0);
      sEdge = results[threadIdx.x].SlaveEdge = sMesh.Edges[slaveInd];
      mEdge = results[threadIdx.x].MasterEdge = mMesh.Edges[masterInd];	    

      results[threadIdx.x].A1 = sMesh.NodeCoordinates[sEdge.x];
      results[threadIdx.x].B1 = sMesh.NodeCoordinates[sEdge.y];
      results[threadIdx.x].C1 = mMesh.NodeCoordinates[mEdge.x];
      results[threadIdx.x].D1 = mMesh.NodeCoordinates[mEdge.y];

      if (tledBVHTraverserGPU_kernels::ComputeEdgeEdgeClosestPointParameters(results[threadIdx.x].Xi.x, results[threadIdx.x].Xi.y, results[threadIdx.x].A1, results[threadIdx.x].B1, results[threadIdx.x].C1, results[threadIdx.x].D1)) {
	numResults[threadIdx.x] = 1;
      }
    }

    tledCUDAHelpers::MakeSortedUnique<t_blockSize, 1, TEdgeEdgeResult,  __OutputOrdering>(results, numResults, __OutputOrdering(), true);    
    tledCUDAHelpers::ParallelFlush(p_edgeList, p_numEdges, results, numResults[0]);
  }

  template <const int t_blockSize, class TFinalProjection, class TInitialProjection, class TMasterMesh, class TSlaveMesh>
  __global__ void NodeFacetNarrowPhaseStage2(TFinalProjection *p_projections, int *p_numProjections, const TInitialProjection *initialProjections, const int numInitial, const TMasterMesh *pc_mMesh, const TSlaveMesh *pc_sMesh) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    const TMasterMesh &mMesh = *pc_mMesh;
    const TSlaveMesh &sMesh = *pc_sMesh;

    __shared__ TFinalProjection results[t_blockSize];
    __shared__ unsigned short numItems[t_blockSize];
    
    tledCudaAssert(t_blockSize == blockDim.x);
    numItems[threadIdx.x] = 0;
    if (tid < numInitial) {
      float3 normals[TMasterMesh::NumberOfFacetVertices];

      results[threadIdx.x].GapValue = initialProjections[tid].Xi.z;
      tledContactSurfaceGPU_kernels::ComputeShapeValues<TMasterMesh>(results[threadIdx.x].ShapeValues, initialProjections[tid].Xi);

      tledCudaAssert(initialProjections[tid].ContactNodeIndices[0] >= 0 && initialProjections[tid].ContactNodeIndices[0] < sMesh.NumberOfNodes);
      results[threadIdx.x].ContactNodeIndices[0] = initialProjections[tid].ContactNodeIndices[0];
      for (int v = 0; v < TMasterMesh::NumberOfFacetVertices; v++) {
	results[threadIdx.x].ContactNodeIndices[v+1] = initialProjections[tid].ContactNodeIndices[v+1];
	normals[v] = tledContactSurfaceGPU_kernels::GetNodeNormal(mMesh, results[threadIdx.x].ContactNodeIndices[v+1]);
	tledCudaAssert(initialProjections[tid].ContactNodeIndices[1+v] >= 0 && initialProjections[tid].ContactNodeIndices[1+v] < mMesh.NumberOfNodes);
      }

      if (tledBVHTraverserGPU_kernels::ProjectOntoSurfaceC0<TFinalProjection, TMasterMesh>(results[threadIdx.x], sMesh.NodeCoordinates[results[threadIdx.x].ContactNodeIndices[0]], mMesh, normals)) {
	if (results[threadIdx.x].GapValue < tledUnstructuredContactManager_kernels::GetCloseDistance()) {
	  numItems[threadIdx.x] += 1;
	}
      }      
    }

    tledCUDAHelpers::CollectResults<t_blockSize, TFinalProjection>(results, numItems, 1);
    tledCUDAHelpers::ParallelFlush(p_projections, p_numProjections, results, numItems[0]);
  }
}
#endif
