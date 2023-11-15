// =========================================================================
// File:       tledTrackingBVHTraverserGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    December 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTrackingBVHTraverserGPU_kernels_CU
#define tledTrackingBVHTraverserGPU_kernels_CU

#include "tledBVHTraverserGPU_kernels.h"
#include "tledCUDAHelpers.h"

#include "tledCUDA_operators.cu"

namespace tledTrackingBVHTraverserGPU_kernels {
  __device__ void _FlushContacts(int2 *p_dst, int *p_counter, const int2 contacts[], const bool haveContact[]) {
    __syncthreads();
    if (haveContact[threadIdx.x]) {
      __shared__ int baseInd;

      tledCudaAssert(threadIdx.x == 0 || haveContact[threadIdx.x-1]);
      if (threadIdx.x == blockDim.x - 1 || !haveContact[threadIdx.x+1]) {
	baseInd = atomicAdd(p_counter, threadIdx.x + 1);
	tledCudaPrintf(3, "block %d: Flushing %d contact candidates.\n", blockIdx.x, threadIdx.x + 1);
      }
      __syncthreads();

      p_dst[baseInd+threadIdx.x] = contacts[threadIdx.x];
    }
  }
  
  template <class TNodeFacetNarrowPhaseResult, const int t_blockSize, const int t_numVerticesPerFacet>
  __global__ void UpdateTrackedNodeFacetContacts(int2 *p_dst, int *p_counter, const TNodeFacetNarrowPhaseResult *nodeFacetContacts, const int numContacts, const int *masterNodeFacetNeighbours, const int2 *masterNodeFacetNeighbourRanges) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ bool haveContact[t_blockSize];
    __shared__ int2 contacts[t_blockSize];
    int threadContactInd = 0;
    
    tledCudaAssert(t_blockSize == blockDim.x);
    haveContact[threadIdx.x] = false;
    __syncthreads();
    for (int c = 0; c < numContacts && !haveContact[threadIdx.x] && tid >= threadContactInd; c++) {
      const TNodeFacetNarrowPhaseResult &contact = nodeFacetContacts[c];
      const int slaveNodeInd = contact.ContactNodeIndices[0];

      for (int v = 0; v < t_numVerticesPerFacet && !haveContact[threadIdx.x] && tid >= threadContactInd; v++) {
	const int centralMasterNode = contact.ContactNodeIndices[1+v];
	const int numFacets = masterNodeFacetNeighbourRanges[centralMasterNode].y - masterNodeFacetNeighbourRanges[centralMasterNode].x;

	if (threadContactInd + numFacets > tid) {
	  const int currentRangeInd = tid - threadContactInd;
	  const int masterFacetInd = masterNodeFacetNeighbours[masterNodeFacetNeighbourRanges[centralMasterNode].x+currentRangeInd];

	  tledCudaAssert(currentRangeInd < numFacets && masterNodeFacetNeighbourRanges[centralMasterNode].x+currentRangeInd < masterNodeFacetNeighbourRanges[centralMasterNode].y);
	  tledCudaAssert(!haveContact[threadIdx.x]);
	  contacts[threadIdx.x].x = slaveNodeInd;
	  contacts[threadIdx.x].y = masterFacetInd;
	  tledCudaPrintf(3, "thread %d-%d: writing node-facet contact %d-%d\n", blockIdx.x, threadIdx.x, contacts[threadIdx.x].x, contacts[threadIdx.x].y);
	  haveContact[threadIdx.x] = true;
	} else {
	  threadContactInd += numFacets;
	}
      }
    }    
    
    _FlushContacts(p_dst, p_counter, contacts, haveContact);
  }

  template <class TEdgeEdgeNarrowPhaseResult, const int t_blockSize>
  __global__ void UpdateTrackedEdgeEdgeContacts(int2 *p_dst, int *p_counter, const TEdgeEdgeNarrowPhaseResult *edgeEdgeContacts, const int numContacts, const int *masterNodeEdgeNeigbours, const int2 *masterNodeEdgeNeigbourRanges, const int *slaveNodeEdgeNeigbours, const int2 *slaveNodeEdgeNeigbourRanges) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ bool haveContact[t_blockSize];
    __shared__ int2 contacts[t_blockSize];
    int threadContactInd = 0;
    
    haveContact[threadIdx.x] = false;
    for (int c = 0; c < numContacts && !haveContact[threadIdx.x] && tid >= threadContactInd; c++) {
      const TEdgeEdgeNarrowPhaseResult &contact = edgeEdgeContacts[c];

      for (int s = 0; s < 2 && !haveContact[threadIdx.x] && tid >= threadContactInd; s++) for (int m = 0; m < 2 && !haveContact[threadIdx.x] && tid >= threadContactInd; m++) {
	  const int masterNode = m == 0? contact.MasterEdge.x : contact.MasterEdge.y;
	  const int slaveNode = s == 0? contact.SlaveEdge.x : contact.SlaveEdge.y;
	  const int numMasterEdges = masterNodeEdgeNeigbourRanges[masterNode].y - masterNodeEdgeNeigbourRanges[masterNode].x;
	  const int numSlaveEdges = slaveNodeEdgeNeigbourRanges[slaveNode].y - slaveNodeEdgeNeigbourRanges[slaveNode].x;

	  if (numMasterEdges*numSlaveEdges + threadContactInd > tid && !haveContact[threadIdx.x]) {
	    const int currentRangeInd = tid - threadContactInd;
	    const int slaveEdgeInd = slaveNodeEdgeNeigbours[slaveNodeEdgeNeigbourRanges[slaveNode].x+currentRangeInd/numMasterEdges];
	    const int masterEdgeInd = masterNodeEdgeNeigbours[masterNodeEdgeNeigbourRanges[masterNode].x+currentRangeInd%numMasterEdges];

	    contacts[threadIdx.x].x = slaveEdgeInd;
	    contacts[threadIdx.x].x = masterEdgeInd;
	    tledCudaPrintf(3, "Writing edge-edge cand: %d - %d\n", slaveEdgeInd, masterEdgeInd);
	    haveContact[threadIdx.x] = true;
	  } else {
	    threadContactInd += numMasterEdges*numSlaveEdges;
	  }
	}      
    }

    _FlushContacts(p_dst, p_counter, contacts, haveContact);
  }
}

#endif
