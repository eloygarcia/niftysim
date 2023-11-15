// =========================================================================
// File:       tledDeformableContactSurfaceGPU_kernels.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableContactSurfaceGPU_kernels_TPP
#define tledDeformableContactSurfaceGPU_kernels_TPP

#include "tledContactSurfaceGPU_kernels.h"
#include "tledDeformableContactSurfaceGPU.h"

#include "tledCUDA_operators.cu"

namespace tledDeformableContactSurfaceGPU_kernels {
  template <class TGPUSurface>
  __device__ bool IsAdjacent(const TGPUSurface &mesh, const int f0Ind, const int f1Ind) {
    const int *f0 = mesh.Facets[f0Ind].NodeIndices;
    const int *f1 = mesh.Facets[f1Ind].NodeIndices;

    for (int const *pc_v0 = f0; pc_v0 < f0 + TGPUSurface::NumberOfFacetVertices; pc_v0++) for (int const *pc_v1 = f1; pc_v1 < f1 + TGPUSurface::NumberOfFacetVertices; pc_v1++) {
	if (*pc_v0 == *pc_v1) return true;
      }

    return false;
  }

  template <class TGPUSurface, const int t_blockSize, const int t_threadsPerNode, const int t_numFacetVtcs>
  __global__ void ComputeNodeNormals(TGPUSurface *p_mesh, const int *nodeList, const int numNodes) {
    const TGPUSurface &mesh = *p_mesh;
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    const int nInd = tid/t_threadsPerNode;
    
    __shared__ float3 normals[t_blockSize];

    normals[threadIdx.x] = make_float3(0.f, 0.f, 0.f);
    tledCudaAssert(nInd >= numNodes || (nodeList[nInd] >= 0 && nodeList[nInd] < mesh.NumberOfNodes));
    if (nInd < numNodes) {
      const int2 bounds = mesh.NodeFacetIndexRanges[nodeList[nInd]];
      const int numFacets = bounds.y - bounds.x;
      const int facetsPerThread = numFacets/t_threadsPerNode + (numFacets%t_threadsPerNode > 0);
      const int subBlockThreadInd = threadIdx.x%t_threadsPerNode;
      const int facetStartIndex = facetsPerThread*subBlockThreadInd;
      const int facetEndIndex = facetStartIndex + facetsPerThread > numFacets? numFacets : facetStartIndex + facetsPerThread;

      tledCudaAssert(bounds.x + facetEndIndex <= bounds.y);
      for (int f = bounds.x + facetStartIndex; f < bounds.x + facetEndIndex; f++) {
	const int *facetDef = mesh.Facets[mesh.NodeFacetList[f]].NodeIndices;

	float3 a, b;
	int vInd;
      
	for (vInd = 0; facetDef[vInd] != nodeList[nInd]; vInd++) {
	  tledCudaAssert(vInd < t_numFacetVtcs);
	}

	a = mesh.NodeCoordinates[facetDef[(vInd+1)%t_numFacetVtcs]] - mesh.NodeCoordinates[nodeList[nInd]];
	b = mesh.NodeCoordinates[facetDef[(vInd+2)%t_numFacetVtcs]] - mesh.NodeCoordinates[nodeList[nInd]];
	
	normals[threadIdx.x] += cross(a, b)*acosf(dot(a, b)/(norm(a)*norm(b)));	
	__syncthreads();
	if (subBlockThreadInd == 0) {
	  for (int i = 0; i < t_threadsPerNode; i++) {
	    tledCudaAssert(threadIdx.x + i < t_blockSize);
	    normals[threadIdx.x] += normals[threadIdx.x+i];
	  }
	  p_mesh->NodeNormals[nodeList[nInd]] = normals[threadIdx.x]/norm(normals[threadIdx.x]);
	}
      }
    }
  }

  template <class TGPUSurface, const int t_numFacetVtcs>
  __global__ void UpdateAllFacetNormals(float3 *p_facetNormalVtxAngles, const TGPUSurface *pc_mesh) {
    const TGPUSurface &mesh = *pc_mesh;
    const int numFacets = mesh.NumberOfFacets;
    const int fInd = blockIdx.x*blockDim.x + threadIdx.x;

    if (fInd < numFacets) {
      const int *facetDef = mesh.Facets[fInd].NodeIndices;

      float3 x, y, n;
      float tAngle;

      if (t_numFacetVtcs == 3) {
	x = mesh.NodeCoordinates[facetDef[1]] - mesh.NodeCoordinates[facetDef[0]];
	y = mesh.NodeCoordinates[facetDef[2]] - mesh.NodeCoordinates[facetDef[0]];
      } else {
	x = mesh.NodeCoordinates[facetDef[2]] - mesh.NodeCoordinates[facetDef[0]];
	y = mesh.NodeCoordinates[facetDef[3]] - mesh.NodeCoordinates[facetDef[1]];	
      }
      n = cross(x, y);

      tAngle = 0.f;
      for (int v = 0; v < t_numFacetVtcs - 1; v++) {
	float alpha;

	x = mesh.NodeCoordinates[facetDef[(v+1)%t_numFacetVtcs]] - mesh.NodeCoordinates[facetDef[v]];
	y = mesh.NodeCoordinates[facetDef[(v+2)%t_numFacetVtcs]] - mesh.NodeCoordinates[facetDef[v]];
	alpha = acosf(dot(x, y)/(norm(x)*norm(y)));
	p_facetNormalVtxAngles[fInd*t_numFacetVtcs+v] = alpha*n;
	tAngle += alpha;	
      }
      tledCudaAssert(tAngle != tAngle || tAngle < tledPi*(t_numFacetVtcs - 2));
      p_facetNormalVtxAngles[(fInd+1)*t_numFacetVtcs-1] = (tledPi*(t_numFacetVtcs - 2) - tAngle)*n;
    }
  }

  template <class TGPUSurface, const int t_numFacetVtcs>
  __global__ void UpdateAllNodeNormals(TGPUSurface *p_mesh, const float3 *facetNormalVtxAngles) {
    const TGPUSurface &mesh = *p_mesh;
    const int numNodes = mesh.NumberOfNodes;
    const int nInd = blockIdx.x*blockDim.x + threadIdx.x;

    if (nInd < numNodes) {
      const int2 facetIndexLookUpBounds = mesh.NodeFacetIndexRanges[nInd];
    
      float3 n = make_float3(0.f, 0.f, 0.f);
      float3 *p_nodeNormals = p_mesh->NodeNormals;
      
      for (int f = facetIndexLookUpBounds.x; f < facetIndexLookUpBounds.y; f++) {
	const int fInd = mesh.NodeFacetList[f];
	const int *facetDef = mesh.Facets[fInd].NodeIndices;
	
	int v = 0;
	
	for (; v < t_numFacetVtcs && facetDef[v] != nInd; v++);
	tledCudaAssert(v < t_numFacetVtcs && facetDef[v] == nInd);
	n += facetNormalVtxAngles[fInd*t_numFacetVtcs+v];
      }
      p_nodeNormals[nInd] = n/norm(n);
    }
  }

  template <class TSurface>
  __device__ float3 ComputeOldNormal(const TSurface &mesh, const int fInd) {
    float3 n;

    n = cross(mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[1]] - mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[0]], mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[2]] - mesh.OldNodeCoordinates[mesh.Facets[fInd].NodeIndices[0]]);    

    return n/norm(n);
  }

  template <class TGPUSurface>
  __global__ void SaveNodeCoordinates(TGPUSurface *p_surface, float3 *p_head, float3 *p_old) {
    const int nInd = blockIdx.x*blockDim.x + threadIdx.x;
    const int numNodes = p_surface->NumberOfNodes;
    
    if (nInd < numNodes) {
      p_head[nInd] = p_surface->NodeCoordinates[nInd];
    }

    if (blockIdx.x == 0 && threadIdx.x == 0) {
      p_surface->OldNodeCoordinates = p_old;
    }
  }
}

#endif
