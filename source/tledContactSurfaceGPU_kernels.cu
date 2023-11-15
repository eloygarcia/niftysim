// =========================================================================
// File:       tledContactSurfaceGPU_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSurfaceGPU_kernels_CU
#define tledContactSurfaceGPU_kernels_CU

namespace tledContactSurfaceGPU_kernels {
  template <class TSurface>
  __device__ void ComputeShapeValues(float *p_out, const float3 &xi) {
    if (TSurface::NumberOfFacetVertices == 3) {
      p_out[0] = 1 - xi.x - xi.y;
      p_out[1] = xi.x;
      p_out[2] = xi.y;
    } 
  }

  __device__ bool ProjectOntoFacet(float3 &r_xi, const float3 &x, const float4 projOp[]) {
    const float facetMargin = 1e-4f;
    const float4 &R0 = projOp[0];
    const float4 &R1 = projOp[1];
    const float4 &R2 = projOp[2];

    r_xi.x = R0.x*x.x + R0.y*x.y + R0.z*x.z + R0.w;
    if (r_xi.x < -facetMargin || r_xi.x > 1 + facetMargin) return false;
    r_xi.y = R1.x*x.x + R1.y*x.y + R1.z*x.z + R1.w;
    if (r_xi.y < -facetMargin || r_xi.y + r_xi.x > 1 + facetMargin) return false;
    else {
      float sum;

      r_xi.x = fmaxf(r_xi.x, 0);
      r_xi.y = fmaxf(r_xi.y, 0);
      if ((sum = r_xi.x + r_xi.y) > 1) {
	r_xi.x /= sum;
	r_xi.y /= sum;
      }

      r_xi.z = R2.x*x.x + R2.y*x.y + R2.z*x.z + R2.w;

      return true;
    } /* if not inside facet else ... */
  }

  template <class TSurface>
  __device__ void ComputeProjectionOperator(float4 *p_projOp, const TSurface &mesh, const int fInd) {
    float3 normal, a, b, facetVtcs[3];
    float inDet, normN;

    for (int i = 0; i < 3; i++) facetVtcs[i] = mesh.NodeCoordinates[mesh.Facets[fInd].NodeIndices[i]];
    a = facetVtcs[1] - facetVtcs[0];
    b = facetVtcs[2] - facetVtcs[0];
    normal = cross(a, b);
    tledCudaAssert(norm(normal) > 1e-6);
      
    inDet = 1/(normN = dot(normal, normal));
    normN = sqrtf(normN);

    p_projOp[0].x = (b.y*normal.z - b.z*normal.y)*inDet;
    p_projOp[0].y = (b.z*normal.x - b.x*normal.z)*inDet;
    p_projOp[0].z = (b.x*normal.y - b.y*normal.x)*inDet;
      
    p_projOp[1].x = (normal.y*a.z - normal.z*a.y)*inDet;
    p_projOp[1].y = (normal.z*a.x - normal.x*a.z)*inDet;
    p_projOp[1].z = (normal.x*a.y - normal.y*a.x)*inDet;
      
    p_projOp[2].x = normal.x*inDet*normN;
    p_projOp[2].y = normal.y*inDet*normN;
    p_projOp[2].z = normal.z*inDet*normN;

    p_projOp[0].w = -(p_projOp[0].x*facetVtcs[0].x + p_projOp[0].y*facetVtcs[0].y + p_projOp[0].z*facetVtcs[0].z);
    p_projOp[1].w = -(p_projOp[1].x*facetVtcs[0].x + p_projOp[1].y*facetVtcs[0].y + p_projOp[1].z*facetVtcs[0].z);
    p_projOp[2].w = -(p_projOp[2].x*facetVtcs[0].x + p_projOp[2].y*facetVtcs[0].y + p_projOp[2].z*facetVtcs[0].z);
  }

  template <class TSurface>
  __device__ const float3& GetNodeNormal(const TSurface &mesh, const int nodeIndex) {
    return mesh.NodeNormals[nodeIndex];
  }
}
#endif
