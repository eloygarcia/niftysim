// =========================================================================
// File:       tledAABB_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledAABB_kernels_CU
#define tledAABB_kernels_CU

#include "tledBV_kernels.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledAABB.h"

namespace tledAABB_kernels {
  inline __device__ bool _DoIntersectAABBAABB(const float2 *aabb0Bnds, const float2 *aabb1Bnds) {
    int c;

    for (c = 0; c < 3 && aabb0Bnds[c].y >= aabb1Bnds[c].x && aabb1Bnds[c].y >= aabb0Bnds[c].x; c++);

    return c == 3;
  }  

  inline __device__ void Translate(float2 *p_bounds, const float3 &t) {
    p_bounds[0].x += t.x, p_bounds[0].y += t.x;
    p_bounds[1].x += t.y, p_bounds[1].y += t.y;
    p_bounds[2].x += t.z, p_bounds[2].y += t.z;
  }

  inline __device__ void ExpandWithNodeList(float2 *p_bounds, const int *nodeListBegin, const int *nodeListEnd, const float3 *nodes) {
    for (int const *pc_n = nodeListBegin; pc_n < nodeListEnd; pc_n++) {
      p_bounds[0].x = fminf(p_bounds[0].x, nodes[*pc_n].x);
      p_bounds[0].y = fmaxf(p_bounds[0].y, nodes[*pc_n].x);

      p_bounds[1].x = fminf(p_bounds[1].x, nodes[*pc_n].y);
      p_bounds[1].y = fmaxf(p_bounds[1].y, nodes[*pc_n].y);

      p_bounds[2].x = fminf(p_bounds[2].x, nodes[*pc_n].z);
      p_bounds[2].y = fmaxf(p_bounds[2].y, nodes[*pc_n].z);
    }
  }

  inline __device__ void ComputeFromNodeList(float2 *p_bounds, const int *nodeListBegin, const int *nodeListEnd, const float3 *nodes) {
    p_bounds[0].x = p_bounds[0].y = nodes[*nodeListBegin].x;
    p_bounds[1].x = p_bounds[1].y = nodes[*nodeListBegin].y;
    p_bounds[2].x = p_bounds[2].y = nodes[*nodeListBegin].z;

    ExpandWithNodeList(p_bounds, nodeListBegin + 1, nodeListEnd, nodes);
  }

  inline __device__ void MergeAABBs(float2 *p_bounds, const float2 *cBounds0, const float2 *cBounds1) {
    for (int c = 0; c < 3; c++) {
      p_bounds[c].x = fminf(cBounds0[c].x, cBounds1[c].x);
      p_bounds[c].y = fmaxf(cBounds0[c].y, cBounds1[c].y);
    }
  }

  inline __device__ void AddMargin(float2 *p_bounds, const float margin) {
    for (int c = 0; c < 3; c++) {
      p_bounds[c].x -= margin;
      p_bounds[c].y += margin;
    }
  }

  inline __device__ void CopyBounds(float2 *p_dst, const float2 *srcBnds) {
    for (int c = 0; c < 3; c++) p_dst[c] = srcBnds[c];
  }

  template <class TSurface>
  __device__ void RefitLeafFromDynamicSurface(float2 *p_bounds, const TSurface &mesh, const int primitiveIndex, const float margin) {
    const int *nodeIndices = mesh.Facets[primitiveIndex].NodeIndices;

    ComputeFromNodeList(p_bounds, nodeIndices, nodeIndices + 3, mesh.NodeCoordinates);
    ExpandWithNodeList(p_bounds, nodeIndices, nodeIndices + 3, mesh.OldNodeCoordinates);
    AddMargin(p_bounds, margin);

#ifdef __CUDA_DEBUG
    for (int const *pc_n = nodeIndices; pc_n < nodeIndices + 3; pc_n++) {
      tledCudaAssert(p_bounds[0].x < mesh.NodeCoordinates[*pc_n].x && p_bounds[0].y > mesh.NodeCoordinates[*pc_n].x);
      tledCudaAssert(p_bounds[0].x < mesh.OldNodeCoordinates[*pc_n].x && p_bounds[0].y > mesh.OldNodeCoordinates[*pc_n].x);

      tledCudaAssert(p_bounds[1].x < mesh.NodeCoordinates[*pc_n].y && p_bounds[1].y > mesh.NodeCoordinates[*pc_n].y);
      tledCudaAssert(p_bounds[1].x < mesh.OldNodeCoordinates[*pc_n].y && p_bounds[1].y > mesh.OldNodeCoordinates[*pc_n].y);

      tledCudaAssert(p_bounds[2].x < mesh.NodeCoordinates[*pc_n].z && p_bounds[2].y > mesh.NodeCoordinates[*pc_n].z);
      tledCudaAssert(p_bounds[2].x < mesh.OldNodeCoordinates[*pc_n].z && p_bounds[2].y > mesh.OldNodeCoordinates[*pc_n].z);
    }
#endif
  }
}

namespace tledBV_kernels {
  template <>
  __device__ bool DoIntersect<tledAABB<2>::GPUBV, tledAABB<2>::GPUBV>(const tledAABB<2>::GPUBV &aabb0, const tledAABB<2>::GPUBV &aabb1) {
    return tledAABB_kernels::_DoIntersectAABBAABB(aabb0.Bounds, aabb1.Bounds);
  }

  template <>
  __device__ bool DoIntersect<tledAABB<4>::GPUBV, tledAABB<4>::GPUBV>(const tledAABB<4>::GPUBV &aabb0, const tledAABB<4>::GPUBV &aabb1) {
    return tledAABB_kernels::_DoIntersectAABBAABB(aabb0.Bounds, aabb1.Bounds);
  }

  template <>
  __device__ void RefitLeaf(tledAABB<2>::GPUBV &r_bv, const tledDeformableContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledAABB_kernels::RefitLeafFromDynamicSurface(r_bv.Bounds, mesh, r_bv.PrimitiveIndex, margin);
  }

  template <>
  __device__ void RefitLeaf(tledAABB<2>::GPUBV &r_bv, const tledMovingRigidContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledAABB_kernels::RefitLeafFromDynamicSurface(r_bv.Bounds, mesh, r_bv.PrimitiveIndex, margin);
  }

  template <>
  __device__ void RefitLeaf(tledAABB<4>::GPUBV &r_bv, const tledDeformableContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledAABB_kernels::RefitLeafFromDynamicSurface(r_bv.Bounds, mesh, r_bv.PrimitiveIndex, margin);
  }

  template <>
  __device__ void RefitLeaf(tledAABB<4>::GPUBV &r_bv, const tledMovingRigidContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledAABB_kernels::RefitLeafFromDynamicSurface(r_bv.Bounds, mesh, r_bv.PrimitiveIndex, margin);
  }

  template <>
  __device__ void RefitInterior(tledAABB<2>::GPUBV *p_bvs, const int bvInd) {
    tledAABB<2>::GPUBV &r_bv = p_bvs[bvInd];

    tledAABB_kernels::MergeAABBs(r_bv.Bounds, p_bvs[r_bv.ChildIndices[0]].Bounds, p_bvs[r_bv.ChildIndices[1]].Bounds);
  }

  template <>
  __device__ void RefitInterior(tledAABB<4>::GPUBV *p_bvs, const int bvInd) {
    tledAABB<4>::GPUBV &r_bv = p_bvs[bvInd];
    float2 *p_dstBnds = r_bv.Bounds;

    tledCudaAssert(r_bv.ChildIndices[0] >= 0 && r_bv.ChildIndices[1] >= 0);
    tledAABB_kernels::MergeAABBs(p_dstBnds, p_bvs[r_bv.ChildIndices[0]].Bounds, p_bvs[r_bv.ChildIndices[1]].Bounds);
    for (int const *pc_c = r_bv.ChildIndices + 2; *pc_c >= 0 && pc_c - r_bv.ChildIndices < 4; pc_c++) {
      tledAABB_kernels::MergeAABBs(p_dstBnds, p_dstBnds, p_bvs[*pc_c].Bounds);
    }
  }

  template <>
  __device__ void Translate(tledAABB<2>::GPUBV &r_bv, const float3 &t) {
    tledAABB_kernels::Translate(r_bv.Bounds, t);
  }

  template <>
  __device__ void Translate(tledAABB<4>::GPUBV &r_bv, const float3 &t) {
    tledAABB_kernels::Translate(r_bv.Bounds, t);
  }
}

#endif
