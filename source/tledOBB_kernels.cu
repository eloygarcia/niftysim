// =========================================================================
// File:       tledOBB_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledOBB_kernels_CU
#define tledOBB_kernels_CU

#include "tledOBB.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledBV_kernels.h"

namespace tledOBB_kernels {
  template <class TOBB>
  __device__ bool DoIntersectOBBOBB(const TOBB &obb0, const TOBB &obb1) {
    tledCudaAssert(false);
    return false;
  }

  template <class TOBB>
  __device__ void MergeOBBs(TOBB &r_bv0, const TOBB &bv1) {
    tledCudaAssert(false);
  }

  __device__ void AddMargin(float3 &r_extents, const float margin) {
    r_extents = r_extents + make_float3(margin, margin, margin);
  } 

  template <class TOBB> 
  __device__ void CopyBounds(TOBB &r_dst, const TOBB &src) {
    r_dst.Extents = src.Extents;
    r_dst.Centroid = src.Centroid;
    for (int c = 0; c < 3; c++) r_dst.Axes[c] = src.Axes[c];
  }

  template <class TOBB>
  __device__ void Translate(TOBB &r_bv, const float3 &t) {
    r_bv.Centroid += t;
  }
}

namespace tledBV_kernels {
  template <>
  __device__ bool DoIntersect<tledOBB<2>::GPUBV, tledOBB<2>::GPUBV>(const tledOBB<2>::GPUBV &obb0, const tledOBB<2>::GPUBV &obb1) {
    return tledOBB_kernels::DoIntersectOBBOBB(obb0, obb1);
  }

  template <>
  __device__ void RefitLeaf<tledOBB<2>::GPUBV, tledDeformableContactSurfaceT3GPU::GPUSurface>(tledOBB<2>::GPUBV &r_bv, const tledDeformableContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledOBB_kernels::AddMargin(r_bv.Extents, margin);
  }

  template <>
  __device__ void RefitLeaf<tledOBB<2>::GPUBV, tledMovingRigidContactSurfaceT3GPU::GPUSurface>(tledOBB<2>::GPUBV &r_bv, const tledMovingRigidContactSurfaceT3GPU::GPUSurface &mesh, const float margin) {
    tledOBB_kernels::AddMargin(r_bv.Extents, margin);
  }

  template <>
  __device__ void RefitInterior<tledOBB<2>::GPUBV>(tledOBB<2>::GPUBV *p_bvs, const int bvInd) {
    tledOBB<2>::GPUBV &r_bv = p_bvs[bvInd];

    tledOBB_kernels::CopyBounds(r_bv, p_bvs[r_bv.ChildIndices[0]]);
    tledOBB_kernels::MergeOBBs(r_bv, p_bvs[r_bv.ChildIndices[1]]);
  }

  template <>
  __device__ void Translate<tledOBB<2>::GPUBV>(tledOBB<2>::GPUBV &r_bv, const float3 &t) {
    tledOBB_kernels::Translate(r_bv, t);
  }

  template <>
  __device__ void Translate<tledOBB<4>::GPUBV>(tledOBB<4>::GPUBV &r_bv, const float3 &t) {
    tledOBB_kernels::Translate(r_bv, t);
  }
}

#endif
