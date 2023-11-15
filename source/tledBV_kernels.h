// =========================================================================
// File:       tledBV_kernels.h
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
#ifndef tledBV_kernels_H
#define tledBV_kernels_H

namespace tledBV_kernels {
  template <class TBV0, class TBV1>
  __device__ bool DoIntersect(const TBV0 &bv0, const TBV1 &bv1);

  template <class TBV>
  __device__ bool IsLeaf(const TBV &bv);

  /** Refits a BVH leaf node */
  template <class TGPUBV, class TGPUMesh>
  __device__ void RefitLeaf(TGPUBV &r_bv, const TGPUMesh &mesh, const float margin);

  /** Refits an interior BVH node without performing any validity checks on bounds of the node's children */
  template <class TGPUBV>
  __device__ void RefitInterior(TGPUBV *p_bvs, const int bvInd);

  /** Applies translation t to a BV */
  template <class TGPUBV>
  __device__ void Translate(TGPUBV &r_bv, const float3 &t);
}

#include "tledBV_kernels.tpp"

#endif
