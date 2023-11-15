// =========================================================================
// File:       tledOBB.cu
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
#ifndef tledOBB_CU
#define tledOBB_CU
#include "tledOBB.h"

template <class TBV>
static void _CopyBoundsGPU(typename TBV::GPUBV &r_gpuOBB, const TBV &bv) {
  r_gpuOBB.Extents = make_float3(bv.Extents[0], bv.Extents[1], bv.Extents[2]);
  r_gpuOBB.Centroid = make_float3(bv.Centroid[0], bv.Centroid[1], bv.Centroid[2]);
  for (int c = 0; c < 3; c++) r_gpuOBB.Axes[c] = make_float3(bv.Axes[c][0], bv.Axes[c][1], bv.Axes[c][2]);

  /* Device code not yet available! */
  tledFatalNotYetImplementedError;
}

template <>
void tledOBB<2>::InitGPU(Superclass::GPUBV &r_dst) {
  Superclass::InitGPU(r_dst);
  _CopyBoundsGPU(static_cast<GPUBV&>(r_dst), *this);  
}

template <>
void tledOBB<4>::InitGPU(Superclass::GPUBV &r_dst) {
  Superclass::InitGPU(r_dst);
  _CopyBoundsGPU(static_cast<GPUBV&>(r_dst), *this);  
}

#endif
