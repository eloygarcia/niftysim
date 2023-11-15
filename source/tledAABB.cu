// =========================================================================
// File:       tledAABB.cu
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
#ifndef tledAABB_CU
#define tledAABB_CU

#include "tledAABB.h"

template <class TBV>
static void _CopyBoundsGPU(typename TBV::GPUBV &r_dst, const TBV &bv) {
  /* Specialisation is ugly but necessary due to severe linker issues in conjuction w/ tledSelfCollisionBV wrapper */
  for (int c = 0; c < 3; c++) {
    r_dst.Bounds[c].x = bv.Bounds[c][0];
    r_dst.Bounds[c].y = bv.Bounds[c][1];
  }
}

template <>
void tledAABB<2>::InitGPU(Superclass::GPUBV &r_dst) {
  Superclass::InitGPU(r_dst);
  _CopyBoundsGPU(static_cast<GPUBV&>(r_dst), *this);
}


template <>
void tledAABB<4>::InitGPU(Superclass::GPUBV &r_dst) {
  Superclass::InitGPU(r_dst);
  _CopyBoundsGPU(static_cast<GPUBV&>(r_dst), *this);
}

#endif
