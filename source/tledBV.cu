// =========================================================================
// File:       tledBV.cu
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
#ifndef tledBV_CU
#define tledBV_CU

#include "tledBV.h"

template <class TBV>
static void _InitGPUImpl(typename TBV::GPUBV &r_dst, const TBV &bv) {
  r_dst.PrimitiveIndex = bv.PrimitiveIndex;
  r_dst.ParentIndex = bv.ParentIndex;
  std::copy(bv.ChildIndices, bv.ChildIndices + TBV::NumberOfChildBVs, r_dst.ChildIndices);
}

template <>
void tledBV<2>::InitGPU(GPUBV &r_dst) {
  _InitGPUImpl(r_dst, *this);
}

template <>
void tledBV<4>::InitGPU(GPUBV &r_dst) {
  _InitGPUImpl(r_dst, *this);
}

#endif
