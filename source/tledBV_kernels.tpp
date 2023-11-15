// =========================================================================
// File:       tledBV_kernels.cu
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
#ifndef tledBV_kernels_CU
#define tledBV_kernels_CU

#include "tledBV.h"

namespace tledBV_kernels {
    template <class TBV>
    __device__ bool IsLeaf(const TBV &bv) {
      return bv.PrimitiveIndex >= 0;
    }
}

#endif
