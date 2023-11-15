// =========================================================================
// File:       tledUnstructuredContactManager_kernels.cu
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
#ifndef tledUnstructuredContactManager_kernels_CU
#define tledUnstructuredContactManager_kernels_CU

namespace tledUnstructuredContactManager_kernels {
  __device__ __constant__ float c_SurfaceCloseDistance;

  __device__ float GetCloseDistance() {
    return c_SurfaceCloseDistance;
  }
}

#endif
