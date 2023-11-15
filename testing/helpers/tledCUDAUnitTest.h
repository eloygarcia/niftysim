// =========================================================================
// File:       tledCUDAUnitTest.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledCUDAUnitTest_H
#define tledCUDAUnitTest_H
#include "tledCUDAHelpers.h"

namespace tledCUDAUnitTest {
  __tled_inline void InitCUDATests(void) {
#ifdef _GPU_
#ifndef _CUDA_5PLUS_SDK
    cudaSetDevice(cutGetMaxGflopsDeviceId());
#else
    cudaSetDevice(gpuGetMaxGflopsDeviceId());
#endif
#endif
  }

  __tled_inline void FinishCUDATests(void) {
#if defined _GPU_ && !defined _CUDA_3MINUS_SDK
    tledCheckCUDAErrors(cudaDeviceReset());
#endif
  }
}
#endif
