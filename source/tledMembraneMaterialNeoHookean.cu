// =========================================================================
// File:       tledMembraneMaterialNeoHookean.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMembraneMaterialNeoHookean_CU
#define tledMembraneMaterialNeoHookean_CU

#include "tledMembraneMaterialNeoHookean.h"

tledShellMaterial::GPUMaterial* tledMembraneMaterialNeoHookean::InitGPU() const {
  GPUMaterial hostMem;
  GPUMaterial *dp_mem;

  hostMem.Mu = m_Mu;
  hostMem.Thickness = GetThickness();
  hostMem.Rho = GetDensity();
  
  cudaMalloc((void**)&dp_mem, sizeof(GPUMaterial));
  cudaMemcpy(dp_mem, &hostMem, sizeof(GPUMaterial), cudaMemcpyHostToDevice);

  return dp_mem;
}
#endif
