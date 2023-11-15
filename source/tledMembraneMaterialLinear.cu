// =========================================================================
// File:       tledMembraneMaterialLinear.cu
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
#ifndef tledMembraneMaterialLinear_CU
#define tledMembraneMaterialLinear_CU

#include "tledMembraneMaterialLinear.h"

void tledMembraneMaterialLinear::InitHostGPURepresentation(tledShellMaterial::GPUMaterial &r_sdst) const {
  GPUMaterial &r_dst = static_cast<GPUMaterial&>(r_sdst);

  Superclass::InitHostGPURepresentation(r_dst);
  r_dst.E = m_E;
  r_dst.Nu = m_Nu;
} 

tledShellMaterial::GPUMaterial* tledMembraneMaterialLinear::InitGPU() const {
  GPUMaterial hostMem;
  GPUMaterial *dp_mem;

  InitHostGPURepresentation(hostMem);
  cudaMalloc((void**)&dp_mem, sizeof(GPUMaterial));
  cudaMemcpy(dp_mem, &hostMem, sizeof(GPUMaterial), cudaMemcpyHostToDevice);

  return dp_mem;
}

#endif
