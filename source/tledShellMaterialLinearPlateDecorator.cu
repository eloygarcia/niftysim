// =========================================================================
// File:       tledShellMaterialLinearPlateDecorator.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    August 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterialLinearPlateDecorator_CU
#define tledShellMaterialLinearPlateDecorator_CU

#include "tledShellMaterialLinearPlateDecorator.h"
#include "tledMembraneMaterialLinear.h"

template <>
tledShellMaterial::GPUMaterial* tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear>::InitGPU() const {
  GPUMaterial hostMem;
  GPUMaterial *dp_mem;
  
  InitHostGPURepresentation(hostMem);
  cudaMalloc((void**)&dp_mem, sizeof(GPUMaterial));
  cudaMemcpy(dp_mem, &hostMem, sizeof(GPUMaterial), cudaMemcpyHostToDevice);

  return dp_mem;
}

#endif
