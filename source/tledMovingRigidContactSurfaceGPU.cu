// =========================================================================
// File:       tledMovingRigidContactSurfaceGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMovingRigidContactSurfaceGPU_CU
#define tledMovingRigidContactSurfaceGPU_CU

#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledCUDAHelpers.h"

#include "tledMovingRigidContactSurfaceGPU_kernels.cu"

void tledMovingRigidContactSurfaceGPU::TranslateNodes(float3 *dp_cds, float3 *dp_oldCds, const float3 *dpc_x0, const int numNodes, const float currentT[], const float oldT[]) {
  using namespace tledCUDAHelpers;

  const int blockSize = 256;
  const int numBlocks = GetNumberOfBlocks(numNodes, blockSize);
  const float3 tCurr3 = ConvertToFloat3(currentT), tOld3 = ConvertToFloat3(oldT);
  
  tledMovingRigidContactSurfaceGPU_kernels::TranslateNodesKernel <<<numBlocks, blockSize>>> (dp_cds, dp_oldCds, dpc_x0, numNodes, tCurr3, tOld3);
}


tledMovingRigidContactSurfaceGPU* tledMovingRigidContactSurfaceGPU::CreateSurface(const std::string &type) {
  tledMovingRigidContactSurfaceGPU *p_surface = NULL;

  if (type == "T3") {
    p_surface = new tledMovingRigidContactSurfaceT3GPU();
  } else if (type == "Q4") {
    p_surface = new tledMovingRigidContactSurfaceQ4GPU();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Surface type \"" << type << "\" not recognised.");
  }

  return p_surface;
}

#endif
