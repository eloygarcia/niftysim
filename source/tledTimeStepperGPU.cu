// =========================================================================
// File:       tledTimeStepperGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledTimeStepperGPU_CU
#define tledTimeStepperGPU_CU

#include "tledTimeStepperGPU.h"
#include "tledDeviceDeclarations.h"
#include "tledCUDAHelpers.h"

#include <algorithm>

void tledTimeStepperGPU::ConvertToHostFloatArray(float *p_outBuffer, float4 *hp_tmpBuffer, const float4 *dp_deviceArray, const int numNodes) {
  float4 const *pc_nv4;
  float *p_out;

  tledCUDAHelpers::CopyFromDevice(hp_tmpBuffer, dp_deviceArray, numNodes);
  for (p_out = p_outBuffer, pc_nv4 = hp_tmpBuffer; pc_nv4 < hp_tmpBuffer + numNodes; pc_nv4++) {
    *(p_out++) = pc_nv4->x;
    *(p_out++) = pc_nv4->y;
    *(p_out++) = pc_nv4->z;
  }
}

void tledTimeStepperGPU::RetrieveSolutionFromDevice() {
  ConvertToHostFloatArray(GetCurrentDisplacements(), mp_HostDeviceCopyBuffer, mdp_CurrentDisplacements, GetNumberOfNodes());
  ConvertToHostFloatArray(GetNextDisplacements(), mp_HostDeviceCopyBuffer, mdp_NextDisplacements, GetNumberOfNodes());
}

void tledTimeStepperGPU::FinishTimeStep() {
  std::iter_swap(&mdp_NextDisplacements, &mdp_CurrentDisplacements);
  tledCheckCUDAErrors(cudaBindTexture(0, txUcurr, mdp_CurrentDisplacements, GetNumberOfNodes()*sizeof(float4)));  
}

void tledTimeStepperGPU::SetCurrentDisplacements(const float U[]) {
  float const *pc_u = U;
  float4 *p_out = mp_HostDeviceCopyBuffer;

  for (; p_out < mp_HostDeviceCopyBuffer + GetNumberOfNodes(); p_out++) {
    p_out->x = *(pc_u++);
    p_out->y = *(pc_u++);
    p_out->z = *(pc_u++);
    p_out->w = 0;
  }

  tledCUDAHelpers::CopyToDevice(mdp_CurrentDisplacements, mp_HostDeviceCopyBuffer, this->GetNumberOfNodes());
}

tledTimeStepperGPU::tledTimeStepperGPU(const int numNodes) : Superclass(numNodes) {
  tledCUDAHelpers::AllocateDeviceMemory(mdp_CurrentDisplacements, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_NextDisplacements, numNodes);
  tledCUDAHelpers::AllocateHostMemory(mp_HostDeviceCopyBuffer, numNodes);

  tledCheckCUDAErrors(cudaMemset(mdp_CurrentDisplacements, 0, sizeof(float4)*GetNumberOfNodes()));
  tledCheckCUDAErrors(cudaBindTexture(0, txUcurr, mdp_CurrentDisplacements, GetNumberOfNodes()*sizeof(float4)));
}

tledTimeStepperGPU::~tledTimeStepperGPU() {
  cudaFree(mdp_NextDisplacements);
  cudaFree(mdp_CurrentDisplacements);
  cudaFreeHost(mp_HostDeviceCopyBuffer);
}

/* Gather CU files of all time stepper modules for inclusion in tledSolverGPU here: */
#include "tledCentralDifferenceTimeStepperGPU.cu"
#include "tledNewmarkTimeStepperGPU.cu"

#endif
