// =========================================================================
// File:       tledNewmarkTimeStepperGPU.cu
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
#ifndef tledNewmarkTimeStepperGPU_CU
#define tledNewmarkTimeStepperGPU_CU

#include "tledNewmarkTimeStepperGPU.h"
#include "tledCUDAHelpers.h"

#include <algorithm>

#include "tledNewmarkTimeStepperGPU_kernels.cu"

tledNewmarkTimeStepperGPU::tledNewmarkTimeStepperGPU(const int numNodes, const float dt, const float alpha, const float M[]) : Superclass(numNodes, dt) {
  tledCUDAHelpers::AllocateDeviceMemory(mdp_PreviousAccelerations, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_CurrentAccelerations, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_PreviousVelocities, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_CurrentVelocities, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_M, numNodes);

  tledCUDAHelpers::AllocateDeviceMemory(_GetOnDeviceRepresentation());
  
  tledCUDAHelpers::CopyToDevice(mdp_M, M, numNodes);
  tledCheckCUDAErrors(cudaBindTexture(0, tledNewmarkTimeStepper_kernels::tx_M, mdp_M, this->GetNumberOfNodes()*sizeof(float)));    
  tledCheckCUDAErrors(cudaMemcpyToSymbol(tledNewmarkTimeStepper_kernels::c_Alpha, &alpha, sizeof(float)));
  tledCheckCUDAErrors(cudaMemcpyToSymbol(tledNewmarkTimeStepper_kernels::c_Dt, &dt, sizeof(float)));

  tledCheckCUDAErrors(cudaMemset(mdp_PreviousVelocities, 0, sizeof(float4)*this->GetNumberOfNodes()));
  tledCheckCUDAErrors(cudaMemset(mdp_PreviousAccelerations, 0, sizeof(float4)*this->GetNumberOfNodes()));

  _UpdateDeviceDataStructure();
}

tledNewmarkTimeStepperGPU::SolutionVariables::SolutionVariables(float4 *p_uCurr, float4 *p_uNext, float4 *p_aCurr, float4 *p_aPrev, float4 *p_vCurr, float4 *p_vPrev) {
  CurrentDisplacements = p_uCurr;
  NextDisplacements = p_uNext;
  CurrentAccelerations = p_aCurr;
  PreviousAccelerations = p_aPrev;
  CurrentVelocities = p_vCurr;
  PreviousVelocities = p_vPrev;
}

void tledNewmarkTimeStepperGPU::_UpdateDeviceDataStructure() {
  SolutionVariables hostBuffer(this->GetOnDeviceCurrentDisplacements(), this->GetOnDeviceNextDisplacements(), mdp_CurrentAccelerations, mdp_PreviousAccelerations, mdp_CurrentVelocities, mdp_PreviousVelocities);

  tledCUDAHelpers::CopyToDevice(_GetOnDeviceRepresentation(), &hostBuffer);
}

void tledNewmarkTimeStepperGPU::RetrieveSolutionFromDevice() {
  Superclass::RetrieveSolutionFromDevice();
  ConvertToHostFloatArray(this->GetCurrentVelocities(), this->GetHostDeviceCopyBuffer(), mdp_CurrentVelocities, this->GetNumberOfNodes());
  ConvertToHostFloatArray(this->GetCurrentAccelerations(), this->GetHostDeviceCopyBuffer(), mdp_CurrentAccelerations, this->GetNumberOfNodes());
}

void tledNewmarkTimeStepperGPU::FinishTimeStep() {
  Superclass::FinishTimeStep();
  std::iter_swap(&mdp_CurrentAccelerations, &mdp_PreviousAccelerations);
  std::iter_swap(&mdp_CurrentVelocities, &mdp_PreviousVelocities);
  _UpdateDeviceDataStructure();
}

tledNewmarkTimeStepperGPU::~tledNewmarkTimeStepperGPU() {
  cudaFree(mdp_M);
  cudaFree(mdp_PreviousAccelerations);
  cudaFree(mdp_PreviousVelocities);
  cudaFree(mdp_CurrentAccelerations);
  cudaFree(mdp_CurrentVelocities);
  cudaFree(_GetOnDeviceRepresentation());
}

#endif
