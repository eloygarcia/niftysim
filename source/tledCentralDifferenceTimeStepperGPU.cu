// =========================================================================
// File:       tledCentralDifferenceTimeStepperGPU.cu
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
#ifndef tledCentralDifferenceTimeStepperGPU_CU
#define tledCentralDifferenceTimeStepperGPU_CU

#include "tledCentralDifferenceTimeStepperGPU.h"
#include "tledCUDAHelpers.h"

#include <algorithm>

#include "tledCentralDifferenceTimeStepperGPU_kernels.cu"

void tledCentralDifferenceTimeStepperGPU::RetrieveSolutionFromDevice() {
  Superclass::RetrieveSolutionFromDevice();
  ConvertToHostFloatArray(GetPreviousDisplacements(), this->GetHostDeviceCopyBuffer(), mdp_PreviousDisplacements, this->GetNumberOfNodes());
}

tledCentralDifferenceTimeStepperGPU::SolutionVariables::SolutionVariables(float4 *dp_UCurr, float4 *dp_UNext, float4 *dp_UPrev) {
  CurrentDisplacements = dp_UCurr;
  NextDisplacements = dp_UNext;
  PreviousDisplacements = dp_UPrev;
}

void tledCentralDifferenceTimeStepperGPU::FinishTimeStep() {
  Superclass::FinishTimeStep();
  std::iter_swap(&mdp_PreviousDisplacements, &this->GetOnDeviceNextDisplacements());

  _UpdateDeviceDataStructure();
}

void tledCentralDifferenceTimeStepperGPU::_UpdateDeviceDataStructure() {
  SolutionVariables hostBuffer(this->GetOnDeviceCurrentDisplacements(), this->GetOnDeviceNextDisplacements(), mdp_PreviousDisplacements);

  tledCUDAHelpers::CopyToDevice(_GetOnDeviceRepresentation(), &hostBuffer);
}

tledCentralDifferenceTimeStepperGPU::tledCentralDifferenceTimeStepperGPU(const int numNodes, const float dt, const float alpha, const float M[]) : Superclass(numNodes) {
  tledCUDAHelpers::AllocateDeviceMemory(mdp_PreviousDisplacements, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_CD, numNodes);
  tledCUDAHelpers::AllocateDeviceMemory(_GetOnDeviceRepresentation());

  for (int node = 0; node < numNodes; node++) {
    this->GetHostDeviceCopyBuffer()[node].x = 1/(M[node]/dt/dt + alpha*M[node]/2/dt);
    this->GetHostDeviceCopyBuffer()[node].y = 2*M[node]*this->GetHostDeviceCopyBuffer()[node].x/dt/dt;
    this->GetHostDeviceCopyBuffer()[node].z = alpha*M[node]*this->GetHostDeviceCopyBuffer()[node].x/2/dt - this->GetHostDeviceCopyBuffer()[node].y/2;
  }  
  tledCUDAHelpers::CopyToDevice(mdp_CD, this->GetHostDeviceCopyBuffer(), numNodes);
  tledCheckCUDAErrors(cudaBindTexture(0, tledCentralDifferenceTimeStepper_kernels::tx_CD, mdp_CD, GetNumberOfNodes()*sizeof(float4)));

  tledCheckCUDAErrors(cudaMemset(mdp_PreviousDisplacements, 0, sizeof(float4)*GetNumberOfNodes()));
  _UpdateDeviceDataStructure();  
}

tledCentralDifferenceTimeStepperGPU::~tledCentralDifferenceTimeStepperGPU() {
  cudaFree(mdp_CD);
  cudaFree(mdp_PreviousDisplacements);
  cudaFree(_GetOnDeviceRepresentation());
}
#endif
