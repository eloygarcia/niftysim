// =========================================================================
// File:       tledBVHGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHGPU_CU
#define tledBVHGPU_CU

template <class TBaseBVH>
tledBVHImplGPU<TBaseBVH>::~tledBVHImplGPU() {
  tledCheckCUDAErrors(cudaFreeHost(mhp_GPUBVs));
  tledCheckCUDAErrors(cudaFree(mdp_GPUBVs));
}

template <class TBaseBVH>
void tledBVHImplGPU<TBaseBVH>::InitDeviceMemory() {
  tledCUDAHelpers::AllocateHostMemory(mhp_GPUBVs, this->GetNumberOfBVs());
  tledCUDAHelpers::AllocateDeviceMemory(mdp_GPUBVs, this->GetNumberOfBVs());

  for (int i = 0; i < this->GetNumberOfBVs(); i++) this->GetBV(i).InitGPU(mhp_GPUBVs[i]);
  
  this->CopyBVsToGPU();  
}

template <class TBaseBVH>
void tledBVHImplGPU<TBaseBVH>::CopyBVsToGPU() {
  tledCUDAHelpers::CopyToDevice(this->mdp_GPUBVs, this->mhp_GPUBVs, this->GetNumberOfBVs());
}

template <class TBaseBVH>
void tledBVHImplGPU<TBaseBVH>::Init(tledBVHCreator &r_bvhBuilder) {
  Superclass::Init(r_bvhBuilder);
  this->InitDeviceMemory();
}

template <class TBaseBVH>
void tledBVHImplGPU<TBaseBVH>::LoadFromXMLPostloadHook() {
  Superclass::LoadFromXMLPostloadHook();
  this->InitDeviceMemory();
}
#endif
