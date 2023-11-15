// =========================================================================
// File:       tledDynamicContactSurfaceGPU.tpp
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

#include "tledDynamicContactSurfaceGPU_kernels.tpp"

template <class TBaseSurface>
void tledDynamicContactSurfaceGPU<TBaseSurface>::AllocateDeviceSurface() {
  tledCUDAHelpers::AllocateDeviceMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetDeviceGPUSurface()));
  tledCUDAHelpers::AllocateHostMemory<GPUSurface>(reinterpret_cast<GPUSurface*&>(this->GetHostGPUSurfacePointer()));
}

template <class TBaseSurface>
void tledDynamicContactSurfaceGPU<TBaseSurface>::InitNodes() {
  Superclass::InitNodes();  

  tledCUDAHelpers::AllocateDeviceMemory(_GetHostSurface().NodeCoordinates0, this->GetNumberOfNodes());
  this->CopyNodeVectorToDevice(this->GetAllOnDeviceNodeCoordinates0(), this->GetAllNodeCoordinates());
  assert(this->GetCoordinateHistoryBufferMultiplier() > 0);
  tledCUDAHelpers::AllocateDeviceMemory(mdp_OldCoordinateBuffer, this->GetCoordinateHistoryBufferMultiplier()*this->GetNumberOfNodes());
  _GetHostSurface().OldNodeCoordinates = mdp_OldCoordinateBuffer;
  this->CopyNodeVectorToDevice(mdp_OldCoordinateBuffer, this->GetAllNodeCoordinates());
}

template <class TBaseSurface>
void tledDynamicContactSurfaceGPU<TBaseSurface>::CopySurfaceToDevice() {  
  tledCUDAHelpers::CopyToDevice(reinterpret_cast<GPUSurface*>(this->GetDeviceGPUSurface()), &_GetHostSurface());
}

template <class TBaseSurface>
void tledDynamicContactSurfaceGPU<TBaseSurface>::LoadFromXMLPostloadHook() {  
}

template <class TBaseSurface>
tledDynamicContactSurfaceGPU<TBaseSurface>::tledDynamicContactSurfaceGPU() {
  mdp_OldCoordinateBuffer = NULL;
}

template <class TBaseSurface>
tledDynamicContactSurfaceGPU<TBaseSurface>::~tledDynamicContactSurfaceGPU() {
  if (mdp_OldCoordinateBuffer != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_OldCoordinateBuffer)); 
    tledCheckCUDAErrors(cudaFree(_GetHostSurface().NodeCoordinates0));
  }
}
