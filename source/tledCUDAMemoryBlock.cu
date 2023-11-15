// =========================================================================
// File:       tledCUDAMemoryBlock.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledCUDAMemoryBlock_CU
#define tledCUDAMemoryBlock_CU

#include "tledCUDAMemoryBlock.h"

#include <cassert>
#include <algorithm>
#include <cstdlib>

template <class TMemory>
int tledCUDAMemoryBlock::GetNumberOfActiveBuffers(const std::vector<TMemory*> &vp_blocks) {
  int numActive = 0;

  for (typename std::vector<TMemory*>::const_iterator icp_m = vp_blocks.begin(); icp_m < vp_blocks.end(); icp_m++) {
    numActive += (*icp_m)->IsActive();
  }

  return numActive;
}

void tledCUDAMemoryBlock::Release() { 
  this->Release(mp_Memory); 
  this->mp_Memory = NULL;
  this->m_Size = 0;
}

void tledCUDAMemoryBlock::SwapMemoryBlock(tledCUDAMemoryBlock &r_otherBlock) { 
  assert(this->IsActive() && r_otherBlock.IsActive());
  std::iter_swap(&mp_Memory, &r_otherBlock.mp_Memory); 
  std::iter_swap(&m_Size, &r_otherBlock.m_Size);
}

template <class TBuffer>
void tledCUDAMemoryBlock::ReleaseAll(std::vector<TBuffer*> &rvp_buffers) {
  for (typename std::vector<TBuffer*>::iterator ip_b = rvp_buffers.begin(); ip_b < rvp_buffers.end(); ip_b++) {
    if ((*ip_b)->IsActive()) (*ip_b)->ToggleActive();
    delete *ip_b;
  }
  rvp_buffers.resize(0);
}

void tledCUDAMemoryBlock::PushMarker(std::stack<int> &r_markerStack, const int currentAllocation) {
#ifndef NDEBUG
  r_markerStack.push(currentAllocation);
#endif
}

void tledCUDAMemoryBlock::CheckMarker(std::stack<int> &r_markerStack, const int currentAllocation) {
#ifndef NDEBUG
  assert(!r_markerStack.empty());
  if (r_markerStack.top() != currentAllocation) {
    tledLogErrorStream(tledHelper::FatalError() << "Have a leak of " << currentAllocation - r_markerStack.top() << " CUDA memory buffers.");
  }
  r_markerStack.pop();
#endif
}

/************************************************************************************************************/
std::vector<tledCUDADeviceMemoryBlock*> tledCUDADeviceMemoryBlock::svp_MemoryBlocks(0);
std::stack<int> tledCUDADeviceMemoryBlock::s_Allocation;

void* tledCUDADeviceMemoryBlock::Allocate(const int sz) {
  void *dp_mem;

  assert(sz > 0);
  tledCheckCUDAErrors(cudaMalloc(&dp_mem, sz));

  return dp_mem;
}

int tledCUDADeviceMemoryBlock::GetNumberOfActiveBuffers() { 
  return tledCUDAMemoryBlock::GetNumberOfActiveBuffers(svp_MemoryBlocks); 
}

void tledCUDADeviceMemoryBlock::Release(void *dp_mem) {
  tledCheckCUDAErrors(cudaFree(dp_mem));
}

void tledCUDADeviceMemoryBlock::ReleaseAll() {
  tledCUDAMemoryBlock::ReleaseAll(svp_MemoryBlocks);
}

/************************************************************************************************************/
std::vector<tledCUDAHostMemoryBlock*> tledCUDAHostMemoryBlock::svp_MemoryBlocks(0);
std::stack<int> tledCUDAHostMemoryBlock::s_Allocation;

void* tledCUDAHostMemoryBlock::Allocate(const int sz) {
  void *hp_mem;

  assert(sz > 0);
  tledCheckCUDAErrors(cudaMallocHost(&hp_mem, sz));

  return hp_mem;
}

void tledCUDAHostMemoryBlock::Release(void *hp_mem) {
  tledCheckCUDAErrors(cudaFreeHost(hp_mem));
}

void tledCUDAHostMemoryBlock::ReleaseAll() {
  tledCUDAMemoryBlock::ReleaseAll(svp_MemoryBlocks);
}

int tledCUDAHostMemoryBlock::GetNumberOfActiveBuffers() { 
  return tledCUDAMemoryBlock::GetNumberOfActiveBuffers(svp_MemoryBlocks); 
}

#endif //tledCUDAMemoryBlock_CU
