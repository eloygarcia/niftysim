// =========================================================================
// File:       tledCUDAMemoryBlock.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <typename T>
T* tledCUDAMemoryBlock::Resize(const int numElements) {
  if (numElements > this->template GetMaxNumberOfItems<T>()) {
    this->Release();
    mp_Memory = this->Allocate(sizeof(T)*numElements);
    this->m_Size = sizeof(T)*numElements;
  }

  this->SetActive();

  return this->template GetBuffer<T>();
}

template <class TBuffer>
TBuffer& tledCUDAMemoryBlock::GetNextFreeBuffer(std::vector<TBuffer*> &rvp_buffers) {
  TBuffer *p_buffer = NULL;

  for (typename std::vector<TBuffer*>::iterator ip_buffer = rvp_buffers.begin(); ip_buffer < rvp_buffers.end(); ip_buffer++) {
    if (!(*ip_buffer)->IsActive()) {
      p_buffer = *ip_buffer;
      break;
    }
  }
  
  if (p_buffer == NULL) {
    rvp_buffers.push_back(new TBuffer);
    p_buffer = rvp_buffers.back();
  }
  p_buffer->SetActive();

  return *p_buffer;
}

template <class TBuffer, typename TValue>
TBuffer& tledCUDAMemoryBlock::GetNextFreeBufferWithSize(std::vector<TBuffer*> &rvp_buffers, const int size) {
  TBuffer *p_buffer = NULL;

  for (typename std::vector<TBuffer*>::iterator ip_buffer = rvp_buffers.begin(); ip_buffer < rvp_buffers.end(); ip_buffer++) {
    if (!(*ip_buffer)->IsActive() && (p_buffer == NULL || (p_buffer->template GetMaxNumberOfItems<TValue>() < size  && (*ip_buffer)->template GetMaxNumberOfItems<TValue>() > size) || (*ip_buffer)->template GetMaxNumberOfItems<TValue>() - size < p_buffer->template GetMaxNumberOfItems<TValue>() - size)) {
      p_buffer = *ip_buffer;
    }
  }
  
  if (p_buffer == NULL) {
    rvp_buffers.push_back(new TBuffer);
    p_buffer = rvp_buffers.back();
  }

  p_buffer->template Resize<TValue>(size);

  return *p_buffer;
}

template <typename TValue>
void tledCUDADeviceMemoryBlock::Grow(const int size) {
  assert(this->IsActive());
  if (size > this->template GetMaxNumberOfItems<TValue>()) {
    tledCUDADeviceMemoryBlock &r_swapBlock = tledCUDADeviceMemoryBlock::template GetNextFreeBufferWithSize<TValue>(size);

    tledCUDAHelpers::CopyDeviceToDevice(r_swapBlock.template GetBuffer<TValue>(), this->template GetBuffer<TValue>(), this->template GetMaxNumberOfItems<TValue>());
    this->SwapMemoryBlock(r_swapBlock);
    r_swapBlock.ToggleActive();
  }
}

template <typename TValue>
void tledCUDADeviceMemoryBlock::CopyFromHost(const tledCUDAHostMemoryBlock &data, const int num) {
  assert(this->GetBufferSize() >= (int)(num*sizeof(TValue)));
  assert(data.GetBufferSize() >= (int)(num*sizeof(TValue)));
  tledCUDAHelpers::CopyToDevice(this->template GetBuffer<TValue>(), data.template GetBuffer<TValue>(), num);
}

template <typename TValue>
void tledCUDAHostMemoryBlock::CopyFromDevice(const tledCUDADeviceMemoryBlock &data, const int num) {
  assert(this->GetBufferSize() >= (int)(num*sizeof(TValue)));
  assert(data.GetBufferSize() >= (int)(num*sizeof(TValue)));
  tledCUDAHelpers::CopyFromDevice(this->template GetBuffer<TValue>(), data.template GetBuffer<TValue>(), num);
}

