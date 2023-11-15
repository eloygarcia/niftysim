// =========================================================================
// File:       tledCUDAMemoryBlock.h
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
#ifndef tledCUDAMemoryBlock_H
#define tledCUDAMemoryBlock_H

#include "tledHelper.h"
#include "tledCUDAHelpers.h"

#ifdef _USE_THRUST_
#include <thrust/device_ptr.h>
#endif

#include <algorithm>
#include <stack>
#include <cassert>

/** 
 * \brief Base class for reusable temporary memory blocks of CUDA memory.
 */
class tledCUDAMemoryBlock {
  /**
   * \name Buffer Management and Access
   * @{
   */
private:
  void *mp_Memory;
  int m_Size;
  bool m_IsActive;

protected:
  virtual void* Allocate(const int sz) = 0;
  virtual void Release(void *p_memory) = 0;

  template <class TBuffer>
  static TBuffer& GetNextFreeBuffer(std::vector<TBuffer*> &rvp_buffers);

  template <class TBuffer, class TValue>
  static TBuffer& GetNextFreeBufferWithSize(std::vector<TBuffer*> &rvp_buffers, const int size);

  template <class TBuffer>
  static int GetNumberOfActiveBuffers(const std::vector<TBuffer*> &vpc_buffers);

  template <class TBuffer>
  static void ReleaseAll(std::vector<TBuffer*> &rvp_buffers);

  void Release(void);

  void SwapMemoryBlock(tledCUDAMemoryBlock &r_otherBlock);

public:
  /** Resizes the internal buffer so it can hold numElements elements of type T. Automatically sets the active flag. */
  template <typename T>
  T* Resize(const int numElements);

  /** Returns the current buffer size in bytes */
  int GetBufferSize(void) const { return m_Size; }

  /** Returns the number of elements of given type currently fitting inside the buffer */
  template <typename T>
  int GetMaxNumberOfItems(void) const { return m_Size/sizeof(T); }  

  template <typename T>
  T* GetBuffer(void) { return reinterpret_cast<T*>(mp_Memory); }
  template <typename T>
  const T* GetBuffer(void) const { return reinterpret_cast<T*>(mp_Memory); }

#ifdef _USE_THRUST_
  template <typename T>
  thrust::device_ptr<T> GetThrustBuffer(void) { return thrust::device_pointer_cast(this->GetBuffer<T>()); }
  template <typename T>
  thrust::device_ptr<const T> GetThrustBuffer(void) const { return thrust::device_pointer_cast(this->GetBuffer<T>()); }
#endif

  /** Indicates if the buffer is currently being used */
  bool IsActive(void) const { return m_IsActive; }
  void SetActive(void) { m_IsActive = true; }
  void ToggleActive(void) { m_IsActive = !IsActive(); }
  /** @} */

  /**
   * \name Debug Interface
   * @{
   */
protected:
  static void PushMarker(std::stack<int> &r_markerStack, const int currentAllocation);

  /** Throws an error if the top element of markerStack does not hold the same value as currentAllocation. Only done if NDEBUG not set. */
  static void CheckMarker(std::stack<int> &r_markerStack, const int currentAllocation);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledCUDAMemoryBlock(void) : mp_Memory(NULL), m_Size(0), m_IsActive(false) {}
  virtual ~tledCUDAMemoryBlock(void) {}
  /** @} */
};

class tledCUDAHostMemoryBlock;

/** CUDA device memory block */
class tledCUDADeviceMemoryBlock : public tledCUDAMemoryBlock {
  /**
   * \name Buffer Management
   * @{
   */
private:
  static std::vector<tledCUDADeviceMemoryBlock*> svp_MemoryBlocks;
  static std::vector<tledCUDADeviceMemoryBlock*> svp_SwapBlocks;

private:
  static void _MarkSwapBufferForRelease(tledCUDADeviceMemoryBlock &r_swapBlock);

protected:
  virtual void* Allocate(const int sz);
  virtual void Release(void *dp_mem);

public:
  static tledCUDADeviceMemoryBlock& GetNextFreeBuffer(void) { return tledCUDAMemoryBlock::GetNextFreeBuffer(svp_MemoryBlocks); }

  /** Returns a buffer of size similar to "size" */
  template <typename TValue>
  static tledCUDADeviceMemoryBlock& GetNextFreeBufferWithSize(const int size) { return tledCUDAMemoryBlock::GetNextFreeBufferWithSize<tledCUDADeviceMemoryBlock, TValue>(svp_MemoryBlocks, size); }

  static int GetNumberOfActiveBuffers(void);

  /** Deallocates all memory blocks */
  static void ReleaseAll(void);

  /** 
   * \brief Resizes an active buffer (only upward) preserving its content. Roughly equivalent to calloc. 
   */
  template <typename TValue>
  void Grow(const int newSize);

  /** 
   * \brief Swaps the back-end storage buffer and its parameters with those of another device memory block.
   *
   *
   * Both blocks must be active.
   */
  void Swap(tledCUDADeviceMemoryBlock &r_otherBlock) { tledCUDAMemoryBlock::SwapMemoryBlock(r_otherBlock); }
  /** @} */

  /**
   * \name Debug Interface
   * @{
   */
private:
  static std::stack<int> s_Allocation;

public:
  /** Saves the current number of active buffers for a subsequent check. */
  static void SaveAllocationCounter(const int offset = 0) { tledCUDAMemoryBlock::PushMarker(s_Allocation, tledCUDADeviceMemoryBlock::GetNumberOfActiveBuffers() + offset); }
  static void CheckAllocationCounter(void) { tledCUDAMemoryBlock::CheckMarker(s_Allocation, tledCUDADeviceMemoryBlock::GetNumberOfActiveBuffers()); }
  /** @} */

  /**
   * \name Copying
   * @{
   */
public:
  template <typename TValue>
  void CopyFromHost(const tledCUDAHostMemoryBlock &data, const int num);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledCUDADeviceMemoryBlock(void) { if (this->GetBufferSize() > 0) tledCUDAMemoryBlock::Release(); }
  /** @} */
};

/** \brief Aligned host memory */
class tledCUDAHostMemoryBlock : public tledCUDAMemoryBlock {
  /**
   * \name Buffer Management
   * @{
   */
private:
  static std::vector<tledCUDAHostMemoryBlock*> svp_MemoryBlocks;

protected:
  virtual void* Allocate(const int sz);
  virtual void Release(void *hp_mem);

public:
  static tledCUDAHostMemoryBlock& GetNextFreeBuffer(void) { return tledCUDAMemoryBlock::GetNextFreeBuffer(svp_MemoryBlocks); }
  
  template <typename TValue>
  static tledCUDAHostMemoryBlock& GetNextFreeBufferWithSize(const int size) { return tledCUDAMemoryBlock::GetNextFreeBufferWithSize<tledCUDAHostMemoryBlock, TValue>(svp_MemoryBlocks, size); }

  static int GetNumberOfActiveBuffers(void);

  /** Deallocates all memory blocks */
  static void ReleaseAll(void);
  /** @} */

  /**
   * \name Debug Interface
   * @{
   */
private:
  static std::stack<int> s_Allocation;

public:
  static void SaveAllocationCounter(const int offset = 0) { tledCUDAMemoryBlock::PushMarker(s_Allocation, tledCUDAHostMemoryBlock::GetNumberOfActiveBuffers() + offset); }
  static void CheckAllocationCounter(void) { tledCUDAMemoryBlock::CheckMarker(s_Allocation, tledCUDAHostMemoryBlock::GetNumberOfActiveBuffers()); }
  /** @} */

  /**
   * \name Copying
   * @{
   */
public:
  template <typename TValue>
  void CopyFromDevice(const tledCUDADeviceMemoryBlock &data, const int num);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledCUDAHostMemoryBlock(void) { if (this->GetBufferSize() > 0) tledCUDAMemoryBlock::Release(); }
  /** @} */
};

#include "tledCUDAMemoryBlock.tpp"
#endif
