// =========================================================================
// File:       tledCUDAHelpers.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    June 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledCUDAHelpers_H
#define tledCUDAHelpers_H

#ifdef _GPU_
#ifdef __CUDA_DEBUG
#include <cassert>
#include <stdio.h>
#endif

#include <vector>

#include <cuda_runtime.h>
#include <vector_types.h>
#ifndef _CUDA_5PLUS_SDK
#include <cutil.h>
#include <cutil_inline_runtime.h>
#else
#include <helper_cuda.h>
#include <helper_math.h>
#endif

#if defined _USE_CUDA_ATOMICS && (defined __CUDA_ARCH__ &&  __CUDA_ARCH__ < 130)
#warning Enabling atomic operations was requested but GPU architecture too old
#undef _USE_CUDA_ATOMICS
#endif

#if defined __CUDA_ARCH__ && __CUDA_ARCH__ < 200
#undef __CUDA_DEBUG
#undef __CUDA_VERBOSITY
#endif

#if __CUDA_VERBOSITY > 0 && !defined __CUDA_DEBUG
#error Verbose mode only availably in debug builds, set __CUDA_DEBUG=1
#endif

#if defined __CUDA_DEBUG && defined __CUDA_VERBOSITY
#define tledCudaPrintf(DBG_LEVEL, fmt, ...) if (__CUDA_VERBOSITY >= DBG_LEVEL) printf(fmt, __VA_ARGS__)
#else
#define tledCudaPrintf(DBG_LEVEL, fmt, ...) 
#endif

#if defined __CUDA_DEBUG && defined __CUDA_VERBOSITY
#define tledCudaBlockPrintf(DBG_LEVEL, fmt, ...) if (threadIdx.x == 0 && __CUDA_VERBOSITY >= DBG_LEVEL) printf(fmt, __VA_ARGS__)
#else
#define tledCudaBlockPrintf(DBG_LEVEL, fmt, ...) 
#endif

#ifdef __CUDA_DEBUG
#define tledCudaAssert(cond) if (!(cond)) { tledCudaPrintf(0, __FILE__":%d: assertion: \"%s\" failed.\n", __LINE__, #cond); } //
#else
#define tledCudaAssert(cond)
#endif

#ifndef _CUDA_5PLUS_SDK
#define tledCheckCUDAErrors(expr) CUDA_SAFE_CALL(expr)
#else
#define tledCheckCUDAErrors(expr) checkCudaErrors(expr)
#endif

#if !defined NDEBUG && !defined _CUDA_3MINUS_SDK
#define tledDeviceSyncDebug tledCheckCUDAErrors(cudaDeviceSynchronize())
#else 
#define tledDeviceSyncDebug
#endif

/**
 * \brief Namespace for CUDA utility functions
 * \ingroup helper
 */
namespace tledCUDAHelpers {
  /** Returns the number of blocks required for processing a given number of items with blocks of a given size. */
  inline __host__ __device__ int GetNumberOfBlocks(const int numberOfProcessItems, const int blockSize);

  /** Allocates on-device memory */
  template <typename TValue>
  __host__ void AllocateDeviceMemory(TValue* (&rdp_dst), const int numEls = 1);

  /** Allocates alligned host memory */
  template <typename TValue>
  __host__ void AllocateHostMemory(TValue* (&rhp_dst), const int numEls = 1);

  /** Copies data to the device */
  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const TValue *hp_data, const int numEls = 1);

  /** Copies data to the device (async) */
  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const TValue *hp_data, const int numEls, const cudaStream_t stream);

  /** Copies data between device locations */
  template <typename TValue>
  __host__ void CopyDeviceToDevice(TValue *dp_dst, const TValue *dp_data, const int numEls = 1);

  /** Copies an entire STL vector to the device */
  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const std::vector<TValue> &data);

  /** Copies data from the device */
  template <typename TValue>
  __host__ void CopyFromDevice(TValue *hp_dst, const TValue *dp_data, const int numEls = 1);

  /** Copies data from the device */
  template <typename TValue>
  __host__ void CopyFromDevice(TValue *hp_dst, const TValue *dp_data, const int numEls, const cudaStream_t stream);

  /** Zeros out an on-device memory block */
  template <typename TValue>
  __host__ void SetZero(TValue *dp_mem, const int numEls = 1);

  /** Zeros out an on-device memory block (asynchronous version) */
  template <typename TValue>
  __host__ void SetZero(TValue *dp_mem, const int numEls, const cudaStream_t stream);

  /** Converts a 3-component array to a CUDA float3 */
  template <typename TValue>
  __host__ __device__ float3 ConvertToFloat3(const TValue *v);

  /** Converts a 3-component array to a CUDA float4 */
  template <typename TValue>
  __host__ __device__ float4 ConvertToFloat4(const TValue *v);

  /** Converts a CUDA float3 or float4 (ignoring 4-th component) to a 3-component float array */
  template <class TFloatN>
  __host__ __device__ float* ConvertFromFloatN(float *p_dst, const TFloatN &f);

  /** 
   * \brief Consolidates the results computed by the threads of one block (max. threadDataBlockSize/thread) in one list starting at p_list[0] and whose total number of elements is stored in p_numItems[0]. 
   *
   * 
   * All threads of the block must enter this function.
   * \param threadDataBlockSize the maximum number of results one individual thread can produce, i.e. the size of the block in p_list that is assigned to one thread.
   */
  template <const int t_blockSize, class TResult, typename TSize>
  __device__ void CollectResults(TResult *p_list, TSize *p_numItems, const int threadDataBlockSize);

  /** Same as CollectResults but also sorts the results. */
  template <const int t_blockSize, const int t_threadDataBlockSize, typename TComparable, class TStrictOrdering>
  __device__ void MergeSort(TComparable *p_list, unsigned short *p_numItems, TStrictOrdering ordering);

  /** 
   * Retains only elements in a list satisfying a given predicate. 
   * \param p_list R/W contiguous list of sorted items.
   * \param r_numItems R/W number of items.
   */
  template <const int t_blockSize, typename TValue, class TPredicate, typename TSize>
  __device__ void KeepIf(TValue *p_list, TSize &r_numItems, TPredicate pred);

  /** 
   * Removes dupes from a <i>sorted</i> list retaining only the first one of a group of items satisfying a given equality predicate. 
   * \param p_list R/W contiguous list of sorted items.
   * \param r_numItems R/W number of items.
   */
  template <const int t_blockSize, typename TComparable, class TEqualPredicate, typename TSize>
  __device__ void Unique(TComparable *p_list, TSize &r_numItems, TEqualPredicate pred);

  /** Same as CollectResults but also sorts the output and removes dupes. If doThreadBlockUnique, input thread blocks are checked internally for dupes. */
  template <const int t_blockSize, const int t_threadDataBlockSize, typename TComparable, class TStrictOrdering>
  __device__ void MakeSortedUnique(TComparable *p_list, unsigned short *p_numItems, TStrictOrdering ordering, const bool doThreadBlockUnique = true);

  /** Returns the number of elements in a sorted list smaller than or equal to "item". Useful for computing offsets. */
  template <typename TComparable, class TStrictOrdering, typename TSize>
  __device__ int GetNumberOfLowerEqual(const TComparable *list, const TComparable &item, const TSize numItems, const TStrictOrdering ordering);

  /** Returns the number of elements in a sorted list strictly smaller than "item". Useful for computing offsets. */
  template <typename TComparable, class TStrictOrdering, typename TSize>
  __device__ int GetNumberOfLower(const TComparable *list, const TComparable &item, const TSize numItems, const TStrictOrdering ordering);

  /** 
   * Computes (in log-time) the number boolean flags set corresponding to threads w/ an index < than the entering one.
   * \param flags, bool array of size t_blockSize.
   */
  template <const int t_blockSize>
  __device__ unsigned short ComputeOffsetFromFlag(const bool *flags);

  /** Parallel copy. All threads must enter, synchronicity at exit not guaranteed. */
  template <typename TItem, typename TSize>
  __device__ void Copy(TItem *p_out, const TItem *in, const TSize numItems);

#ifdef _USE_CUDA_ATOMICS
  /** 
   * \brief Writes results from shared memory to a global list, uses atomics to avoid concurrency issues.
   *
   * 
   * All threads of the block must enter this function.
   */
  template <class TResult>
  __device__ void ParallelFlush(TResult *gp_list, int *p_numGItems, const TResult *localList, const int numItems);
#endif
}

/**
 * \name Former gpuUtils functions.
 * @{
 */
/** Copy V3 data to a V4 variable */
template <class T3, class T4> void gpuV3ToV4(const T3* A, const int n, T4* B);

/** Copy V4 data to a V3 variable (x,y,z components only) */
template <class T4, class T3> void gpuV4ToV3(const T4* A, const int n, T3* B);
/** @} */

#ifdef __CUDACC__
#include "tledCUDAHelpers.tpp"
#endif

#endif
#endif
