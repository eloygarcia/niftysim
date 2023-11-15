// =========================================================================
// File:       tledCUDAHelpers.tpp
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

namespace tledCUDAHelpers {
  template <typename TValue>
  __host__ void AllocateDeviceMemory(TValue* (&rdp_dst), const int numEls) {
    tledDeviceSyncDebug;
    tledCheckCUDAErrors(cudaMalloc((void**)&rdp_dst, sizeof(TValue)*numEls));
  }

  template <typename TValue>
  __host__ void AllocateHostMemory(TValue* (&rhp_dst), const int numEls) {
    tledDeviceSyncDebug;
    tledCheckCUDAErrors(cudaMallocHost((void**)&rhp_dst, sizeof(TValue)*numEls));
  }

  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const TValue *hpc_data, const int numEls) {
    tledCheckCUDAErrors(cudaMemcpy(dp_dst, hpc_data, numEls*sizeof(TValue), cudaMemcpyHostToDevice));
  }

  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const TValue *hpc_data, const int numEls, const cudaStream_t stream) {
    tledCheckCUDAErrors(cudaMemcpyAsync(dp_dst, hpc_data, numEls*sizeof(TValue), cudaMemcpyHostToDevice, stream));
  }

  template <typename TValue>
  __host__ void CopyDeviceToDevice(TValue *dp_dst, const TValue *dpc_data, const int numEls) {
    tledCheckCUDAErrors(cudaMemcpy(dp_dst, dpc_data, numEls*sizeof(TValue), cudaMemcpyDeviceToDevice));
  }

  template <typename TValue>
  __host__ void CopyToDevice(TValue *dp_dst, const std::vector<TValue> &data) {
    CopyToDevice(dp_dst, &data.front(), data.size());
  }

  template <typename TValue>
  __host__ void CopyFromDevice(TValue *hp_dst, const TValue *dpc_data, const int numEls, const cudaStream_t stream) {
    tledCheckCUDAErrors(cudaMemcpyAsync(hp_dst, dpc_data, numEls*sizeof(TValue), cudaMemcpyDeviceToHost, stream));
  }

  template <typename TValue>
  __host__ void CopyFromDevice(TValue *hp_dst, const TValue *dpc_data, const int numEls) {
    tledCheckCUDAErrors(cudaMemcpy(hp_dst, dpc_data, numEls*sizeof(TValue), cudaMemcpyDeviceToHost));
  }

  template <typename TValue>
  __host__ void SetZero(TValue *dp_mem, const int numEls) {
    tledCheckCUDAErrors(cudaMemset(dp_mem, 0, sizeof(TValue)*numEls));
  }

  template <typename TValue>
  __host__ void SetZero(TValue *dp_mem, const int numEls, const cudaStream_t stream) {
    tledCheckCUDAErrors(cudaMemsetAsync(dp_mem, 0, sizeof(TValue)*numEls, stream));
  }

  template <typename TValue>
  __host__ __device__ float3 ConvertToFloat3(const TValue *v) {
    return make_float3(float(v[0]), float(v[1]), float(v[2]));
  }

  template <typename TValue>
  __host__ __device__ float4 ConvertToFloat4(const TValue *v) {
    return make_float4(float(v[0]), float(v[1]), float(v[2]), 0.f);
  }

  template <class TFloatN>
  __host__ __device__ float* ConvertFromFloatN(float *p_dst, const TFloatN &f) {
    p_dst[0] = f.x;
    p_dst[1] = f.y;
    p_dst[2] = f.z;

    return p_dst;
  }

  template <const int t_blockSize, typename TValue, class TListPredicate, typename TSize>
  __device__ void _KeepIf(TValue *p_list, TSize &r_numItems, TListPredicate pred) {
    __shared__ TSize offBlock;

    tledCudaAssert(t_blockSize == blockDim.x);
    if (threadIdx.x == 0) offBlock = 0;
    __syncthreads();   
    for (int k = 0; k < r_numItems; k += t_blockSize) {
      __shared__ bool doKeep[t_blockSize];
      __shared__ TValue copyBuffer[t_blockSize];
      TValue &item = copyBuffer[threadIdx.x];
      TSize lOff;

      if (k + threadIdx.x < r_numItems) {
	item = p_list[k+threadIdx.x];
	doKeep[threadIdx.x] = pred(k + threadIdx.x, p_list);
      } else doKeep[threadIdx.x] = false;

      lOff = ComputeOffsetFromFlag<t_blockSize>(doKeep);
      if (doKeep[threadIdx.x]) p_list[offBlock+lOff] = item;
      __syncthreads();
      if (threadIdx.x == t_blockSize - 1) offBlock += lOff + doKeep[threadIdx.x];
      __syncthreads();
    }

    __syncthreads();
    if (threadIdx.x == 0) r_numItems = offBlock;
    __syncthreads();
  }  

  template <typename TSize, class TResult>
  class _CollectKeepIfPred {
  private:
    const int mc_MaxNumThreadBlockItems;
    const TSize *mpc_NumThreadBlockItems;

  public:
    __device__ bool operator()(const int itemInd, const TResult *items) const {
      const int bInd = (itemInd)/mc_MaxNumThreadBlockItems;
      const int bItemInd = (itemInd)%mc_MaxNumThreadBlockItems;

      return mpc_NumThreadBlockItems[bInd] > bItemInd;
    }

  public:
    __device__ _CollectKeepIfPred(const TSize *numThreadBlkItems, const int maxNumThreadBlkItms) : mc_MaxNumThreadBlockItems(maxNumThreadBlkItms), mpc_NumThreadBlockItems(numThreadBlkItems) {}
  };

  template <const int t_blockSize, class TResult, typename TSize>
  __device__ void CollectResults(TResult *p_list, TSize *p_numItems, const int threadDataBlockSize) {
    typedef _CollectKeepIfPred<TSize, TResult> __Pred;

    __shared__ TSize numItems;

    if (threadIdx.x == 0) numItems = threadDataBlockSize*t_blockSize;
    __syncthreads();
    _KeepIf<t_blockSize, TResult, __Pred>(p_list, numItems, __Pred(p_numItems, threadDataBlockSize));
    if (threadIdx.x == 0) p_numItems[0] = numItems;
    __syncthreads();
  }

  template <const int t_blockSize>
  __device__ unsigned short ComputeOffsetFromFlag(const bool *flags) {
    __shared__ unsigned short numItems[t_blockSize];

    tledCudaAssert(t_blockSize == blockDim.x);
    numItems[threadIdx.x] = flags[threadIdx.x];
    for (int m = 2; m <= t_blockSize; m *= 2) {  
      unsigned short oldV;

      __syncthreads();
      if (threadIdx.x >= m/2) {
	oldV = numItems[threadIdx.x-m/2];
      } else oldV = 0;

      __syncthreads();
      numItems[threadIdx.x] += oldV;
    }

    return numItems[threadIdx.x] - static_cast<unsigned short>(flags[threadIdx.x]);
  }

  template <typename TComparable, class TStrictOrdering, typename TSize>
  __device__ int GetNumberOfLower(const TComparable *list, const TComparable &item, const TSize numItems, const TStrictOrdering ordering) {
    if (numItems == 0) return 0;
    else if (ordering(list[numItems-1], item)) return numItems;
    else {
      TSize highInd = numItems - 1, lowInd = 0, dstInd = 0;

      while (highInd - lowInd > 1) {	    
	dstInd = (highInd + lowInd)/2;
	if (!ordering(list[dstInd], item)) {
	  highInd = dstInd;
	} else {
	  lowInd = dstInd;
	}
      }	    
    
      if (ordering(list[lowInd], item)) return highInd;
      else return lowInd;
    }
  }

  template <typename TComparable, class TStrictOrdering, typename TSize>
  __device__ int GetNumberOfLowerEqual(const TComparable *list, const TComparable &item, const TSize numItems, const TStrictOrdering ordering) {
    if (numItems == 0) return 0;
    else if (!ordering(item, list[numItems-1])) return numItems;
    else if (ordering(item, list[0])) return 0;
    else {
      TSize highInd = numItems - 1, lowInd = 0, dstInd = 0;

      while (highInd - lowInd > 1) {	    
	dstInd = (highInd + lowInd)/2;
	if (ordering(item, list[dstInd])) {
	  highInd = dstInd;
	} else {
	  lowInd = dstInd;
	}
      }

      return lowInd + 1;    
    }
  }

  template <const int t_blockSize, const int t_threadDataBlockSize, typename TComparable, class TStrictOrdering>
  __device__ void MergeSort(TComparable *p_list, unsigned short *p_numItems, TStrictOrdering ordering) {
    __shared__ TComparable tmpBuffer[t_blockSize*t_threadDataBlockSize];

    tledCudaAssert(p_numItems[threadIdx.x] <= t_threadDataBlockSize);
    tledCudaAssert(t_blockSize == blockDim.x);

    for (int i = 0; i < p_numItems[threadIdx.x]; i++) {
      const TComparable item = p_list[t_threadDataBlockSize*threadIdx.x+i];

      for (int j = 0; j < i; j++) {
	if (ordering(item, p_list[t_threadDataBlockSize*threadIdx.x+j])) {
	  for (int k = i; k > j; k--) {
	    p_list[t_threadDataBlockSize*threadIdx.x+k] = p_list[t_threadDataBlockSize*threadIdx.x+k-1];
	  }
	  p_list[t_threadDataBlockSize*threadIdx.x+j] = item;		
	  break;
	}
      }
    }
    
    for (int m = 2; m <= t_blockSize; m *= 2) {
      const int bInd = (threadIdx.x/m)*m;
      const int oInd = (threadIdx.x/m)*m + m/2;
      const int cpyInd = threadIdx.x - bInd;

      __syncthreads();
      {
	const int numIn = p_numItems[oInd];
	const int numOut = p_numItems[bInd];

	for (int k = 0; k < numOut; k += m) if (cpyInd + k < numOut) {
	    const TComparable item = p_list[bInd*t_threadDataBlockSize+cpyInd+k];
	    const int dstInd = GetNumberOfLower(p_list + oInd*t_threadDataBlockSize, item, numIn, ordering);

	    tledCudaAssert(dstInd >= 0 && bInd*t_threadDataBlockSize + dstInd + k + cpyInd < bInd*t_threadDataBlockSize + m*t_threadDataBlockSize);
	    tmpBuffer[bInd*t_threadDataBlockSize+dstInd+k+cpyInd] = item;
	  }

	for (int k = 0; k < numIn; k += m) if (cpyInd + k < numIn) {
	    const TComparable item = p_list[oInd*t_threadDataBlockSize+cpyInd+k];
	    const int dstInd = GetNumberOfLowerEqual(p_list + bInd*t_threadDataBlockSize, item, numOut, ordering);

	    tledCudaAssert(dstInd >= 0 && bInd*t_threadDataBlockSize + dstInd + k + cpyInd < bInd*t_threadDataBlockSize + m*t_threadDataBlockSize);
	    tmpBuffer[bInd*t_threadDataBlockSize+dstInd+k+cpyInd] = item;
	  }
      
	__syncthreads();
	if (threadIdx.x == bInd) p_numItems[bInd] = numIn + numOut;
	for (int k = 0; k < t_threadDataBlockSize*t_blockSize; k += t_blockSize) {
	  p_list[k+threadIdx.x] = tmpBuffer[k+threadIdx.x];
	}
      }
    }

    __syncthreads();
  }

  template <class TComparable, class TEqualPred>
  class _UniqueKeepIfPred {
  private:
    TEqualPred m_Equal;

  public:
    __device__ bool operator()(const int itemInd, const TComparable *items) {
      if (itemInd == 0) return true;
      else return !m_Equal(items[itemInd-1], items[itemInd]);
    }

  public:
    __device__ _UniqueKeepIfPred(TEqualPred pred) : m_Equal(pred) {}
  };

  template <const int t_blockSize, typename TComparable, class TEqualPredicate, typename TSize>
  __device__ void Unique(TComparable *p_list, TSize &r_numItems, TEqualPredicate equal) {
    typedef _UniqueKeepIfPred<TComparable, TEqualPredicate> __Pred;

    _KeepIf<t_blockSize, TComparable, __Pred>(p_list, r_numItems, __Pred(equal));
  }

  template <typename TValue, class TUnaryPred>
  struct _UnaryPredListPredWrapper {
  private:
    TUnaryPred &mr_pred;
    
  public:
    __device__ bool operator()(const int itemInd, const TValue *list) const {
      return mr_pred(list[itemInd]);
    }

  public:
    __device__ _UnaryPredListPredWrapper(TUnaryPred &r_pred) : mr_pred(r_pred) {}
  };

  template <typename TValue, class TOrderingPred>
  class _NotOrderingPred {
  private:
    TOrderingPred &mr_Ordering;

  public:
    __device__ bool operator()(const TValue &a, const TValue &b) {
      return !mr_Ordering(a, b);
    }

  public:
    __device__ _NotOrderingPred(TOrderingPred &r_pred) : mr_Ordering(r_pred) {}
  };

  template <const int t_blockSize, typename TValue, class TPredicate, typename TSize>
  __device__ void KeepIf(TValue *p_list, TSize &r_numItems, TPredicate pred) {
    return _KeepIf<t_blockSize>(p_list, r_numItems, _UnaryPredListPredWrapper<TValue, TPredicate>(pred));
  }

  template <const int t_blockSize, const int t_threadDataBlockSize, typename TComparable, class TStrictOrdering>
  __device__ void MakeSortedUnique(TComparable *p_list, unsigned short *p_numItems, TStrictOrdering ordering, const bool doThreadBlockUnique) {
    MergeSort<t_blockSize, t_threadDataBlockSize, TComparable, TStrictOrdering>(p_list, p_numItems, ordering);
    Unique<t_blockSize, TComparable, _NotOrderingPred<TComparable, TStrictOrdering> >(p_list, p_numItems[0], _NotOrderingPred<TComparable, TStrictOrdering>(ordering));
  }
  
  template <typename TItem, typename TSize>
  __device__ void Copy(TItem *p_out, const TItem *in, const TSize numItems) {
    for (int j = 0; j < numItems; j += blockDim.x) {
      if (j + threadIdx.x < numItems) p_out[threadIdx.x+j] = in[j+threadIdx.x];
    }
  }

#ifdef _USE_CUDA_ATOMICS
  template <class TResult>
  __device__ void ParallelFlush(TResult *gp_list, int *p_numGItems, const TResult *localList, const int numItems) {
    if (numItems > 0) {
      __shared__ int baseInd;

      if (threadIdx.x == 0) {
	baseInd = atomicAdd(p_numGItems, numItems);
      }
      __syncthreads();

      Copy<TResult, int>(gp_list + baseInd, localList, numItems);
    }
  }
#endif

  inline __host__ __device__ int GetNumberOfBlocks(const int numberOfProcessItems, const int blockSize) {
    return numberOfProcessItems/blockSize + (numberOfProcessItems%blockSize > 0);
  }
}

template <class T3, class T4> void gpuV3ToV4(const T3* A, const int n, T4* B) {
   for (int i = 0; i < n; i++) {
      B[i].x = A[i].x;
      B[i].y = A[i].y;
      B[i].z = A[i].z;
      B[i].w = 0;
   }
}

template <class T4, class T3> void gpuV4ToV3(const T4* A, const int n, T3* B) {
   for (int i = 0; i < n; i++) {
      B[i].x = A[i].x;
      B[i].y = A[i].y;
      B[i].z = A[i].z;
   }
}

