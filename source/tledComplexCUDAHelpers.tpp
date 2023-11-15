// =========================================================================
// File:       tledComplexCUDAHelpers.tpp
// Purpose:    Utility
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

namespace tledComplexCUDAHelpers_kernels {
  template <const int t_blockSize, typename TComparable, class TStrictOrdering>
  __global__ void SortInit(TComparable *p_outData, int2 *p_blockBounds, TStrictOrdering ordering, const TComparable *inData, const int numIn) {
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    
    __shared__ unsigned short numItems[t_blockSize];
    __shared__ TComparable blockItems[t_blockSize];

    tledCudaAssert(t_blockSize == blockDim.x);
    numItems[threadIdx.x] = 0;
    if (tid < numIn) {      
      numItems[threadIdx.x] += 1;
      blockItems[threadIdx.x] = inData[tid];
    }
    __syncthreads();
    
    tledCUDAHelpers::MergeSort<t_blockSize, 1, TComparable, TStrictOrdering>(blockItems, numItems, ordering);    

    if (threadIdx.x == 0) {
      p_blockBounds[blockIdx.x].x = blockIdx.x*blockDim.x;
      p_blockBounds[blockIdx.x].y = p_blockBounds[blockIdx.x].x + numItems[0];
    }

    tledCudaAssert(numItems[0] <= t_blockSize);
    if (threadIdx.x < numItems[0]) {
      p_outData[tid] = blockItems[threadIdx.x];  
    }
  }

  template <const int t_blockSize, typename TComparable>
  __device__ void _CopyBlock(TComparable *p_outData, int2 *p_outBounds, const TComparable *inData, const int2 blockBounds) {
    p_outBounds[blockIdx.x] = blockBounds;
    tledCUDAHelpers::Copy<TComparable, int>(p_outData + blockBounds.x, inData + blockBounds.x, blockBounds.y - blockBounds.x);
  }

  template <const int t_blockSize, typename TComparable, class TStrictOrdering>
  __global__ void UniqueBlockWise(TComparable *p_data, int2 *p_bounds, TStrictOrdering ordering) {
    using namespace tledCUDAHelpers;

    const int2 inBounds = p_bounds[blockIdx.x];

    __shared__ int numItems;

    if (threadIdx.x == 0) numItems = inBounds.y - inBounds.x;
    Unique<t_blockSize, TComparable, _NotOrderingPred<TComparable, TStrictOrdering> >(p_data + inBounds.x, numItems, _NotOrderingPred<TComparable, TStrictOrdering>(ordering));
    if (threadIdx.x == 0) {
      p_bounds[blockIdx.x].y = numItems + inBounds.x;
    }
  }

  template <const int t_blockSize, typename TComparable, class TStrictOrdering>
  __global__ void UniqueSubBlockWise(TComparable *p_outData, int2 *p_outBounds, TStrictOrdering ordering, const TComparable *inData, const int2 *inBounds, const int numBlocks) {
    using namespace tledCUDAHelpers;

    const int2 inBlkBounds = inBounds[0];
    const int numIn = inBlkBounds.y - inBlkBounds.x;
    const int numSubBlk = numIn/numBlocks + (numIn%numBlocks > 0);

    int2 inSubBlkBounds;

    tledCudaAssert(blockDim.x == t_blockSize);
    inSubBlkBounds.x = inBlkBounds.x + blockIdx.x*numSubBlk;
    if (inSubBlkBounds.x >= inBlkBounds.y) {
      p_outBounds[blockIdx.x].y = p_outBounds[blockIdx.x].x = inBlkBounds.y;
    } else {
      __shared__ int offBlock;

      inSubBlkBounds.y = inSubBlkBounds.x + numSubBlk;
      if (inSubBlkBounds.y > inBlkBounds.y) inSubBlkBounds.y = inBlkBounds.y;
      __syncthreads();

      if (threadIdx.x == 0) offBlock = inSubBlkBounds.x;
      tledCudaAssert(inSubBlkBounds.y - inSubBlkBounds.x <= blockDim.x);
      {
	__shared__ bool isUnique[t_blockSize];
	int lOff;

	if (inSubBlkBounds.x + threadIdx.x < inSubBlkBounds.y) {
	  if (inSubBlkBounds.x + threadIdx.x == 0) isUnique[threadIdx.x] = true;
	  else {
	    isUnique[threadIdx.x] = ordering(inData[inSubBlkBounds.x+threadIdx.x-1], inData[inSubBlkBounds.x+threadIdx.x]);
	  }
	} else isUnique[threadIdx.x] = false;
	lOff = tledCUDAHelpers::ComputeOffsetFromFlag<t_blockSize>(isUnique);
	tledCudaAssert(lOff <= threadIdx.x);
	if (isUnique[threadIdx.x]) p_outData[offBlock+lOff] = inData[inSubBlkBounds.x+threadIdx.x];
	__syncthreads();
	if (threadIdx.x == t_blockSize - 1) offBlock += lOff + isUnique[threadIdx.x];
	__syncthreads();
      }

      if (threadIdx.x == 0) {
	p_outBounds[blockIdx.x].x = inSubBlkBounds.x;
	p_outBounds[blockIdx.x].y = offBlock;
      }
    }

    __syncthreads();
    tledCudaAssert(p_outBounds[blockIdx.x].x <= p_outBounds[blockIdx.x].y);
    tledCudaAssert(p_outBounds[blockIdx.x].y > 0);
  }

  template <typename TComparable>
  __global__ void MergeUnique(TComparable *p_outData, int *p_numItems, const TComparable *inData, const int2 *inBounds, const int numBlocks) {
    const int2 inBlkBounds = inBounds[blockIdx.x];
    const int numIn = inBlkBounds.y - inBlkBounds.x;
    
    int2 outBlkBounds;

    tledCudaAssert(inBlkBounds.y >= inBlkBounds.x);
    outBlkBounds.x = 0;
    for (int b = 0; b < blockIdx.x; b++) {
      tledCudaAssert(inBounds[b].y >= inBounds[b].x && inBounds[b].x >= 0);
      outBlkBounds.x += inBounds[b].y - inBounds[b].x;
    }
    outBlkBounds.y = outBlkBounds.x + numIn;

    tledCudaAssert(inBlkBounds.y - inBlkBounds.x <= blockDim.x);
    if (inBlkBounds.x + threadIdx.x < inBlkBounds.y) {
      p_outData[outBlkBounds.x+threadIdx.x] = inData[inBlkBounds.x+threadIdx.x];
    }

    tledCudaAssert(blockIdx.x == 0 || outBlkBounds.x > 0);
    if (threadIdx.x == 0 && numBlocks == blockIdx.x + 1) *p_numItems = outBlkBounds.y;
  }

  template <typename TComparable, class TStrictOrdering>
  __global__ void MergeBlocks(TComparable *p_outData, int2 *p_outBounds, TStrictOrdering ordering, const TComparable *inData, const int2 *inBounds, const int blockMultiplier, const int numBlocks) {
    using namespace tledCUDAHelpers;
    
    const int outBlockInd = blockIdx.x/blockMultiplier;
    const int subBlockInd = blockIdx.x%blockMultiplier;
    const int inBlockInd0 = 2*outBlockInd, inBlockInd1 = inBlockInd0 + 1;

    tledCudaAssert(p_outData != inData && p_outBounds != inBounds);
    if (inBlockInd1 < numBlocks) {
      const int2 src0Bounds = inBounds[inBlockInd0], src1Bounds = inBounds[inBlockInd1];
      const int numSrc0 = src0Bounds.y - src0Bounds.x, numSrc1 = src1Bounds.y - src1Bounds.x;
      const int subSrc0BlkSize = numSrc0/blockMultiplier + (numSrc0%blockMultiplier > 0), subSrc1BlkSize = numSrc1/blockMultiplier + (numSrc1%blockMultiplier > 0);

      int2 subSrc0Bounds, subSrc1Bounds, searchSrc0Bounds = src0Bounds, searchSrc1Bounds = src1Bounds;

      tledCudaAssert(numSrc0 > 0 && numSrc1 > 0);
      if (threadIdx.x == 0) {
	p_outBounds[outBlockInd].x = src0Bounds.x;
	p_outBounds[outBlockInd].y = numSrc0 + numSrc1 + src0Bounds.x;
      }
      tledCudaAssert(src0Bounds.y <= src1Bounds.x);

      subSrc0Bounds.x = src0Bounds.x + subBlockInd*subSrc0BlkSize;
      subSrc1Bounds.x = src1Bounds.x + subBlockInd*subSrc1BlkSize;

      subSrc0Bounds.y = src0Bounds.x + (subBlockInd + 1)*subSrc0BlkSize;
      if (subSrc0Bounds.y > src0Bounds.y) subSrc0Bounds.y = src0Bounds.y;
      subSrc1Bounds.y = src1Bounds.x + (subBlockInd + 1)*subSrc1BlkSize;
      if (subSrc1Bounds.y > src1Bounds.y) subSrc1Bounds.y = src1Bounds.y;

#ifdef __CUDA_DEBUG
      for (int j = subSrc0Bounds.x; j < subSrc0Bounds.y; j += blockDim.x) if (j + threadIdx.x < subSrc0Bounds.y) {
	  if (j + threadIdx.x > src0Bounds.x) {
	    tledCudaAssert(j + threadIdx.x - 1 >= src0Bounds.x && j + threadIdx.x < src0Bounds.y);
	    tledCudaAssert(ordering(inData[j+threadIdx.x-1], inData[j+threadIdx.x]) || !ordering(inData[j+threadIdx.x], inData[j+threadIdx.x-1]));
	  }
	}

      for (int j = subSrc0Bounds.x; j < subSrc0Bounds.y; j += blockDim.x) if (j + threadIdx.x < subSrc0Bounds.y) {
	  if (j + threadIdx.x > src1Bounds.x) {
	    tledCudaAssert(j + threadIdx.x - 1 >= src1Bounds.x && j + threadIdx.x < src1Bounds.y);
	    tledCudaAssert(ordering(inData[j+threadIdx.x-1], inData[j+threadIdx.x]) || !ordering(inData[j+threadIdx.x], inData[j+threadIdx.x-1]));
	  }
	}
#endif

      tledCudaAssert(subSrc0Bounds.y - subSrc0Bounds.x <= blockDim.x);
      if (subSrc0Bounds.x + threadIdx.x < subSrc0Bounds.y) {
	const TComparable item = inData[subSrc0Bounds.x+threadIdx.x];
	const int relInd = subSrc0Bounds.x + threadIdx.x - src0Bounds.x;
	const int offItem = GetNumberOfLower(inData + searchSrc1Bounds.x, item, searchSrc1Bounds.y - searchSrc1Bounds.x, ordering);

	tledCudaAssert(offItem == 0 || ordering(inData[src1Bounds.x+offItem-1], item) || !ordering(item, inData[src1Bounds.x+offItem-1]));
	tledCudaAssert(relInd >= 0 && offItem >= 0);
	tledCudaAssert(relInd + offItem + src0Bounds.x < src1Bounds.y);
	p_outData[src0Bounds.x+relInd+offItem] = item;
	searchSrc1Bounds.x += offItem;
      }

      tledCudaAssert(subSrc1Bounds.y - subSrc1Bounds.x <= blockDim.x);
      if (subSrc1Bounds.x + threadIdx.x < subSrc1Bounds.y) {
	const TComparable item = inData[subSrc1Bounds.x+threadIdx.x];
	const int relInd = subSrc1Bounds.x + threadIdx.x - src1Bounds.x;
	const int offItem = GetNumberOfLowerEqual(inData + searchSrc0Bounds.x, item, searchSrc0Bounds.y - searchSrc0Bounds.x, ordering);

	tledCudaAssert(offItem == 0 || ordering(inData[src0Bounds.x+offItem-1], item) || !ordering(item, inData[src0Bounds.x+offItem-1]));
	tledCudaAssert(relInd >= 0 && offItem >= 0);
	tledCudaAssert(relInd + offItem + src0Bounds.x < src1Bounds.y);
	p_outData[src0Bounds.x+relInd+offItem] = item;
	searchSrc1Bounds.x += offItem;
      }
    } /* if have input */
  }  

  template <const int t_blockSize, typename TComparable>
  __global__ void CopyRemainderBlock(TComparable *p_outData, int2 *p_outBounds, const TComparable *inData, const int2 *inBounds, const int blockInd) {
    const int2 blockBounds = inBounds[blockInd];

    p_outBounds[blockInd/2] = blockBounds;
    tledCUDAHelpers::Copy<TComparable, int>(p_outData + blockBounds.x, inData + blockBounds.x, blockBounds.y - blockBounds.x);
  }
  
  template <const int t_blockSize, typename TComparable, class TStrictOrdering>
  __global__ void SortSmallBlock(TComparable *p_data, int *p_numItems, TStrictOrdering ordering, const bool makeUnique) {
    const int tid = threadIdx.x + blockDim.x*blockIdx.x;

    __shared__ unsigned short numItems[t_blockSize];
    __shared__ TComparable items[t_blockSize];

    numItems[threadIdx.x] = 0;
    if (tid < *p_numItems) {
      items[threadIdx.x] = p_data[tid];
      numItems[threadIdx.x] += 1;
    }
    __syncthreads();

    if (makeUnique) {
      tledCUDAHelpers::MakeSortedUnique<t_blockSize, 1, TComparable, TStrictOrdering>(items, numItems, ordering, false);
    } else {
      tledCUDAHelpers::MergeSort<t_blockSize, 1, TComparable, TStrictOrdering>(items, numItems, ordering);
    }

    tledCudaAssert(numItems[0] <= t_blockSize);
    tledCUDAHelpers::Copy<TComparable, unsigned short>(p_data, items, numItems[0]);
    if (threadIdx.x == 0) *p_numItems = numItems[0];
  }  
}

namespace tledComplexCUDAHelpers {
  template <typename TComparable, class TStrictOrdering>
  __host__ void _MergeSortAndUnique(tledCUDADeviceMemoryBlock &r_data, int *dp_numItems, TStrictOrdering ordering, const bool makeUnique, const cudaStream_t stream) {    
    const int blockSize = 128;

    int numIn = -1;

    tledCUDAHelpers::CopyFromDevice(&numIn, dp_numItems, 1, stream);    
    tledLogDebugStream(tledHelper::Info() << "Sorting " << numIn << " items" << (makeUnique? " and performing \"unique\"." : "."));
    tledCUDADeviceMemoryBlock::SaveAllocationCounter();
    if (numIn <= blockSize) {
      const int numBlocks = 1;
      
      tledComplexCUDAHelpers_kernels::SortSmallBlock<blockSize, TComparable, TStrictOrdering> <<<numBlocks, blockSize, 0, stream>>> (r_data.GetBuffer<TComparable>(), dp_numItems, ordering, makeUnique);
    } else {
      const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(numIn, blockSize);

      tledCUDADeviceMemoryBlock &r_bounds0 = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(numBlocks);  
      tledCUDADeviceMemoryBlock &r_bounds1 = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(numBlocks/2 + numBlocks%2);  
      tledCUDADeviceMemoryBlock &r_mergeBuffer = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<TComparable>(numIn);
      tledCUDADeviceMemoryBlock *pp_dataBuffers[2], *pp_boundsBuffers[2];

      assert(numBlocks > 1);
#ifndef NDEBUG
      tledCheckCUDAErrors(cudaMemset(r_mergeBuffer.GetBuffer<TComparable>(), -1, sizeof(TComparable)*numIn));
#endif 
      tledComplexCUDAHelpers_kernels::SortInit<blockSize, TComparable, TStrictOrdering> <<<numBlocks, blockSize, 0, stream>>> (r_mergeBuffer.GetBuffer<TComparable>(), r_bounds0.GetBuffer<int2>(), ordering, r_data.GetBuffer<TComparable>(), numIn);
      tledDeviceSyncDebug;

      pp_dataBuffers[0] = &r_mergeBuffer;
      pp_dataBuffers[1] = &r_data;
      pp_boundsBuffers[0] = &r_bounds0;
      pp_boundsBuffers[1] = &r_bounds1;

      if (makeUnique) {
	tledComplexCUDAHelpers_kernels::UniqueBlockWise<blockSize, TComparable, TStrictOrdering> <<<numBlocks, blockSize, 0, stream>>> (pp_dataBuffers[0]->GetBuffer<TComparable>(), pp_boundsBuffers[0]->GetBuffer<int2>(), ordering);
      }

      for (int numRemBlocks = numBlocks, blockMult = 1; numRemBlocks > 1; blockMult *= 2) {
#ifndef NDEBUG
	tledCheckCUDAErrors(cudaMemset(pp_dataBuffers[1]->GetBuffer<TComparable>(), -1, sizeof(TComparable)*numIn));
#endif 
	tledComplexCUDAHelpers_kernels::MergeBlocks<TComparable, TStrictOrdering> <<<(numRemBlocks/2)*blockMult, blockSize, 0, stream>>> (pp_dataBuffers[1]->GetBuffer<TComparable>(), pp_boundsBuffers[1]->GetBuffer<int2>(), ordering, pp_dataBuffers[0]->GetBuffer<TComparable>(), pp_boundsBuffers[0]->GetBuffer<int2>(), blockMult, numRemBlocks);
	tledDeviceSyncDebug;
	if (numRemBlocks%2 == 1) {
	  assert(numRemBlocks > 2);
	  tledComplexCUDAHelpers_kernels::CopyRemainderBlock<blockSize, TComparable> <<<1, blockSize, 0, stream>>> (pp_dataBuffers[1]->GetBuffer<TComparable>(), pp_boundsBuffers[1]->GetBuffer<int2>(), pp_dataBuffers[0]->GetBuffer<TComparable>(), pp_boundsBuffers[0]->GetBuffer<int2>(), numRemBlocks - 1);
	  tledDeviceSyncDebug;
	}

	std::iter_swap(pp_boundsBuffers, pp_boundsBuffers + 1);
	std::iter_swap(pp_dataBuffers, pp_dataBuffers + 1);

	numRemBlocks = numRemBlocks/2 + numRemBlocks%2;	
      }

      if (makeUnique) {
	tledComplexCUDAHelpers_kernels::UniqueSubBlockWise<blockSize, TComparable, TStrictOrdering> <<<numBlocks, blockSize, 0, stream>>> (pp_dataBuffers[1]->GetBuffer<TComparable>(), pp_boundsBuffers[1]->GetBuffer<int2>(), ordering, pp_dataBuffers[0]->GetBuffer<TComparable>(), pp_boundsBuffers[0]->GetBuffer<int2>(), numBlocks);
	std::iter_swap(pp_boundsBuffers, pp_boundsBuffers + 1);
	std::iter_swap(pp_dataBuffers, pp_dataBuffers + 1);	  
	tledComplexCUDAHelpers_kernels::MergeUnique<TComparable> <<<numBlocks, blockSize, 0, stream>>> (pp_dataBuffers[1]->GetBuffer<TComparable>(), dp_numItems, pp_dataBuffers[0]->GetBuffer<TComparable>(), pp_boundsBuffers[0]->GetBuffer<int2>(), numBlocks);
	std::iter_swap(pp_dataBuffers, pp_dataBuffers + 1);
      }

      if (pp_dataBuffers[0] != &r_data) {	
	r_data.Swap(*pp_dataBuffers[0]);	
      }

      r_bounds1.ToggleActive();
      r_bounds0.ToggleActive();
      r_mergeBuffer.ToggleActive();
      assert(!(r_mergeBuffer.IsActive() || r_bounds0.IsActive() || r_bounds1.IsActive()) && r_data.IsActive());
    }
    tledCUDADeviceMemoryBlock::CheckAllocationCounter();
  }  

  template <typename TComparable, class TStrictOrdering>
  __host__ void MakeSortedUnique(tledCUDADeviceMemoryBlock &r_data, int *dp_numItems, TStrictOrdering ordering, const cudaStream_t stream) {
    _MergeSortAndUnique<TComparable, TStrictOrdering>(r_data, dp_numItems, ordering, true, stream);
  }

  template <typename TComparable, class TStrictOrdering>
  __host__ void MergeSort(tledCUDADeviceMemoryBlock &r_data, TStrictOrdering ordering, const int *dpc_numItems, const cudaStream_t stream) {
    _MergeSortAndUnique<TComparable, TStrictOrdering>(r_data, const_cast<int*>(dpc_numItems), ordering, false, stream);
  }
}
