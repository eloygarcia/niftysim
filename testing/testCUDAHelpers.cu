// =========================================================================
// File:       testCUDAHelpers.cu
// Purpose:    CUDA helper module unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledCUDAHelpers.h"
#include "tledHelper.h"

#include <cstdlib>
#include <iterator>

#include "tledCUDAHelpers.cu"

using namespace std;
using namespace tledUnitTest;

struct _IntOrdering {
  __device__ bool operator()(const int a, const int b) const {
    return a < b;
  }
};

struct _IntEquality {
  __device__ bool operator()(const int a, const int b) const {
    return a == b;
  }
};

#ifdef  __USE_CUDA_ATOMICS
template <const int t_threadDataBlockSize, const int t_blockSize>
__global__ void _TestCollectKernel(int *p_dataOut, int *p_numOut, const int *inData, const int numIn) {
  const int baseInd = (blockDim.x*blockIdx.x + threadIdx.x)*t_threadDataBlockSize;

  __shared__ int items[t_threadDataBlockSize*t_blockSize];
  __shared__ unsigned short numItems[t_blockSize];

  numItems[threadIdx.x] = 0;
  for (int c = 0; c < t_threadDataBlockSize; c++) {
    if (baseInd + c < numIn && inData[baseInd+c]%2 == 0) {
      items[threadIdx.x*t_threadDataBlockSize+numItems[threadIdx.x]] = inData[baseInd+c];
      numItems[threadIdx.x] += 1;
    } 
  }
  
  tledCUDAHelpers::CollectResults<t_blockSize, int>(items, numItems, t_threadDataBlockSize);
  tledCUDAHelpers::ParallelFlush<int>(p_dataOut, p_numOut, items, numItems[0]);
}

template <const int t_threadDataBlockSize, const int t_blockSize>
static void _TestCollect(void) {
  const int numItems = 2617;  
  const int numBlks = tledCUDAHelpers::GetNumberOfBlocks(numItems, t_blockSize);

  int *hp_data, *dp_outData, *dp_data, *dp_numResults, h_numResults;
  std::vector<int> dataRef;
  int numOdd = 0, numEven = 0;
  
  tledCUDAHelpers::AllocateHostMemory(hp_data, numItems);
  dataRef.reserve(numItems/2);
  for (int c = 0; c < numItems; c++) {
    hp_data[c] = rand()%123;

    if (hp_data[c]%2 == 0) {
      dataRef.push_back(hp_data[c]);
      numEven += 1;
    } else numOdd += 1;
  }

  tledCUDAHelpers::AllocateDeviceMemory(dp_numResults);
  tledCUDAHelpers::SetZero(dp_numResults);
  tledCUDAHelpers::AllocateDeviceMemory(dp_data, numItems);
  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_outData, numItems);

  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, numItems);

  _TestCollectKernel<t_threadDataBlockSize, t_blockSize> <<<numBlks, t_blockSize>>> (dp_outData, dp_numResults, dp_data, numItems);

  tledCUDAHelpers::CopyFromDevice(&h_numResults, dp_numResults);
  tledCUDAHelpers::CopyFromDevice(hp_data, dp_outData, h_numResults);
  tledDeviceSyncDebug;

  tledUnitTestAssert(h_numResults == numEven);
  std::sort(dataRef.begin(), dataRef.end());
  std::sort(hp_data, hp_data + h_numResults);
  tledUnitTestAssert(std::equal(dataRef.begin(), dataRef.end(), hp_data));

  tledCheckCUDAErrors(cudaFree(dp_data));
  tledCheckCUDAErrors(cudaFree(dp_outData));
  tledCheckCUDAErrors(cudaFree(dp_numResults));
  tledCheckCUDAErrors(cudaFreeHost(hp_data));
}
#else
template <const int t_threadDataBlockSize, const int t_blockSize>
static void _TestCollect(void) {
}
#endif

template <const int t_threadDataBlockSize, const int t_blockSize>
__global__ void _TestSortedUniqueKernel(int *p_outData, int *p_numOut, const int *inData, const int numIn) {
  __shared__ int s_items[t_threadDataBlockSize*t_blockSize];
  __shared__ unsigned short s_numItems[t_blockSize];

  tledCudaAssert(numIn == blockDim.x*t_threadDataBlockSize);
  s_numItems[threadIdx.x] = t_threadDataBlockSize;
  for (int i = 0; i < t_threadDataBlockSize; i++) s_items[t_threadDataBlockSize*threadIdx.x+i] = inData[t_threadDataBlockSize*threadIdx.x+i];

  tledCUDAHelpers::MakeSortedUnique<t_blockSize, t_threadDataBlockSize, int, _IntOrdering>(s_items, s_numItems, _IntOrdering());
  if (threadIdx.x == 0) {
    *p_numOut = s_numItems[0];
    for (int i = 0; i < s_numItems[0]; i++) {
      p_outData[i] = s_items[i];
    }
  }
  __syncthreads();
}

template <const int t_threadDataBlockSize, const int t_blockSize>
static void _TestSortedUnique(void) {
  const int numItems = t_threadDataBlockSize*t_blockSize;
  
  int *dp_data, *hp_data, *dp_outData, *dp_numItems, h_numItems;
  std::vector<int> refData;

  tledCUDAHelpers::AllocateHostMemory(hp_data, numItems);
  for (int i = 0; i < numItems; i++) {
    hp_data[i] = rand()%(numItems/2 - 1);
  }

  refData.insert(refData.end(), hp_data, hp_data + numItems);
  refData = tledHelper::MakeSortedUnique(refData);

  tledCUDAHelpers::AllocateDeviceMemory(dp_data, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_outData, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_numItems);
  tledCUDAHelpers::SetZero(dp_numItems);
  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, numItems);

  _TestSortedUniqueKernel<t_threadDataBlockSize, t_blockSize> <<<1, t_blockSize>>> (dp_outData, dp_numItems, dp_data, numItems);
  tledDeviceSyncDebug;

  tledCUDAHelpers::CopyFromDevice(&h_numItems, dp_numItems);
  std::cout << "Copying " << h_numItems << " (of " << refData.size() << " expected) from device\n";
  tledCUDAHelpers::CopyFromDevice(hp_data, dp_outData, h_numItems);  
  tledDeviceSyncDebug;

  std::copy(hp_data, hp_data + h_numItems, std::ostream_iterator<int>(cout << "\nResult:\n", " "));
  std::copy(refData.begin(), refData.end(), std::ostream_iterator<int>(cout << "\nRef:   \n", " "));
  tledUnitTestAssert(h_numItems == (int)refData.size());
  tledUnitTestAssert(std::equal(refData.begin(), refData.end(), hp_data));

  tledCheckCUDAErrors(cudaFree(dp_data));
  tledCheckCUDAErrors(cudaFree(dp_outData));
  tledCheckCUDAErrors(cudaFree(dp_numItems));
  tledCheckCUDAErrors(cudaFreeHost(hp_data));
}

template <const int t_threadDataBlockSize, const int t_blockSize>
__global__ void _TestUniqueKernel(int *p_outData, int *p_numOut, const int *inData, const int numIn) {
  __shared__ unsigned short numItems;
  __shared__ int items[t_threadDataBlockSize*t_blockSize];

  tledCudaAssert(t_threadDataBlockSize*t_blockSize >= numIn);
  for (int j = 0; j < t_threadDataBlockSize*t_blockSize; j += t_threadDataBlockSize) {
    __syncthreads();
    if (threadIdx.x + j < numIn) {
      items[threadIdx.x+j] = inData[threadIdx.x+j];
    } else if (threadIdx.x+j < t_threadDataBlockSize*t_blockSize) items[threadIdx.x+j] = 1234567;
  }
  if (threadIdx.x == 0) numItems = numIn;
  tledCUDAHelpers::Unique<t_blockSize>(items, numItems, _IntEquality());
  for (int j = 0; j < numItems; j += t_threadDataBlockSize) {
    if (threadIdx.x + j < numItems) {
      p_outData[threadIdx.x+j] = items[threadIdx.x+j];
    }
  }
  if (threadIdx.x == 0) *p_numOut = numItems;
}

template <const int t_threadDataBlockSize, const int t_blockSize>
static void _TestUnique(void) {
  const int maxNumItems = t_threadDataBlockSize*t_blockSize;
  
  int *dp_data, *hp_data, *dp_outData, *dp_numItems, h_numItems;
  std::vector<int> refData;

  h_numItems = rand()%maxNumItems + 1;
  tledCUDAHelpers::AllocateHostMemory(hp_data, h_numItems);
  for (int i = 0; i < h_numItems; i++) {
    hp_data[i] = rand()%h_numItems;
  }
  std::sort(hp_data, hp_data + h_numItems);
  refData.insert(refData.end(), hp_data, hp_data + h_numItems);
  refData.erase(std::unique(refData.begin(), refData.end()), refData.end());

  tledCUDAHelpers::AllocateDeviceMemory(dp_data, h_numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_outData, h_numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_numItems);
  tledCUDAHelpers::SetZero(dp_numItems);
  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, h_numItems);

  _TestUniqueKernel<t_threadDataBlockSize, t_blockSize> <<<1, t_blockSize>>> (dp_outData, dp_numItems, dp_data, h_numItems);
  tledDeviceSyncDebug;

  tledCUDAHelpers::CopyFromDevice(&h_numItems, dp_numItems);
  std::cout << "Copying " << h_numItems << " (of " << refData.size() << " expected) from device\n";
  tledCUDAHelpers::CopyFromDevice(hp_data, dp_outData, h_numItems);  

  tledUnitTestAssert(h_numItems == (int)refData.size());
  tledUnitTestAssert(std::equal(refData.begin(), refData.end(), hp_data));

  tledCheckCUDAErrors(cudaFree(dp_data));
  tledCheckCUDAErrors(cudaFree(dp_outData));
  tledCheckCUDAErrors(cudaFree(dp_numItems));
  tledCheckCUDAErrors(cudaFreeHost(hp_data));
}

template <const int t_threadBlockSize, const int t_blockSize>
__global__ void _TestMergeSortKernel(int *p_outData, int *p_numOut, const int *inData, const int numIn) {
  const int baseInd = (blockDim.x*blockIdx.x + threadIdx.x)*t_threadBlockSize;

  __shared__ int s_items[t_threadBlockSize*t_blockSize];
  __shared__ unsigned short s_numItems[t_blockSize];

  tledCudaAssert(numIn == blockDim.x*t_threadBlockSize);
  s_numItems[threadIdx.x] = t_threadBlockSize;
  for (int i = 0; i < t_threadBlockSize; i++) s_items[t_threadBlockSize*threadIdx.x+i] = inData[baseInd+i];

  __syncthreads();
  tledCUDAHelpers::MergeSort<t_blockSize, t_threadBlockSize, int>(s_items, s_numItems, _IntOrdering());
  if (threadIdx.x == 0) {
    *p_numOut = s_numItems[0];
    for (int i = 0; i < s_numItems[0]; i++) {
      p_outData[i] = s_items[i];
    }
  }
}

template <const int t_threadBlockSize, const int t_blockSize>
static void _TestMergeSort(void) {
  const int numItems = t_threadBlockSize*t_blockSize;
  
  int *dp_data, *hp_data, *dp_outData, *dp_numItems, h_numItems;
  std::vector<int> refData;

  tledCUDAHelpers::AllocateHostMemory(hp_data, numItems);
  for (int i = 0; i < numItems; i++) {
    hp_data[i] = rand()%(numItems - 1);
  }
  refData.insert(refData.end(), hp_data, hp_data + numItems);
  std::sort(refData.begin(), refData.end());

  tledCUDAHelpers::AllocateDeviceMemory(dp_data, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_outData, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_numItems);
  tledCUDAHelpers::SetZero(dp_numItems);
  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, numItems);

  _TestMergeSortKernel<t_threadBlockSize, t_blockSize> <<<1, t_blockSize>>> (dp_outData, dp_numItems, dp_data, numItems);
  tledDeviceSyncDebug;

  tledCUDAHelpers::CopyFromDevice(&h_numItems, dp_numItems);
  std::cout << "Copying " << h_numItems << " (of " << refData.size() << " expected) from device\n";
  tledCUDAHelpers::CopyFromDevice(hp_data, dp_outData, h_numItems);  
  tledDeviceSyncDebug;

  tledUnitTestAssert(h_numItems == (int)refData.size());
  tledUnitTestAssert(std::equal(refData.begin(), refData.end(), hp_data));

  tledCheckCUDAErrors(cudaFree(dp_data));
  tledCheckCUDAErrors(cudaFree(dp_outData));
  tledCheckCUDAErrors(cudaFree(dp_numItems));
  tledCheckCUDAErrors(cudaFreeHost(hp_data));
}

struct _KeepIfTestStruct {
  int payload;
  bool doKeep;
};

struct _KeepIfPredicate {
  __device__ bool operator()(const _KeepIfTestStruct &item) const {
    return item.doKeep;
  }
};

template <const int t_blockSize>
__global__ void _TestKeepIfKernel(_KeepIfTestStruct *p_data, int *p_numItems) {
  tledCUDAHelpers::KeepIf<t_blockSize, _KeepIfTestStruct, _KeepIfPredicate>(p_data, *p_numItems, _KeepIfPredicate());
}

template <const int t_threadBlockSize, const int t_blockSize>
static void _TestKeepIf(void) {
  const int numItems = t_threadBlockSize*t_blockSize;
  
  _KeepIfTestStruct *dp_data, *hp_data;
  int *dp_numItems, h_numItems;
  std::vector<_KeepIfTestStruct> refData;

  tledCUDAHelpers::AllocateHostMemory(hp_data, numItems);
  for (int i = 0; i < numItems; i++) {
    hp_data[i].payload = rand()%(numItems - 1);
    hp_data[i].doKeep = bool(rand()%2);
    if (hp_data[i].doKeep) refData.push_back(hp_data[i]);
  }

  tledCUDAHelpers::AllocateDeviceMemory(dp_data, numItems);
  tledCUDAHelpers::AllocateDeviceMemory(dp_numItems);
  tledCUDAHelpers::CopyToDevice(dp_data, hp_data, numItems);
  tledCUDAHelpers::CopyToDevice(dp_numItems, &numItems);

  _TestKeepIfKernel<t_blockSize> <<<1, t_blockSize>>> (dp_data, dp_numItems);
  tledDeviceSyncDebug;

  tledCUDAHelpers::CopyFromDevice(&h_numItems, dp_numItems);
  std::cout << "Copying " << h_numItems << " (of " << refData.size() << " expected) from device\n";
  tledCUDAHelpers::CopyFromDevice(hp_data, dp_data, h_numItems);  
  tledDeviceSyncDebug;

  tledUnitTestAssert(h_numItems == (int)refData.size());
  for (int i = 0; i < h_numItems; i++) {
    /* Test implies preservation of ordering, but this is not really guaranteed to remain the case! */
    tledCudaAssert(hp_data[i].doKeep);
    tledCudaAssert(hp_data[i].payload == refData[i].payload);
  }

  tledCheckCUDAErrors(cudaFree(dp_data));
  tledCheckCUDAErrors(cudaFree(dp_numItems));
  tledCheckCUDAErrors(cudaFreeHost(hp_data));
}

int main(void) {
  InitUnitTest();
  tledCUDAUnitTest::InitCUDATests();

  _TestCollect<2, 32>();
  _TestCollect<8, 32>();

  _TestCollect<3, 128>();
  _TestCollect<6, 128>();

  _TestCollect<2, 256>();
  _TestCollect<4, 256>();

  _TestKeepIf<6, 32>();
  _TestKeepIf<2, 128>();
  _TestKeepIf<3, 256>();

  _TestUnique<3, 32>();
  _TestUnique<1, 128>();
  _TestUnique<5, 256>();
  
  _TestMergeSort<3, 32>();
  _TestMergeSort<7, 32>();

  _TestMergeSort<9, 64>();
  _TestMergeSort<12, 64>();

  _TestMergeSort<3, 256>();
  _TestMergeSort<6, 256>();

  _TestSortedUnique<3, 32>();
  _TestSortedUnique<15, 32>();

  _TestSortedUnique<3, 64>();
  _TestSortedUnique<5, 64>();

  _TestSortedUnique<7, 128>();
  _TestSortedUnique<2, 128>();

  _TestSortedUnique<4, 256>();
  _TestSortedUnique<6, 256>();

  tledCUDAUnitTest::FinishCUDATests();
  tledUnitTestPrintSuccess;  

  return EXIT_SUCCESS;
}
