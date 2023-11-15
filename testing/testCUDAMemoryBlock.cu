// =========================================================================
// File:       testCUDAMemoryBlock.cu
// Purpose:    tledCUDAMemoryBlock unit test
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
#ifndef _CUDA_3MINUS_SDK
#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledCUDAMemoryBlock.h"

using namespace std;

struct _TestStruct {
  float3 v0, v1;
  char a, b, c;
  int d, e, f;
};

template <typename TValue>
struct _TestValueGenerator {
  TValue operator()() {
    return TValue(rand());
  }
};

template <class TCUDAMemoryBlock, typename TValue>
static void _TestSize() {
  TCUDAMemoryBlock block;
  int sz = 0;

  block.template Resize<TValue>(sz);
  tledUnitTestAssert(block.GetBufferSize() == 0 && block.template GetMaxNumberOfItems<TValue>() == 0);

  sz = (rand()%1024) + 1;
  block.template Resize<TValue>(sz);
  tledUnitTestAssert(block.GetBufferSize() == (int)(sz*sizeof(TValue)) && block.template GetMaxNumberOfItems<TValue>() == sz);
  tledUnitTestAssert(block.IsActive());
  block.ToggleActive();

  block.template Resize<TValue>(sz/2);
  tledUnitTestAssert(block.GetBufferSize() == (int)(sz*sizeof(TValue)) && block.template GetMaxNumberOfItems<TValue>() == sz);  
  tledUnitTestAssert(block.IsActive());
  block.ToggleActive();

  block.template Resize<char>(sz);
  tledUnitTestAssert(block.GetBufferSize() == (int)(sz*sizeof(TValue)) && block.template GetMaxNumberOfItems<TValue>() == sz);  
  tledUnitTestAssert(block.IsActive());
  block.ToggleActive();
}

template <typename TValue>
static void _TestGrow() {
  const int initSize = 13;
  const int growSize = 3*initSize/2;

  vector<TValue> refVals(growSize), testVals(growSize);
  tledCUDADeviceMemoryBlock block;

  generate(refVals.begin(), refVals.end(), _TestValueGenerator<TValue>());

  block.Resize<TValue>(initSize);
  tledCUDAHelpers::CopyToDevice(block.GetBuffer<TValue>(), &refVals.front(), initSize);
  block.Grow<TValue>(growSize);
  tledCUDAHelpers::CopyToDevice(block.GetBuffer<TValue>() + initSize, &refVals.front() + initSize, growSize - initSize);
  tledCUDAHelpers::CopyFromDevice(&testVals.front(), block.GetBuffer<TValue>(), growSize);
  tledCheckCUDAErrors(cudaDeviceSynchronize());
  tledUnitTestAssert(equal(testVals.begin(), testVals.end(), refVals.begin()));
}

template <class TCUDAMemoryBlock, typename TValue>
static void _TestManagement() {
  TCUDAMemoryBlock *pp_blocks[10];
  int sizes[10];
  
  TCUDAMemoryBlock::ReleaseAll();
  tledUnitTestAssert(TCUDAMemoryBlock::GetNumberOfActiveBuffers() == 0);
  for (int i = 0; i < 10; i++) {
    int sz = -1;

    pp_blocks[i] = &TCUDAMemoryBlock::template GetNextFreeBufferWithSize<TValue>((sz = (rand()%1320) + 1));
    tledUnitTestAssert(pp_blocks[i]->template GetMaxNumberOfItems<TValue>() == sz);
    sizes[i] = sz;
  }
  tledUnitTestAssert(TCUDAMemoryBlock::GetNumberOfActiveBuffers() == 10);

  for (int i = 0; i < 10; i += 2) {
    pp_blocks[i]->ToggleActive();
    tledUnitTestAssert(!pp_blocks[i]->IsActive());
  }
  tledUnitTestAssert(TCUDAMemoryBlock::GetNumberOfActiveBuffers() == 5);

  for (int i = 1; i < 10; i += 2) {
    tledUnitTestAssert(pp_blocks[i]->IsActive());
    tledUnitTestAssert(pp_blocks[i]->template GetMaxNumberOfItems<TValue>() == sizes[i]);
  }
  TCUDAMemoryBlock::ReleaseAll();
  tledUnitTestAssert(TCUDAMemoryBlock::GetNumberOfActiveBuffers() == 0);
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();
  tledCUDAUnitTest::InitCUDATests();

  _TestSize<tledCUDADeviceMemoryBlock, char>();
  _TestSize<tledCUDADeviceMemoryBlock, int>();
  _TestSize<tledCUDADeviceMemoryBlock, _TestStruct>();

  _TestSize<tledCUDAHostMemoryBlock, char>();
  _TestSize<tledCUDAHostMemoryBlock, int>();
  _TestSize<tledCUDAHostMemoryBlock, _TestStruct>();

  _TestGrow<int>();
  _TestGrow<char>();
  _TestGrow<double>();

  _TestManagement<tledCUDADeviceMemoryBlock, char>();
  _TestManagement<tledCUDADeviceMemoryBlock, int>();
  _TestManagement<tledCUDADeviceMemoryBlock, _TestStruct>();

  _TestManagement<tledCUDAHostMemoryBlock, char>();
  _TestManagement<tledCUDAHostMemoryBlock, int>();
  _TestManagement<tledCUDAHostMemoryBlock, _TestStruct>();

  tledCUDAUnitTest::FinishCUDATests();
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
#else
#include <iostream>

int main(int argc, char *argv[]) {
  std::cerr << "Warning: test omitted due to old CUDA toolkit version.\n";

  return EXIT_SUCCESS;
} /* main */
#endif
