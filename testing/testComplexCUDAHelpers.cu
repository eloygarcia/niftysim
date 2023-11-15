// =========================================================================
// File:       testComplexCUDAHelpers.cpp
// Purpose:    tledComplexCUDAHelpers unit test
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
#ifndef _GPU_
#error "CUDA is not available"
#endif

#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledCUDAMemoryBlock.h"
#include "tledComplexCUDAHelpers.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstdlib>

class _NumberGenerator {
private:
  const int mc_BlockSize, mc_RedundancyFactor;

public:
  int operator()() const {
    return std::rand()%(mc_BlockSize/mc_RedundancyFactor);
  }

public:
  _NumberGenerator(const int blockSize, const int redFact = 2) : mc_BlockSize(blockSize), mc_RedundancyFactor(redFact) {}
};

struct _IntLower {
  __device__ bool operator()(const int a, const int b) const {
    return a < b;
  }
};

static void _TestSort(const int dataSize, const bool makeUnique) {
  const int numTests = 40;

  tledCUDADeviceMemoryBlock &r_devMem = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int>(dataSize);
  std::vector<int> hostMem, refMem;
  int *dp_counter, hostCounter;
  
  tledCUDAHelpers::AllocateDeviceMemory(dp_counter);  
  assert(r_devMem.IsActive());
  for (int t = 0; t < numTests; t++) {
    hostMem.resize(dataSize);
    std::generate(hostMem.begin(), hostMem.end(), _NumberGenerator(dataSize));
    assert(int(hostMem.size()) == dataSize && r_devMem.IsActive());
    
    refMem = hostMem;
    std::sort(refMem.begin(), refMem.end());
    std::copy(hostMem.begin(), hostMem.end(), std::ostream_iterator<int>(std::cout << "in(unsorted) = ", " "));    
    std::copy(refMem.begin(), refMem.end(), std::ostream_iterator<int>(std::cout << "\nin(sorted) = ", " "));
    std::cout << std::endl;
    if (makeUnique) {
      refMem = tledHelper::MakeSortedUnique(hostMem);
    } else {
      refMem = hostMem;
      std::sort(refMem.begin(), refMem.end());
    }

    tledCUDAHelpers::CopyToDevice(r_devMem.GetBuffer<int>(), &hostMem.front(), hostMem.size());
    tledCUDAHelpers::CopyToDevice(dp_counter, &dataSize);    
    if (makeUnique) {
      tledComplexCUDAHelpers::MakeSortedUnique<int, _IntLower>(r_devMem, dp_counter, _IntLower(), 0);
    } else {
      tledComplexCUDAHelpers::MergeSort<int, _IntLower>(r_devMem, _IntLower(), dp_counter, 0);
    }
    tledDeviceSyncDebug;

    tledCUDAHelpers::CopyFromDevice(&hostCounter, dp_counter);
    tledCUDAHelpers::CopyFromDevice(&hostMem.front(), r_devMem.GetBuffer<int>(), hostCounter);    
    hostMem.erase(hostMem.begin() + hostCounter, hostMem.end());

    std::cout << "found/expected/initial: " << hostCounter << "/" << refMem.size() << "/" << dataSize << std::endl;
    std::cout << "ref/gpu\n";
    for (std::vector<int>::const_iterator ic_r = refMem.begin(), ic_t = hostMem.begin(); ic_r < refMem.end(); ic_r++, ic_t++) {
      if (*ic_t != *ic_r) std::cout << "!!!!!! ";
      std::cout << ic_r - refMem.begin() << ": " << *ic_r << " / " << *ic_t << std::endl;
    }		  

    tledUnitTestAssert(hostCounter == int(refMem.size()));        
    tledUnitTestAssert(std::equal(hostMem.begin(), hostMem.end(), refMem.begin()));
  }	
  r_devMem.ToggleActive();
  assert(tledCUDADeviceMemoryBlock::GetNumberOfActiveBuffers() == 0);
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();
  tledCUDAUnitTest::InitCUDATests();
  tledCUDADeviceMemoryBlock::SaveAllocationCounter();

  _TestSort(128, false);
  _TestSort(133, false);
  _TestSort(256, false);
  _TestSort(381, false);
  _TestSort(1027, false);  
  _TestSort(4125, false);

  _TestSort(128, true);
  _TestSort(132, true);
  _TestSort(255, true);
  _TestSort(381, true);
  _TestSort(1024, true);
  _TestSort(1981, true);
  _TestSort(4095, true);

  tledCUDADeviceMemoryBlock::CheckAllocationCounter();
  tledCUDAUnitTest::FinishCUDATests();
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
