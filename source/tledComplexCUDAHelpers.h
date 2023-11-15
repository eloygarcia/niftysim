// =========================================================================
// File:       tledComplexCUDAHelpers.h
// Purpose:    
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
#ifndef tledComplexCUDAHelpers_H
#define tledComplexCUDAHelpers_H
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#include "tledCUDAMemoryBlock.h"

/**
 * \brief Namespace for complex CUDA utility functions with multiple module dependencies.
 * \ingroup helper
 */
namespace tledComplexCUDAHelpers {
  /**
   * \brief Streamable sort and unique 
   */
  template <typename TComparable, class TStrictOrdering>
  __host__ void MakeSortedUnique(tledCUDADeviceMemoryBlock &r_data, int *dp_numItems, TStrictOrdering ordering, const cudaStream_t stream);

  /**
   * \brief Streamable merge sort
   */
  template <typename TComparable, class TStrictOrdering>
  __host__ void MergeSort(tledCUDADeviceMemoryBlock &r_data, TStrictOrdering ordering, const int *dpc_numItems, const cudaStream_t stream);  
}

#ifdef __CUDACC__
#include "tledComplexCUDAHelpers.tpp"
#endif

#endif /* _GPU_ */
#endif
