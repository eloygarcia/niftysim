// =========================================================================
// File:       tledBVH.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVH_CU
#define tledBVH_CU

#include "tledBVH.h"
#include "tledCUDAHelpers.h"

void AllocateStaticIndexList(int* &rdp_list, const std::vector<int> &list) {
  tledCUDAHelpers::AllocateDeviceMemory(rdp_list, list.size());
  tledCUDAHelpers::CopyToDevice(rdp_list, &list.front(), list.size());
}

template <>
void AllocateGPUBVs(tledAABB<2>::GPUBV* &rdp_bvs, tledAABB<2>::GPUBV* &rhp_bvs, const std::vector<tledAABB<2> > &hostBVs) {
}

template <>
void CopyGPUBVs(tledAABB<2>::GPUBV *rdp_bvs, const tledAABB<2>::GPUBV hpc_bvs[], const int numBVs) {
  tledCUDAHelpers::CopyToDevice(rdp_bvs, hpc_bvs, numBVs);
}
#endif
