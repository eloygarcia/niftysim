// =========================================================================
// File:       tledSelfCollisionBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHTraverserGPU_CU
#define tledSelfCollisionBVHTraverserGPU_CU


template <class TBVH, class TAPI>
const std::vector<int>& tledSelfCollisionBVHTraverserImplGPU<TBVH, TAPI>::GetHostStartBVs() const {
  return this->GetMasterBVH().GetSelfCollisionCandidates();
}


#endif
