// =========================================================================
// File:       tledShellMaterial_kernels.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterial_kernels_H
#define tledShellMaterial_kernels_H

namespace tledShellMaterial_kernels {
  template <class TMaterial>
  __device__ void ComputeStress(float *p_stress, const float strain[], const typename TMaterial::GPUMaterial *gpc_mat);
}

#endif
