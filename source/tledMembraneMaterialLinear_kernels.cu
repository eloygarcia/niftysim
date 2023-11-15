// =========================================================================
// File:       tledMembraneMaterialLinear_kernels.cu
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
#ifndef tledMembraneMaterialLinear_kernels_CU
#define tledMembraneMaterialLinear_kernels_CU

#include "tledMembraneMaterialLinear.h"
#include "tledShellMaterial_kernels.h"

namespace tledShellMaterial_kernels {
  template <>
  __device__ void ComputeStress<tledMembraneMaterialLinear>(float *p_stress, const float strain[], const tledMembraneMaterialLinear::GPUMaterial *pc_mat) {        
    p_stress[0] = strain[0] + pc_mat->Nu*strain[1];
    p_stress[1] = pc_mat->Nu*strain[0] + strain[1];
    p_stress[2] = (1 - pc_mat->Nu)/2*strain[2];

    for (int cInd = 0; cInd < 3; cInd++) p_stress[cInd] *= pc_mat->E/(1 - pc_mat->Nu*pc_mat->Nu);
  }
}
#endif
