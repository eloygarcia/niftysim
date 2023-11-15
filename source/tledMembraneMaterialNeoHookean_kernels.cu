// =========================================================================
// File:       tledMembraneMaterialNeoHookean_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMembraneMaterialNeoHookean_kernels_CU
#define tledMembraneMaterialNeoHookean_kernels_CU
#include "tledMatrixFunctions.h"

namespace tledShellMaterial_kernels {
  template <>
  __device__ void ComputeStress<tledMembraneMaterialNeoHookean>(float *p_stress, const float strain[], const tledMembraneMaterialNeoHookean::GPUMaterial *pc_mat) {        
    const float *C0 = strain;
    const float *Cn = strain + 4;
    const float mu = pc_mat->Mu;

    float psk[2*2], invCn[2*2];

    MatInverse22(psk, C0);
    MatInverse22(invCn, Cn);

    MatSubtract(psk, MatMultScalar(invCn, 2, 2, MatDet22(C0)/MatDet22(Cn), invCn), 2, 2, psk);
    MatMultScalar(psk, 2, 2, mu, p_stress);
  }
}
#endif
