// =========================================================================
// File:       tledShellMaterialLinearPlateDecorator_kernels.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   Cuda
// Created:    August 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMaterialLinearPlateDecorator_kernels_CU
#define tledShellMaterialLinearPlateDecorator_kernels_CU

namespace tledShellMaterial_kernels {
  template <class TMembraneMaterial>
  __device__ void _ComputeStress(float *p_stress, const float strain[], const typename tledShellMaterialLinearPlateDecorator<TMembraneMaterial>::GPUMaterial *pc_mat) {        
    ComputeStress<TMembraneMaterial>(p_stress, strain, pc_mat);

    {
      const float *curvs = strain + TMembraneMaterial::NumberOfStrainComponents;
      const float a = pc_mat->BendingE*pc_mat->Thickness*pc_mat->Thickness*pc_mat->Thickness/(12*(1 - pc_mat->BendingNu*pc_mat->BendingNu));
    
      float *p_stressOut = p_stress + TMembraneMaterial::NumberOfStressComponents;

      p_stressOut[0] = a*(curvs[0] + curvs[1]*pc_mat->BendingNu);
      p_stressOut[1] = a*(curvs[0]*pc_mat->BendingNu + curvs[1]);
      p_stressOut[2] = a*curvs[2]*(1 - pc_mat->BendingNu)/2;
    }
  }

  template <>
  __device__ void ComputeStress<tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear> >(float *p_stress, const float strain[], const tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear>::GPUMaterial *pc_mat) {        
    _ComputeStress<tledMembraneMaterialLinear>(p_stress, strain, pc_mat);
  }
}


#endif
