// =========================================================================
// File:       tledRigidMotionBVHUpdaterGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledRigidMotionBVHUpdaterGPU_kernels.tpp"

template <class TBVH>
void tledRigidMotionBVHUpdaterGPU<TBVH>::_TranslateBVH(const float tInc[]) {
  const int blockSize = 256;
  const int numBlocks = tledCUDAHelpers::GetNumberOfBlocks(this->GetBVH().GetNumberOfBVs(), blockSize);

  tledRigidMotionBVHUpdaterGPU_kernels::TranslateBVHKernel <<<numBlocks, blockSize>>> (this->GetBVH().GetOnDeviceBVs(), this->GetBVH().GetNumberOfBVs(), tledCUDAHelpers::ConvertToFloat3(tInc));
}

template <class TBVH>
void tledRigidMotionBVHUpdaterGPU<TBVH>::UpdateBVH() {
  const ContactMesh &mesh = this->GetBVH().GetMesh();

  if (mesh.HasRotation()) {
    tledFatalNotYetImplementedError;
  } else {
    float tInc[3], tCurr[3];

    tledVectorArithmetic::Sub(tInc, mesh.GetTranslation(tCurr), this->GetRootTranslation());
    this->SetRootTranslation(tCurr);
    _TranslateBVH(tInc);
  } 
}
