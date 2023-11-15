// =========================================================================
// File:       tledMovingRigidContactSurfaceGPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseSurface> 
void tledMovingRigidContactSurfaceImplGPU<TBaseSurface>::TranslateNodes() {
  float tCurr[3], tOld[3];

  tledMovingRigidContactSurfaceGPU::TranslateNodes(this->GetAllOnDeviceNodeCoordinates(), this->GetOnDeviceOldNodeCoordinateBuffer(), this->GetAllOnDeviceNodeCoordinates0(), this->GetNumberOfNodes(), this->GetTranslation(tCurr), this->GetTranslation(tOld, std::max(0, this->GetCurrentStep() - this->GetHistoryLength())));
}

template <class TBaseSurface> 
void tledMovingRigidContactSurfaceImplGPU<TBaseSurface>::Update() {
  Superclass::Update();
  if (this->HasRotation()) {
    tledFatalNotYetImplementedError;
  } else if (tledVectorArithmetic::Norm(this->GetTotalTranslation()) > 0.f) {
    this->TranslateNodes();
  }
}
