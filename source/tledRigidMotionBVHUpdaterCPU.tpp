// =========================================================================
// File:       tledRigidMotionBVHUpdaterCPU.tpp
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

template <class TBVH>
void tledRigidMotionBVHUpdaterCPU<TBVH>::UpdateBVH() {
  const ContactMesh &mesh = this->GetBVH().GetMesh();

  if (!mesh.HasRotation() && mesh.GetHistoryLength() < mesh.GetCurrentStep() && mesh.GetCurrentStep()%100 != 0) {
    float tInc[3], tCurr[3];

    tledVectorArithmetic::Sub(tInc, mesh.GetTranslation(tCurr), this->GetRootTranslation());
    this->SetRootTranslation(tCurr);
    this->GetBVH().TranslateSubTree(0, tInc);
  } else this->GetBVH().UpdateTopDownRecursive(0);
}
