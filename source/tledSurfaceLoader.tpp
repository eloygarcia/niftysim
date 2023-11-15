// =========================================================================
// File:       tledSurfaceLoader.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
template <class TSurface>
void tledSurfaceLoader<TSurface>::ApplyTransforms() {
  using namespace tledVectorArithmetic;

  if (this->GetRotationX() == this->GetRotationX()) {
    float R[3*3];

    this->AssembleRotationMatrix(R);
    for (float *p_nodeCd = this->GetOutput().GetAllNodeCoordinates(); p_nodeCd < this->GetOutput().GetAllNodeCoordinates() + 3*this->GetOutput().GetNumberOfNodes(); p_nodeCd += 3) {
      float xCnt[3];
      
      Sub(xCnt, p_nodeCd, this->GetCentreOfRotation());
      MatMultAB(R, 3, 3, xCnt, 3, 1, p_nodeCd);
      Add(p_nodeCd, p_nodeCd, this->GetCentreOfRotation());
    }
  }

  for (float *p_nodeCd = this->GetOutput().GetAllNodeCoordinates(); p_nodeCd < this->GetOutput().GetAllNodeCoordinates() + 3*this->GetOutput().GetNumberOfNodes(); p_nodeCd += 3) {
    Add(p_nodeCd, ScalarMul(p_nodeCd, this->GetScaleFactor()), this->GetTranslation());
  }  
}

template <class TSurface>
void tledSurfaceLoader<TSurface>::Read() {
  this->ReadFile();
  this->ApplyTransforms();
}

template <class TSurface>
void tledSurfaceLoaderSurfaceCreatorAdapter<TSurface>::Create() {
  this->GetLoader().SetOutputMesh(this->GetOutput());
  this->GetLoader().Read(); 
}
