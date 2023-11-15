// =========================================================================
// File:       tledMeshLoader.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledVectorArithmetic.h"
#include "tledMeshLoader.h"

using namespace std;

void tledMeshLoader::SetMeshType(const char meshType[]) { 
  m_MeshType = meshType; 
  if (m_MeshType == "T4" || m_MeshType == "T4ANP") m_NumElementNodes = 4;
  else if (m_MeshType == "H8") m_NumElementNodes = 8;
  else {
    tledLogErrorStream(tledHelper::FatalError() << "Element type " << meshType << " not supported.");
  }
}

void tledMeshLoader::ApplyTransforms() {
  using namespace tledVectorArithmetic;

  float *p_nodesStart = this->GetOutput().GetAllNodeCds();
  float *p_nodesEnd = p_nodesStart + 3*this->GetOutput().GetNumNodes();

  if (this->GetRotationX() == this->GetRotationX()) {
    float R[3*3];

    this->AssembleRotationMatrix(R);
    for (float *p_nodeCd = p_nodesStart; p_nodeCd < p_nodesEnd; p_nodeCd += 3) {
      float xCnt[3];
      
      Sub(xCnt, p_nodeCd, this->GetCentreOfRotation());
      MatMultAB(R, 3, 3, xCnt, 3, 1, p_nodeCd);
      Add(p_nodeCd, p_nodeCd, this->GetCentreOfRotation());
    }
  }

  for (float *p_nodeCd = p_nodesStart; p_nodeCd < p_nodesEnd; p_nodeCd += 3) {
    Add(p_nodeCd, ScalarMul(p_nodeCd, this->GetScaleFactor()), this->GetTranslation());
  }
}

void tledMeshLoader::Read() {
  ReadFile();
  ApplyTransforms();
}
