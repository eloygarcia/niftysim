// =========================================================================
// File:       tledElementMembraneNonLinear.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledVectorArithmetic.h"
#include "tledMatrixFunctions.h"
#include "tledElementMembraneNonLinear.h"

//TEST
#include "tledMembraneMaterialNeoHookean.h"

using namespace tledVectorArithmetic;

template <>
void tledElementMembraneNonLinear<3>::_ComputeCauchyGreenTensor(float *p_tensor, const float J[][3]) {
  for (int r = 0; r < 2; r++) for (int c = 0; c < 2; c++) p_tensor[2*r+c] = Dot(J[r], J[c]);
  assert(*(p_tensor + 1) == *(p_tensor + 2));
}

template <>
void tledElementMembraneNonLinear<3>::InitialiseElement(const Surface &surf, const int facetIndex) {
  Superclass::InitialiseElement(surf, facetIndex);

  Sub(m_InitialJacobian[0], surf.GetNodeCoordinates(this->GetFacetVertexIndices()[1]), surf.GetNodeCoordinates(this->GetFacetVertexIndices()[0]));
  Sub(m_InitialJacobian[1], surf.GetNodeCoordinates(this->GetFacetVertexIndices()[2]), surf.GetNodeCoordinates(this->GetFacetVertexIndices()[0]));

  {
    float tmp[3];
	
    this->SetArea(Norm(Cross(tmp, m_InitialJacobian[0], m_InitialJacobian[1]))/2);
  }

  _ComputeCauchyGreenTensor(&m_InitialCauchyGreenTensor[0][0], m_InitialJacobian);
}

template <>
float* tledElementMembraneNonLinear<3>::ComputeStrain(float *p_dst, const float U[]) {
  std::copy(&m_InitialCauchyGreenTensor[0][0], &m_InitialCauchyGreenTensor[0][0] + 4, p_dst);
  for (int vInd = 0; vInd < 2; vInd++) Sub(m_J[vInd], Add(m_J[vInd], m_InitialJacobian[vInd], U + 3*this->GetFacetVertexIndices()[vInd+1]), U + 3*this->GetFacetVertexIndices()[0]);
  _ComputeCauchyGreenTensor(p_dst + 4, m_J);
  std::copy(p_dst, p_dst + 8, m_LastStrain);

  return p_dst;
}

template <>
void tledElementMembraneNonLinear<3>::ComputeForces(float *p_dst, const float psk[]) {
  for (int nInd = 0; nInd < 3; nInd++) {
    static const float dN[][2] = {{-1, -1}, {1, 0}, {0, 1}};

    for (int dofInd = 0; dofInd < 3; dofInd++) {
      float tmpF;

      tmpF = 0;
      for (int c = 0; c < 2; c++) tmpF += (this->m_J[0][dofInd]*psk[c] + this->m_J[1][dofInd]*psk[2+c])*dN[nInd][c];
      *(p_dst + 3*this->GetFacetVertexIndices()[nInd] + dofInd) += tmpF*this->GetArea();
    }
  }      																	
}

template <>
float tledElementMembraneNonLinear<3>::ComputeCurrentThickness(const float U[]) const {
  float II_C0, II_CN;

  II_C0 = MatDet22(m_LastStrain), II_CN = MatDet22(m_LastStrain + 4);
  
  return GetMaterial().GetThickness()*std::sqrt(II_C0/II_CN);
}
