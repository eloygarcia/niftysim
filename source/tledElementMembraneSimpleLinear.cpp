// =========================================================================
// File:       tledElementMembraneSimpleLinear.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include <algorithm>
#include <iterator>
#include <iostream>

#include "tledVectorArithmetic.h"
#include "tledMatrixFunctions.h"
#include "tledElementMembraneSimpleLinear.h"
#include "tledMembraneMaterialLinear.h"

template <>
void tledElementMembraneSimpleLinear<3>::InitialiseElement(const Surface &surf, const int facetIndex) {
  using namespace tledVectorArithmetic;

  float J[2][3], phi[2][2];

  Superclass::InitialiseElement(surf, facetIndex);

  for (int v = 0; v < 3; v++) std::copy(surf.GetNodeCoordinates(this->GetFacetVertexIndices()[v]), surf.GetNodeCoordinates(this->GetFacetVertexIndices()[v]) + 3, m_X0[v]);

  Sub(J[0], m_X0[1], m_X0[0]);
  Sub(J[1], m_X0[2], m_X0[0]);

  ScalarDiv(m_ElementBasis[0], J[0], Norm(J[0]));  
  surf.ComputeFacetNormal(m_ElementBasis[2], facetIndex);  
  this->SetArea(Norm(m_ElementBasis[2]));
  ScalarDiv(m_ElementBasis[2], this->GetArea());
  ScalarDiv(m_ElementBasis[1], Norm(Cross(m_ElementBasis[1], m_ElementBasis[2], m_ElementBasis[0])));

  this->SetArea(this->GetArea()/2);
  assert(this->GetArea() > 0);

  for (int r = 0; r < 2; r++) for (int c = 0; c < 2; c++) phi[r][c] = Dot(m_ElementBasis[r], J[c]);
  assert(std::fabs(phi[1][0]) < 1e-3*std::fabs(phi[0][0]));
  MatInverse22(&m_PhiInvT[0][0], &phi[0][0]);
  std::iter_swap(&m_PhiInvT[0][0] + 1, &m_PhiInvT[0][0] + 1*2);
}

template <>
float* tledElementMembraneSimpleLinear<3>::ComputeStrain(float *p_dst, const float U[]) {
  using namespace tledVectorArithmetic;

  float uElem[2*2];

  // for (int vInd = 0; vInd < 2; vInd++) {
  //   Add(m_ElementBasis[vInd], m_X0[vInd+1], U + 3*m_FacetVertexIndices[vInd+1]);
  //   Sub(m_ElementBasis[vInd], Sub(m_ElementBasis[vInd], m_ElementBasis[vInd], m_X0[0]), U + 3*m_FacetVertexIndices[0]);
  // }
  // ScalarDiv(m_ElementBasis[0], Norm(m_ElementBasis[0]));
  // Cross(m_ElementBasis[2], m_ElementBasis[0], m_ElementBasis[1]);
  // ScalarDiv(m_ElementBasis[2], Norm(m_ElementBasis[2]));
  // Cross(m_ElementBasis[1], m_ElementBasis[2], m_ElementBasis[0]);
  //for (float (*p_e)[3] = m_ElementBasis; p_e < m_ElementBasis + 3; p_e++) ScalarDiv(*p_e, Norm(*p_e));

  for (int vInd = 0; vInd < 2; vInd++) {
    float uLocal[3];

    Sub(uLocal, U + 3*this->GetFacetVertexIndices()[vInd+1], U + 3*this->GetFacetVertexIndices()[0]);

    uElem[2*vInd] = Dot(m_ElementBasis[0], uLocal);
    uElem[2*vInd+1] = Dot(m_ElementBasis[1], uLocal);
    assert(!std::isnan(uElem[2*vInd+1] + uElem[2*vInd]));
  }

  p_dst[0] = m_PhiInvT[0][0]*uElem[0] + m_PhiInvT[0][1]*uElem[2];
  p_dst[1] = m_PhiInvT[1][0]*uElem[1] + m_PhiInvT[1][1]*uElem[3];
  p_dst[2] = m_PhiInvT[1][0]*uElem[0] + m_PhiInvT[0][0]*uElem[1];
  p_dst[2] += m_PhiInvT[1][1]*uElem[2] + m_PhiInvT[0][1]*uElem[3];

  return p_dst;
}

template <>
void tledElementMembraneSimpleLinear<3>::ComputeForces(float *p_dst, const float planeStress[]) {
  for (int vInd = 0; vInd < 3; vInd++) {
    float planeF[2];

    switch (vInd) {
    case 0:
      planeF[0] = -(m_PhiInvT[0][0] + m_PhiInvT[0][1])*planeStress[0] - (m_PhiInvT[1][0] + m_PhiInvT[1][1])*planeStress[2];
      planeF[1] = -(m_PhiInvT[1][0] + m_PhiInvT[1][1])*planeStress[1] - (m_PhiInvT[0][0] + m_PhiInvT[0][1])*planeStress[2];
      break;

    case 1:
    case 2:
      planeF[0] = m_PhiInvT[0][vInd-1]*planeStress[0] + m_PhiInvT[1][vInd-1]*planeStress[2];
      planeF[1] = m_PhiInvT[1][vInd-1]*planeStress[1] + m_PhiInvT[0][vInd-1]*planeStress[2];
      break;
    }
    
    for (int cInd = 0; cInd < 3; cInd++) {
      p_dst[3*this->GetFacetVertexIndices()[vInd]+cInd] += m_ElementBasis[0][cInd]*planeF[0] + m_ElementBasis[1][cInd]*planeF[1];
    }
  }
}

template <>
float tledElementMembraneSimpleLinear<3>::ComputeCurrentThickness(const float U[]) const {
  using namespace tledVectorArithmetic;

  const tledMembraneMaterialLinear &mat = static_cast<const tledMembraneMaterialLinear&>(GetMaterial());

  float R, currX[3][3];

  for (int v = 0; v < 3; v++) Add(currX[v], m_X0[v], U + 3*this->GetFacetVertexIndices()[v]);
  Sub(currX[1], currX[1], currX[0]);
  Sub(currX[2], currX[2], currX[0]);
  R = Norm(Cross(currX[0], currX[1], currX[2]))/(2*GetArea());

  return mat.GetThickness()*(1 - mat.GetNu()*(R - 1)/(1 - mat.GetNu()));
}
