// =========================================================================
// File:       tledTractionConstraint.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledTractionConstraint.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

#include <limits>
#include <algorithm>
#include <cassert>

void tledTractionConstraint::SetFaceTractions(const std::vector<float> &ts) { 
  if (ts.size() < 3) {
    tledFatalError("Traction vector cannot have size < 3.");
  }

  if (ts.size() == 3) {
    m_FaceTractions.clear();
    m_FaceTractions.reserve(this->GetNumberOfFaces()*3);
    for (int f = 0; f < this->GetNumberOfFaces(); f++) m_FaceTractions.insert(m_FaceTractions.end(), ts.begin(), ts.end());
  } else {
    if (int(ts.size()) != this->GetNumberOfFaces()*3) {
      tledFatalError("Traction vector must have size = 3 x number of faces.");
    }

    m_FaceTractions = ts; 
  }
}

tledTractionConstraint::tledTractionConstraint(tledSolver &r_solver, const int faceType, const std::vector<int> &faces, const std::vector<float> &tractions, const loadShape loadShape) {
  this->SetSolver(r_solver);
  this->SetFaces(faceType, faces);
  this->SetFaceTractions(tractions);
  this->SetLoadShape(loadShape);
}

void tledTractionConstraint::SetFaces(const int type, const std::vector<int> &faces) {
  Superclass::SetFaces(type, faces);
  m_CurrentNodeTractions.clear();
  m_CurrentNodeTractions.insert(m_CurrentNodeTractions.end(), this->GetSolver().GetMesh()->GetNumNodes()*3, std::numeric_limits<float>::quiet_NaN());
}

template <class TSurface>
void tledTractionConstraint::_UpdateNodeForces(const int step, const double dt, const double T) {
  using namespace tledVectorArithmetic;

  const double TR = dt*(step+1)/T;
  const float amp = this->ComputeAmplitude(TR, this->GetLoadShape())/TSurface::Facet::NumberOfVertices;

  std::fill(m_CurrentNodeTractions.begin(), m_CurrentNodeTractions.end(), 0.f);
  for (int f = 0; f < this->GetNumberOfFaces(); f++) {
    const int *nInds = static_cast<const TSurface&>(this->GetSurface()).GetFacet(f).NodeIndices;
    const float area = static_cast<const TSurface&>(this->GetSurface()).ComputeFacetArea(f);

    float fNodeTrac[3];

    ScalarMul(fNodeTrac, this->GetFacetPeakTraction(f), -amp*area);
    for (int const *pc_n = nInds; pc_n < nInds + TSurface::Facet::NumberOfVertices; pc_n++) {
      Add(&m_CurrentNodeTractions[*pc_n*3], &m_CurrentNodeTractions[*pc_n*3], fNodeTrac);
    }
  }
}

std::vector<float>* tledTractionConstraint::GetForceVal(int dof, int step, double dt, double T) {
  if (this->GetLastUpdateStep() < step) {
    this->UpdateGeometry(step);
    
    if (this->GetFaceType() == 1) {
      _UpdateNodeForces<TriangleConstraintSurface>(step, dt, T);
    } else {
      _UpdateNodeForces<QuadConstraintSurface>(step, dt, T);
    }
  }
       
  return this->ConvertToConstraintDOFForces(m_CurrentNodeTractions.begin(), m_CurrentNodeTractions.end(), dof);
}
