// =========================================================================
// File:       tledPressureConstraint.cpp
// Purpose:    Pressure constraint class. Used for applying uniform surface
//             pressure loads.
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledPressureConstraint.h"
#include "tledMatrixFunctions.h"

using namespace std;

tledPressureConstraint::tledPressureConstraint() {
}

tledPressureConstraint::tledPressureConstraint(vector<int> faces, int facetype, tledSolver* solver, float mag, enum loadShape ls) {
  this->SetLoadShape(ls);
  this->SetMag(mag);
  this->SetSolver(*solver);
  this->SetFaces(facetype, faces);
}

tledPressureConstraint::tledPressureConstraint(tledSolver &r_solver, const int faceType, const std::vector<int> &faces, const float mag, const loadShape loadShape) {
  this->SetLoadShape(loadShape);
  this->SetMag(mag);
  this->SetSolver(r_solver);
  this->SetFaces(faceType, faces);
}

tledPressureConstraint::~tledPressureConstraint() {
}

void tledPressureConstraint::SetFaces(const int type, const std::vector<int> &faces) {
  Superclass::SetFaces(type, faces);
  m_CurrentNodePressureForces.clear();
  m_CurrentNodePressureForces.insert(m_CurrentNodePressureForces.end(), this->GetSolver().GetMesh()->GetNumNodes()*3, std::numeric_limits<float>::quiet_NaN());
}

template <class TSurface>
void tledPressureConstraint::_UpdateNodePressureForces(const int step, const double dt, const double T) {
  using namespace tledVectorArithmetic;

  const double TR = dt*(step+1)/T;
  const float amp = this->ComputeAmplitude(TR, this->GetLoadShape())/TSurface::Facet::NumberOfVertices;

  std::fill(m_CurrentNodePressureForces.begin(), m_CurrentNodePressureForces.end(), 0.f);
  for (int f = 0; f < this->GetNumberOfFaces(); f++) {
    const int *nInds = static_cast<const TSurface&>(this->GetSurface()).GetFacet(f).NodeIndices;
    
    float fNodePressure[3];
    
    /* Note normal is not normalised -> ||n||*2*area */
    ScalarMul(static_cast<const TSurface&>(this->GetSurface()).ComputeFacetNormal(fNodePressure, f), amp*this->GetMag()/2);
    for (int const *pc_n = nInds; pc_n < nInds + TSurface::Facet::NumberOfVertices; pc_n++) {
      Sub(&m_CurrentNodePressureForces[*pc_n*3], &m_CurrentNodePressureForces[*pc_n*3], fNodePressure);
    }
  }
}

std::vector<float>* tledPressureConstraint::GetForceVal(int dof, int step, double dt, double T) {
  if (this->GetLastUpdateStep() < step) {
    this->UpdateGeometry(step);
    
    if (this->GetFaceType() == 1) {
      _UpdateNodePressureForces<TriangleConstraintSurface>(step, dt, T);
    } else {
      _UpdateNodePressureForces<QuadConstraintSurface>(step, dt, T);
    }
  }
   
  return this->ConvertToConstraintDOFForces(m_CurrentNodePressureForces.begin(), m_CurrentNodePressureForces.end(), dof);
}
