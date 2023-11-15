// =========================================================================
// File:       tledSurfaceConstraint.h
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

#include "tledSurfaceConstraint.h"

#include <algorithm>
#include <limits>

void tledSurfaceConstraint::ConstraintSurface::Update(const float X0[], const float U[]) {
  for (std::vector<int>::const_iterator ic_n = m_NodeIndices.begin(); ic_n < m_NodeIndices.end(); ic_n++) {
    tledVectorArithmetic::Add(&m_NodeCoordinates[*ic_n*3], &X0[*ic_n*3], &U[*ic_n*3]);
  }
}

void tledSurfaceConstraint::ConstraintSurface::Init(const std::vector<int> &facetDefs, const int numMeshNodes) {
  m_NodeCoordinates.clear();
  m_NodeCoordinates.insert(m_NodeCoordinates.end(), 3*numMeshNodes, std::numeric_limits<float>::quiet_NaN());
  m_NodeIndices = tledHelper::MakeSortedUnique(facetDefs);
}

/***************************************************************************************************************/

template <class TSurface>
static void _InitFacetDefinitions(TSurface &r_surface, const std::vector<int> &facetDefs) {
  const int numFVtcs = TSurface::Facet::NumberOfVertices;

  r_surface.SetNumberOfFacets(facetDefs.size()/numFVtcs);
  for (int f = 0; f < r_surface.GetNumberOfFacets(); f++) {
    std::copy(facetDefs.begin() + numFVtcs*f, facetDefs.begin() + numFVtcs*(f + 1), r_surface.GetFacet(f).NodeIndices);
  }
}

void tledSurfaceConstraint::TriangleConstraintSurface::Init(const std::vector<int> &facetDefs, const int numMeshNodes) {
  ConstraintSurface::Init(facetDefs, numMeshNodes);
  _InitFacetDefinitions(*this, facetDefs);
  this->SetNodeVector(this->GetNodeVector(), this->GetNodeVectorSize()/3);
}

void tledSurfaceConstraint::QuadConstraintSurface::Init(const std::vector<int> &facetDefs, const int numMeshNodes) {
  ConstraintSurface::Init(facetDefs, numMeshNodes);
  _InitFacetDefinitions(*this, facetDefs);
  this->SetNodeVector(this->GetNodeVector(), this->GetNodeVectorSize()/3);
}

/***************************************************************************************************************/

tledSurfaceConstraint::tledSurfaceConstraint(void) : mp_CurrentConfSolver(NULL), mp_Surface(NULL), m_LastUpdateStep(-1) {}

tledSurfaceConstraint::~tledSurfaceConstraint(void) {
  if (mp_Surface != NULL) delete mp_Surface;
}

int tledSurfaceConstraint::GetNumberOfFaceVertices(void) const { 
  return this->GetSurface().GetNumberOfFacetVertices(); 
}

int tledSurfaceConstraint::GetNumberOfFaces(void) const { 
  return this->GetSurface().GetNumberOfFacets(); 
}

int tledSurfaceConstraint::GetFaceType(void) const { 
  return this->GetNumberOfFaceVertices() == 3? 1 : 0; 
}

std::vector<int>* tledSurfaceConstraint::GetForceInd(int dof) { 
  return &this->GetSurface().GetNodeIndices(); 
}

std::vector<int>& tledSurfaceConstraint::GetInd() { 
  return this->GetSurface().GetNodeIndices(); 
}

const std::vector<int>& tledSurfaceConstraint::GetInd() const { 
  return this->GetSurface().GetNodeIndices(); 
}

tledSurfaceConstraint::ConstraintSurface* tledSurfaceConstraint::CreateSurface(const int type, const std::vector<int> &faces) {
  if (type == 1) {
    return new TriangleConstraintSurface();
  } else if (type == 0) {
    return new QuadConstraintSurface();
  } else {
    tledLogErrorStream(tledHelper::FatalError() << "Invalid face type (in {0, 1}): " << type);
  }

  return NULL;
}

void tledSurfaceConstraint::SetFaces(const int type, const std::vector<int> &faces) {  
  if (!this->HaveSolver()) {
    tledFatalError("SetFaces called before setting of solver!");
  }

  mp_Surface = this->CreateSurface(type, faces);
  mp_Surface->Init(faces, this->GetSolver().GetMesh()->GetNumNodes());

  m_ConstraintForces.reserve(this->GetInd().size());
}

void tledSurfaceConstraint::UpdateGeometry(const int step) {
  if (step > m_LastUpdateStep) {
    m_LastUpdateStep = step;
    this->GetSurface().Update(this->GetSolver().GetMesh()->GetAllNodeCds(), this->GetSolver().GetAllDisps());
  }
}

std::vector<float>* tledSurfaceConstraint::ConvertToConstraintDOFForces(const std::vector<float>::const_iterator allForcesBegin, const std::vector<float>::const_iterator allForcesEnd, const int dof) {
  assert(allForcesEnd - allForcesBegin >= this->GetSurface().GetNumberOfNodes());

  m_ConstraintForces.clear();
  for (std::vector<int>::const_iterator ic_n = this->GetSurface().GetNodeIndices().begin(); ic_n < this->GetSurface().GetNodeIndices().end(); ic_n++) {
    assert(allForcesBegin + *ic_n*3 < allForcesEnd);
    m_ConstraintForces.push_back(*(allForcesBegin + *ic_n*3 + dof));
  }
  assert(m_ConstraintForces.size() == this->GetSurface().GetNodeIndices().size());

  return &m_ConstraintForces;
}
