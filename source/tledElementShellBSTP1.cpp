// =========================================================================
// File:       tledElementShellBSTP1.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledElementShellBSTP1.h"
#include "tledElementShellBSTP1Edge.h"
#include "tledVectorArithmetic.h"
#include "tledShellMaterialLinearPlateDecorator.h"
#include "tledMembraneMaterialLinear.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>

void tledElementShellBSTP1::InitialiseElement(const tledSurfaceTopology<Surface> &surfTopo, const int facetIndex) {
  using namespace tledVectorArithmetic;

  Superclass::InitialiseElement(surfTopo.GetSurface(), facetIndex);

  {
    float area;

    mc_X0 = surfTopo.GetSurface().GetAllNodeCoordinates();
    ComputeElementBasis(m_ElementBasis, area, mc_X0);
    this->SetArea(area);
  }

  {
    float J[2][2], JInv[2][2];
    
    for (int xic = 0; xic < 2; xic++) {
      for (int c = 0; c < 2; c++) J[c][xic] = Dot(surfTopo.GetSurface().GetNodeCoordinates(this->GetFacetVertexIndices()[xic+1]), m_ElementBasis[c]);
    }

    for (int c = 0; c < 2; c++) {
      for (int xic = 0; xic < 2; xic++) J[c][xic] -= Dot(surfTopo.GetSurface().GetNodeCoordinates(this->GetFacetVertexIndices()[0]), m_ElementBasis[c]);
    } 
    
    MatInverse22(&JInv[0][0], &J[0][0]);
    
    for (int n = 0; n < 3; n++) {
      static const int r[][2] = {{-1, -1}, {1, 0}, {0, 1}};
      
      for (int c = 0; c < 2; c++) m_LinearShapeGradients[n][c] = JInv[0][c]*r[n][0] + JInv[1][c]*r[n][1];
    }
  } 

  for (int e = 0; e < 3; e++) {
    const int globalEdgeIndex = surfTopo.GetFacetEdges()[facetIndex][e];

    m_Edges[e].InitPatchTopologyAndGeometry(surfTopo, globalEdgeIndex);
    m_Edges[e].RegisterPrimaryTriangle(*this, (e + 2)%3);
  }
}

void tledElementShellBSTP1::InitialiseEdgeElements(std::vector<tledElementShellBSTP1> &r_elements, const tledSurfaceTopology<Surface> &surfTopo, const int facetIndex, const std::vector<int> &reverseIndexMap) {
  for (int e = 0; e < 3; e++) {
    const int globalEdgeIndex = surfTopo.GetFacetEdges()[facetIndex][e];
    const int neighbourFacetIndex = surfTopo.GetEdgeNeighbours()[globalEdgeIndex].first == facetIndex? surfTopo.GetEdgeNeighbours()[globalEdgeIndex].second : surfTopo.GetEdgeNeighbours()[globalEdgeIndex].first;

    assert(globalEdgeIndex >= 0 && globalEdgeIndex < (int)surfTopo.GetEdges().size());
    assert(neighbourFacetIndex == -1 || neighbourFacetIndex < surfTopo.GetSurface().GetNumberOfFacets());
    if (neighbourFacetIndex >= 0) {
      const int neighbourElSetIndex = reverseIndexMap[neighbourFacetIndex];

      if (neighbourFacetIndex >= 0) {
	int ne;

	for (ne = 0; surfTopo.GetFacetEdges()[neighbourFacetIndex][ne] != globalEdgeIndex && ne < 3; ne++);
	assert(ne < 3 && neighbourElSetIndex >= 0 && neighbourElSetIndex < (int)r_elements.size());
	r_elements[neighbourElSetIndex].RegisterSecondaryTriangle(ne, e, *this);
      }
    }
  }  
}

void tledElementShellBSTP1::RegisterSecondaryTriangle(const int localEdgeIndex, const int neighbourLocalEdgeIndex, const tledElementShellBSTP1 &neighbour) {
  assert(localEdgeIndex >= 0 && localEdgeIndex < 3 && neighbourLocalEdgeIndex >= 0 && neighbourLocalEdgeIndex < 3);
  m_Edges[localEdgeIndex].RegisterSecondaryTriangle(neighbour, (neighbourLocalEdgeIndex + 2)%3);
}

void tledElementShellBSTP1::FinaliseInitialisation(void) {
  float eb[3][3], area;

  ComputeElementBasis(eb, area, mc_X0);
  for (tledElementShellBSTP1Edge *p_edge = m_Edges; p_edge < m_Edges + 3; p_edge++) p_edge->InitTensors(eb);
}

void tledElementShellBSTP1::_ComputeCurrentConfiguration(float (*p_deformationGrad)[3], float (*p_nodePos)[3], const float U[]) const {
  using namespace tledVectorArithmetic;

  std::fill(&p_deformationGrad[0][0], &p_deformationGrad[0][0] + 2*3, 0.0f);
  for (int n = 0; n < 3; n++) {
    float dxCurr[3];

    Add(p_nodePos[n], mc_X0 + 3*this->GetFacetVertexIndices()[n], U + 3*this->GetFacetVertexIndices()[n]);
    for (int c = 0; c < 2; c++) {
      Add(p_deformationGrad[c], p_deformationGrad[c], ScalarMul(dxCurr, p_nodePos[n], m_LinearShapeGradients[n][c]));
    }
  }  
}

float* tledElementShellBSTP1::ComputeStrain(float *p_dst, const float U[]) {
  using namespace tledVectorArithmetic;

  float centralDeformationGradient[2][3], xs[3][3];  

  _ComputeCurrentConfiguration(centralDeformationGradient, xs, U);
  ScalarMul(m_CurrentNormal, (m_ThicknessRatio = 1/Norm(Cross(m_CurrentNormal, centralDeformationGradient[0], centralDeformationGradient[1]))));
  assert(this->GetThicknessRatio() > 0);

#ifndef NDEBUG
  {
    float e0[3], e1[3];
    float nRef[3];

    ScalarDiv(nRef, Norm(Cross(nRef, Sub(e0, xs[1], xs[0]), Sub(e1, xs[2], xs[0]))));
    assert(std::fabs(1 - Dot(m_CurrentNormal, nRef)) < 1e-2);
  }
#endif

  ScalarMul(Cross(m_ContravariantBasis[0], centralDeformationGradient[1], m_CurrentNormal), this->GetThicknessRatio());
  ScalarMul(Cross(m_ContravariantBasis[1], m_CurrentNormal, centralDeformationGradient[0]), this->GetThicknessRatio());

  std::fill(p_dst, p_dst + 6, 0.0f);
  std::fill(&m_Hs[0][0], &m_Hs[0][0] + 3*3, 0.0f);
  for (tledElementShellBSTP1Edge *p_edge = m_Edges; p_edge < m_Edges + 3; p_edge++) {
    if (p_edge->IsComplete()) {
      p_edge->ComputeDeformationGradient(xs, U);      
    } else {
      std::copy(centralDeformationGradient[0], centralDeformationGradient[0] + 3, p_edge->GetDeformationGradient(0));
      std::copy(centralDeformationGradient[1], centralDeformationGradient[1] + 3, p_edge->GetDeformationGradient(1));
    }

    if (p_edge->IsClamped()) p_edge->ComputeClampedOrSymmetryH(m_Hs, U, this->GetThicknessRatio());
    else p_edge->ComputeH(m_Hs);

    for (int c = 0; c < 2; c++) p_dst[c] += Dot(p_edge->GetDeformationGradient(c), p_edge->GetDeformationGradient(c));
    p_dst[2] += Dot(p_edge->GetDeformationGradient(1), p_edge->GetDeformationGradient(0));
  } 

  for (int c = 0; c < 2; c++) p_dst[c] = 0.5f*(p_dst[c]/3 - 1);
  p_dst[2] /= 3;
  for (int c = 0; c < 3; c++) p_dst[3+c] = Dot(GetCurrentNormal(), m_Hs[c]);

  return p_dst;
}

void tledElementShellBSTP1::ComputeForces(float *p_dst, const float stress[]) {
  using namespace tledVectorArithmetic;

  for (tledElementShellBSTP1Edge const *pc_edge = m_Edges; pc_edge < m_Edges + 3; pc_edge++) {
    pc_edge->ComputeMembraneForces(p_dst, stress);
    pc_edge->ComputeBendingForces(p_dst, stress + 3, m_Hs);
  }
}

void tledElementShellBSTP1::ComputeElementForces(float *p_dst, const float U[]) { 
  ComputeElementForcesTemplate<6, 6>(p_dst, U); 
}

float tledElementShellBSTP1::ComputeCurrentThickness(const float U[]) const { 
  return this->GetThicknessRatio()*GetMaterial().GetThickness(); 
}

void tledElementShellBSTP1::ClampNodes(const std::vector<bool> &clampedNodeMask) {
  for (int e = 0; e < 3; e++) {
    if (!m_Edges[e].IsComplete() && clampedNodeMask[m_Edges[e].m_EdgePatchNodeIndices[0]] && clampedNodeMask[m_Edges[e].m_EdgePatchNodeIndices[1]]) m_Edges[e].ClampEdge();
  } 
}

tledElementShellBSTP1::tledElementShellBSTP1() {
#ifndef NDEBUG
  std::fill(&m_LinearShapeGradients[0][0], &m_LinearShapeGradients[0][0] + 3*2, std::numeric_limits<float>::quiet_NaN());
  std::fill(m_CurrentNormal, m_CurrentNormal + 3, std::numeric_limits<float>::quiet_NaN());
  std::fill(&m_ContravariantBasis[0][0], &m_ContravariantBasis[0][0] + 2*3, std::numeric_limits<float>::quiet_NaN());
  std::fill(&m_Hs[0][0], &m_Hs[0][0] + 3*3, std::numeric_limits<float>::quiet_NaN());
#endif
}
