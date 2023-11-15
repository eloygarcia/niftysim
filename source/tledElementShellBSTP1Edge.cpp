// =========================================================================
// File:       tledElementShellBSTP1Edge.cpp
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

#include "tledElementShellBSTP1Edge.h"
#include "tledElementShellBSTP1.h"
#include "tledVectorArithmetic.h"

#include <cassert>

const float tledElementShellBSTP1Edge::s_PatchShapeGradients[][4][2] = {{{0.5, -0.5}, {-0.5, 0.5}, {-0.5, -0.5}, {0.5, 0.5}},
									{{-0.5, -1}, {0.5, 1}, {0.5, 0}, {-0.5, 0}},
									{{-1, -0.5}, {1, 0.5}, {0, 0.5}, {0, -0.5}}};

void tledElementShellBSTP1Edge::ComputeDeformationGradient(const float triangleXs[][3], const float U[]) {
  using namespace tledVectorArithmetic;

  assert(IsComplete());
  
  {
    float neighbourX[3];

    Add(neighbourX, mc_X0s + 3*m_EdgePatchNodeIndices[3], U + 3*m_EdgePatchNodeIndices[3]);
    for (int c = 0; c < 2; c++) ScalarMul(m_DeformationGradient[c], neighbourX, m_PatchShapeGradients[3][c]);
  }

  for (int n = 0; n < 3; n++) {
    for (int c = 0; c < 2; c++) {    
      float tmp[3];

      Add(m_DeformationGradient[c], m_DeformationGradient[c], ScalarMul(tmp, triangleXs[(n+m_ElementEdgeIndex+1)%3], m_PatchShapeGradients[n][c]));	
    }
  } 
}
 
void tledElementShellBSTP1Edge::ComputeH(float (*p_hs)[3]) {
  using namespace tledVectorArithmetic;

  const float *linearGradient = mpc_Owner->GetLinearGradient(m_ElementEdgeIndex);

  float hLoc[3];

  for (int c = 0; c < 2; c++) Add(p_hs[c], p_hs[c], ScalarMul(hLoc, m_DeformationGradient[c], 2*linearGradient[c]));
  Add(p_hs[2], p_hs[2], ScalarMul(hLoc, m_DeformationGradient[0], 2*linearGradient[1]));
  Add(p_hs[2], p_hs[2], ScalarMul(hLoc, m_DeformationGradient[1], 2*linearGradient[0]));
}

void tledElementShellBSTP1Edge::ComputeClampedOrSymmetryH(float (*p_hs)[3], const float U[], const float tr) {
  using namespace tledVectorArithmetic;

  const float *linearGradient = mpc_Owner->GetLinearGradient(m_ElementEdgeIndex);

  float phiS[3], phiN[3], tmp[3];

  assert(!std::isnan(m_L0));
  Sub(phiS, Add(phiS, Sub(phiS, mc_X0s + 3*m_EdgePatchNodeIndices[1],  mc_X0s + 3*m_EdgePatchNodeIndices[0]), U + 3*m_EdgePatchNodeIndices[1]), U + 3*m_EdgePatchNodeIndices[0]);
  ScalarMul(phiS, 2/m_L0);  
  ScalarDiv(phiN, m_PhiN0, Norm(phiS)*tr/4);

  Add(p_hs[0], p_hs[0], ScalarMul(tmp, phiS, -linearGradient[0]*m_EdgeNormal[1]));
  Add(p_hs[0], p_hs[0], ScalarMul(tmp, phiN, linearGradient[0]*m_EdgeNormal[0]));

  Add(p_hs[1], p_hs[1], ScalarMul(tmp, phiS, linearGradient[1]*m_EdgeNormal[0]));
  Add(p_hs[1], p_hs[1], ScalarMul(tmp, phiN, linearGradient[1]*m_EdgeNormal[1]));

  Add(p_hs[2], p_hs[2], ScalarMul(tmp, phiS, linearGradient[0]*m_EdgeNormal[0] - linearGradient[1]*m_EdgeNormal[1]));
  Add(p_hs[2], p_hs[2], ScalarMul(tmp, phiN, linearGradient[1]*m_EdgeNormal[0] + linearGradient[0]*m_EdgeNormal[1]));
}

void tledElementShellBSTP1Edge::ClampEdge(void) {
  using namespace tledVectorArithmetic;

  float n3[3], e[3], tmp[3];

  Cross(n3, Sub(e, mc_X0s + 3*m_EdgePatchNodeIndices[1], mc_X0s + 3*m_EdgePatchNodeIndices[0]), mpc_Owner->GetElementAxis(2));
  m_L0 = Norm(e);  
  if (Dot(n3, Sub(tmp, mc_X0s + 3*m_EdgePatchNodeIndices[2], mc_X0s + 3*m_EdgePatchNodeIndices[0])) > 0) ScalarMul(n3, -1);
  for (int c = 0; c < 2; c++) m_EdgeNormal[c] = Dot(n3, mpc_Owner->GetElementAxis(c));
  for (int c = 0; c < 2; c++) m_EdgeNormal[c] /= Norm(m_EdgeNormal, 2);

  ScalarDiv(m_PhiN0, n3, Norm(n3));
}

void tledElementShellBSTP1Edge::ComputeMembraneForces(float *p_fU, const float sigma[]) const {
  using namespace tledVectorArithmetic;

  const float area = mpc_Owner->GetArea();

  float tmpF[6], tmp[3];
  float sigmaEdge[3];    

  for (int c = 0; c < 3; c++) sigmaEdge[c] = sigma[c]/3;

  ScalarMul(tmpF, m_DeformationGradient[0], sigmaEdge[0]);
  Add(tmpF, tmpF, ScalarMul(tmp, m_DeformationGradient[1], sigmaEdge[2]));
    
  ScalarMul(tmpF + 3, m_DeformationGradient[1], sigmaEdge[1]);
  Add(tmpF + 3, tmpF + 3, ScalarMul(tmp, m_DeformationGradient[0], sigmaEdge[2]));

  for (int n = 0; n < m_NumPatchNodes; n++) {
    for (int c = 0; c < 2; c++) {
      Add(p_fU + 3*m_EdgePatchNodeIndices[n], p_fU + 3*m_EdgePatchNodeIndices[n], ScalarMul(tmp, tmpF + 3*c, m_PatchShapeGradients[n][c]*area*mpc_Owner->GetMaterial().GetThickness()));
      assert(!std::isnan(*(p_fU + 3*m_EdgePatchNodeIndices[n])));  
    }
  }
}

void tledElementShellBSTP1Edge::ComputeBendingForces(float *p_fU, const float M[], const float hs[][3]) const {
  using namespace tledVectorArithmetic;

  const float *linearGradient = mpc_Owner->GetLinearGradient(m_ElementEdgeIndex);
  const float area = mpc_Owner->GetArea();

  {
    float gradientRho[3], fGlob[3];

    for (int c = 0; c < 3; c++) {
      float rho[2];

      rho[0] = Dot(mpc_Owner->GetContravariantAxis(0), hs[c]);
      rho[1] = Dot(mpc_Owner->GetContravariantAxis(1), hs[c]);

      gradientRho[c] = ::Dot(linearGradient, 2, rho);
    }
    
    assert(!std::isnan(gradientRho[0] + gradientRho[1] + gradientRho[2]));
    
    Add(p_fU + 3*m_EdgePatchNodeIndices[2], p_fU + 3*m_EdgePatchNodeIndices[2], ScalarMul(fGlob, mpc_Owner->GetCurrentNormal(), -2*area*Dot(gradientRho, M)));
  }

  if (IsClamped()) {
    float tmpHM, f[3];

    tmpHM = 2*area/m_L0*(-M[0]*linearGradient[0]*m_EdgeNormal[1] + M[1]*linearGradient[1]*m_EdgeNormal[0] + M[2]*(linearGradient[0]*m_EdgeNormal[0] - linearGradient[1]*m_EdgeNormal[1]));
    Add(p_fU + 3*m_EdgePatchNodeIndices[0], p_fU + 3*m_EdgePatchNodeIndices[0], ScalarMul(f, mpc_Owner->GetCurrentNormal(), -tmpHM));
    Add(p_fU + 3*m_EdgePatchNodeIndices[1], p_fU + 3*m_EdgePatchNodeIndices[1], ScalarMul(f, mpc_Owner->GetCurrentNormal(), tmpHM));
  } else {
    float tmp[2];

    tmp[0] = 2*area*(linearGradient[0]*M[0] + linearGradient[1]*M[2]);
    tmp[1] = 2*area*(linearGradient[1]*M[1] + linearGradient[0]*M[2]);

    for (int n = 0; n < m_NumPatchNodes; n++) {
      float f, fGlob[3];

      f = ::Dot(m_PatchShapeGradients[n], 2, tmp);
      Add(p_fU + 3*m_EdgePatchNodeIndices[n], p_fU + 3*m_EdgePatchNodeIndices[n], ScalarMul(fGlob, mpc_Owner->GetCurrentNormal(), f));
      assert(!std::isnan(f));  
    }    
  } 
}

void tledElementShellBSTP1Edge::RegisterPrimaryTriangle(const tledElementShellBSTP1 &te, const int localEdgeIndex) {
  assert(localEdgeIndex < 3 && localEdgeIndex >= 0);

  m_NumPatchNodes += 1;  
  mpc_Owner = &te;
  m_ElementEdgeIndex = localEdgeIndex;

  m_EdgePatchNodeIndices[2] = te.GetFacetVertexIndices()[m_ElementEdgeIndex];
  m_EdgePatchNodeIndices[0] = mpc_Owner->GetFacetVertexIndices()[(m_ElementEdgeIndex+1)%3];
  m_EdgePatchNodeIndices[1] = mpc_Owner->GetFacetVertexIndices()[(m_ElementEdgeIndex+2)%3];

  assert(m_EdgePatchNodeIndices[0] != m_EdgePatchNodeIndices[2] && m_EdgePatchNodeIndices[1] != m_EdgePatchNodeIndices[2]);
  assert(m_EdgePatchNodeIndices[0] == te.GetFacetVertexIndices()[(localEdgeIndex+1)%3] || m_EdgePatchNodeIndices[0] == te.GetFacetVertexIndices()[(localEdgeIndex+2)%3]);
  assert(m_EdgePatchNodeIndices[1] == te.GetFacetVertexIndices()[(localEdgeIndex+1)%3] || m_EdgePatchNodeIndices[1] == te.GetFacetVertexIndices()[(localEdgeIndex+2)%3]);
}

void tledElementShellBSTP1Edge::RegisterSecondaryTriangle(const tledElementShellBSTP1 &te, const int localEdgeIndex) {
  assert(localEdgeIndex < 3 && localEdgeIndex >= 0);

  m_NumPatchNodes += 1;
  m_EdgePatchNodeIndices[3] = te.GetFacetVertexIndices()[localEdgeIndex];
  assert(m_EdgePatchNodeIndices[2] != m_EdgePatchNodeIndices[3] && m_EdgePatchNodeIndices[0] != m_EdgePatchNodeIndices[3] && m_EdgePatchNodeIndices[1] != m_EdgePatchNodeIndices[3]);
}

void tledElementShellBSTP1Edge::InitTensors(const float elementBasis[][3]) {
  using namespace tledVectorArithmetic;

  if (m_NumPatchNodes == 4) {
    float dx_dXi[2][3], J[2][2], JInv[2][2];

    for (int c = 0; c < 2; c++) {
      for (int xc = 0; xc < 3; xc++) {
	dx_dXi[c][xc] = 0;
	for (int n = 0; n < m_NumPatchNodes; n++) {
	  dx_dXi[c][xc] += mc_X0s[3*m_EdgePatchNodeIndices[n]+xc]*s_PatchShapeGradients[m_ElementEdgeIndex][n][c];
	}
      }
    }

    assert(Norm(dx_dXi[0]) > 0 && Norm(dx_dXi[1]) > 0);

    for (int r = 0; r < 2; r++) for (int c = 0; c < 2; c++) {
	J[r][c] = Dot(elementBasis[c], dx_dXi[r]);
	assert(!std::isnan(J[r][c]));
      }
    MatInverse22(&JInv[0][0], &J[0][0]);

    for (int n = 0; n < 4; n++) for (int c = 0; c < 2; c++) { 
	m_PatchShapeGradients[n][c] = JInv[c][0]*s_PatchShapeGradients[m_ElementEdgeIndex][n][0] + JInv[c][1]*s_PatchShapeGradients[m_ElementEdgeIndex][n][1];
	assert(!std::isnan(m_PatchShapeGradients[n][c]));
      }  
  } else {
    assert(m_EdgePatchNodeIndices[2] == mpc_Owner->GetFacetVertexIndices()[m_ElementEdgeIndex]);
    std::copy(mpc_Owner->GetLinearGradient(m_ElementEdgeIndex), mpc_Owner->GetLinearGradient(m_ElementEdgeIndex) + 2, m_PatchShapeGradients[2]);

    for (int n = 0; n < 2; n++) {
      const int oVInd = (m_ElementEdgeIndex + n + 1)%3;

      assert(m_EdgePatchNodeIndices[n] == mpc_Owner->GetFacetVertexIndices()[oVInd]);
      std::copy(mpc_Owner->GetLinearGradient(oVInd), mpc_Owner->GetLinearGradient(oVInd) + 2, m_PatchShapeGradients[n]);
    }
  }
}

void tledElementShellBSTP1Edge::InitPatchTopologyAndGeometry(const tledSurfaceTopology<Surface> &surfTopo, const int edgeIndex) {
  assert(surfTopo.GetEdges()[edgeIndex].first < surfTopo.GetEdges()[edgeIndex].second);
  mc_X0s = surfTopo.GetSurface().GetAllNodeCoordinates();
}

tledElementShellBSTP1Edge::tledElementShellBSTP1Edge(void) {
  m_NumPatchNodes = 2;
  m_L0 = std::numeric_limits<float>::quiet_NaN();

#ifndef NDEBUG
  mc_X0s = NULL;
  m_ElementEdgeIndex = -1000;
  std::fill(m_EdgePatchNodeIndices, m_EdgePatchNodeIndices + 4, -1000);
  std::fill(&m_DeformationGradient[0][0], &m_DeformationGradient[0][0] + 2*3, std::numeric_limits<float>::quiet_NaN());
  std::fill(&m_PatchShapeGradients[0][0], &m_PatchShapeGradients[0][0] + 2*4, std::numeric_limits<float>::quiet_NaN());
  std::fill(m_EdgeNormal, m_EdgeNormal + 2, std::numeric_limits<float>::quiet_NaN());
#endif
}
