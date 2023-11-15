// =========================================================================
// File:       tledContactSurface.tpp
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
  
template <const int t_numFacetVertices, class TInterface>
void tledContactSurfaceImpl<t_numFacetVertices, TInterface>::SetNumberOfNodes(const int numNodes) {
  m_NodeCoordinates.resize(numNodes*3);
  this->SetNodeVector(&m_NodeCoordinates.front(), numNodes);
  TInterface::SetNumberOfNodes(numNodes);
}

template <const int t_numFacetVertices, class TInterface>
void tledContactSurfaceImpl<t_numFacetVertices, TInterface>::ComputeDiameters() {
  using namespace tledVectorArithmetic;

  this->m_MaxH = -(this->m_MinH = std::numeric_limits<float>::max());
  for (int fInd = 0; fInd < this->GetNumberOfFacets(); fInd++) {
    const Facet &facet = this->GetFacet(fInd);
    const float area = this->ComputeFacetArea(fInd);

    float minH, maxH, currH, edgeV[3];

    maxH = -(minH = std::numeric_limits<float>::max());
    for (int eInd = 0; eInd < Facet::NumberOfVertices; eInd++) {
      const std::pair<int, int> &edge = this->GetEdge(facet.EdgeIndices[eInd]);

      if ((currH = area/Norm(Sub(edgeV, this->GetNodeCoordinates(edge.first), this->GetNodeCoordinates(edge.second)))) < minH) minH = currH;
      maxH = std::max(currH, maxH);
    }

    this->m_MinH = std::min(minH, this->m_MinH);
    this->m_MaxH = std::max(maxH, this->m_MaxH);
  }
}

template <const int t_numFacetVertices, class TInterface>
float* tledContactSurfaceImpl<t_numFacetVertices, TInterface>::ComputeFacetProjectionOperator(float *p_R, const int facetInd) const {
  using namespace tledVectorArithmetic;

  const int *facet = this->GetFacet(facetInd).NodeIndices;

  if (t_numFacetVertices == 3) {
    const float *v0 = this->GetNodeCoordinates(facet[0]);
    const float *v1 = this->GetNodeCoordinates(facet[1]);
    const float *v2 = this->GetNodeCoordinates(facet[2]);

    float a[3], b[3], n[3], in_det;

    /*
     * Moller-Trumbore facet projection operator, modified for normalised normals
     */
    Sub(a, v1, v0);
    Sub(b, v2, v0);
    Cross(n, a, b);
  
    in_det = 1/Dot(n, n);
    p_R[0] = (b[1]*n[2] - b[2]*n[1])*in_det;
    p_R[1] = (b[2]*n[0] - b[0]*n[2])*in_det;
    p_R[2] = (b[0]*n[1] - b[1]*n[0])*in_det;
    p_R[4] = (n[1]*a[2] - n[2]*a[1])*in_det;
    p_R[5] = (n[2]*a[0] - n[0]*a[2])*in_det;
    p_R[6] = (n[0]*a[1] - n[1]*a[0])*in_det;
    in_det = std::sqrt(in_det);
    p_R[8] = (n[0])*in_det;
    p_R[9] = (n[1])*in_det;
    p_R[10] = (n[2])*in_det;
  
    p_R[3] = -(p_R[0]*v0[0] + p_R[1]*v0[1] + p_R[2]*v0[2]);
    p_R[7] = -(p_R[4]*v0[0] + p_R[5]*v0[1] + p_R[6]*v0[2]);
    p_R[11] = -(p_R[8]*v0[0] + p_R[9]*v0[1] + p_R[10]*v0[2]);
  } else {
    tledFatalNotYetImplementedError;
  }

  return p_R;
}
