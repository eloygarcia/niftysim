// =========================================================================
// File:       tledSurface.tpp
// Purpose:    General surface mesh representation (implementation)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TFacet, class TInterface>
float* tledSurfaceImpl<TFacet, TInterface>::ComputeFacetNormal(float *p_normalDst, const int facet[], const float nodes[]) {
  if (Facet::NumberOfVertices == 3) {    
    float x[3], y[3];

    tledSurfaceImpl::ComputeFacetBasis(x, y, p_normalDst, facet, nodes);
  } else {
    using namespace tledVectorArithmetic;

    float ac[3], bd[3];

    /* Formula from old pressure constraint code, keep for consistency! */
    Sub(ac, nodes + 3*facet[2], nodes + 3*facet[0]);
    Sub(bd, nodes + 3*facet[3], nodes + 3*facet[1]);
    Cross(p_normalDst, ac, bd);
  }

  return p_normalDst;
}

template <class TFacet, class TInterface>
void tledSurfaceImpl<TFacet, TInterface>::ComputeFacetBasis(float *p_dxDXi, float *p_dxDEta, float *p_normalDst, const int facet[], const float nodes[]) {
  using namespace tledVectorArithmetic;

  if (Facet::NumberOfVertices == 3) {    
    Sub(p_dxDXi, nodes + 3*facet[1], nodes + 3*facet[0]), Sub(p_dxDEta, nodes + 3*facet[2], nodes + 3*facet[0]);
  } else if (Facet::NumberOfVertices == 4) {
    ScalarMul(Sub(p_dxDXi, Sub(p_dxDXi, Add(p_dxDXi, nodes + 3*facet[1], nodes + 3*facet[2]), nodes + 3*facet[0]), nodes + 3*facet[3]), 0.5);
    ScalarMul(Sub(p_dxDEta, Sub(p_dxDEta, Add(p_dxDEta, nodes + 3*facet[2], nodes + 3*facet[3]), nodes + 3*facet[0]), nodes + 3*facet[1]), 0.5);
  } else {
    tledFatalError("Unsupported facet type.");
  }

  Cross(p_normalDst, p_dxDXi, p_dxDEta);
}

template <class TFacet, class TInterface>
float* tledSurfaceImpl<TFacet, TInterface>::ComputeCentroid(float *p_centroid, const int facetInd) const {
  const int *facetVtxInds = this->GetFacet(facetInd).NodeIndices;

  std::copy(this->GetNodeCoordinates(*facetVtxInds), this->GetNodeCoordinates(*facetVtxInds) + 3, p_centroid);
  for (int const *pc_nInd = facetVtxInds + 1; pc_nInd < facetVtxInds + Facet::NumberOfVertices; pc_nInd++) {
    tledVectorArithmetic::Add(p_centroid, p_centroid, this->GetNodeCoordinates(*pc_nInd));
  }

  return tledVectorArithmetic::ScalarDiv(p_centroid, (float)Facet::NumberOfVertices);
}

template <class TFacet, class TInterface>
float* tledSurfaceImpl<TFacet, TInterface>::ComputeNormalisedFacetNormal(float *p_normalDst, const int facetInd) const {
  return tledVectorArithmetic::ScalarDiv(p_normalDst, tledVectorArithmetic::Norm(ComputeFacetNormal(p_normalDst, facetInd)));
}

template <class TFacet, class TInterface>
void tledSurfaceImpl<TFacet, TInterface>::SetNodeVector(const float nodes[], const int numNodes) {
  mpc_NodeCoordinates = nodes;
  m_NumNodes = numNodes;
}

template <class TFacet, class TInterface>
template <typename TConstIterator>
std::vector<int> tledSurfaceImpl<TFacet, TInterface>::CompileNodeListFromFacetList(const TConstIterator pIndsBegin, const TConstIterator pIndsEnd) const {
  std::vector<int> surfNodeInds;  

  surfNodeInds.reserve(Facet::NumberOfVertices*(pIndsEnd - pIndsBegin));
  for (TConstIterator ic_primInd = pIndsBegin; ic_primInd < pIndsEnd; ic_primInd++) {
    const int *facet = GetFacet(*ic_primInd).NodeIndices;

    assert(*ic_primInd >= 0 && *ic_primInd < GetNumberOfFacets());
    surfNodeInds.insert(surfNodeInds.end(), facet, facet + Facet::NumberOfVertices);
  }
  surfNodeInds = tledHelper::MakeSortedUnique(surfNodeInds);

  return surfNodeInds;
}

template <class TFacet, class TInterface>
float tledSurfaceImpl<TFacet, TInterface>::ComputeFacetArea(const int facetInd) const {
  float n[3];

  return tledVectorArithmetic::Norm(this->ComputeFacetNormal(n, facetInd))/2;
}

template <class TFacet, class TInterface>
void tledSurfaceImpl<TFacet, TInterface>::SetNumberOfFacets(const int numFacets) { 
  m_Facets.resize(numFacets); 
  TInterface::SetNumberOfFacets(numFacets);
}
