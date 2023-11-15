// =========================================================================
// File:       tledRigidContactSurfaceCPU.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBaseSurface>
float* tledRigidContactSurfaceImplCPU<TBaseSurface>::ComputeC0Normal(float *p_dstNormal, const int facet[], const float shapeValues[]) const {  
  using namespace tledVectorArithmetic;

  float tmp[3];
  int tInd;

  ScalarMul(p_dstNormal, GetNodeNormal(facet[0]), shapeValues[0]);
  for (tInd = 1; tInd < Facet::NumberOfVertices; tInd++) Add(p_dstNormal, p_dstNormal, ScalarMul(tmp, GetNodeNormal(facet[tInd]), shapeValues[tInd]));
  
  return ScalarDiv(p_dstNormal, Norm(p_dstNormal));
} 

template <class TBaseSurface>
void tledRigidContactSurfaceImplCPU<TBaseSurface>::SetNumberOfFacets(const int numFacets) {
  Superclass::SetNumberOfFacets(numFacets);
  m_FacetData.resize(this->GetNumberOfFacets());  
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplCPU<TBaseSurface>::InitNormals() {
  assert(this->GetAllFacets().size() == m_FacetData.size());
  for (int facetInd = 0; facetInd < this->GetNumberOfFacets(); facetInd++) this->ComputeNormalisedFacetNormal(m_FacetData[facetInd].Normal, facetInd);

  assert((int)this->m_NodeNormals.size() == 3*this->GetNumberOfNodes());
  std::fill(this->m_NodeNormals.begin(), this->m_NodeNormals.end(), 0.0f);
  this->ComputeNodeNormals(&m_NodeNormals.front());
   
  for (int facetInd = 0; facetInd < this->GetNumberOfFacets(); facetInd++) {
    assert(std::fabs(tledVectorArithmetic::Norm(m_FacetData[facetInd].Normal) - 1) < 1e-3f);
    this->ComputeFacetProjectionOperator(m_FacetData[facetInd].ProjectionOperator, facetInd);
  }
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplCPU<TBaseSurface>::SetAllNodeNormals(const float normals[]) {
  m_NodeNormals.resize(3*this->GetNumberOfNodes());
  std::copy(normals, normals + 3*this->GetNumberOfNodes(), m_NodeNormals.begin());
}

template <class TBaseSurface>
void tledRigidContactSurfaceImplCPU<TBaseSurface>::SetNumberOfNodes(const int numNodes) {
  Superclass::SetNumberOfNodes(numNodes);
  m_NodeNormals.resize(3*numNodes);
#ifndef NDEBUG
  std::fill(m_NodeNormals.begin(), m_NodeNormals.end(), std::numeric_limits<float>::quiet_NaN());
#endif
}
