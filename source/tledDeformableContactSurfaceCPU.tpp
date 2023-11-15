// =========================================================================
// File:       tledDeformableContactSurfaceCPU.tpp
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
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetNormalisedOldFacetNormalCached(const int facetInd) {
  using namespace tledVectorArithmetic;

  CachedNormal &r_dst = this->GetOldFacetNormalCache(facetInd);

  assert(facetInd >= 0 && facetInd < this->GetNumberOfFacets());
  if (this->GetSaveCount() > r_dst.UpdateCounter) {
    ScalarDiv(r_dst.Normal, Norm(this->ComputeFacetNormal(r_dst.Normal, this->GetFacet(facetInd).NodeIndices, this->GetAllOldNodeCoordinates())));
    r_dst.UpdateCounter = this->GetSaveCount();
  }

  return r_dst.Normal;
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetNormalisedOldFacetNormal(const int facetInd) const {
  assert(GetOldFacetNormalCache(facetInd).UpdateCounter == this->GetSaveCount());

  return this->GetOldFacetNormalCache(facetInd).Normal;
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetFacetNormalCached(const int facetInd) {
  CachedNormal &r_dst = this->GetFacetNormalCache(facetInd);

  if (r_dst.UpdateCounter < this->GetUpdateCount()) {
    this->ComputeFacetNormal(r_dst.Normal, facetInd);
    r_dst.UpdateCounter = r_dst.UpdateCounter;
  } 

  return r_dst.Normal;
}

template <class TBaseSurface>
float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::ComputeNormalisedFacetNormalCached(float *p_n, const int facetInd) {
  using namespace tledVectorArithmetic;

  return ScalarDiv(p_n, GetFacetNormalCached(facetInd), Norm(GetFacetNormalCached(facetInd)));
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetFacetProjectionOperatorCached(const int facetIndex) {
  CachedProjectionOperator &r_dst = this->GetProjectionOperatorCache(facetIndex);
  
  if (this->GetUpdateCount() > r_dst.UpdateCounter) {
    Superclass::ComputeFacetProjectionOperator(r_dst.Operator, facetIndex);
    r_dst.UpdateCounter = this->GetUpdateCount();
  }

  return r_dst.Operator;
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetFacetProjectionOperator(const int facetIndex) const {
  assert(this->GetProjectionOperatorCache(facetIndex).UpdateCounter == this->GetUpdateCount());

  return this->GetProjectionOperatorCache(facetIndex).Operator;
}

template <class TBaseSurface>
float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::ComputeC0NormalCached(float *p_dstNormal, const int facet[], const float shapeValues[]) {  
  using namespace tledVectorArithmetic;

  float tmp[3];
  int tInd;

  ScalarMul(p_dstNormal, this->GetNodeNormalCached(facet[0]), shapeValues[0]);
  for (tInd = 1; tInd < Facet::NumberOfVertices; tInd++) Add(p_dstNormal, p_dstNormal, ScalarMul(tmp, this->GetNodeNormalCached(facet[tInd]), shapeValues[tInd]));
  
  return ScalarDiv(p_dstNormal, Norm(p_dstNormal));
} 

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::Update(const float uNexts[]) {
  float* p_nPos;
  std::vector<float>::const_iterator ic_nPos0;
  std::vector<int>::const_iterator ic_gnInd;
  
  assert(this->GetNumberOfNodes() == (int)this->GetSurface2VolumeNodeMap().size());
  assert(this->GetAllNodeCoordinates0().size() == 3*this->GetSurface2VolumeNodeMap().size());
  for (ic_gnInd = this->GetSurface2VolumeNodeMap().begin(), ic_nPos0 = this->GetAllNodeCoordinates0().begin(), p_nPos = this->GetAllNodeCoordinates(); ic_gnInd < this->GetSurface2VolumeNodeMap().end(); ic_gnInd++, ic_nPos0 += 3, p_nPos += 3) {
    tledVectorArithmetic::Add(p_nPos, &(*ic_nPos0), uNexts + 3*(*ic_gnInd));
  }
  
  Superclass::Update(uNexts);
} /* Update */

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::Init() {  
  Superclass::Init();

  m_FacetData.resize(this->GetNumberOfFacets());
  this->SetAllOldNodeCoordinates(this->GetAllNodeCoordinates());
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetNodeNormal(const int nInd) const {
  assert(this->GetNodeNormalCache(nInd).UpdateCounter == this->GetUpdateCount());
  return this->GetNodeNormalCache(nInd).Normal;
}

template <class TBaseSurface>
const float* tledDeformableContactSurfaceImplCPU<TBaseSurface>::GetNodeNormalCached(const int nInd) {
  using namespace tledVectorArithmetic;

  CachedNormal &r_n = this->GetNodeNormalCache(nInd);  
  float *p_n = r_n.Normal;

  if (r_n.UpdateCounter < this->GetUpdateCount()) {
    const int *facetIndicesBegin = this->GetNodeFacetIndices(nInd);
    const int *facetIndicesEnd = facetIndicesBegin + this->GetNumberOfNodeFacets(nInd);

    float tmpNorm[3];

    std::fill(tmpNorm, tmpNorm + 3, 0.0f);
    for (int const *pc_fInd = facetIndicesBegin; pc_fInd < facetIndicesEnd; pc_fInd++) {
      const int numFacetVertices = Facet::NumberOfVertices;
      const int *facet = this->GetFacet(*pc_fInd).NodeIndices;    
      const int lNodeInd = std::find(facet, facet + numFacetVertices, nInd) - facet;

      float x[3], y[3], sn[3], angle;        

      assert(lNodeInd < numFacetVertices && facet[lNodeInd] == nInd);
      Sub(x, this->GetNodeCoordinates(facet[(lNodeInd+1)%numFacetVertices]), this->GetNodeCoordinates(facet[lNodeInd]));
      Sub(y, this->GetNodeCoordinates(facet[(lNodeInd-1+numFacetVertices)%numFacetVertices]), this->GetNodeCoordinates(facet[lNodeInd]));

      angle = ComputeAngleFast(x, y);
      assert(angle != angle || angle + 1 == angle || (angle > 0 && angle < tledPi));

      Add(tmpNorm, tmpNorm, ScalarMul(sn, this->GetFacetNormalCached(*pc_fInd), angle));	
    }

    ScalarDiv(p_n, tmpNorm, Norm(tmpNorm));    
    r_n.UpdateCounter = this->GetUpdateCount();
  }

  return p_n;
} 

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::SetNumberOfFacets(const int numFacets) {
  Superclass::SetNumberOfFacets(numFacets);
  m_FacetData.resize(numFacets);
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::SetAllOldNodeCoordinates(const float *nodeCds) {
  if ((int)m_OldNodeCoordinates.size() != this->GetCoordinateHistorySize()*this->GetNumberOfNodes()) m_OldNodeCoordinates.resize(this->GetCoordinateHistorySize()*3*this->GetNumberOfNodes());
  std::copy(nodeCds, nodeCds + this->GetNumberOfNodes()*3, m_OldNodeCoordinates.begin());
#ifndef NDEBUG
  std::fill(m_OldNodeCoordinates.begin() + this->GetNumberOfNodes()*3, m_OldNodeCoordinates.end(), std::numeric_limits<float>::quiet_NaN());
#endif
  mpc_OldNodeCoordinates = &m_OldNodeCoordinates.front();
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::Save() {
  assert(3*this->GetNumberOfNodes()*this->GetCoordinateHistorySize() == (int)m_OldNodeCoordinates.size());  
  assert(this->GetCoordinateHistorySize() > 0);

  Superclass::Save();
  std::copy(this->GetAllNodeCoordinates(), this->GetAllNodeCoordinates() + 3*this->GetNumberOfNodes(), m_OldNodeCoordinates.begin() + (this->GetSaveCount()%this->GetCoordinateHistorySize())*3*this->GetNumberOfNodes());  
  if (this->GetSaveCount() >= this->GetCoordinateHistorySize()) {
    mpc_OldNodeCoordinates = &m_OldNodeCoordinates.front() + ((this->GetSaveCount() + 1)%this->GetCoordinateHistorySize())*3*this->GetNumberOfNodes();
  }
}

template <class TBaseSurface>
void tledDeformableContactSurfaceImplCPU<TBaseSurface>::LoadFromXMLPostloadHook(void) {
  this->SetAllOldNodeCoordinates(this->GetAllNodeCoordinates());
}
