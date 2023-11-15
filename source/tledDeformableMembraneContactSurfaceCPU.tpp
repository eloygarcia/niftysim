// =========================================================================
// File:       tledDeformableMembraneContactSurfaceCPU.tpp
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
const float* tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::GetFacetNormalCached(const int facetIndex) {
  if (facetIndex < this->GetMembraneFacetBaseIndex() + this->GetNumberOfMembraneElements()) return Superclass::GetFacetNormalCached(facetIndex);    
  else {
    typename Superclass::CachedNormal &r_nData = this->GetFacetNormalCache(facetIndex);    
    
    if (r_nData.UpdateCounter < this->GetUpdateCount()) {
      const float *n = Superclass::GetFacetNormalCached(facetIndex - this->GetNumberOfMembraneElements());

      tledVectorArithmetic::ScalarMul(r_nData.Normal, n, -1);
      r_nData.UpdateCounter = this->GetUpdateCount();
    }

    return r_nData.Normal;
  }
}

template <class TBaseSurface>
const float* tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::GetNodeNormalCached(const int nodeIndex) {
  if (nodeIndex < this->GetMembraneNodeBaseIndex() + this->GetNumberOfMembraneNodes()) return Superclass::GetNodeNormalCached(nodeIndex);    
  else {
    typename Superclass::CachedNormal &r_n = this->GetNodeNormalCache(nodeIndex);    
    
    if (r_n.UpdateCounter < this->GetUpdateCount()) {
      const float *n = Superclass::GetNodeNormalCached(nodeIndex - this->GetNumberOfMembraneNodes());

      tledVectorArithmetic::ScalarMul(r_n.Normal, n, -1);
      r_n.UpdateCounter = this->GetUpdateCount();
    }

    return r_n.Normal;
  }
}

template <class TBaseSurface>
void tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::Update(const float u[]) {
  Superclass::Update(u);
  this->OffsetNodes();
}

template <class TBaseSurface>
void tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::UpdateMembraneFacetNormal(const int membraneFacetIndex, const float n[]) {
  {
    typename Superclass::CachedNormal &r_nData = this->GetFacetNormalCache(membraneFacetIndex + this->GetMembraneFacetBaseIndex());

    r_nData.UpdateCounter = this->GetUpdateCount() + 1;
    std::copy(n, n + 3, r_nData.Normal);
  }

  {
    typename Superclass::CachedNormal &r_nData = this->GetFacetNormalCache(membraneFacetIndex + this->GetMembraneFacetBaseIndex() + this->GetNumberOfMembraneElements());

    r_nData.UpdateCounter = this->GetUpdateCount() + 1;
    tledVectorArithmetic::ScalarMul(r_nData.Normal, n, -1);
  }
}

template <class TBaseSurface>
void tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::OffsetNodes() {
  using namespace tledVectorArithmetic;

  for (int fInd = this->GetMembraneFacetBaseIndex(); fInd < this->GetMembraneFacetBaseIndex() + this->GetNumberOfMembraneElements(); fInd++) {
    const float off = this->GetFacetThickness(fInd)/2;
    const float *un = this->GetFacetNormalCached(fInd);

    float d[3];
    
    ScalarMul(d, un, off/Norm(un));
    for (int const *pc_n = this->GetFacet(fInd).NodeIndices; pc_n < this->GetFacet(fInd).NodeIndices + Facet::NumberOfVertices; pc_n++) {
      float tmpD[3];

      assert(*pc_n >= this->GetMembraneNodeBaseIndex() && *pc_n < this->GetMembraneNodeBaseIndex() + this->GetNumberOfMembraneNodes());
      ScalarDiv(tmpD, d, (float)this->GetNumberOfNodeFacets(*pc_n));
      Add(this->GetAllNodeCoordinates() + 3*(*pc_n), this->GetNodeCoordinates(*pc_n), tmpD);
      Sub(this->GetAllNodeCoordinates() + 3*(*pc_n + this->GetNumberOfMembraneNodes()), this->GetNodeCoordinates(*pc_n + this->GetNumberOfMembraneNodes()), tmpD);
    }
  }
}

template <class TBaseSurface>
float tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::GetFacetThickness(const int facetIndex) const {
  if (facetIndex < this->GetMembraneFacetBaseIndex()) return 0;
  else if (facetIndex - this->GetMembraneFacetBaseIndex() < this->GetNumberOfMembraneElements()) return this->GetMembraneElementThickness(facetIndex - this->GetMembraneFacetBaseIndex());
  else {
    assert(facetIndex < this->GetNumberOfFacets());

    return this->GetMembraneElementThickness(facetIndex - this->GetMembraneFacetBaseIndex() - this->GetNumberOfMembraneElements());
  }
}

template <class TBaseSurface>
void tledDeformableMembraneContactSurfaceImplCPU<TBaseSurface>::Init() {
  Superclass::Init();
  this->OffsetNodes();
  this->SetAllOldNodeCoordinates(this->GetAllNodeCoordinates());
}
