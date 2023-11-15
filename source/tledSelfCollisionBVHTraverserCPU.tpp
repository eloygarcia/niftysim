// =========================================================================
// File:       tledSelfCollisionBVHTraverserCPU.tpp
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

template <class TBVH, class TAPI>
void tledSelfCollisionBVHTraverserImplCPU<TBVH, TAPI>::AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd) {
  typedef typename SlaveMesh::Facet __SlaveFacet;
  typedef typename MasterMesh::Facet __MasterFacet;

  const __SlaveFacet &sFacet = this->GetSlaveMesh().GetFacet(sFacetInd);
  const __MasterFacet &mFacet = this->GetMasterMesh().GetFacet(mFacetInd);

  for (int const *pc_sn = sFacet.NodeIndices; pc_sn < sFacet.NodeIndices + __SlaveFacet::NumberOfVertices; pc_sn++) {
    for (int const *pc_mn = mFacet.NodeIndices; pc_mn < mFacet.NodeIndices + __MasterFacet::NumberOfVertices; pc_mn++) {
      if (*pc_mn == *pc_sn) return;
    }
  }

  Superclass::AddNarrowPhaseTests(sFacetInd, mFacetInd);
}
