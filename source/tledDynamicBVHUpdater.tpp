// =========================================================================
// File:       tledDynamicBVHUpdater.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef NDEBUG
template <class TBVH>
void tledDynamicBVHUpdaterImpl<TBVH>::CheckContainmentRecursive(std::vector<int> &r_nodeInds, const int bvInd, const float trans[]) const {
  std::vector<int> allNodeInds;

  if (this->GetBVH().IsLeaf(bvInd)) {
    const Facet &facet = this->GetMesh().GetFacet(this->GetBVH().GetBV(bvInd).PrimitiveIndex);

    allNodeInds.insert(allNodeInds.end(), facet.NodeIndices, facet.NodeIndices + Facet::NumberOfVertices);
  } else {
    const int *children = this->GetBVH().GetBV(bvInd).ChildIndices;

    for (int const *pc_c = children; pc_c < children + BoundingVolume::NumberOfChildBVs; pc_c++) {
      this->CheckContainmentRecursive(allNodeInds, *pc_c, trans);
    }
    allNodeInds = tledHelper::MakeSortedUnique(allNodeInds);
  }
  
  for (std::vector<int>::const_iterator ic_n = allNodeInds.begin(); ic_n < allNodeInds.end(); ic_n++) {
    if (trans == NULL) {
      assert(this->GetBVH().IsPointInside(this->GetMesh().GetNodeCoordinates(*ic_n), bvInd));
    } else {
      using namespace tledVectorArithmetic;

      float tmp[3];

      assert(this->GetBVH().IsPointInside(Add(tmp, this->GetMesh().GetNodeCoordinates(*ic_n), trans), bvInd));
    }
  }

  r_nodeInds.insert(r_nodeInds.end(), allNodeInds.begin(), allNodeInds.end());
}

template <class TBVH>
void tledDynamicBVHUpdaterImpl<TBVH>::CheckContainmentRecursive(const int bvInd, const float trans[]) const {
  std::vector<int> nodeInds;

  this->CheckContainmentRecursive(nodeInds, bvInd, trans);
}

template <class TBVH>
void tledDynamicBVHUpdaterImpl<TBVH>::CheckContainmentRecursive(const int bvInd) const {
  float trans[] = {0, 0, 0};  

  this->CheckContainmentRecursive(bvInd, trans);
}

#endif
