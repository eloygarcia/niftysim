// =========================================================================
// File:       tledDynamicBVH.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::TransformSubTree(const int rootBVInd, const float m[][3], const float cor[], const float t[]) {
  BoundingVolume &r_bv = this->GetBV(rootBVInd);

  assert(BoundingVolume::CanRotate);
  r_bv.Translate(t);
  r_bv.Rotate(m, cor);  

  if (r_bv.PrimitiveIndex < 0) {
    for (int const *pc_childInd = r_bv.ChildIndices; pc_childInd < r_bv.ChildIndices + Superclass::BVHOrder; pc_childInd++) if (BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
	this->TransformSubTree(*pc_childInd, m, cor, t);
      }
  }
}

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::TranslateSubTree(const int rootBVInd, const float t[]) {
  BoundingVolume &r_bv = this->GetBV(rootBVInd);

  r_bv.Translate(t);
  if (r_bv.PrimitiveIndex < 0) {
    for (int const *pc_childInd = r_bv.ChildIndices; pc_childInd < r_bv.ChildIndices + Superclass::BVHOrder; pc_childInd++) if (BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
	this->TranslateSubTree(*pc_childInd, t);      
      }
  }
}

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::Init(tledBVHCreator &r_bvhBuilder) {
  Superclass::Init(r_bvhBuilder);
  this->GetUpdater().Init();
}

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::Update() {
  m_UpdateCounter += 1;
  tledLogDebugStream(tledHelper::Info() << "Performing BVH-update " << m_UpdateCounter);
  this->GetUpdater().UpdateBVH();
}

#ifndef NDEBUG
template <class TBVH>
static std::vector<int> _SubtreeNodeList(const TBVH &bvh, const int bv) {
  if (bvh.IsLeaf(bv)) {
    const int *facet = bvh.GetMesh().GetFacet(bvh.GetBV(bv).PrimitiveIndex).NodeIndices; 

    return std::vector<int>(facet, facet + TBVH::ContactMesh::Facet::NumberOfVertices);    
  } else {
    std::vector<int> nl = _SubtreeNodeList(bvh, bvh.GetBV(bv).ChildIndices[0]);

    for (int const *pc_c = bvh.GetBV(bv).ChildIndices + 1; pc_c < bvh.GetBV(bv).ChildIndices + TBVH::BoundingVolume::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
	std::vector<int> tmp = _SubtreeNodeList(bvh, *pc_c);

	nl.insert(nl.end(), tmp.begin(), tmp.end());
      }

    return nl;
  }
}

template <class TBVH>
static bool _CheckContainment(const TBVH &bvh, const int bvInd) {
  const typename TBVH::BoundingVolume &bv = bvh.GetBV(bvInd);

  std::vector<int> nodes = _SubtreeNodeList(bvh, bvInd);

  for (std::vector<int>::const_iterator ic_n = nodes.begin(); ic_n < nodes.end(); ic_n++) {
    for (int c = 0; c < 3; c++) {
      if (bvh.GetMesh().GetNodeCoordinates(*ic_n)[c] < bv.Bounds[c][0]
	  || bvh.GetMesh().GetNodeCoordinates(*ic_n)[c] > bv.Bounds[c][1]
	  || bvh.GetMesh().GetOldNodeCoordinates(*ic_n)[c] < bv.Bounds[c][0]
	  || bvh.GetMesh().GetOldNodeCoordinates(*ic_n)[c] > bv.Bounds[c][1]) {
	return false;
      }
    }
  }
  
  return true;
}
#endif

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::UpdateTopDownRecursive(const int bvInd) {
  const BoundingVolume &bv = this->GetBV(bvInd);  

  assert(bvInd >= 0 && bvInd < this->GetNumberOfBVs());

  if (bv.PrimitiveIndex < 0) {
    for (int const *pc_childInd = bv.ChildIndices; pc_childInd < bv.ChildIndices + TBV::NumberOfChildBVs; pc_childInd++) if (BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
	this->UpdateTopDownRecursive(*pc_childInd);
      }

    this->RefitInteriorBV(bvInd);
  } else this->RefitLeafBV(bvInd);
} /* UpdateTopDownRecursive */

template <class TSurface, class TBV, class TAPI>
void tledDynamicBVHImpl<TSurface, TBV, TAPI>::ComputeBoundsFromNodes(BoundingVolume &r_bv, const int *nodeListStart, const int *nodeListEnd) const {
  r_bv.ComputeFromNodeList(nodeListStart, nodeListEnd, this->GetMesh().GetAllNodeCoordinates());
  r_bv.ExpandWithNodeList(nodeListStart, nodeListEnd, this->GetMesh().GetAllOldNodeCoordinates());
  r_bv.AddMargin(this->GetMargin());
}
