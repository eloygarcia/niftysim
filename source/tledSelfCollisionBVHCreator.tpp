// =========================================================================
// File:       tledSelfCollisionBVHCreator.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    August 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::BVPrimitiveSet::Insert(const Facet &facet) {
  m_NodeSet.insert(facet.NodeIndices, facet.NodeIndices + Facet::NumberOfVertices);
}

template <class TBVH>
bool tledSelfCollisionBVHCreator<TBVH>::BVPrimitiveSet::IsAdjacent(const Facet &facet) const {
  for (int const *pc_nInd = facet.NodeIndices; pc_nInd < facet.NodeIndices + Facet::NumberOfVertices; pc_nInd++) {
    if (m_NodeSet.find(*pc_nInd) != m_NodeSet.end()) return true;
  }

  return false;
}

/********************************************************************************************************************************************/

template <class TBVH>
int tledSelfCollisionBVHCreator<TBVH>::_ReinsertBVsRecursive(const std::vector<BoundingVolume> &src, const int bvInd, const int parentInd) {
  const int bvhOrder = tledSelfCollisionBVHCreator::BVHOrder;
  const int dstBVInd = this->GetOutput().GetBVs().size();

  assert(bvInd < (int)src.size() && bvInd >= 0);
  assert(bvInd < (int)(src.size() - 1) || dstBVInd == 0);
  assert(this->GetOutput().GetUpdateCounter() == 0);
  this->GetOutput().GetBVs().push_back(src[bvInd]);
  this->GetOutput().GetBVs().back().ParentIndex = parentInd;
#ifndef NDEBUG
  this->GetOutput().GetBVs().back().UpdateCounter = 0;
#endif

  if (this->GetOutput().GetBVs().back().ChildIndices[0] >= 0) {
    assert(this->GetOutput().GetBVs().back().PrimitiveIndex == -1);
    for (int childInd = 0; childInd < bvhOrder; childInd++) {
      if (bvhOrder == 2 || this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd] >= 0) {    
	assert(this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd] >= 0 && this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd] < (int)src.size());
	this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd] = _ReinsertBVsRecursive(src, this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd], dstBVInd);
      } else {
	this->GetOutput().GetBVs()[dstBVInd].ChildIndices[childInd] = -1;
      }
    }
  } else {
    assert(this->GetOutput().GetBVs().back().PrimitiveIndex < this->GetMesh().GetNumberOfFacets());
#ifndef NDEBUG
    for (int childInd = 0; childInd < bvhOrder; childInd++) assert(this->GetOutput().GetBVs().back().ChildIndices[childInd] == -1);
#endif
    assert(dstBVInd >= 0 && dstBVInd < (int)this->GetOutput().GetBVs().size());
    assert(this->GetOutput().GetBVs()[dstBVInd].PrimitiveIndex >= 0 && this->GetOutput().GetBVs()[dstBVInd].PrimitiveIndex < this->GetMesh().GetNumberOfFacets());
    this->GetOutput().GetLeafBVIndices().push_back(dstBVInd);
  }

  return dstBVInd;
}

template <class TBVH>
std::vector<int> tledSelfCollisionBVHCreator<TBVH>::FindClusters(std::vector<int>::iterator i_pIndsBegin, std::vector<int>::iterator i_pIndsEnd) const {
  std::vector<int> newPIndList;
  std::vector<int>::iterator i_currPInd, i_nextPInd;
  std::vector<int> clusterStartInds;
  BVPrimitiveSet clusterSet;  
  
  newPIndList.push_back(*i_pIndsBegin);
  clusterStartInds.push_back(0);
  clusterSet.Insert(this->GetMesh().GetFacet(*i_pIndsBegin));
  for (i_currPInd = i_pIndsBegin + 1; i_currPInd < i_pIndsEnd; i_currPInd++) {
    assert(clusterStartInds.back() < i_pIndsEnd - i_pIndsBegin);
    for (i_nextPInd = i_currPInd; i_nextPInd < i_pIndsEnd && !clusterSet.IsAdjacent(this->GetMesh().GetFacet(*i_nextPInd)); i_nextPInd++);

    if (i_nextPInd < i_pIndsEnd) {
      std::iter_swap(i_currPInd, i_nextPInd);
      assert(*i_currPInd >= 0 && *i_currPInd < this->GetMesh().GetNumberOfFacets());
      assert(*i_nextPInd >= 0 && *i_nextPInd < this->GetMesh().GetNumberOfFacets());
    } else {
      /* Have a new cluster */
      clusterSet = BVPrimitiveSet();
      clusterStartInds.push_back(newPIndList.size());
    }
    newPIndList.push_back(*i_currPInd);    
    clusterSet.Insert(this->GetMesh().GetFacet(*i_currPInd));
  }

  assert(i_pIndsEnd - i_pIndsBegin == (int)newPIndList.size() && newPIndList.size() == tledHelper::MakeSortedUnique(newPIndList).size());
  std::copy(newPIndList.begin(), newPIndList.end(), i_pIndsBegin);

  return clusterStartInds;
}

namespace tledSelfCollisionBVHCreator_internal {
  template <class TBV>
  void InitialiseChildPrimitiveSets(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const TBV &bv, const std::vector<float> &primitiveCentroids) {
    tledFatalError("BV geometry not supported");
  }

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledAABB<2> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledAABB<2> > &bv, const std::vector<float> &primitiveCentroids);

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledAABB<4> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledAABB<4> > &bv, const std::vector<float> &primitiveCentroids);

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledOBB<2> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledOBB<2> > &bv, const std::vector<float> &primitiveCentroids);

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledOBB<4> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledOBB<4> > &bv, const std::vector<float> &primitiveCentroids);
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::_InitialiseChildPrimitiveSets(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const BoundingVolume &bv) {
  tledSelfCollisionBVHCreator_internal::InitialiseChildPrimitiveSets<BoundingVolume>(pIndsBegin, pIndsEnd, bv, this->GetPrimitiveCentroids());  
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::_ComputeClusterVolinoAxis(float *p_vAxis, const std::vector<int>::const_iterator pIndsBegin, const std::vector<int>::const_iterator pIndsEnd) {
  using namespace tledVectorArithmetic;

  std::vector<int>::const_iterator ic_pInd;

  assert(pIndsEnd > pIndsBegin);
  std::copy(this->GetPrimitiveBV(*pIndsBegin).VolinoAxis, this->GetPrimitiveBV(*pIndsBegin).VolinoAxis + 3, p_vAxis);
  for (ic_pInd = pIndsBegin + 1; ic_pInd < pIndsEnd; ic_pInd++) {
    Add(p_vAxis, p_vAxis, this->GetPrimitiveBVs()[*ic_pInd].VolinoAxis);
  }
  if (Norm(p_vAxis) > 0.f) {
    ScalarDiv(p_vAxis, Norm(p_vAxis));
  } 
}

namespace tledSelfCollisionBVHCreator_internal {
#ifndef NDEBUG
  template <class TBV>
  void CheckBV(const TBV &bv) {
    tledFatalError("BV type not supported");
  }

  template <>
  void CheckBV<tledAABB<2> >(const tledAABB<2> &bv);
  template <>
  void CheckBV<tledAABB<4> >(const tledAABB<4> &bv);

  template <>
  void CheckBV<tledOBB<2> >(const tledOBB<2> &bv);
  template <>
  void CheckBV<tledOBB<4> >(const tledOBB<4> &bv);
#else
  template <class TBV>
  void CheckBV(const TBV &bv) {
  }
#endif
}

template <class TBVH>
std::vector<int>::iterator tledSelfCollisionBVHCreator<TBVH>::_FindOptimumAdjacentPrimitive(std::vector<int>::iterator pIndsBegin, std::vector<int>::iterator pIndsEnd, const BVPrimitiveSet &cluster, const BoundingVolume &clusterBV) {
  using namespace tledVectorArithmetic;

  const float cAxisMag = Norm(clusterBV.VolinoAxis);

  std::vector<int>::iterator i_minVolPInd, i_currPInd;
  std::vector<int> adjacentPrims;
  float minVol, vol;
  BoundingVolume tmpVol;

  minVol = std::numeric_limits<float>::max();
  i_minVolPInd = pIndsEnd;
  for (i_currPInd = pIndsBegin; i_currPInd < pIndsEnd; i_currPInd++) if (cluster.IsAdjacent(this->GetMesh().GetFacet(*i_currPInd))) {
      assert(*i_currPInd >= 0 && *i_currPInd < this->GetMesh().GetNumberOfFacets() && int(this->GetPrimitiveBVs().size()) >= *i_currPInd);

#ifndef NDEBUG
      {
	float tmp[3], n[3];

	assert(Norm(Sub(tmp, this->GetMesh().ComputeNormalisedFacetNormal(n, *i_currPInd), this->GetPrimitiveBVs()[*i_currPInd].VolinoAxis)) < 1e-4);
	for (int v = 0; v < TBVH::ContactMesh::Facet::NumberOfVertices; v++) {
	  assert(this->GetPrimitiveBVs()[*i_currPInd].IsInside(this->GetMesh().GetNodeCoordinates(this->GetMesh().GetFacet(*i_currPInd).NodeIndices[v])));
	}

	tledSelfCollisionBVHCreator_internal::CheckBV<typename BoundingVolume::BaseBV>(clusterBV);
	tledSelfCollisionBVHCreator_internal::CheckBV<typename BoundingVolume::BaseBV>(this->GetPrimitiveBVs()[*i_currPInd]);
      }
#endif

      BoundingVolume::CopyBoundsFromBV(tmpVol, clusterBV);
      BoundingVolume::Merge(tmpVol, this->GetPrimitiveBVs()[*i_currPInd]);
      if (cAxisMag > 0.f) {
	if ((vol = tmpVol.ComputeVolume()*(1 + ComputeAngleNormalised(clusterBV.VolinoAxis, this->GetPrimitiveBVs()[*i_currPInd].VolinoAxis)/(tledPi/2))) < minVol) {
	  i_minVolPInd = i_currPInd;
	  minVol = vol;
	}
      } else {
	if ((vol = tmpVol.ComputeVolume()) < minVol) {
	  i_minVolPInd = i_currPInd;
	  minVol = vol;
	}
      }
      assert(vol == vol);
    }

  return i_minVolPInd;
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::_SplitBV2(std::vector<int>::iterator (*ppi_childIndexBounds)[2], const BoundingVolume &bv, const std::vector<int>::iterator pIndsBegin, const std::vector<int>::iterator pIndsEnd, const int child1Offset) {
  BoundingVolume child0, child1;

  assert(CheckConnectivity(std::vector<int>(pIndsBegin, pIndsEnd), this->GetMesh()));
  _InitialiseChildPrimitiveSets(pIndsBegin, pIndsEnd, bv);

  {
    std::vector<int>::iterator i_cluster0End, i_cluster1Start, i_minPInd;
    BVPrimitiveSet childSets[2] = {BVPrimitiveSet(this->GetMesh().GetFacet(*pIndsBegin)), BVPrimitiveSet(this->GetMesh().GetFacet(*(pIndsEnd - 1)))};

    child0 = this->GetPrimitiveBV(*pIndsBegin), child1 = this->GetPrimitiveBV(*(pIndsEnd - 1));

    std::copy(this->GetPrimitiveBV(*pIndsBegin).VolinoAxis, this->GetPrimitiveBV(*pIndsBegin).VolinoAxis + 3, child0.VolinoAxis);
    std::copy(this->GetPrimitiveBV(*(pIndsEnd - 1)).VolinoAxis, this->GetPrimitiveBV(*(pIndsEnd - 1)).VolinoAxis + 3, child1.VolinoAxis);

    for (i_cluster0End = pIndsBegin + 1, i_cluster1Start = pIndsEnd - 2; i_cluster0End <= i_cluster1Start;) {
      if ((i_minPInd = _FindOptimumAdjacentPrimitive(i_cluster0End, i_cluster1Start + 1, childSets[0], child0)) < i_cluster1Start + 1) {
	BoundingVolume::Merge(child0, this->GetPrimitiveBV(*i_minPInd));
	childSets[0].Insert(this->GetMesh().GetFacet(*i_minPInd));
	_ComputeClusterVolinoAxis(child0.VolinoAxis, pIndsBegin, i_cluster0End);
	
	std::iter_swap(i_minPInd, i_cluster0End);
	i_cluster0End++;	  
	assert(CheckConnectivity(std::vector<int>(pIndsBegin, i_cluster0End), this->GetMesh()));
      } 
      tledSelfCollisionBVHCreator_internal::CheckBV<typename BoundingVolume::BaseBV>(child0);

      if (i_cluster0End > i_cluster1Start) break;

      if ((i_minPInd = _FindOptimumAdjacentPrimitive(i_cluster0End, i_cluster1Start + 1, childSets[1], child1)) < i_cluster1Start + 1) {
	BoundingVolume::Merge(child1, this->GetPrimitiveBV(*i_minPInd));
	childSets[1].Insert(this->GetMesh().GetFacet(*i_minPInd));
	_ComputeClusterVolinoAxis(child1.VolinoAxis, i_cluster1Start + 1, pIndsEnd);

	std::iter_swap(i_minPInd, i_cluster1Start);
	i_cluster1Start--;
	assert(CheckConnectivity(std::vector<int>(i_cluster1Start + 1, pIndsEnd), this->GetMesh()));
      }
      tledSelfCollisionBVHCreator_internal::CheckBV<typename BoundingVolume::BaseBV>(child1);

#ifndef NDEBUG
      {
	std::vector<int>::iterator i_p;

	for (i_p = i_cluster0End; i_p < i_cluster1Start + 1; i_p++) {
	  if (childSets[0].IsAdjacent(this->GetMesh().GetFacet(*i_p)) || childSets[1].IsAdjacent(this->GetMesh().GetFacet(*i_p))) {
	    break;
	  }
	}
	assert(i_p < i_cluster1Start + 1 || i_cluster0End >= i_cluster1Start + 1);
      }
#endif
    }

    ppi_childIndexBounds[0][0] = pIndsBegin;
    ppi_childIndexBounds[0][1] = i_cluster0End;

    ppi_childIndexBounds[child1Offset][0] = i_cluster0End;
    ppi_childIndexBounds[child1Offset][1] = pIndsEnd;
  }

  if (child1Offset > 1) {
    if (ppi_childIndexBounds[0][1] - ppi_childIndexBounds[0][0] > 1) _SplitBV2(ppi_childIndexBounds, child0, ppi_childIndexBounds[0][0], ppi_childIndexBounds[0][1], child1Offset/2);
    else {
      for (int i = 1; i < child1Offset; i++) ppi_childIndexBounds[i][0] = ppi_childIndexBounds[i][1] = ppi_childIndexBounds[0][1];
    }

    if (ppi_childIndexBounds[child1Offset][1] - ppi_childIndexBounds[child1Offset][0] > 1) _SplitBV2(ppi_childIndexBounds + child1Offset, child1, ppi_childIndexBounds[child1Offset][0], ppi_childIndexBounds[child1Offset][1], child1Offset/2);
    else {
      for (int i = 1; i < child1Offset; i++) ppi_childIndexBounds[child1Offset+i][0] = ppi_childIndexBounds[child1Offset+i][1] = ppi_childIndexBounds[child1Offset][1];
    }
  }
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::SplitBV(std::vector<int>::iterator ppi_childIndexBounds[][2], const int bvInd, const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd) {
  const int bvhOrder = tledSelfCollisionBVHCreator::BVHOrder;

  assert(pIndsEnd - pIndsBegin >= 1);
  assert(CheckConnectivity(std::vector<int>(pIndsBegin, pIndsEnd), this->GetMesh()));
  _SplitBV2(ppi_childIndexBounds, this->GetOutput().GetBV(bvInd), pIndsBegin, pIndsEnd, bvhOrder/2);  
  for (int i = 0; i < bvhOrder; i++) {
    if (ppi_childIndexBounds[i][1] - ppi_childIndexBounds[i][0] == 0) {
      for (int k = i + 1; k < bvhOrder; k++) {
	ppi_childIndexBounds[k-1][1] = ppi_childIndexBounds[k][1];
	ppi_childIndexBounds[k-1][0] = ppi_childIndexBounds[k][0];
      }
      ppi_childIndexBounds[bvhOrder-1][0] = ppi_childIndexBounds[bvhOrder-1][1] = pIndsEnd;
    }
  }
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::InitBVHCreationData() {
  Superclass::InitBVHCreationData();

  m_PrimitiveBVs.resize(this->GetMesh().GetNumberOfFacets());
  for (typename std::vector<BoundingVolume>::iterator i_pBV = m_PrimitiveBVs.begin(); i_pBV < m_PrimitiveBVs.end(); i_pBV++) {
    i_pBV->PrimitiveIndex = i_pBV - m_PrimitiveBVs.begin();
    i_pBV->ComputeFromNodeList(this->GetMesh().GetFacet(i_pBV->PrimitiveIndex).NodeIndices, this->GetMesh().GetFacet(i_pBV->PrimitiveIndex).NodeIndices + TBVH::ContactMesh::Facet::NumberOfVertices, this->GetMesh().GetAllNodeCoordinates());
    i_pBV->AddMargin(this->GetOutput().GetMargin());

#ifndef NDEBUG
    tledSelfCollisionBVHCreator_internal::CheckBV<typename BoundingVolume::BaseBV>(*i_pBV);
    for (int const *pc_n = this->GetMesh().GetFacet(i_pBV->PrimitiveIndex).NodeIndices; pc_n < this->GetMesh().GetFacet(i_pBV->PrimitiveIndex).NodeIndices + TBVH::ContactMesh::Facet::NumberOfVertices; pc_n++) {
      assert(i_pBV->IsInside(this->GetMesh().GetNodeCoordinates(*pc_n)));
    }
#endif
  }

  for (int pInd = 0; pInd < this->GetMesh().GetNumberOfFacets(); pInd++) this->GetMesh().ComputeNormalisedFacetNormal(this->GetPrimitiveBV(pInd).VolinoAxis, pInd);
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::GenerateMain() {
  this->InitialiseClusters();
  tledLogDebugStream(tledHelper::Info() << "Clusters initialised");

  /* 
   * Top-down, recursive splitting of cluster BVs
   */
  this->GetOutput().GetGeometryClusterSubtreeRootIndices() = this->GetOutput().GetLeafBVIndices();
  this->GetOutput().GetLeafBVIndices().clear();

  for (std::vector<int>::const_iterator ic_clusterBVInd = this->GetOutput().GetGeometryClusterSubtreeRootIndices().begin(); ic_clusterBVInd < this->GetOutput().GetGeometryClusterSubtreeRootIndices().end(); ic_clusterBVInd++) {
    BoundingVolume &r_bv = this->GetOutput().GetBVs()[*ic_clusterBVInd];

    if (r_bv.PrimitiveIndex + 1 >= (int)m_ClusterStartInds.size()) {
      tledLogDebugStream(tledHelper::Info() << "Descending in last cluster");
      this->InitialiseTopDownRecursive(*ic_clusterBVInd, this->GetBVPrimitiveIndices().begin() + m_ClusterStartInds[r_bv.PrimitiveIndex], this->GetBVPrimitiveIndices().end());
      tledLogDebugStream(tledHelper::Info() << "Cluster top-down processing finished.");
    } else {
      tledLogDebugStream(tledHelper::Info() << "Descending in cluster " << int(ic_clusterBVInd - this->GetOutput().GetGeometryClusterSubtreeRootIndices().begin()));
      this->InitialiseTopDownRecursive(*ic_clusterBVInd, this->GetBVPrimitiveIndices().begin() + m_ClusterStartInds[r_bv.PrimitiveIndex], this->GetBVPrimitiveIndices().begin() + m_ClusterStartInds[r_bv.PrimitiveIndex+1]);
      tledLogDebugStream(tledHelper::Info() << "Cluster " << int(ic_clusterBVInd - this->GetOutput().GetGeometryClusterSubtreeRootIndices().begin()) << " processed.");
    }
    assert(this->GetOutput().GetBV(*ic_clusterBVInd).PrimitiveIndex == -1 || std::find(this->GetOutput().GetLeafBVIndices().begin(), this->GetOutput().GetLeafBVIndices().end(), *ic_clusterBVInd) != this->GetOutput().GetLeafBVIndices().end());
    assert((int)tledHelper::MakeSortedUnique(this->GetBVPrimitiveIndices()).size() == this->GetMesh().GetNumberOfFacets());
  }    
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::FreeBVHCreationData() {
  Superclass::FreeBVHCreationData();
  m_PrimitiveBVs.resize(0);
  m_ClusterStartInds.resize(0);
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::_CollapseSubTree(const int subtreeRoot, const int subTreeOrder) {  
  const int bvhOrder = tledSelfCollisionBVHCreator::BVHOrder;

  BoundingVolume &r_bv = this->GetOutput().GetBV(subtreeRoot);
  int oldChildren[bvhOrder/2];
  int *p_cIndDst = r_bv.ChildIndices;

  assert(subTreeOrder <= bvhOrder/2 && bvhOrder > 2);
  std::copy(r_bv.ChildIndices, r_bv.ChildIndices + subTreeOrder, oldChildren);
  assert(r_bv.ChildIndices[subTreeOrder] == -1);
  for (int const *pc_c = oldChildren; pc_c < oldChildren + subTreeOrder; pc_c++) if (*pc_c >= 0) {
      BoundingVolume &r_child = this->GetOutput().GetBV(*pc_c);

      if (r_child.PrimitiveIndex >= 0) {
	assert(r_child.ChildIndices[0] == -1 && r_child.ChildIndices[bvhOrder-1] == -1);
	*(p_cIndDst++) = *pc_c;
      } else {
	assert(r_child.ChildIndices[subTreeOrder] == -1 && r_child.ChildIndices[bvhOrder-1] == -1);
	for (int const *pc_gc = r_child.ChildIndices; pc_gc < r_child.ChildIndices + subTreeOrder; pc_gc++) if (*pc_gc >= 0) {
	    *(p_cIndDst++) = *pc_gc;
	    this->GetOutput().GetBV(*pc_gc).ParentIndex = subtreeRoot;
	    _CollapseSubTree(*pc_gc, subTreeOrder);
	  }
	std::fill(r_child.ChildIndices, r_child.ChildIndices + subTreeOrder, -1);
	r_child.ParentIndex = -1;
      }
    }  
}

template <class TBVH>
int tledSelfCollisionBVHCreator<TBVH>::InsertParentBinary(const int child0Ind, const int child1Ind) {
  const int parentInd = this->GetOutput().GetBVs().size();
  const int bvhOrder = tledSelfCollisionBVHCreator::BVHOrder;

  BoundingVolume &r_child0 = this->GetOutput().GetBV(child0Ind);
  BoundingVolume &r_child1 = this->GetOutput().GetBV(child1Ind);
  BoundingVolume parent;
  
  parent.PrimitiveIndex = parent.ParentIndex = -1;
  std::fill(parent.ChildIndices + 2, parent.ChildIndices + bvhOrder, -1);

  if (r_child1.ComputeVolume() > r_child0.ComputeVolume()) parent.ChildIndices[0] = child1Ind, parent.ChildIndices[1] = child0Ind;
  else parent.ChildIndices[0] = child0Ind, parent.ChildIndices[1] = child1Ind;
  r_child0.ParentIndex = r_child1.ParentIndex = parentInd;

  this->GetOutput().ComputeBoundsFromChildren(parent);
  this->GetOutput().GetBVs().push_back(parent);

  return parentInd;
}

template <class TBVH>
void tledSelfCollisionBVHCreator<TBVH>::InitialiseClusters() {
  const int bvhOrder = tledSelfCollisionBVHCreator::BVHOrder;

  std::vector<int> nodeList;
  std::vector<int> activeBVInds;
  std::vector<int>::const_iterator ic_clusterStartInd;

  m_ClusterStartInds = this->FindClusters(this->GetBVPrimitiveIndices().begin(), this->GetBVPrimitiveIndices().end());    
  for (ic_clusterStartInd = m_ClusterStartInds.begin(); ic_clusterStartInd < m_ClusterStartInds.end(); ic_clusterStartInd++) {
    BoundingVolume newBV;

    if (ic_clusterStartInd + 1 < m_ClusterStartInds.end()) {	
      nodeList = this->GetMesh().CompileNodeListFromFacetList(this->GetBVPrimitiveIndices().begin() + *ic_clusterStartInd, this->GetBVPrimitiveIndices().begin() + *(ic_clusterStartInd + 1));
    } else {
      nodeList = this->GetMesh().CompileNodeListFromFacetList(this->GetBVPrimitiveIndices().begin() + *ic_clusterStartInd, this->GetBVPrimitiveIndices().end());
    }

    newBV.ComputeFromNodeList(&nodeList.front(), &nodeList.back() + 1, this->GetMesh().GetAllNodeCoordinates());	
    newBV.AddMargin(this->GetOutput().GetMargin());
    newBV.ParentIndex = -1;
    newBV.PrimitiveIndex = ic_clusterStartInd - m_ClusterStartInds.begin();
    std::fill(newBV.ChildIndices, newBV.ChildIndices + bvhOrder, -1);

#ifndef NDEBUG
    for (std::vector<int>::const_iterator ic_n = nodeList.begin(); ic_n < nodeList.end(); ic_n++) {
      assert(newBV.IsInside(this->GetMesh().GetNodeCoordinates(*ic_n)));
    }
#endif

    this->GetOutput().GetBVs().push_back(newBV);
    activeBVInds.push_back(this->GetOutput().GetNumberOfBVs() - 1);
  }

  this->PreBottomUpHook(activeBVInds);

  while (activeBVInds.size() > 1) {
    float minVol, vol;
    std::vector<int>::iterator i_bvInd0, i_bvInd1;
    int minBVInd0, minBVInd1;  

    /*
     * Find pair of BVs that yields the smallest volume, group them:
     * remove from active BV list, add new parent node
     */
    minVol = std::numeric_limits<float>::max();
    for (i_bvInd0 = activeBVInds.begin(); i_bvInd0 < activeBVInds.end() - 1; i_bvInd0++) for (i_bvInd1 = i_bvInd0 + 1; i_bvInd1 < activeBVInds.end(); i_bvInd1++) {
	BoundingVolume tmpParent;

	std::fill(tmpParent.ChildIndices + 2, tmpParent.ChildIndices + BoundingVolume::NumberOfChildBVs, -1);
	tmpParent.ChildIndices[0] = *i_bvInd0, tmpParent.ChildIndices[1] = *i_bvInd1;
	this->GetOutput().ComputeBoundsFromChildren(tmpParent);
	if ((vol = tmpParent.ComputeVolume()) < minVol) {
	  minVol = vol;
	  minBVInd0 = *i_bvInd0, minBVInd1 = *i_bvInd1;
	}
      }
        
    activeBVInds.erase(std::find(activeBVInds.begin(), activeBVInds.end(), minBVInd0));
    activeBVInds.erase(std::find(activeBVInds.begin(), activeBVInds.end(), minBVInd1));

    activeBVInds.push_back(this->InsertParentBinary(minBVInd0, minBVInd1));        

    assert(this->GetOutput().GetBV(activeBVInds.back()).ChildIndices[0] >= 0 && this->GetOutput().GetBV(activeBVInds.back()).ChildIndices[1] >= 0);
    assert(this->GetOutput().GetBV(activeBVInds.back()).ChildIndices[0] < (int)this->GetOutput().GetBVs().size() && this->GetOutput().GetBV(activeBVInds.back()).ChildIndices[1] < (int)this->GetOutput().GetBVs().size());
  }

  if (bvhOrder > 2) {
    std::vector<BoundingVolume> initialBVs;
    
    /* Collapse */
    for (int m = 2; m < bvhOrder; m *= 2) {
      _CollapseSubTree(this->GetOutput().GetNumberOfBVs() - 1, m);
    }
    initialBVs.swap(this->GetOutput().GetBVs());
    this->GetOutput().GetBVs().reserve(bvhOrder*this->GetMesh().GetNumberOfFacets()/(bvhOrder - 1));
    _ReinsertBVsRecursive(initialBVs, initialBVs.size() - 1, -1);    
  } else {
    std::vector<BoundingVolume> initialBVs;
    
    initialBVs.swap(this->GetOutput().GetBVs());
    this->GetOutput().GetBVs().reserve(2*this->GetMesh().GetNumberOfFacets());
    _ReinsertBVsRecursive(initialBVs, initialBVs.size() - 1, -1);
    assert(this->GetOutput().GetBVs().size() == initialBVs.size());
  }

  assert(m_ClusterStartInds.size() == this->GetOutput().GetLeafBVIndices().size());

  /*
   * From this point on downwards the tree only consists of BVs bounding connected primitives.
   * Store this information.
   */
  this->GetOutput().GetNonAdjacentGeometryNodes().clear();
  for (int bvInd = 0; bvInd < this->GetOutput().GetNumberOfBVs(); bvInd++) {      
    BoundingVolume &r_bv = this->GetOutput().GetBVs()[bvInd];

    if (r_bv.PrimitiveIndex < 0) {
      r_bv.SetDisconnectedGeometryFlag();
      this->GetOutput().GetNonAdjacentGeometryNodes().push_back(bvInd);
      assert(std::find(this->GetOutput().GetLeafBVIndices().begin(), this->GetOutput().GetLeafBVIndices().end(), bvInd) == this->GetOutput().GetLeafBVIndices().end());
    } 
  }
}
