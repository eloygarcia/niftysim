// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdater.tpp
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
/* Check that all paths from the root end in an update node, i.e. no leafs can be reached w/o encountering a UN */
template <class TBVH>
bool tledGeometricSelfCollisionBVHUpdater<TBVH>::_CheckPathHasUpdateNode(const int bvInd) const {
  if (IsUpdateNode(bvInd)) return true;
  else {
    if (this->GetBVH().IsLeaf(bvInd)) return false;
    for (int const *pc_childInd = this->GetBVH().GetBV(bvInd).ChildIndices; pc_childInd < this->GetBVH().GetBV(bvInd).ChildIndices + BoundingVolume::NumberOfChildBVs; pc_childInd++) if (BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
	if (!_CheckPathHasUpdateNode(*pc_childInd)) return false;
      }
    
    return true;
  }
}
#endif

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdater<TBVH>::Init() {
  using namespace tledVectorArithmetic;

  static const int bvOrder = BoundingVolume::NumberOfChildBVs;

  std::vector<int> updateNodeInds;
  float minPairVolAngle;

  if (m_MaxNumUN == -1) this->m_MaxNumUN = std::max(this->GetMesh().GetNumberOfFacets()/4, bvOrder);
  if (m_MinNumUN == -1) this->m_MinNumUN = bvOrder;
  if (!(m_MaxConeAngle == m_MaxConeAngle)) this->m_MaxConeAngle = BoundingVolume::GetVolinoThresholdAngle()/2;

  /*
   * While update node count > target num.
   * - Remove leafs from update node set
   * - Find two nodes with smallest merged volino angle
   * - Merge
   */
  assert(m_MaxNumUN >= 2 && m_MaxConeAngle > 0 && m_MaxConeAngle < BoundingVolume::GetVolinoThresholdAngle());

  updateNodeInds = this->GetBVH().GetLeafBVIndices();
  std::sort(updateNodeInds.begin(), updateNodeInds.end());

  do {
    const float maxAcceptAngle = 3*BoundingVolume::GetVolinoThresholdAngle()/4 + m_MaxConeAngle/4;

    int minParentInd = -1;
    int remIdcs[bvOrder], numRem = -1;

#ifndef NDEBUG
    std::fill(remIdcs, remIdcs + bvOrder, -1);
#endif
    minPairVolAngle = maxAcceptAngle;
    for (std::vector<int>::const_iterator ic_un = updateNodeInds.begin(); ic_un < updateNodeInds.end(); ic_un++) {      
      if (!this->GetBVH().GetBV(this->GetBVH().GetBV(*ic_un).ParentIndex).HasDisconnectedGeometry() && (this->GetBVH().GetBV(this->GetBVH().GetBV(*ic_un).ParentIndex).VolinoAngle < minPairVolAngle || this->GetBVH().IsLeaf(*ic_un))) {
	const int candBVInd = this->GetBVH().GetBV(*ic_un).ParentIndex;
	const BoundingVolume &candBV = this->GetBVH().GetBV(candBVInd);

	int candRemIdcs[bvOrder];
	int numFound = 0, numUnused = 0;

	for (int const *pc_c = candBV.ChildIndices; pc_c < candBV.ChildIndices + bvOrder; pc_c++) {
	  if (bvOrder == 2 || *pc_c >= 0) {
	    std::vector<int>::const_iterator ic_c = std::lower_bound(updateNodeInds.begin(), updateNodeInds.end(), *pc_c);

	    if (ic_c != updateNodeInds.end() && *ic_c == *pc_c) {
	      candRemIdcs[numFound++] = ic_c - updateNodeInds.begin();
	    } else break;
	  } else numUnused += 1;
	}

	if (numFound == bvOrder - numUnused) {
	  std::copy(candRemIdcs, candRemIdcs + numFound, remIdcs);
	  numRem = numFound;
	  minParentInd = candBVInd;
	  if (this->GetBVH().IsLeaf(*ic_un)) {
	    minPairVolAngle = 0;
	    break;
	  } else {
	    minPairVolAngle = candBV.VolinoAngle;
	  }
	}
      }
    }
    
    if (minPairVolAngle >= maxAcceptAngle) break;

    assert(minPairVolAngle >= 0);
    assert(minParentInd >= 0 && minParentInd < this->GetBVH().GetNumberOfBVs());
    
    {
      const int parentInsertInd = std::lower_bound(updateNodeInds.begin(), updateNodeInds.end(), minParentInd) - updateNodeInds.begin();
#ifndef NDEBUG
      const std::vector<int> oldUN = updateNodeInds;
#endif

      std::vector<int>::iterator i_out;
      std::vector<int>::const_iterator ic_in;

      /* By construction of the BVH: we know that always parent index < child indices */           
      std::sort(remIdcs, remIdcs + numRem);
      
      assert(minParentInd < updateNodeInds[remIdcs[0]]);
      assert(parentInsertInd <= remIdcs[0]);
      assert(remIdcs[0] >= 0);

      std::copy_backward(updateNodeInds.begin() + parentInsertInd, updateNodeInds.begin() + remIdcs[0], updateNodeInds.begin() + remIdcs[0] + 1);	
      updateNodeInds[parentInsertInd] = minParentInd;
      i_out = updateNodeInds.begin() + remIdcs[0] + 1;
      ic_in = updateNodeInds.begin() + remIdcs[0] + 1;
      
      for (int const *pc_remInd = remIdcs + 1; pc_remInd < remIdcs + numRem; pc_remInd++) {
	const int numCopied = std::vector<int>::const_iterator(updateNodeInds.begin() + *pc_remInd) - ic_in;

	assert(numCopied >= 0 && *pc_remInd >= 0);
	i_out = std::copy(ic_in, std::vector<int>::const_iterator(updateNodeInds.begin() + *pc_remInd), i_out);
	ic_in += numCopied + 1;
      }
      
      std::copy(ic_in, std::vector<int>::const_iterator(updateNodeInds.end()), i_out);
      updateNodeInds.resize(updateNodeInds.size() - numRem + 1);

#ifndef NDEBUG
      for (std::vector<int>::const_iterator ic_un = updateNodeInds.begin(); ic_un < updateNodeInds.end() - 1; ic_un++) assert(*ic_un < *(ic_un + 1)); 
      for (std::vector<int>::const_iterator ic_un = oldUN.begin(); ic_un < oldUN.end(); ic_un++) {
	int i;

	for (i = 0; i < numRem && *ic_un != oldUN[remIdcs[i]]; i++); 
	assert((std::find(updateNodeInds.begin(), updateNodeInds.end(), *ic_un) == updateNodeInds.end()) != (i == numRem));
      }
#endif
    }

    assert(updateNodeInds.size() == tledHelper::MakeSortedUnique(updateNodeInds).size());
  } while ((int)updateNodeInds.size() > m_MaxNumUN || (minPairVolAngle < m_MaxConeAngle && (int)updateNodeInds.size() > m_MinNumUN));
  assert(int(updateNodeInds.size()) <= m_MaxNumUN || minPairVolAngle >= m_MaxConeAngle);

#ifndef NDEBUG
  if (int(updateNodeInds.size()) >= m_MaxNumUN) {
    tledLogErrorStream(tledHelper::Warning() << "Number of UN exceeds requested " << m_MaxNumUN << "; is " << int(updateNodeInds.size()));
  }
#endif

  GetUpdateNodes().reserve(updateNodeInds.size());
  GetUpdateNodeIndices().reserve(updateNodeInds.size()*3);
  for (std::vector<int>::const_iterator ic_updateNodeInd = updateNodeInds.begin(); ic_updateNodeInd < updateNodeInds.end(); ic_updateNodeInd++) {
    std::vector<int> bvhNodePInds, bvhNodeNInds;
    UpdateNodeInfo bvhNodeInfo;

    bvhNodeInfo.BVIndex = *ic_updateNodeInd;
#ifdef NO_LEAF_UN
    assert(!this->GetBVH().IsLeaf(bvhNodeInfo.BVIndex));
#endif

    this->GetBVH().CompilePrimitiveListRecursive(bvhNodePInds, bvhNodeInfo.BVIndex);
    bvhNodeNInds = this->GetMesh().CompileNodeListFromFacetList(bvhNodePInds.begin(), bvhNodePInds.end());
    bvhNodeInfo.NodeStartIndex = GetUpdateNodeIndices().size();
    GetUpdateNodeIndices().insert(GetUpdateNodeIndices().end(), bvhNodeNInds.begin(), bvhNodeNInds.end());
    bvhNodeInfo.NodeEndIndex = GetUpdateNodeIndices().size();
    bvhNodeInfo.LastUpdateConeAngle = this->GetBVH().GetBV(bvhNodeInfo.BVIndex).VolinoAngle;
    bvhNodeInfo.Status = (bvhNodeInfo.LastUpdateConeAngle > BoundingVolume::GetVolinoThresholdAngle()? ACTIVE : 0);

    assert(!this->GetBVH().GetBV(bvhNodeInfo.BVIndex).HasDisconnectedGeometry() || (bvhNodeInfo.Status & ACTIVE) != 0);

#ifdef __USE_OBB
    bvhNodeInfo.RigidUpdateCounter = 0;
#endif

    GetUpdateNodes().push_back(bvhNodeInfo);
  } /* for update BVH nodes */

  assert(_CheckPathHasUpdateNode(0));
  
  InitialiseUpdateNodePositionBuffers();
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdater<TBVH>::InitialiseUpdateNodePositionBuffers() {
  using namespace tledVectorArithmetic;

  this->GetUpdateNodePreviousPositions().resize(this->GetUpdateNodeIndices().size()*3);
  this->GetUpdateNodeCurrentPositions().resize(this->GetUpdateNodePreviousPositions().size());
  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++) {
    std::fill(i_updateNode->CurrentCentroid, i_updateNode->CurrentCentroid + 3, 0.0f);
    for (std::vector<int>::const_iterator ic_nInd = this->GetUpdateNodeIndices().begin() + i_updateNode->NodeStartIndex; ic_nInd < this->GetUpdateNodeIndices().begin() + i_updateNode->NodeEndIndex; ic_nInd++) {
      Add(i_updateNode->CurrentCentroid, i_updateNode->CurrentCentroid, this->GetMesh().GetNodeCoordinates(*ic_nInd));
    }
    ScalarDiv(i_updateNode->CurrentCentroid, float(i_updateNode->NodeEndIndex - i_updateNode->NodeStartIndex));
    
    for (int localNodeInd = i_updateNode->NodeStartIndex; localNodeInd < i_updateNode->NodeEndIndex; localNodeInd++) {
      Sub(&this->GetUpdateNodeCurrentPositions()[3*localNodeInd], this->GetMesh().GetNodeCoordinates(GetUpdateNodeIndices()[localNodeInd]), i_updateNode->CurrentCentroid);    
    }

    this->ComputeBoundedNodeCentroid(*i_updateNode);    
    this->StoreUpdateNodePositions(*i_updateNode);    
  }
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdater<TBVH>::StoreUpdateNodePositions(UpdateNodeInfo &r_updateInfo) {
  assert(r_updateInfo.NodeStartIndex < r_updateInfo.NodeEndIndex);
  assert(3*r_updateInfo.NodeEndIndex <= (int)m_UpdateNodeLastUpdateNodePositions.size());
  std::copy(m_UpdateNodeCurrPositions.begin() + 3*r_updateInfo.NodeStartIndex, m_UpdateNodeCurrPositions.begin() + 3*r_updateInfo.NodeEndIndex, m_UpdateNodeLastUpdateNodePositions.begin() + 3*r_updateInfo.NodeStartIndex);
  std::copy(r_updateInfo.CurrentCentroid, r_updateInfo.CurrentCentroid + 3, r_updateInfo.LastUpdateCentroid);
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdater<TBVH>::ComputeBoundedNodeCentroid(UpdateNodeInfo &r_updateInfo) {
  using namespace tledVectorArithmetic;

  const float *nodes = this->GetMesh().GetAllNodeCoordinates();
  const int nodesStartInd = r_updateInfo.NodeStartIndex;
  const int nodesEndInd = r_updateInfo.NodeEndIndex;

  std::fill(r_updateInfo.CurrentCentroid, r_updateInfo.CurrentCentroid + 3, 0.0f);
  for (std::vector<int>::const_iterator ic_nodeInd = GetUpdateNodeIndices().begin() + nodesStartInd; ic_nodeInd < GetUpdateNodeIndices().begin() + nodesEndInd; ic_nodeInd++) {
    assert(Norm(nodes + 3*(*ic_nodeInd)) == Norm(nodes + 3*(*ic_nodeInd)));
    Add(r_updateInfo.CurrentCentroid, r_updateInfo.CurrentCentroid, nodes + 3*(*ic_nodeInd));
  }  
  ScalarDiv(r_updateInfo.CurrentCentroid, float(nodesEndInd - nodesStartInd));
}

#ifndef NDEBUG
template <class TBVH>
bool tledGeometricSelfCollisionBVHUpdater<TBVH>::IsUpdateNode(const int bvInd) const {
  for (typename std::vector<UpdateNodeInfo>::const_iterator ic_un = GetUpdateNodes().begin(); ic_un < GetUpdateNodes().end(); ic_un++) if (ic_un->BVIndex == bvInd) return true;
  return false;
}
#endif

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdater<TBVH>::_SetDefaults() {
  m_MaxNumUN = m_MinNumUN = -1;
  m_MaxConeAngle = std::numeric_limits<float>::quiet_NaN();
}
