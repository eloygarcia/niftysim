// =========================================================================
// File:       tledSelfCollisionBVH.tpp
// Purpose:    
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

#ifndef NDEBUG
template <class TSurface>
static bool _CheckConnectivity(const std::vector<int> &primIndList, const TSurface &surface) {
  std::vector<int> uniqueInds;
  std::vector<int>::const_iterator ic_fInd;
  
  uniqueInds = tledHelper::MakeSortedUnique(primIndList);
  if (uniqueInds.size() == 1) return true;
  for (ic_fInd = uniqueInds.begin(); ic_fInd < uniqueInds.end(); ic_fInd++) {
    const int *tVtxInds = surface.GetFacet(*ic_fInd).NodeIndices;

    std::vector<int> otherPrimInds;
    std::vector<int> otherNodeInds;

    otherPrimInds.reserve(uniqueInds.size() - 1);
    std::copy(std::vector<int>::const_iterator(uniqueInds.begin()), ic_fInd, back_inserter(otherPrimInds));
    std::copy(ic_fInd + 1, std::vector<int>::const_iterator(uniqueInds.end()), back_inserter(otherPrimInds));
    assert(otherPrimInds.size() == uniqueInds.size() - 1);

    otherNodeInds = surface.CompileNodeListFromPrimitiveList(otherPrimInds.begin(), otherPrimInds.end());
    if (find(otherNodeInds.begin(), otherNodeInds.end(), tVtxInds[0]) == otherNodeInds.end()
	&& find(otherNodeInds.begin(), otherNodeInds.end(), tVtxInds[1]) == otherNodeInds.end()
	&& find(otherNodeInds.begin(), otherNodeInds.end(), tVtxInds[2]) == otherNodeInds.end()) {
      std::cerr << "Didn't find any of (", std::copy(tVtxInds, tVtxInds + 3, std::ostream_iterator<int>(std::cerr, " ")), std::cerr << ") in:\n";
      std::copy(otherNodeInds.begin(), otherNodeInds.end(), std::ostream_iterator<int>(std::cerr, " ")), std::cerr << std::endl;
      break;
    }
  }

  return ic_fInd == uniqueInds.end();
} /* _CheckConnectivity */
#endif

template <class TSlaveFacet, class TMasterFacet>
void _AddNarrowPhaseTests(std::vector<int> &r_nodeList, std::vector<int> &r_edgeList, std::vector<std::vector<int> > &r_nodeFacetIndLists, std::vector<std::vector<int> > &r_edgeMasterLists, const TSlaveFacet &sFacet, const TMasterFacet &mFacet, const int mFacetInd) {
  for (int const *pc_nodeInd = sFacet.NodeIndices; pc_nodeInd < sFacet.NodeIndices + TSlaveFacet::NumberOfVertices; pc_nodeInd++) {
    if (r_nodeFacetIndLists[*pc_nodeInd].size() == 0) r_nodeList.push_back(*pc_nodeInd);
    r_nodeFacetIndLists[*pc_nodeInd].push_back(mFacetInd);
    assert(r_nodeFacetIndLists[*pc_nodeInd].size() == 0 || find(r_nodeList.begin(), r_nodeList.end(), *pc_nodeInd) != r_nodeList.end());
  }

  for (int const *pc_edgeInd = sFacet.EdgeIndices; pc_edgeInd < sFacet.EdgeIndices + TSlaveFacet::NumberOfVertices; pc_edgeInd++) {
    if (r_edgeMasterLists[*pc_edgeInd].size() == 0) r_edgeList.push_back(*pc_edgeInd);	
    for (int const *pc_masterInd = mFacet.EdgeIndices; pc_masterInd < mFacet.EdgeIndices + TMasterFacet::NumberOfVertices; pc_masterInd++) {
      r_edgeMasterLists[*pc_edgeInd].push_back(*pc_masterInd);
    }
  }  
}

template <class TSurface, class TBV, class TAPI>
tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::tledSelfCollisionBVHImpl(ContactMesh &r_mesh) : Superclass(r_mesh) {
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::RefitBVCommon(const int bvInd) {
  Superclass::RefitBVCommon(bvInd);

  assert(this->GetUpdateCounter() == 0 || this->DoesNeedUpdate(bvInd));
  this->GetBV(bvInd).UpdateCounter = this->GetUpdateCounter();  
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::TranslateSubTree(const int rootBVInd, const float t[]) {
  Superclass::TranslateSubTree(rootBVInd, t);
  this->GetBV(rootBVInd).UpdateCounter = this->GetUpdateCounter();
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::TransformSubTree(const int rootBVInd, const float m[][3], const float cor[], const float t[]) {
  BoundingVolume &r_bv = this->GetBV(rootBVInd);
  float oldAxis[3];

  std::copy(r_bv.VolinoAxis, r_bv.VolinoAxis + 3, oldAxis);
  tledVectorArithmetic::MatMul(r_bv.VolinoAxis, m, oldAxis);  

  Superclass::TransformSubTree(rootBVInd, m, cor, t);
  r_bv.UpdateCounter = this->GetUpdateCounter();
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::ComputeSurfaceConePairwise(float &r_angle, float *p_axis, const float angle0, const float axis0[], const float angle1, const float axis1[]) {
  using namespace tledVectorArithmetic;

  float axisNorm;

  assert(Norm(axis0) != Norm(axis0) || std::fabs(Norm(axis0) - 1) < 1e-3);
  assert(Norm(axis1) != Norm(axis1) || std::fabs(Norm(axis1) - 1) < 1e-3);

  Add(p_axis, axis0, axis1);
  axisNorm = Norm(p_axis);
  if (axisNorm < 1e-3f) {
    /* Angle between Vol. axes ~Pi */
    r_angle = 2*BoundingVolume::GetVolinoThresholdAngle();
  } else {
    assert(axisNorm != axisNorm || axisNorm > 0);
    ScalarDiv(p_axis, axisNorm);
    assert(ComputeAngleNormalised(p_axis, axis0) < tledPi && ComputeAngleNormalised(p_axis, axis1) < tledPi);

    r_angle = ComputeAngleNormalised(axis0, axis1)*1.001f + std::max(angle0, angle1);
  }
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::RefitInteriorBV(const int bvInd) {
  using namespace tledVectorArithmetic;

  BoundingVolume &r_bv = this->GetBV(bvInd);

  if (!r_bv.HasDisconnectedGeometry()) {
#ifndef NDEBUG
    static const float volAngleChkVal = 1234;
    static const float volAxisChkVals[] = {234, 567, 890};
#endif

    int const *pc_childInd;

    for (pc_childInd = r_bv.ChildIndices; pc_childInd < r_bv.ChildIndices + BoundingVolume::NumberOfChildBVs; pc_childInd++) if (BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
	assert(!this->GetBV(*pc_childInd).HasDisconnectedGeometry());
	if (this->GetBV(*pc_childInd).VolinoAngle > BoundingVolume::GetVolinoThresholdAngle()) {
	  r_bv.VolinoAngle = 2*BoundingVolume::GetVolinoThresholdAngle();
	  break;
	}
      }

    if (pc_childInd == r_bv.ChildIndices + BoundingVolume::NumberOfChildBVs) {
      if (BoundingVolume::NumberOfChildBVs == 2) {
	const BoundingVolume &child0 = this->GetBV(r_bv.ChildIndices[0]);
	const BoundingVolume &child1 = this->GetBV(r_bv.ChildIndices[1]);

	assert(BVHOrder == 2);

	if (tledSelfCollisionBVHImpl::IsConeContainedInCone(child1.VolinoAxis, child1.VolinoAngle, child0.VolinoAxis, child0.VolinoAngle)) {
	  std::copy(child1.VolinoAxis, child1.VolinoAxis + 3, r_bv.VolinoAxis);
	  r_bv.VolinoAngle = child1.VolinoAngle;
	} else if (tledSelfCollisionBVHImpl::IsConeContainedInCone(child0.VolinoAxis, child0.VolinoAngle, child1.VolinoAxis, child1.VolinoAngle)) {
	  std::copy(child0.VolinoAxis, child0.VolinoAxis + 3, r_bv.VolinoAxis);
	  r_bv.VolinoAngle = child0.VolinoAngle;
	} else {
	  assert(&r_bv >= &this->GetBVs().front() && &r_bv <= &this->GetBVs().back());
	  assert(&child0 >= &this->GetBVs().front() && &child0 <= &this->GetBVs().back());
	  assert(&child1 >= &this->GetBVs().front() && &child1 <= &this->GetBVs().back());
	  
	  assert(!std::equal(child0.VolinoAxis, child0.VolinoAxis + 3, volAxisChkVals));
	  assert(!std::equal(child1.VolinoAxis, child1.VolinoAxis + 3, volAxisChkVals));
	  assert(child0.VolinoAngle != volAngleChkVal);
	  assert(child1.VolinoAngle != volAngleChkVal);
	  
	  assert(Norm(child0.VolinoAxis) != Norm(child0.VolinoAxis) || fabsf(Norm(child0.VolinoAxis) - 1) < 1e-3);
	  assert(Norm(child1.VolinoAxis) != Norm(child1.VolinoAxis) || fabsf(Norm(child1.VolinoAxis) - 1) < 1e-3);
	  
	  this->ComputeSurfaceConePairwise(r_bv.VolinoAngle, r_bv.VolinoAxis, child0.VolinoAngle, child0.VolinoAxis, child1.VolinoAngle, child1.VolinoAxis);
	  
	  assert(r_bv.VolinoAngle >= 0 && r_bv.VolinoAngle < 2.5f*tledPi);
	  assert(r_bv.VolinoAngle >= BoundingVolume::GetConeMinAngle());
	  assert(r_bv.VolinoAngle >= child0.VolinoAngle && r_bv.VolinoAngle >= child1.VolinoAngle);
	}

	assert(r_bv.VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle()
	       || (std::equal(r_bv.VolinoAxis, r_bv.VolinoAxis + 3, child0.VolinoAxis) && child0.VolinoAngle == r_bv.VolinoAngle)
	       || tledSelfCollisionBVHImpl::IsConeContainedInCone(r_bv.VolinoAxis, r_bv.VolinoAngle, child0.VolinoAxis, child0.VolinoAngle)
	       || ComputeAngleNormalised(r_bv.VolinoAxis, child0.VolinoAxis) < tledPi/180.0f);
	assert(r_bv.VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle()
	       || (std::equal(r_bv.VolinoAxis, r_bv.VolinoAxis + 3, child1.VolinoAxis) && child1.VolinoAngle == r_bv.VolinoAngle)
	       || tledSelfCollisionBVHImpl::IsConeContainedInCone(r_bv.VolinoAxis, r_bv.VolinoAngle, child1.VolinoAxis, child1.VolinoAngle)
	       || ComputeAngleNormalised(r_bv.VolinoAxis, child1.VolinoAxis) < tledPi/180.0f);
	
	r_bv.SubtreeMinH = std::min(child0.SubtreeMinH, child1.SubtreeMinH);
	assert(r_bv.SubtreeMinH == r_bv.SubtreeMinH);
      } else {
	assert(r_bv.ChildIndices[0] >= 0);

	{
	  const BoundingVolume &child0 = this->GetBV(r_bv.ChildIndices[0]);

	  r_bv.VolinoAngle = child0.VolinoAngle;
	  std::copy(child0.VolinoAxis, child0.VolinoAxis + 3, r_bv.VolinoAxis);
	  r_bv.SubtreeMinH = child0.SubtreeMinH;
	}

	for (int const *pc_c = r_bv.ChildIndices + 1; pc_c < r_bv.ChildIndices + BoundingVolume::NumberOfChildBVs; pc_c++) if (*pc_c >= 0) {
	    const BoundingVolume &childX = this->GetBV(*pc_c);

	    r_bv.SubtreeMinH = std::min(r_bv.SubtreeMinH, childX.SubtreeMinH);
	    assert(childX.VolinoAngle < BoundingVolume::GetVolinoThresholdAngle());
	    if (r_bv.VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle()) {
	      r_bv.VolinoAngle = 2*BoundingVolume::GetVolinoThresholdAngle();
	      break;
	    } else if (childX.VolinoAngle > r_bv.VolinoAngle) {
	      if (!tledSelfCollisionBVHImpl::IsConeContainedInCone(childX.VolinoAxis, childX.VolinoAngle, r_bv.VolinoAxis, r_bv.VolinoAngle)) {
		float tmpAxis[3];

		std::copy(r_bv.VolinoAxis, r_bv.VolinoAxis + 3, tmpAxis);	      
		this->ComputeSurfaceConePairwise(r_bv.VolinoAngle, r_bv.VolinoAxis, r_bv.VolinoAngle, tmpAxis, childX.VolinoAngle, childX.VolinoAxis);
	      } else {
		std::copy(childX.VolinoAxis, childX.VolinoAxis + 3, r_bv.VolinoAxis);
		r_bv.VolinoAngle = childX.VolinoAngle;
	      }
	    } else if (!tledSelfCollisionBVHImpl::IsConeContainedInCone(r_bv.VolinoAxis, r_bv.VolinoAngle, childX.VolinoAxis, childX.VolinoAngle)) {
	      float tmpAxis[3];

	      std::copy(r_bv.VolinoAxis, r_bv.VolinoAxis + 3, tmpAxis);	      
	      this->ComputeSurfaceConePairwise(r_bv.VolinoAngle, r_bv.VolinoAxis, r_bv.VolinoAngle, tmpAxis, childX.VolinoAngle, childX.VolinoAxis);
	    }
	    assert(r_bv.VolinoAngle >= BoundingVolume::GetVolinoThresholdAngle() || std::fabs(1.f - Norm(r_bv.VolinoAxis)) < 1e-3f);
	  } /* for valid child indices */
      }
    } /* if have to compute new cone */ 
	
    if (r_bv.VolinoAngle > BoundingVolume::GetVolinoThresholdAngle()) {
#ifndef NDEBUG
      r_bv.SubtreeMinH = volAngleChkVal;
      std::copy(volAxisChkVals, volAxisChkVals + 3, r_bv.VolinoAxis);
#endif
      m_NewSelfCollisionCandidateNodes.push_back(bvInd);
    } 
  } /* if cone criterion applicable */
    
  assert((r_bv.VolinoAngle > BoundingVolume::GetVolinoThresholdAngle() || r_bv.HasDisconnectedGeometry()) || std::find(m_NewSelfCollisionCandidateNodes.begin(), m_NewSelfCollisionCandidateNodes.end(), bvInd) == m_NewSelfCollisionCandidateNodes.end());
  assert((r_bv.VolinoAngle <= BoundingVolume::GetVolinoThresholdAngle() || r_bv.HasDisconnectedGeometry()) || std::find(m_NewSelfCollisionCandidateNodes.begin(), m_NewSelfCollisionCandidateNodes.end(), bvInd) < m_NewSelfCollisionCandidateNodes.end());

  /* Do this last, as BV will be marked as updated. Otherwise, race conditions will occur in multi-threaded updating! */
  Superclass::RefitInteriorBV(bvInd);
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::RefitLeafBV(const int bvInd) {
  using namespace tledVectorArithmetic;

  BoundingVolume &r_bv = this->GetBV(bvInd);

  /* acosf very inaccurate for small angles -> have to add some slack */
  r_bv.VolinoAngle = BoundingVolume::GetConeMinAngle();
  this->GetMesh().ComputeNormalisedFacetNormalCached(r_bv.VolinoAxis, r_bv.PrimitiveIndex);
  assert(std::fabs(Norm(r_bv.VolinoAxis) - 1) < 1e-3f);

  /**
   * Approximate diameter: c*sqrt(area), factor c = 0.95 for safety       
   */
  r_bv.SubtreeMinH = std::sqrt(Norm(this->GetMesh().GetFacetNormalCached(r_bv.PrimitiveIndex))/2)*0.95f;

  Superclass::RefitLeafBV(bvInd);
} 

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::Update() {  
  m_SelfCollisionCandidateNodes = m_NewSelfCollisionCandidateNodes;
  m_NewSelfCollisionCandidateNodes.clear();

  Superclass::Update();

  m_SelfCollisionCandidateNodes.insert(m_SelfCollisionCandidateNodes.end(), m_NewSelfCollisionCandidateNodes.begin(), m_NewSelfCollisionCandidateNodes.end());
  this->m_SelfCollisionCandidateNodes.insert(m_SelfCollisionCandidateNodes.end(), this->GetNonAdjacentGeometryNodes().begin(), this->GetNonAdjacentGeometryNodes().end());
  this->m_SelfCollisionCandidateNodes = tledHelper::MakeSortedUnique(m_SelfCollisionCandidateNodes); 
} /* Update */
   
#ifndef NDEBUG
template <class TSurface, class TBV, class TAPI>
bool tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::_CheckVolinoComprehensive(const int bvInd) const {
  using namespace tledVectorArithmetic;

  const BoundingVolume &bv = this->GetBV(bvInd);

  bool isBounded = true;

  if (bv.PrimitiveIndex < 0 && bv.VolinoAngle < BoundingVolume::GetVolinoThresholdAngle()) {
    std::vector<int> primitiveInds;
    float n[3];

    this->CompilePrimitiveListRecursive(primitiveInds, bvInd);
    for (std::vector<int>::const_iterator ic_pInd = primitiveInds.begin(); ic_pInd < primitiveInds.end(); ic_pInd++) {
      if (ComputeAngleNormalised(this->GetMesh().ComputeNormalisedFacetNormal(n, *ic_pInd), bv.VolinoAxis) > bv.VolinoAngle/2) {
	tledLogErrorStream(tledHelper::Warning() << "Surface-cone assumption violation:"
			   << "\nBV: " << bvInd
			   << "\naxis = " << bv.VolinoAxis[0] << " " << bv.VolinoAxis[1] << " " << bv.VolinoAxis[2]
			   << ", n = " << n[0] << " " << n[1] << " " << n[2]
			   << "\nalpha = " << bv.VolinoAngle
			   << "\nangle(a, n) = " << ComputeAngleNormalised(this->GetMesh().ComputeNormalisedFacetNormal(n, *ic_pInd), bv.VolinoAxis)
			   << "\naxis(c0) = " << this->GetBV(bv.ChildIndices[0]).VolinoAxis[0] << " " << this->GetBV(bv.ChildIndices[0]).VolinoAxis[1] << " " << this->GetBV(bv.ChildIndices[0]).VolinoAxis[2] 
			   << "\naxis(c1) = " << this->GetBV(bv.ChildIndices[1]).VolinoAxis[0] << " " << this->GetBV(bv.ChildIndices[1]).VolinoAxis[1] << " " << this->GetBV(bv.ChildIndices[1]).VolinoAxis[2]);	
      }

      isBounded &= ComputeAngleNormalised(this->GetMesh().ComputeNormalisedFacetNormal(n, *ic_pInd), bv.VolinoAxis) <= bv.VolinoAngle/2;
    }	
  }
  
  return isBounded;
}
#endif

template <class TSurface, class TBV, class TAPI>
bool tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::IsConeContainedInCone(const float axis[], const float angle, const float testAxis[], const float testAngle) {
  const float diffAngle = tledVectorArithmetic::ComputeAngleNormalised(axis, testAxis);

  return diffAngle + testAngle/2 <= angle/2;
}

#ifndef NDEBUG
template <class TSurface, class TBV, class TAPI>
bool tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::_CheckIntersectionRecursive(const int subtree0Ind, const int subtree1Ind) const {
  if (this->DoIntersect(subtree0Ind, subtree1Ind)) return true;
  else {
    const BoundingVolume &aabb0 = this->GetBV(subtree0Ind);
    const BoundingVolume &aabb1 = this->GetBV(subtree1Ind);

    if (aabb0.PrimitiveIndex < 0) {
      for (int const *pc_childInd = aabb0.ChildIndices; pc_childInd < aabb0.ChildIndices + 2; pc_childInd++) {
	if (_CheckIntersectionRecursive(*pc_childInd, subtree1Ind)) {
	  return true;
	}
      }
    } 

    if (aabb1.PrimitiveIndex < 0) {
      for (int const *pc_childInd = aabb1.ChildIndices; pc_childInd < aabb1.ChildIndices + 2; pc_childInd++) {
	if (_CheckIntersectionRecursive(*pc_childInd, subtree0Ind)) {
	  return true;
	}
      }
    }         
  }

  return false;
}

template <class TSurface, class TBV, class TAPI>
bool tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::_CheckCompletenessRecursive(const int subtreeInd) const {
  const BoundingVolume &aabb = this->GetBV(subtreeInd);

  if (aabb.PrimitiveIndex >= 0) {
    for (int vInd = 0; vInd < Facet::NumberOfVertices; vInd++) for (int cInd = 0; cInd < 3; cInd++) {
	if (!(this->GetMesh().GetNodeCoordinates(aabb.PrimitiveIndex, vInd)[cInd] > aabb.Bounds[cInd][0]
	      && this->GetMesh().GetNodeCoordinates(aabb.PrimitiveIndex, vInd)[cInd] < aabb.Bounds[cInd][1])) return false;
      }
  } else {
    for (int const *pc_childInd = aabb.ChildIndices; pc_childInd < aabb.ChildIndices + 2; pc_childInd++) {
      const BoundingVolume &child = this->GetBV(*pc_childInd);

      _CheckCompletenessRecursive(*pc_childInd);
      for (int cInd = 0; cInd < 3; cInd++) {
	if (!(child.Bounds[cInd][0] < child.Bounds[cInd][1] && child.Bounds[cInd][0] >= aabb.Bounds[cInd][0] && child.Bounds[cInd][1] <= aabb.Bounds[cInd][1])) {
	  return false;
	}
      }
    }
  }

  return true;
}
#endif

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::UpdateBottomUpRecursive(const int bvInd) {
  assert(!this->DoesNeedUpdate(bvInd));    
  if (this->GetBV(bvInd).ParentIndex >= 0) {
    const int parentInd = this->GetBV(bvInd).ParentIndex;
    const BoundingVolume &parent = this->GetBV(parentInd);

    for (int const *pc_childInd = parent.ChildIndices; pc_childInd < parent.ChildIndices + TBV::NumberOfChildBVs; pc_childInd++) if (TBV::NumberOfChildBVs == 2 || *pc_childInd >= 0) {
      if (this->DoesNeedUpdate(*pc_childInd)) return;
    }

    this->RefitInteriorBV(parentInd);
    this->UpdateBottomUpRecursive(parentInd);
    assert(!this->DoesNeedUpdate(parentInd));
  }    
} /* UpdateBottomUpRecursive */

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::Init(tledBVHCreator &r_bvhBuilder) {
  Superclass::Init(r_bvhBuilder);
  m_SelfCollisionCandidateNodes = m_NewSelfCollisionCandidateNodes;
  m_NewSelfCollisionCandidateNodes.clear();
}

template <class TSurface, class TBV, class TAPI>
void tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::AddSelfCollisionCandidate(const int nodeInd) {
  assert(nodeInd < this->GetNumberOfBVs());
  m_NewSelfCollisionCandidateNodes.push_back(nodeInd);
}

template <class TSurface, class TBV, class TAPI>
XMLNode tledSelfCollisionBVHImpl<TSurface, TBV, TAPI>::ExportToXML() const {
  tledSelfCollisionBVHXMLExporter<tledSelfCollisionBVHImpl> exporter;
  XMLNode root;

  exporter.SetInput(*this);
  exporter.Export();
  root = exporter.GetRootNode();

  return root;
}
