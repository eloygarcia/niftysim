// =========================================================================
// File:       tledBVHTraverserCPU.tpp
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

struct tledBVHTraverserCPU::NarrowPhaseOrdering {
  bool operator()(const std::pair<int, int> &p0, const std::pair<int, int> &p1) const {
    return p0.first < p1.first || (p0.first == p1.first && p0.second < p1.second);
  }
};

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::AddNodeFacetPairs(const int *nodesBegin, const int *nodesEnd, const int masterFacetInd) {
  for (int const *pc_nodeInd = nodesBegin; pc_nodeInd < nodesEnd; pc_nodeInd++) {
    this->GetNodeFacetNarrowPhasePairs().push_back(std::pair<int, int>(*pc_nodeInd, masterFacetInd));
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::AddEdgeEdgePairs(const int *slaveEdgesBegin, const int *slaveEdgesEnd, const int *masterEdgesBegin, const int *masterEdgesEnd) {
  for (int const *pc_edgeInd = slaveEdgesBegin; pc_edgeInd < slaveEdgesEnd; pc_edgeInd++) {
    for (int const *pc_masterInd = masterEdgesBegin; pc_masterInd < masterEdgesEnd; pc_masterInd++) {
      this->GetEdgeEdgeNarrowPhasePairs().push_back(std::pair<int, int>(*pc_edgeInd, *pc_masterInd));
    }
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::AddNarrowPhaseTests(const int sFacetInd, const int mFacetInd) {
  typedef typename SlaveMesh::Facet __SlaveFacet;
  typedef typename MasterMesh::Facet __MasterFacet;

  const __SlaveFacet &sFacet = this->GetSlaveMesh().GetFacet(sFacetInd);
  const __MasterFacet &mFacet = this->GetMasterMesh().GetFacet(mFacetInd);

  this->AddNodeFacetPairs(sFacet.NodeIndices, sFacet.NodeIndices + __SlaveFacet::NumberOfVertices, mFacetInd);
  this->AddEdgeEdgePairs(sFacet.EdgeIndices, sFacet.EdgeIndices + __SlaveFacet::NumberOfVertices, mFacet.EdgeIndices, mFacet.EdgeIndices + __MasterFacet::NumberOfVertices);
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::AddNarrowPhaseTestsSwitching(const int sFacetInd, const int mFacetInd) {
  if (this->DoMaster()) {
    typedef typename SlaveMesh::Facet __SlaveFacet;
    typedef typename MasterMesh::Facet __MasterFacet;

    const __SlaveFacet &sFacet = this->GetSlaveMesh().GetFacet(sFacetInd);
    const __MasterFacet &mFacet = this->GetMasterMesh().GetFacet(mFacetInd);

    assert(sFacetInd >= 0 && sFacetInd < this->GetSlaveMesh().GetNumberOfFacets());
    assert(mFacetInd >= 0 && mFacetInd < this->GetMasterMesh().GetNumberOfFacets());
    this->AddNodeFacetPairs(sFacet.NodeIndices, sFacet.NodeIndices + __SlaveFacet::NumberOfVertices, mFacetInd);
    this->AddEdgeEdgePairs(sFacet.EdgeIndices, sFacet.EdgeIndices + __SlaveFacet::NumberOfVertices, mFacet.EdgeIndices, mFacet.EdgeIndices + __MasterFacet::NumberOfVertices);
  } else {
    typedef typename SlaveMesh::Facet __SlaveFacet;
    typedef typename MasterMesh::Facet __MasterFacet;

    const __SlaveFacet &mFacet = this->GetSlaveMesh().GetFacet(sFacetInd);
    const __MasterFacet &sFacet = this->GetMasterMesh().GetFacet(mFacetInd);

    assert(mFacetInd >= 0 && mFacetInd < this->GetMasterMesh().GetNumberOfFacets());
    assert(sFacetInd >= 0 && sFacetInd < this->GetSlaveMesh().GetNumberOfFacets());
    this->AddNodeFacetPairs(sFacet.NodeIndices, sFacet.NodeIndices + __MasterFacet::NumberOfVertices, sFacetInd);
    this->AddEdgeEdgePairs(sFacet.EdgeIndices, sFacet.EdgeIndices + __MasterFacet::NumberOfVertices, mFacet.EdgeIndices, mFacet.EdgeIndices + __SlaveFacet::NumberOfVertices);
  }
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::DoBroadPhaseDetectionRecursive(const int masterStartInd, const int slaveStartInd) {
  const typename MasterBVH::BoundingVolume &aabb0 = this->GetMasterBVH().GetBV(masterStartInd);
  const typename SlaveBVH::BoundingVolume &aabb1 = this->GetSlaveBVH().GetBV(slaveStartInd);

  if (aabb0.DoesIntersect(aabb1)) {
    if (aabb0.PrimitiveIndex >= 0 && aabb1.PrimitiveIndex >= 0) {
      this->AddNarrowPhaseTests(aabb1.PrimitiveIndex, aabb0.PrimitiveIndex);
    } else {
      if ((aabb0.ComputeVolume() > aabb1.ComputeVolume() && aabb0.PrimitiveIndex < 0) || aabb1.PrimitiveIndex >= 0) {
	assert(aabb0.PrimitiveIndex < 0);
	for (int const *pc_childInd = aabb0.ChildIndices; pc_childInd < aabb0.ChildIndices + MasterBVH::BoundingVolume::NumberOfChildBVs; pc_childInd++) if (MasterBVH::BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0) {	  
	    assert(*pc_childInd >= 0 && *pc_childInd < this->GetMasterBVH().GetNumberOfBVs());
	    this->DoBroadPhaseDetectionRecursive(*pc_childInd, slaveStartInd);
	  }
      } else {
	assert(aabb1.PrimitiveIndex < 0);
	for (int const *pc_childInd = aabb1.ChildIndices; pc_childInd < aabb1.ChildIndices + SlaveBVH::BoundingVolume::NumberOfChildBVs; pc_childInd++) if (SlaveBVH::BoundingVolume::NumberOfChildBVs == 2 || *pc_childInd >= 0)  {	  
	    assert(*pc_childInd >= 0 && *pc_childInd < this->GetSlaveBVH().GetNumberOfBVs());
	    this->DoBroadPhaseDetectionRecursive(masterStartInd, *pc_childInd);
	  }
      }      
    } /* if leafs else .. */
  } /* if intersection */
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::Init(tledUnstructuredContactManager &r_manager) {
  this->GetNodeProjectionBuffer().resize(std::max(this->GetMasterMesh().GetNumberOfNodes(), std::max(this->GetSlaveMesh().GetNumberOfNodes(), (int)this->GetNodeProjectionBuffer().size())));
  this->GetEdgeProjectionBuffer().resize(std::max(this->GetMasterMesh().GetNumberOfEdges(), std::max(this->GetSlaveMesh().GetNumberOfEdges(), (int)this->GetEdgeProjectionBuffer().size())));
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::ProjectOntoFacetC0Iterative(float *p_xi, float *p_n, const float x[], const float facetNodes[][3], const float nodeNormals[][3]) const {  
  using namespace tledVectorArithmetic;

  static const int maxProjIts = 10;

  float shapeVals[MasterMesh::Facet::NumberOfVertices], currItN[3], itXi[3], a[3], b[3], r[3], tmp[3];
  int it;

  /* In plane coordinate computation needs modifying if to be extended to quads. */
  Sub(a, facetNodes[1], facetNodes[0]);
  Sub(b, facetNodes[2], facetNodes[0]);
  Sub(r, x, facetNodes[0]);
  for (it = 0; it < maxProjIts; it++) {
    this->GetMasterMesh().ComputeShapeValues(shapeVals, p_xi[0], p_xi[1]);

    std::fill(currItN, currItN + 3, 0.f);
    for (int vInd = 0; vInd < MasterMesh::Facet::NumberOfVertices; vInd++) Add(currItN, currItN, ScalarMul(tmp, nodeNormals[vInd], shapeVals[vInd]));
    ScalarDiv(currItN, Norm(currItN));

    SolveEquation3x3(itXi[0], itXi[1], itXi[2], a, b, currItN, r);  
    if (itXi[0] > -1e-4f && itXi[1] > -1e-4f && itXi[0] + itXi[1] < 1 + 1e-4f) {
      float sum;

      itXi[0] = std::max(itXi[0], 0.0f);
      itXi[1] = std::max(itXi[1], 0.0f);
      if ((sum = itXi[0] + itXi[1]) > 1) {
	itXi[0] /= sum;
	itXi[1] /= sum;
      }

      if (Norm(Sub(tmp, itXi, p_xi)) <= 1e-5*Norm(p_xi)) {
	std::copy(currItN, currItN + 3, p_n);

	return true;
      } else std::copy(itXi, itXi + 3, p_xi);
    } else {
      std::copy(currItN, currItN + 3, p_n);
      break;
    } 
  }

  tledLogDebugStream(tledHelper::Warning() << "C0 projection failed after " << it << " iterations.");    

  return false;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::DoBehindFacetTest(const float oldSlaveNodePosition[], const float oldFacetNormal[], const float oldFacetV0[]) const {
  using namespace tledVectorArithmetic;

  float tmp[3];

  assert(std::fabs(Norm(oldFacetNormal) - 1.0f) < 1e-3);
  return Dot(Sub(tmp, oldSlaveNodePosition, oldFacetV0), oldFacetNormal) >= -this->GetNarrowPhaseMaxDistance();  
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::DoNodeFacetInitialProjection(float *p_xiInitial, const float slaveNodePosition[], const float facetProjectionOperator[], const float currentBestProjection[]) const {
  return MasterMesh::ProjectOntoFacet(p_xiInitial, slaveNodePosition, facetProjectionOperator) && std::fabs(p_xiInitial[2]) < 1.25*this->GetNarrowPhaseMaxDistance() 
    && !(currentBestProjection[2] < 0 && p_xiInitial[2] > 0) && std::fabs(p_xiInitial[2])/1.5 < std::fabs(currentBestProjection[2]);  
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::IsBetterNodeProjection(const float xi0[], const float xi1[]) const {
  assert(xi0[0] >= 0 && xi0[1] >= 0 && xi0[0] + xi0[1] <= 1 + 1e-6);
      
  return xi0[2] < xi1[2] && xi0[2] < this->GetNarrowPhaseMaxDistance();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::IsBetterEdgeProjection(const float xi0[], const float xi1[]) const {
  return xi0[0] < xi1[0] && xi0[0] < this->GetNarrowPhaseMaxDistance();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::DoEdgeBehindEdgeTest(const float slave0[], const float slave1[], const float master0[], const float master1[], const float masterNormal0[], const float masterNormal1[]) const {
  using namespace tledVectorArithmetic;
  
  float diff[3];

  return (Dot(Sub(diff, slave0, master0), masterNormal0) < 0
	  || Dot(Sub(diff, slave1, master0), masterNormal0) < 0
	  || Dot(Sub(diff, slave0, master1), masterNormal1) < 0
	  || Dot(Sub(diff, slave1, master1), masterNormal1) < 0);
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::ComputeEdgeEdgeClosestPointParameters(float &r_r, float &r_q, const float A[], const float B[], const float C[], const float D[]) {
  using namespace tledVectorArithmetic;

  const float safetyMargin = 1e-2f;

  /*
   * Based on "Distance between lines" from http://geomalgorithms.com/a07-_distance.html, in turn derived from work by (Teller 2000, Eberly 2001)
   */
  float u[3], v[3];
  float a, b, c, d, e, CA[3], denom;
  
  Sub(u, B, A);
  Sub(v, D, C);
  Sub(CA, A, C);
  
  a = Dot(u, u);
  b = Dot(u, v);
  c = Dot(v, v);
  d = Dot(u, CA);
  e = Dot(v, CA);

  denom = a*c - b*b;
  assert(denom >= -1e-4f*(a + c));
  denom = std::max(denom, 0.0f);
  if (tledHelper::IsNumericallyZero(denom, 128*(a + c))) {
    /* 
     * Lines are (almost) parallel. Use r = 0.5 for best response results.
     * Could use optimisation to get value close to centre for both segments... 
     */
    r_r = 0.5f;
    r_q = (d + a/2)/b;
    if (r_q < -safetyMargin) {
      r_q = 0.0f;
      r_r = -d/a;
    } else if (r_q > 1 + safetyMargin) {
      r_q = 1.0f;
      r_r = (c - e)/b;
    }
  } else { 
    r_r = (b*e - c*d)/denom;
    r_q = (a*e - b*d)/denom;
  }

  if (r_r >= -safetyMargin && r_r <= 1.0f + safetyMargin && r_q >= -safetyMargin && r_q <= 1.0f + safetyMargin) {
    const float thr = (a + c)/1024;

    float ab[3], cd[3], d[3];

    /* Clamping and additional check of orthogonality to see if we have sufficient trust in result */
    assert((a != a || a >= 0) && (c != c || c >= 0));
    Interpolate(ab, A, B, r_r);
    Interpolate(cd, C, D, r_q);
    Sub(d, ab, cd);

    r_r = std::max(0.0f, std::min(r_r, 1.0f));
    r_q = std::max(0.0f, std::min(r_q, 1.0f));
    if (!(std::fabs(Dot(u, d)) < thr && std::fabs(Dot(v, d)) < thr)) return false;
    else return true;
  } else return false;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
template <const bool t_isMovingSlave, const bool t_isMovingMaster>
bool tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::ComputeEdgeEdgePenetrationDepth(float &r_g, float *p_n, const float slave0T0[], const float slave1T0[], const float slave0T1[], const float slave1T1[], const float slaveN0[], const float slaveN1[], const float r, const float master0T0[], const float master1T0[], const float master0T1[], const float master1T1[], const float masterN0[], const float masterN1[], const float q) {
  using namespace tledVectorArithmetic;
 
  const float safetyMargin = this->GetNarrowPhaseMaxDistance();

  float slaveX0[3], masterX0[3], d0[3];

  Interpolate(slaveX0, slave0T0, slave1T0, r);
  Interpolate(masterX0, master0T0, master1T0, q);
  ScalarDiv(p_n, Norm(Interpolate(p_n, masterN0, masterN1, q)));
  if (Dot(p_n, Sub(d0, slaveX0, masterX0)) >= -safetyMargin) {
    float *p_slaveX1, *p_masterX1, d1[3], slaveX1[3], masterX1[3];

    if (t_isMovingSlave) {
      p_slaveX1 = Interpolate(slaveX1, slave0T1, slave1T1, r);
    } else {
      p_slaveX1 = slaveX0;
    }

    if (t_isMovingMaster) {
      p_masterX1 = Interpolate(masterX1, master0T1, master1T1, q);    
    } else {
      p_masterX1 = masterX0;
    }

    r_g = Dot(p_n, Sub(d1, p_slaveX1, p_masterX1));
    if (r_g < this->GetNarrowPhaseMaxDistance()) {
      float slaveN[3];

      Interpolate(slaveN, slaveN0, slaveN1, r);

      return Dot(p_n, slaveN) < 0 && std::fabs(r_g) > Norm(d1)/2;
    }
  } 

  return false;
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::RunNarrowPhase() {
  assert(this->GetContactNodeBuffer().size() == 0 && this->GetContactEdgeBuffer().size() == 0);
  this->ProcessNodeFacetNarrowPhaseItems(this->GetNodeFacetNarrowPhasePairs().begin(), this->GetNodeFacetNarrowPhasePairs().end());
  this->GetNodeFacetNarrowPhasePairs().clear();
  this->ProcessEdgeEdgeNarrowPhaseItems(this->GetEdgeEdgeNarrowPhasePairs().begin(), this->GetEdgeEdgeNarrowPhasePairs().end());
  this->GetEdgeEdgeNarrowPhasePairs().clear();
}

template <class TMasterBVH, class TSlaveBVH, class TAPI>
void tledBVHTraverserImplCPU<TMasterBVH, TSlaveBVH, TAPI>::FindCollisions() {
  this->GetContactNodeBuffer().clear();
  this->GetContactEdgeBuffer().clear();

  assert(this->GetNodeFacetNarrowPhasePairs().size() == 0 && this->GetEdgeEdgeNarrowPhasePairs().size() == 0);
  Superclass::FindCollisions();
}
