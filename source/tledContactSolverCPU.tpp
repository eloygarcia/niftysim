// =========================================================================
// File:       tledContactSolverCPU.tpp
// Purpose:    Contact solver interface and general purpose functions
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

template <class TContactMesh, class TContactSolverInterface>
tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::~tledContactSolverImplCPU() {
  if (mp_BVHTraverser != NULL) delete mp_BVHTraverser;
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::InitBVHTraverser() {
  this->GetBVHTraverser().SetNarrowPhaseMaxDistance(this->GetNodeCloseDistance());
  this->GetBVHTraverser().Init(this->GetManager());
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::Init() {
  Superclass::Init();
  this->SetContactForceBufferSize(this->GetMesh().GetNumberOfNodes(), this->DoFriction());  
  mp_BVHTraverser = this->InstantiateBVHTraverser();
  this->InitBVHTraverser();
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_ApplyForceToNode(const int nodeIndex, const float f[]) {
  assert(this->GetMesh().GetNumberOfNodes() > nodeIndex && nodeIndex >= 0 && nodeIndex < this->GetContactForceBufferSize());  
  assert(this->GetNormalResponse(nodeIndex).MaxProjection < 0.f);
  tledVectorArithmetic::Add(this->GetNormalResponse(nodeIndex).AccumulatedResponse, this->GetNormalResponse(nodeIndex).AccumulatedResponse, f);    
  this->MarkActiveNode(nodeIndex);  
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyForceToSlaveNode(NodeFacetConstraintItem &r_constraint, const float df[]) {
  assert(tledVectorArithmetic::Dot(r_constraint.Normal, df) <= 0);
  std::copy(df, df + 3, r_constraint.ContactForces[0]);
  _ApplyForceToNode(r_constraint.ContactNodeIndices[0], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyForceToMasterFacetNode(NodeFacetConstraintItem &r_constraint, const float df[], const int vertexIndex) {
  assert(tledVectorArithmetic::Dot(r_constraint.Normal, df) >= 0 || r_constraint.ShapeValues[vertexIndex] < 0);
  std::copy(df, df + 3, r_constraint.ContactForces[1+vertexIndex]);
  _ApplyForceToNode(r_constraint.ContactNodeIndices[1+vertexIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyForceToSlaveEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex) {
  assert(tledVectorArithmetic::Dot(r_constraint.Normal, df) <= 0 || r_constraint.SlaveShapeValues[nodeIndex] < 0);
  std::copy(df, df + 3, r_constraint.ContactForces[nodeIndex]);
  _ApplyForceToNode(r_constraint.SlaveNodeIndices[nodeIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyForceToMasterEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex) {
  assert(tledVectorArithmetic::Dot(r_constraint.Normal, df) >= 0 || r_constraint.MasterShapeValues[nodeIndex] < 0);
  std::copy(df, df + 3, r_constraint.ContactForces[2+nodeIndex]);
  _ApplyForceToNode(r_constraint.MasterNodeIndices[nodeIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyContactForces(float *p_f, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints) {
  m_ForceResetIndices = tledHelper::MakeSortedUnique(m_ForceResetIndices);
  _RunMaxResponseConsolidationAndApply<t_doMaster, t_doSlave>(p_f, this->GetAllNormalResponses(), nodeFacetConstraints, edgeEdgeConstraints);  
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_ApplyFrictionToNode(const int nodeIndex, const float f[]) {
  assert(this->GetMesh().GetNumberOfNodes() > nodeIndex && nodeIndex >= 0);  
  assert(std::find(m_ForceResetIndices.begin(), m_ForceResetIndices.end(), nodeIndex) < m_ForceResetIndices.end());
  assert(this->GetFrictionResponse(nodeIndex).MaxProjection < 0.f);
  tledVectorArithmetic::Add(this->GetFrictionResponse(nodeIndex).AccumulatedResponse, this->GetFrictionResponse(nodeIndex).AccumulatedResponse, f);  
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyFrictionToSlaveNode(NodeFacetConstraintItem &r_constraint, const float df[]) {
  std::copy(df, df + 3, r_constraint.ContactForces[0]);
  _ApplyFrictionToNode(r_constraint.ContactNodeIndices[0], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyFrictionToMasterFacetNode(NodeFacetConstraintItem &r_constraint, const float df[], const int vertexIndex) {
  std::copy(df, df + 3, r_constraint.ContactForces[1+vertexIndex]);
  _ApplyFrictionToNode(r_constraint.ContactNodeIndices[1+vertexIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyFrictionToSlaveEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex) {
  std::copy(df, df + 3, r_constraint.ContactForces[nodeIndex]);
  _ApplyFrictionToNode(r_constraint.SlaveNodeIndices[nodeIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyFrictionToMasterEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex) {
  std::copy(df, df + 3, r_constraint.ContactForces[2+nodeIndex]);
  _ApplyFrictionToNode(r_constraint.MasterNodeIndices[nodeIndex], df);
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ApplyFrictionForces(float *p_f, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints) {
  _RunMaxResponseConsolidationAndApply<t_doMaster, t_doSlave>(p_f, this->GetAllFrictionResponses(), nodeFacetConstraints, edgeEdgeConstraints);  
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_RunMaxResponseConsolidationAndApply(float *p_f, ContactResponse *p_responses, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints) {
  using namespace tledVectorArithmetic;

  for (std::vector<int>::const_iterator ic_n = m_ForceResetIndices.begin(); ic_n < m_ForceResetIndices.end(); ic_n++) {    
    ContactResponse &r_resp = p_responses[*ic_n];
    float accMag = Norm(r_resp.AccumulatedResponse);

    if (accMag > 0) ScalarDiv(r_resp.AccumulatedResponse, accMag);
    assert(r_resp.MaxProjection < 0);
    assert(Norm(r_resp.AccumulatedResponse) == 0.f || std::fabs(1.f - Norm(r_resp.AccumulatedResponse)) < 1e-4f);
  }

  for (typename std::vector<NodeFacetConstraintItem>::const_iterator ic_c = nodeFacetConstraints.begin(); ic_c < nodeFacetConstraints.end(); ic_c++) {
    const int startInd = t_doSlave? 0 : 1;
    const int endInd = t_doMaster? 1 + Facet::NumberOfVertices : 1;

    for (int n = startInd; n < endInd; n++) {
      ContactResponse &r_resp = p_responses[ic_c->ContactNodeIndices[n]];
      float &r_maxProj = r_resp.MaxProjection;

      assert(std::find(m_ForceResetIndices.begin(), m_ForceResetIndices.end(), ic_c->ContactNodeIndices[n]) != m_ForceResetIndices.end());
      r_maxProj = std::max(r_maxProj, Dot(r_resp.AccumulatedResponse, ic_c->ContactForces[n]));
    }
  }

  for (typename std::vector<EdgeEdgeConstraintItem>::const_iterator ic_c = edgeEdgeConstraints.begin(); ic_c < edgeEdgeConstraints.end(); ic_c++) {
    if (t_doSlave) {
      for (int n = 0; n < 2; n++) {
	ContactResponse &r_resp = p_responses[ic_c->SlaveNodeIndices[n]];
	float &r_maxProj = r_resp.MaxProjection;

	assert(std::find(m_ForceResetIndices.begin(), m_ForceResetIndices.end(), ic_c->SlaveNodeIndices[n]) != m_ForceResetIndices.end());
	r_maxProj = std::max(r_maxProj, Dot(r_resp.AccumulatedResponse, ic_c->ContactForces[n]));
      }
    }

    if (t_doMaster) {
      for (int n = 0; n < 2; n++) {
	ContactResponse &r_resp = p_responses[ic_c->MasterNodeIndices[n]];
	float &r_maxProj = r_resp.MaxProjection;

	assert(std::find(m_ForceResetIndices.begin(), m_ForceResetIndices.end(), ic_c->MasterNodeIndices[n]) != m_ForceResetIndices.end());
	r_maxProj = std::max(r_maxProj, Dot(r_resp.AccumulatedResponse, ic_c->ContactForces[n+2]));
      }
    }
  }

  for (std::vector<int>::const_iterator ic_n = m_ForceResetIndices.begin(); ic_n < m_ForceResetIndices.end(); ic_n++) {
    ContactResponse &r_resp = p_responses[*ic_n];
    float *p_dst = p_f + 3*this->GetMesh().MapSurface2VolumeNode(*ic_n);

    assert(r_resp.MaxProjection >= 0);
    ScalarMul(r_resp.AccumulatedResponse, r_resp.MaxProjection);
    tledLogDebugStream(tledHelper::Info() << "Applying consolidated contact force to node " << *ic_n << "/" << this->GetMesh().MapSurface2VolumeNode(*ic_n) << ": " << r_resp.AccumulatedResponse[0] << ", " << r_resp.AccumulatedResponse[1] << ", " << r_resp.AccumulatedResponse[2]);
    Add(p_dst, p_dst, r_resp.AccumulatedResponse);
  }
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_ApplyResponsesOnly(float *p_f, ContactResponse *p_responses, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints) {
  using namespace tledVectorArithmetic;

  for (std::vector<int>::const_iterator ic_i = m_ForceResetIndices.begin(); ic_i < m_ForceResetIndices.end(); ic_i++) {
    ContactResponse &r_resp = p_responses[*ic_i];

    Add(p_f + 3*this->GetMesh().MapSurface2VolumeNode(*ic_i), p_f + 3*this->GetMesh().MapSurface2VolumeNode(*ic_i), r_resp.AccumulatedResponse);
  }
}

template <class TContactMesh, class TContactSolverInterface>
template <const int t_numSegNodes>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeMass(const int segment[], const float shapeVals[]) const {
  float mass;

  mass = this->GetMesh().GetSurfaceNodeMass(segment[0])*shapeVals[0];
  for (int nInd = 1; nInd < t_numSegNodes; nInd++) mass += this->GetMesh().GetSurfaceNodeMass(segment[nInd])*shapeVals[nInd];

  return mass;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeSlaveBeta(const EdgeEdgeConstraintItem &ci) const {
  if (t_doMaster != t_doSlave) return t_doSlave? 1 : 0;
  else {
    const float mMaster = this->ComputeMass<2>(ci.MasterNodeIndices, ci.MasterShapeValues);
    const float mSlave = this->ComputeMass<2>(ci.SlaveNodeIndices, ci.SlaveShapeValues);

    assert(mMaster > 0 && mSlave > 0);

    return mMaster/(mMaster + mSlave);
  }
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeEdgeEdgeResponse(EdgeEdgeConstraintItem &r_ci) {
  using namespace tledVectorArithmetic;

  if (r_ci.GapValue < 0) {
    const float slaveBeta = ComputeSlaveBeta<t_doMaster, t_doSlave>(r_ci);

    float lambda, shapeSqrSum, df[3];

#ifndef NDEBUG
    if (t_doSlave && t_doMaster) {
      tledLogDebugStream(tledHelper::Info() << "Detected edge-edge collision (" << r_ci.MasterNodeIndices[0] << ", " << r_ci.MasterNodeIndices[1] 
			 << ") - (" << r_ci.SlaveNodeIndices[0] << ", " << r_ci.SlaveNodeIndices[1] 
			 << "): r = " << r_ci.SlaveShapeValues[1] << ", q = " << r_ci.MasterShapeValues[1] << ", gap = " << r_ci.GapValue);
    } else if (t_doSlave) {
      tledLogDebugStream(tledHelper::Info() << "Detected edge-edge collision (" << r_ci.SlaveNodeIndices[0] << ", " << r_ci.SlaveNodeIndices[1] 
			 << "): r = " << r_ci.SlaveShapeValues[1] << ", q = " << r_ci.MasterShapeValues[1] << ", gap = " << r_ci.GapValue);
    } else {
      tledLogDebugStream(tledHelper::Info() << "Detected edge-edge collision (" << r_ci.MasterNodeIndices[0] << ", " << r_ci.MasterNodeIndices[1] 
			 << "): r = " << r_ci.SlaveShapeValues[1] << ", q = " << r_ci.MasterShapeValues[1] << ", gap = " << r_ci.GapValue);
    }
#endif

    lambda = r_ci.GapValue/(this->GetDt()*this->GetDt());

    if (t_doSlave) {
      shapeSqrSum = tledContactSolverImplCPU::ComputeShapeSquareSum(r_ci.SlaveShapeValues, 2);
      for (int nInd = 0; nInd < 2; nInd++) {
	assert(r_ci.SlaveNodeIndices[nInd] >= 0 && r_ci.SlaveNodeIndices[nInd] < this->GetMesh().GetNumberOfNodes());
	assert(this->GetContactForceBufferSize() >= this->GetMesh().GetNumberOfNodes());
	ScalarMul(df, r_ci.Normal, slaveBeta*lambda*this->GetMesh().GetSurfaceNodeMass(r_ci.SlaveNodeIndices[nInd])*r_ci.SlaveShapeValues[nInd]/shapeSqrSum);
	this->ApplyForceToSlaveEdgeNode(r_ci, df, nInd);
	tledLogDebugStream(tledHelper::Info() << "Applying pred. force (edge-edge) to " << r_ci.SlaveNodeIndices[nInd] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.SlaveNodeIndices[nInd]) << " (slave): " << df[0] << ", " << df[1] << ", " << df[2]);
      }
    }

    if (t_doMaster) {
      shapeSqrSum = -tledContactSolverImplCPU::ComputeShapeSquareSum(r_ci.MasterShapeValues, 2);
      for (int nInd = 0; nInd < 2; nInd++) {
	assert(r_ci.MasterNodeIndices[nInd] >= 0 && r_ci.MasterNodeIndices[nInd] < this->GetMesh().GetNumberOfNodes());
	assert(this->GetContactForceBufferSize() >= this->GetMesh().GetNumberOfNodes());
	ScalarMul(df, r_ci.Normal, (1.0f - slaveBeta)*lambda*this->GetMesh().GetSurfaceNodeMass(r_ci.MasterNodeIndices[nInd])*r_ci.MasterShapeValues[nInd]/shapeSqrSum);
	this->ApplyForceToMasterEdgeNode(r_ci, df, nInd);
	tledLogDebugStream(tledHelper::Info() << "Applying pred. force (edge-edge) to " << r_ci.MasterNodeIndices[nInd] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.MasterNodeIndices[nInd]) << " (master): " << df[0] << ", " << df[1] << ", " << df[2]);
      }
    }
    
    r_ci.GapRateWeight = std::numeric_limits<float>::quiet_NaN();
    
    return true;
  } else return false;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeSlaveBeta(const NodeFacetConstraintItem &ci) const {
  if (t_doMaster != t_doSlave) return t_doSlave? 1 : 0;
  else {
    const float mMaster = this->ComputeMass<Facet::NumberOfVertices>(ci.ContactNodeIndices + 1, ci.ShapeValues);
    const float mSlave = this->GetMesh().GetSurfaceNodeMass(ci.ContactNodeIndices[0]);

    assert(mMaster > 0 && mSlave > 0);

    return mMaster/(mMaster + mSlave);
  }
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeNodeFacetResponse(NodeFacetConstraintItem &r_ci) {
  using namespace tledVectorArithmetic;

  if (r_ci.GapValue < 0) {
    const float slaveBeta = ComputeSlaveBeta<t_doMaster, t_doSlave>(r_ci);

    float lambda, df[3];

    lambda = r_ci.GapValue/(this->GetDt()*this->GetDt());

    if (t_doSlave) {
      ScalarMul(df, r_ci.Normal, slaveBeta*lambda*this->GetMesh().GetSurfaceNodeMass(r_ci.ContactNodeIndices[0]));
      this->ApplyForceToSlaveNode(r_ci, df);
      tledLogDebugStream(tledHelper::Info() << "Applying pred. force (node-facet) to " << r_ci.ContactNodeIndices[0] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.ContactNodeIndices[0]) << " (slave): " << df[0] << ", " << df[1] << ", " << df[2]);
    }

    if (t_doMaster) {
      lambda /= this->ComputeShapeSquareSum(r_ci.ShapeValues, Facet::NumberOfVertices);
      for (int tInd = 0; tInd < Facet::NumberOfVertices; tInd++) {
	const int nodeInd = r_ci.ContactNodeIndices[tInd+1];
      
	ScalarMul(df, r_ci.Normal, -(1.0f - slaveBeta)*this->GetMesh().GetSurfaceNodeMass(nodeInd)*r_ci.ShapeValues[tInd]*lambda);
	this->ApplyForceToMasterFacetNode(r_ci, df, tInd);
	tledLogDebugStream(tledHelper::Info() << "Applying pred. force (node-facet) to " << nodeInd << "/" << this->GetMesh().MapSurface2VolumeNode(nodeInd) << " (master): " << df[0] << ", " << df[1] << ", " << df[2]);
      }  
    }

    r_ci.GapRateWeight = std::numeric_limits<float>::quiet_NaN();

    return true;
  } else return false;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeNodeFacetFrictionResponse(std::vector<NodeFacetConstraintItem> &r_constraints, const float uNexts[], const float uCurrs[], const float frictionCoefficient) {
  using namespace tledVectorArithmetic;

  for (typename std::vector<NodeFacetConstraintItem>::iterator i_c = r_constraints.begin(); i_c < r_constraints.end(); i_c++) {
    float tmp[3], df[3], vRel[3], lambda, magRel, magN;

    /*
     * - Compute friction force based on local v. and consolidated force:
     * -- Compute rel. node-facet velocity / edge-edge v.
     * -- Compute required "stopping force"
     * -- Compare to Coulomb threshold force
     * -- Adjust slip as required
     */      
    if (t_doMaster && t_doSlave) this->ComputeRelativeNodeFacetVelocityBoth(vRel, *i_c, uNexts, uCurrs);
    else if (t_doMaster) this->ComputeRelativeNodeFacetVelocitySlave(vRel, *i_c, uNexts, uCurrs);
    else this->ComputeRelativeNodeFacetVelocityMaster(vRel, *i_c, uNexts, uCurrs);

    Sub(vRel, vRel, ScalarMul(tmp, i_c->Normal, Dot(vRel, i_c->Normal)));

    magRel = Norm(vRel);
    lambda = _ComputeStoppingNodeFacetForce<t_doMaster, t_doSlave>(magRel, *i_c);
    if (t_doSlave) magN = Norm(i_c->ContactForces[0]);
    else {
      float fN[3], fNodeN[3];

      std::fill(fN, fN + 3, 0.f);
      for (int vInd = 0; vInd < Facet::NumberOfVertices; vInd++) {
	ScalarMul(fNodeN, i_c->ContactForces[1+vInd], i_c->ShapeValues[vInd]);
	Add(fN, fN, fNodeN);
      }
      magN = Norm(fN);
    }
    if (i_c->GapRateWeight == i_c->GapRateWeight) magN /= i_c->GapRateWeight;

    if (tledHelper::IsNumericallyZero(magRel, frictionCoefficient*magN)) continue;
    ScalarMul(df, vRel, std::min(frictionCoefficient*magN, lambda)/magRel);    

    if (t_doSlave) {
      assert(i_c->ContactNodeIndices[0] >= 0 && i_c->ContactNodeIndices[0] < this->GetMesh().GetNumberOfNodes());
      tledLogDebugStream(tledHelper::Info() << "Applying friction force (node, slave) " << df[0] << ", " << df[1] << ", " << df[2] << " to " << this->GetMesh().MapSurface2VolumeNode(i_c->ContactNodeIndices[0]));
      this->ApplyFrictionToSlaveNode(*i_c, df);
    }
      
    if (t_doMaster) {
      ScalarDiv(df, -this->ComputeShapeSquareSum(i_c->ShapeValues, Facet::NumberOfVertices));
      for (int vInd = 0; vInd < Facet::NumberOfVertices; vInd++) {
	ScalarMul(tmp, df, i_c->ShapeValues[vInd]);
	tledLogDebugStream(tledHelper::Info() << "Applying friction force (facet, master) " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << " to " << this->GetMesh().MapSurface2VolumeNode(i_c->ContactNodeIndices[1+vInd]));
	this->ApplyFrictionToMasterFacetNode(*i_c, tmp, vInd);
      }
    }
  }
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeEdgeEdgeFrictionResponse(std::vector<EdgeEdgeConstraintItem> &r_constraints, const float uNexts[], const float uCurrs[], const float frictionCoefficient) {
  using namespace tledVectorArithmetic;

  for (typename std::vector<EdgeEdgeConstraintItem>::iterator i_c = r_constraints.begin(); i_c < r_constraints.end(); i_c++) {
    float tmp[3], df[3], vRel[3], lambda, magRel, magN;

    /*
     * - Compute tangential speed by removing normal component from relative velocity
     * - Compute required "stopping force"
     * - Compare to Coulomb threshold force
     * - Adjust slip as required
     */      
    if (t_doMaster && t_doSlave) this->ComputeRelativeEdgeEdgeVelocityBoth(vRel, *i_c, uNexts, uCurrs);
    else if (t_doMaster) this->ComputeRelativeEdgeEdgeVelocitySlave(vRel, *i_c, uNexts, uCurrs);
    else this->ComputeRelativeEdgeEdgeVelocityMaster(vRel, *i_c, uNexts, uCurrs);

    Sub(vRel, vRel, ScalarMul(tmp, i_c->Normal, Dot(vRel, i_c->Normal)));

    magRel = Norm(vRel);
    lambda = _ComputeStoppingEdgeEdgeForce<t_doMaster, t_doSlave>(magRel, *i_c);
    if (t_doSlave) {
      float fNormal[3];

      std::fill(fNormal, fNormal + 3, 0.f);
      for (int i = 0; i < 2; i++) Add(fNormal, fNormal, ScalarMul(tmp, i_c->ContactForces[i], i_c->SlaveShapeValues[i]));
      magN = Norm(fNormal);
    } else {
      float fNormal[3];

      std::fill(fNormal, fNormal + 3, 0.f);
      for (int i = 0; i < 2; i++) Add(fNormal, fNormal, ScalarMul(tmp, i_c->ContactForces[2+i], i_c->MasterShapeValues[i]));
      magN = Norm(fNormal);
    }
    if (i_c->GapRateWeight == i_c->GapRateWeight) magN /= i_c->GapRateWeight;

    if (tledHelper::IsNumericallyZero(magRel, frictionCoefficient*128*magN)) continue;
    ScalarMul(df, vRel, std::min(frictionCoefficient*magN, lambda)/magRel);    

    if (t_doSlave) {
      const float gammaBase = this->ComputeShapeSquareSum(i_c->SlaveShapeValues, 2);

      for (int i = 0; i < 2; i++) {
	ScalarMul(tmp, df, i_c->SlaveShapeValues[i]/gammaBase);
	this->ApplyFrictionToSlaveEdgeNode(*i_c, tmp, i);
	tledLogDebugStream(tledHelper::Info() << "Applying friction force (edge, slave) " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << " to " << this->GetMesh().MapSurface2VolumeNode(i_c->SlaveNodeIndices[i]));
      }
    }
      
    if (t_doMaster) {
      const float gammaBase = this->ComputeShapeSquareSum(i_c->MasterShapeValues, 2);

      for (int i = 0; i < 2; i++) {
	ScalarMul(tmp, df, -i_c->MasterShapeValues[i]/gammaBase);
	this->ApplyFrictionToMasterEdgeNode(*i_c, tmp, i);
	tledLogDebugStream(tledHelper::Info() << "Applying friction force (edge, master) " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << " to " << this->GetMesh().MapSurface2VolumeNode(i_c->MasterNodeIndices[i]));
      }
    }
  }
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeRelativeNodeFacetVelocity(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  using namespace tledVectorArithmetic;

  if (t_doSlave) Sub(p_v, uNexts + 3*this->GetMesh().MapSurface2VolumeNode(ci.ContactNodeIndices[0]), uCurrs + 3*this->GetMesh().MapSurface2VolumeNode(ci.ContactNodeIndices[0]));
  else std::fill(p_v, p_v + 3, 0.0f);
  
  if (t_doMaster) for (int m = 0; m < Facet::NumberOfVertices; m++) {
      float vN[3];

      ScalarMul(Sub(vN, uNexts + 3*this->GetMesh().MapSurface2VolumeNode(ci.ContactNodeIndices[1+m]), uCurrs + 3*this->GetMesh().MapSurface2VolumeNode(ci.ContactNodeIndices[1+m])), ci.ShapeValues[m]);
      Sub(p_v, p_v, vN);
    }
  ScalarDiv(p_v, this->GetDt());
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeRelativeEdgeEdgeVelocity(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const {
  using namespace tledVectorArithmetic;

  float vN[3];

  std::fill(p_v, p_v + 3, 0.0f);
  if (t_doSlave) for (int i = 0; i < 2; i++) {
      ScalarMul(Sub(vN, uNexts + 3*this->GetMesh().MapSurface2VolumeNode(ci.SlaveNodeIndices[i]), uCurrs + 3*this->GetMesh().MapSurface2VolumeNode(ci.SlaveNodeIndices[i])), ci.SlaveShapeValues[i]);
      Add(p_v, p_v, vN);
    }
  
  if (t_doMaster) for (int i = 0; i < 2; i++) {
      ScalarMul(Sub(vN, uNexts + 3*this->GetMesh().MapSurface2VolumeNode(ci.MasterNodeIndices[i]), uCurrs + 3*this->GetMesh().MapSurface2VolumeNode(ci.MasterNodeIndices[i])), ci.MasterShapeValues[i]);
      Sub(p_v, p_v, vN);
    }
  ScalarDiv(p_v, this->GetDt());
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_ComputeStoppingNodeFacetForce(const float vMag, const NodeFacetConstraintItem &ci) const {
  float lambda = 0;

  if (t_doMaster) {
    for (int vInd = 0; vInd < Facet::NumberOfVertices; vInd++) {
      lambda += ci.ShapeValues[vInd]*ci.ShapeValues[vInd]/this->GetMesh().GetSurfaceNodeMass(ci.ContactNodeIndices[1+vInd]);
    }
    lambda /= this->ComputeShapeSquareSum(ci.ShapeValues, Facet::NumberOfVertices);
  } 
  
  if (t_doSlave) {
    lambda += 1/this->GetMesh().GetSurfaceNodeMass(ci.ContactNodeIndices[0]);
  }
      
  lambda = std::fabs(vMag/(this->GetDt()*lambda));  

  return lambda;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::_ComputeStoppingEdgeEdgeForce(const float vMag, const EdgeEdgeConstraintItem &ci) const {
  float lambda = 0.f;

  if (t_doMaster) {
    float masterDenom = 0.f;

    for (int i = 0; i < 2; i++) {
      masterDenom += ci.MasterShapeValues[i]*ci.MasterShapeValues[i]/this->GetMesh().GetSurfaceNodeMass(ci.MasterNodeIndices[i]);
    }
    lambda += masterDenom/this->ComputeShapeSquareSum(ci.MasterShapeValues, 2);
  } 

  if (t_doSlave) {
    float slaveDenom = 0.f;
    
    for (int i = 0; i < 2; i++) {
      slaveDenom += ci.SlaveShapeValues[i]*ci.SlaveShapeValues[i]/this->GetMesh().GetSurfaceNodeMass(ci.SlaveNodeIndices[i]);
    }
    lambda += slaveDenom/this->ComputeShapeSquareSum(ci.SlaveShapeValues, 2);
  }
  
  lambda = std::fabs(vMag/(this->GetDt()*lambda));

  return lambda;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeNodeFacetRateResponse(NodeFacetConstraintItem &r_ci, const float uNexts[], const float uCurrs[]) {
  using namespace tledVectorArithmetic;

  if (r_ci.GapValue < this->GetNodeCloseDistance()) {
    float vRel[3], vMag;

    if (t_doMaster && t_doSlave) this->ComputeRelativeNodeFacetVelocityBoth(vRel, r_ci, uNexts, uCurrs);
    else if (t_doMaster) this->ComputeRelativeNodeFacetVelocitySlave(vRel, r_ci, uNexts, uCurrs);
    else this->ComputeRelativeNodeFacetVelocityMaster(vRel, r_ci, uNexts, uCurrs);

    vMag = Dot(r_ci.Normal, vRel);
    if (vMag >= 0) return false;
    else {
      float df[3];
      
      r_ci.GapRateWeight = 1 - r_ci.GapValue/this->GetNodeCloseDistance();      
      ScalarMul(df, r_ci.Normal, -_ComputeStoppingNodeFacetForce<t_doMaster, t_doSlave>(vMag, r_ci)*r_ci.GapRateWeight);
      assert(Dot(df, r_ci.Normal) < 1e-4f*Norm(df));
      if (t_doSlave) {	
	this->ApplyForceToSlaveNode(r_ci, df);
	tledLogDebugStream(tledHelper::Info() << "Applying node-facet rate response (slave) to " << r_ci.ContactNodeIndices[0] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.ContactNodeIndices[0]) << ": " << df[0] << ", " << df[1] << ", " << df[2]);
      }
      
      if (t_doMaster) {
	float tmp[3];
	
	ScalarDiv(df, -this->ComputeShapeSquareSum(r_ci.ShapeValues, Facet::NumberOfVertices));
	for (int vInd = 0; vInd < Facet::NumberOfVertices; vInd++) {
	  ScalarMul(tmp, df, r_ci.ShapeValues[vInd]);
	  this->ApplyForceToMasterFacetNode(r_ci, tmp, vInd);
	  tledLogDebugStream(tledHelper::Info() << "Applying node-facet rate response (master) to " << r_ci.ContactNodeIndices[1+vInd] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.ContactNodeIndices[1+vInd]) << ": " << tmp[0] << ", " << tmp[1] << ", " << tmp[2]);
	}
      }      

      return true;
    }
  } else return false;
}

template <class TContactMesh, class TContactSolverInterface>
template <const bool t_doMaster, const bool t_doSlave>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeEdgeEdgeRateResponse(EdgeEdgeConstraintItem &r_ci, const float uNexts[], const float uCurrs[]) {
  using namespace tledVectorArithmetic;

  assert(r_ci.GapValue >= 0.f);
  if (r_ci.GapValue < this->GetNodeCloseDistance()) {
    float vRel[3], vMag;

    if (t_doMaster && t_doSlave) this->ComputeRelativeEdgeEdgeVelocityBoth(vRel, r_ci, uNexts, uCurrs);
    else if (t_doMaster) this->ComputeRelativeEdgeEdgeVelocitySlave(vRel, r_ci, uNexts, uCurrs);
    else this->ComputeRelativeEdgeEdgeVelocityMaster(vRel, r_ci, uNexts, uCurrs);

    vMag = Dot(r_ci.Normal, vRel);
    if (vMag >= 0) return false;
    else {
      float df[3], ndf[3];
     
      r_ci.GapRateWeight = 1 - r_ci.GapValue/this->GetNodeCloseDistance();
      ScalarMul(df, r_ci.Normal, -_ComputeStoppingEdgeEdgeForce<t_doMaster, t_doSlave>(vMag, r_ci)*r_ci.GapRateWeight);
      assert(Dot(df, r_ci.Normal) < 1e-4f*Norm(df));
      if (t_doSlave) {
	const float gammaBase = this->ComputeShapeSquareSum(r_ci.SlaveShapeValues, 2);

	for (int i = 0; i < 2; i++) {
	  ScalarMul(ndf, df, r_ci.SlaveShapeValues[i]/gammaBase);
	  this->ApplyForceToSlaveEdgeNode(r_ci, ndf, i);
	  tledLogDebugStream(tledHelper::Info() << "Applying edge-edge rate response (slave) to " << r_ci.SlaveNodeIndices[i] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.SlaveNodeIndices[i]) << ": " << ndf[0] << ", " << ndf[1] << ", " << ndf[2]);
	  assert(r_ci.SlaveShapeValues[i]*Norm(ndf) <= 1.001f*Norm(df));
	}
      } 

      if (t_doMaster) {
	const float gammaBase = this->ComputeShapeSquareSum(r_ci.MasterShapeValues, 2);

	for (int i = 0; i < 2; i++) {
	  ScalarMul(ndf, df, -r_ci.MasterShapeValues[i]/gammaBase);
	  this->ApplyForceToMasterEdgeNode(r_ci, ndf, i);
	  tledLogDebugStream(tledHelper::Info() << "Applying edge-edge rate response (master) to " << r_ci.MasterNodeIndices[i] << "/" << this->GetMesh().MapSurface2VolumeNode(r_ci.MasterNodeIndices[i]) << ": " << ndf[0] << ", " << ndf[1] << ", " << ndf[2]);
	  assert(r_ci.MasterShapeValues[i]*Norm(ndf) <= 1.001f*Norm(df));
	}
      }
    }    

    return true;
  }

  return false;
}

template <class TContactMesh, class TContactSolverInterface>
void tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ResetContactForces() {
  if (!this->DoFriction()) {
    for (std::vector<int>::const_iterator ic_n = m_ForceResetIndices.begin(); ic_n < m_ForceResetIndices.end(); ic_n++) {
      this->ResetNormalResponse(*ic_n);
    }
  } else {
    for (std::vector<int>::const_iterator ic_n = m_ForceResetIndices.begin(); ic_n < m_ForceResetIndices.end(); ic_n++) {
      this->ResetNormalResponse(*ic_n);
      this->ResetFrictionResponse(*ic_n);
    }
  }
  m_ForceResetIndices.clear();

#ifndef NDEBUG
  for (int n = 0; n < this->GetMesh().GetNumberOfNodes(); n++) {    
    assert(tledVectorArithmetic::Norm(this->GetNormalResponse(n).AccumulatedResponse) == 0.f);
    assert(this->GetNormalResponse(n).MaxProjection < 0);
    if (this->DoFriction()) {
      assert(tledVectorArithmetic::Norm(this->GetFrictionResponse(n).AccumulatedResponse) == 0.f);
      assert(this->GetFrictionResponse(n).MaxProjection < 0);
    }
  }
#endif
}

template <class TContactMesh, class TContactSolverInterface>
float tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeShapeSquareSum(const float shapeValues[], const int numNodes) const {
  float shapeSqrSum;
  float const *pc_shapeVal;

  shapeSqrSum = shapeValues[0]*shapeValues[0];
  for (pc_shapeVal = shapeValues + 1; pc_shapeVal < shapeValues + numNodes; pc_shapeVal++) shapeSqrSum += (*pc_shapeVal)*(*pc_shapeVal);

  return shapeSqrSum;
}

template <class TContactMesh, class TContactSolverInterface>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::FindCollisions(const float uNexts[], const float uCurrs[]) {
  this->GetBVHTraverser().SetDoMaster(this->DoMaster());  
  this->GetBVHTraverser().FindCollisions();

  return this->GetBVHTraverser().GetSlaveEdgeIndices().size() + this->GetBVHTraverser().GetSlaveNodeIndices().size() > 0;
}

template <class TContactMesh, class TContactSolverInterface>
bool tledContactSolverImplCPU<TContactMesh, TContactSolverInterface>::ComputeContactResponses(float *p_f, const float uNexts[], const float uCurrs[]) {
  this->ToggleDoMaster();
  if (this->FindCollisions(uNexts, uCurrs)) {
    const bool haveContacts = this->ComputeContactForces(p_f, uNexts, uCurrs);

    if (m_ForceResetIndices.size() > 0) this->ResetContactForces();

    return haveContacts;
  } else return false;
}
