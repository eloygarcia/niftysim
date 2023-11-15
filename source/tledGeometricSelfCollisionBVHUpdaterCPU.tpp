// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdaterCPU.tpp
// Purpose:
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
//
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
//
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBVH>
int tledGeometricSelfCollisionBVHUpdaterCPU<TBVH>::ComputeUpdateStatus(UpdateNodeInfo &r_updateInfo) {
  using namespace tledVectorArithmetic;

  const float *nodes = this->GetMesh().GetAllNodeCoordinates();
  const float uNonRigidThrsh = this->GetBVH().GetBV(r_updateInfo.BVIndex).SubtreeMinH/2*std::tan((BoundingVolume::GetVolinoThresholdAngle() - r_updateInfo.LastUpdateConeAngle)/2)/(1 + 1e-3f);
  const int nodesStartInd = r_updateInfo.NodeStartIndex;
  const int nodesEndInd = r_updateInfo.NodeEndIndex;

  float uNonRigidMax, uNonTrans;
  float m[9];

  uNonRigidMax = std::numeric_limits<float>::quiet_NaN();
  assert(!this->GetBVH().GetBV(r_updateInfo.BVIndex).HasDisconnectedGeometry());
  assert(this->GetBVH().GetBV(r_updateInfo.BVIndex).SubtreeMinH > 0);
  assert(this->GetBVH().GetBV(r_updateInfo.BVIndex).VolinoAngle > 0 && this->GetBVH().GetBV(r_updateInfo.BVIndex).VolinoAngle < BoundingVolume::GetVolinoThresholdAngle());

  this->ComputeBoundedNodeCentroid(r_updateInfo);
  Sub(r_updateInfo.Translation, r_updateInfo.CurrentCentroid, r_updateInfo.LastUpdateCentroid);

  assert(this->m_UpdateNodeCurrPositions.size() == this->m_UpdateNodeLastUpdateNodePositions.size());
  for (int localNodeInd = nodesStartInd; localNodeInd < nodesEndInd; localNodeInd++) {
    assert(localNodeInd < (int)this->GetUpdateNodeIndices().size());
    assert(this->GetUpdateNodeIndices()[localNodeInd] < (int)this->GetMesh().GetNumberOfNodes());
    Sub(&this->m_UpdateNodeCurrPositions[3*localNodeInd], nodes + 3*this->GetUpdateNodeIndices()[localNodeInd], r_updateInfo.CurrentCentroid);
  }

  uNonTrans = 0;
  for (float const *pc_currNode = &this->GetUpdateNodeCurrentPositions()[3*nodesStartInd], *pc_oldNode = &this->GetUpdateNodePreviousPositions()[3*nodesStartInd]; pc_currNode < &this->GetUpdateNodeCurrentPositions().front() + 3*nodesEndInd; pc_currNode += 3, pc_oldNode += 3) {
    float norm, diff[3];

    if (uNonTrans < (norm = Norm(Sub(diff, pc_currNode, pc_oldNode)))) uNonTrans = norm;
  }

  r_updateInfo.MaxNonRigidDisplacement = uNonTrans;
  if (Norm(r_updateInfo.Translation) + uNonTrans < this->GetBVH().GetBVMaxDisplacement() && uNonTrans < uNonRigidThrsh) {
    return Superclass::INACTIVE;
  }
  if (uNonTrans < this->GetBVH().GetBVMaxDisplacement() && uNonTrans < uNonRigidThrsh) {
    return Superclass::TRANSLATION;
  }

  std::fill(m, m + 9, 0.0f);
  for (int localNodeInd = nodesStartInd; localNodeInd < nodesEndInd; localNodeInd++) {
    const float *currPt = &this->GetUpdateNodeCurrentPositions()[3*localNodeInd];
    const float *oldPt = &this->GetUpdateNodePreviousPositions()[3*localNodeInd];

    for (int rInd = 0; rInd < 3; rInd++) for (int cInd = 0; cInd < 3; cInd++) m[cInd*3+rInd] += oldPt[rInd]*currPt[cInd];
  }

  {
    float u[9], v[9], s[3];

    SVD3x3(u, s, v, m);
    MatMultAB(u, 3, 3, MatTranspose(v, 3, 3), 3, 3, &r_updateInfo.Rotation[0][0]);
  }

  uNonRigidMax = 0;
  for (float const *pc_currNode = &this->GetUpdateNodeCurrentPositions()[3*nodesStartInd], *pc_oldNode = &this->GetUpdateNodePreviousPositions()[3*nodesStartInd]; pc_currNode < &this->GetUpdateNodeCurrentPositions()[3*nodesEndInd]; pc_currNode += 3, pc_oldNode += 3) {
    float uNR[3], normUNR;

    /*
     * Non-rigid max is set s.t. all non-rigid displacement exceeding BV margin leads to full update. Change this if
     * a better update criterion leading to fewer BV updates in self-collsion detection becomes available.
     */
    if ((normUNR = Norm(Sub(uNR, MatMul(uNR, r_updateInfo.Rotation, pc_oldNode), pc_currNode))) > uNonRigidMax) uNonRigidMax = normUNR;
  }

  r_updateInfo.MaxNonRigidDisplacement = uNonRigidMax;
  if (uNonRigidMax < this->GetBVH().GetBVMaxDisplacement() && uNonRigidMax < uNonRigidThrsh) {
    if (BoundingVolume::CanRotate && r_updateInfo.RigidUpdateCounter > 3) return 0;
    else return Superclass::RIGID_MOTION;
  } else {
    if (uNonRigidMax < uNonRigidThrsh) return 0;
    return Superclass::ACTIVE;
  }
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterCPU<TBVH>::UpdateBVH() {
  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++) {
    this->UpdateSubtreeUpdateStatus(*i_updateNode);
  }

  for (typename std::vector<UpdateNodeInfo>::iterator i_updateNode = this->GetUpdateNodes().begin(); i_updateNode < this->GetUpdateNodes().end(); i_updateNode++) {
    this->PerformSubtreeUpdate(*i_updateNode);
  }
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterCPU<TBVH>::UpdateSubtreeUpdateStatus(UpdateNodeInfo &r_updateNode) {
  if ((r_updateNode.Status & (Superclass::ACTIVE|Superclass::NONADJACENT)) == 0) {
    assert(!this->GetBVH().GetBV(r_updateNode.BVIndex).HasDisconnectedGeometry());
    assert(r_updateNode.LastUpdateConeAngle < BoundingVolume::GetVolinoThresholdAngle());
    switch (this->ComputeUpdateStatus(r_updateNode)) {
    case Superclass::ACTIVE: {
      tledLogDebugStream(tledHelper::Info() << "Immediate update required: subtree " << r_updateNode.BVIndex);
      r_updateNode.Status &= ~(Superclass::INACTIVE|Superclass::RIGID_MOTION|Superclass::TRANSLATION);
      r_updateNode.Status |= Superclass::ACTIVE;
      break;
    }

    case Superclass::INACTIVE: {
      tledLogDebugStream(tledHelper::Info() << "Deactivating subtree " << r_updateNode.BVIndex);
      r_updateNode.Status &= ~(Superclass::RIGID_MOTION|Superclass::TRANSLATION|Superclass::ACTIVE);
      r_updateNode.Status |= Superclass::INACTIVE;
      break;
    }

    case Superclass::RIGID_MOTION: {
      /*
       * There is no need to update the subtree, as self-collisions within it are impossible, but it cannot
       * be deactivated either, hence only check that INACTIVE flag isn't set, and set RIGID_MOTION flag.
       */
      tledLogDebugStream(tledHelper::Info() << "Rigid motion detected in subtree " << r_updateNode.BVIndex);
      r_updateNode.Status &= ~(Superclass::TRANSLATION|Superclass::INACTIVE|Superclass::ACTIVE);
      r_updateNode.Status |= Superclass::RIGID_MOTION;
      break;
    }

    case Superclass::TRANSLATION: {
      tledLogDebugStream(tledHelper::Info() << "Translation detected in subtree " << r_updateNode.BVIndex);
      r_updateNode.Status |= Superclass::TRANSLATION;
      r_updateNode.Status &= ~(Superclass::RIGID_MOTION|Superclass::ACTIVE|Superclass::INACTIVE);
      break;
    }

    default:
      r_updateNode.Status = 0;
    } /* if self-collision or non-adjacent node else ... */
  } else {
    assert(this->GetBVH().GetBV(r_updateNode.BVIndex).HasDisconnectedGeometry() || r_updateNode.LastUpdateConeAngle >= BoundingVolume::GetVolinoThresholdAngle());
  }
}

template <class TBVH>
void tledGeometricSelfCollisionBVHUpdaterCPU<TBVH>::PerformSubtreeUpdate(UpdateNodeInfo &r_updateNode) {
  using namespace tledVectorArithmetic;

  const int bvInd = r_updateNode.BVIndex;

  /*
   * Only update bottom up if update-level node too needs updating. Hence, bottom-up update is done on-demand too (see tledDynamic::UpdateBottomUp).
   */
  if ((r_updateNode.Status & Superclass::RIGID_MOTION) && BoundingVolume::CanRotate) {
    this->GetBVH().TransformSubTree(r_updateNode.BVIndex, r_updateNode.Rotation, r_updateNode.LastUpdateCentroid, r_updateNode.Translation);
    for (std::vector<float>::iterator i_oldPos = this->m_UpdateNodeCurrPositions.begin() + 3*r_updateNode.NodeStartIndex; i_oldPos < this->m_UpdateNodeCurrPositions.begin() + 3*r_updateNode.NodeEndIndex; i_oldPos += 3) {
      float tmp[3];

      std::copy(i_oldPos, i_oldPos + 3, tmp);
      MatMul(&*i_oldPos, r_updateNode.Rotation, tmp);
    }
    Add(r_updateNode.LastUpdateCentroid, r_updateNode.LastUpdateCentroid, r_updateNode.Translation);
    if (BoundingVolume::CanRotate) r_updateNode.RigidUpdateCounter += 1;
  } else if (r_updateNode.Status & Superclass::TRANSLATION) {
    this->GetBVH().TranslateSubTree(r_updateNode.BVIndex, r_updateNode.Translation);
    Add(r_updateNode.LastUpdateCentroid, r_updateNode.LastUpdateCentroid, r_updateNode.Translation);
  } else if ((r_updateNode.Status & Superclass::INACTIVE) == Superclass::INACTIVE) {
    this->GetBVH().GetBV(bvInd).UpdateCounter = this->GetBVH().GetUpdateCounter();
  } else {
    assert((r_updateNode.Status & (Superclass::ACTIVE|Superclass::NONADJACENT)) || r_updateNode.Status == 0 || ((r_updateNode.Status & Superclass::RIGID_MOTION) && !BoundingVolume::CanRotate));
    this->GetBVH().UpdateTopDownRecursive(bvInd);
    if (this->GetBVH().GetBV(bvInd).VolinoAngle < BoundingVolume::GetVolinoThresholdAngle()) r_updateNode.Status &= ~Superclass::ACTIVE;
    else r_updateNode.Status |= Superclass::ACTIVE;

    this->StoreUpdateNodePositions(r_updateNode);
    r_updateNode.LastUpdateConeAngle = this->GetBVH().GetBV(bvInd).VolinoAngle;
    assert(this->GetBVH().GetBV(bvInd).VolinoAngle < BoundingVolume::GetVolinoThresholdAngle() || (r_updateNode.Status & (Superclass::ACTIVE|Superclass::NONADJACENT)) != 0);

    if (BoundingVolume::CanRotate) r_updateNode.RigidUpdateCounter = 0;
  }

  this->GetBVH().UpdateBottomUpRecursive(bvInd);

#ifndef NDEBUG
  for (std::vector<int>::const_iterator ic_nodeInd = this->GetUpdateNodeIndices().begin() + r_updateNode.NodeStartIndex; ic_nodeInd < this->GetUpdateNodeIndices().begin() + r_updateNode.NodeEndIndex; ic_nodeInd++) {
    assert(this->GetBVH().IsPointInside(this->GetMesh().GetNodeCoordinates(*ic_nodeInd), bvInd));
  }
#endif

  assert(!this->GetBVH().DoesNeedUpdate(bvInd));
}
