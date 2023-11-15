// =========================================================================
// File:       tledContactSolverCPU.h
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

#ifndef tledContactSolverCPU_H
#define tledContactSolverCPU_H

#include "tledContactSurface.h"
#include "tledContactSolver.h"
#include "tledContactManager.h"
#include "tledBVHTraverserCPU.h"
#include "tledHelper.h"

#include <vector>
#include <cmath>
#include <cassert>

/**
 * \brief Interface for CPU contact solvers 
 * \ingroup contact
 */
class tledContactSolverCPU {
  /**
   * \name Constraints
   * @{
   */
public:
  /** Base class for contact constraints */
  template <const int t_numNodes>
  struct ContactConstraint {
    /** Unconsolidated contact force. */
    float ContactForces[t_numNodes][3];

    /** Predictor gap */
    float GapValue;

    /** Target-surface normal */
    float Normal[3];

    /** Distance-based weight used with gap-rate constraints. Set NaN if not used (i.e. in penetration constraints). */
    float GapRateWeight;
  };

  /** Node-facet contact constraints */
  template <const int t_numFacetVertices>
  struct NodeFacetConstraintItem : public ContactConstraint<1 + t_numFacetVertices> {
    int ContactNodeIndices[t_numFacetVertices+1];
    /** Master-facet shape function values. */
    float ShapeValues[t_numFacetVertices];
  };

  /** Edge-edge contact constraints */
  struct EdgeEdgeConstraintItem : public ContactConstraint<4> {
    int SlaveNodeIndices[2], MasterNodeIndices[2];
    float SlaveShapeValues[2], MasterShapeValues[2];
  };

public:
  /** CPU Contact response force calculation. Includes contact search */
  virtual bool ComputeContactResponses(float *p_f, const float uNexts[], const float uCurrs[]) = 0;
  /** @} */

  /**
   * \name Contact Forces
   * @{
   */
protected:
  struct ContactResponse {
    /** Sum of responses (mean direction) */
    float AccumulatedResponse[3];
    /** Maximum force in mean (over all applied contact forces) direction. Used only in consolidation. */
    float MaxProjection;
  };

private:
  std::vector<ContactResponse> m_NormalForces;
  std::vector<ContactResponse> m_FrictionForces;

protected:
  /** Resizes the contact force buffer and resets all values to 0 (use sparingly) */
  void SetContactForceBufferSize(const int maxNumSlaveNodes, const bool doFriction);
  int GetContactForceBufferSize(void) const { return int(m_NormalForces.size()); }

  ContactResponse* GetAllNormalResponses(void) { return &m_NormalForces.front(); }
  ContactResponse& GetNormalResponse(const int nodeIndex) { return m_NormalForces[nodeIndex]; }

  static void ResetNormalResponse(ContactResponse &r_response);
  void ResetNormalResponse(const int nodeIndex) { tledContactSolverCPU::ResetNormalResponse(this->GetNormalResponse(nodeIndex)); }

  ContactResponse* GetAllFrictionResponses(void) { return &m_FrictionForces.front(); }
  ContactResponse& GetFrictionResponse(const int nodeIndex) { return m_FrictionForces[nodeIndex]; }

  static void ResetFrictionResponse(ContactResponse &r_response);
  void ResetFrictionResponse(const int nodeIndex) { tledContactSolverCPU::ResetFrictionResponse(this->GetFrictionResponse(nodeIndex)); }

  /** Effective forces acting target-surface normal direction. For internal use only: only valid after consolidation and before reset of contact data structures. */
  const float* GetEffectiveNodeNormalForces(const int nodeIndex) const { return m_NormalForces[nodeIndex].AccumulatedResponse; }

  /** Effective friction (tangential) forces. Behaviour undefined prior to call to ApplyFrictionForces, reset with all other contact force buffers at end of contact response computation. */
  const float* GetEffectiveNodeFrictionForces(const int nodeIndex) const { return m_FrictionForces[nodeIndex].AccumulatedResponse; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledContactSolverCPU(void) {}
  /** @} */
};

/**
 * \brief Base class for CPU contact solver implementations, provides some members generally useful in solving contact problems
 */
template <class TContactMesh, class TContactSolverInterface>
class tledContactSolverImplCPU : public tledContactSolverImpl<TContactMesh, TContactSolverInterface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSolverImpl<TContactMesh, TContactSolverInterface> Superclass;
  typedef TContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  /** @} */

  /**
   * \name Response Calculation
   * @{
   */
public:
  typedef tledContactSolverCPU::NodeFacetConstraintItem<Facet::NumberOfVertices> NodeFacetConstraintItem;
  typedef tledContactSolverCPU::EdgeEdgeConstraintItem EdgeEdgeConstraintItem;

protected:
  typedef tledContactSolverCPU::ContactResponse ContactResponse;

private:
  std::vector<int> m_ForceResetIndices;

private:
  template <const bool t_doMaster, const bool t_doSlave>
  void _RunMaxResponseConsolidationAndApply(float *p_f, ContactResponse *p_responses, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints);
  void _ApplyResponsesOnly(float *p_f, ContactResponse *p_responses, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints);

  template <const bool t_doMaster, const bool t_doSlave>
  float _ComputeStoppingNodeFacetForce(const float vMag, const NodeFacetConstraintItem &ci) const;

  template <const bool t_doMaster, const bool t_doSlave>
  float _ComputeStoppingEdgeEdgeForce(const float vMag, const EdgeEdgeConstraintItem &ci) const;

  /** \internal Applies a contact (normal) force to the given node of the surface. Automatically marks the node as active. */
  void _ApplyForceToNode(const int nodeIndex, const float f[]);

  /** \internal Applies a friction force to the given node of the surface. */
  void _ApplyFrictionToNode(const int nodeIndex, const float f[]);

protected:
  void ResetContactForces(void);

  /** Returns list of nodes with active contact constraints (note only valid within ComputeContactForces). */
  const int* GetAllActiveContactNodeIndices(void) const { return &m_ForceResetIndices.front(); }
  /** Number of nodes with active contact constraints (note only valid within ComputeContactForces). */
  int GetNumberOfActiveContactNodeIndices(void) const { return int(m_ForceResetIndices.size()); }

  /** Marks an entry in the contact force vector as active */
  void MarkActiveNode(const int nodeIndex) { m_ForceResetIndices.push_back(nodeIndex); }

  template <const int t_numSegNodes>
  float ComputeMass(const int segment[], const float shapeVals[]) const;

  template <const bool t_doMaster, const bool t_doSlave>
  float ComputeSlaveBeta(const EdgeEdgeConstraintItem &ci) const;
  template <const bool t_doMaster, const bool t_doSlave>
  float ComputeSlaveBeta(const NodeFacetConstraintItem &ci) const;

  float ComputeShapeSquareSum(const float shapeValues[], const int numNodes) const;

  template <const bool t_doMaster, const bool t_doSlave>
  void ComputeRelativeNodeFacetVelocity(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;

  /** Computes the relative velocity when only the slave mesh depends on the current displacement (rigid-deformable: use when DoMaster == true) */
  virtual void ComputeRelativeNodeFacetVelocityMaster(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeNodeFacetVelocity<false, true>(p_v, ci, uNexts, uCurrs); }
  /** Computes the relative velocity when only the master mesh depends on the current displacement */
  virtual void ComputeRelativeNodeFacetVelocitySlave(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeNodeFacetVelocity<true, false>(p_v, ci, uNexts, uCurrs); }
  /** Computes the relative velocity when both the master and slave surfaces are deformable */
  virtual void ComputeRelativeNodeFacetVelocityBoth(float *p_v, const NodeFacetConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeNodeFacetVelocity<true, true>(p_v, ci, uNexts, uCurrs); }

  template <const bool t_doMaster, const bool t_doSlave>
  void ComputeRelativeEdgeEdgeVelocity(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const;

  /** Computes the relative velocity when only the slave mesh depends on the current displacement (rigid-deformable: use when DoMaster == true) */
  virtual void ComputeRelativeEdgeEdgeVelocityMaster(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeEdgeEdgeVelocity<false, true>(p_v, ci, uNexts, uCurrs); }
  /** Computes the relative velocity when only the master mesh depends on the current displacement */
  virtual void ComputeRelativeEdgeEdgeVelocitySlave(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeEdgeEdgeVelocity<true, false>(p_v, ci, uNexts, uCurrs); }
  /** Computes the relative velocity when both the master and slave surfaces are deformable */
  virtual void ComputeRelativeEdgeEdgeVelocityBoth(float *p_v, const EdgeEdgeConstraintItem &ci, const float uNexts[], const float uCurrs[]) const { ComputeRelativeEdgeEdgeVelocity<true, true>(p_v, ci, uNexts, uCurrs); }  

  template <const bool t_doMaster, const bool t_doSlave>
  bool ComputeNodeFacetResponse(NodeFacetConstraintItem &r_constraint);
  template <const bool t_doMaster, const bool t_doSlave>
  bool ComputeEdgeEdgeResponse(EdgeEdgeConstraintItem &r_constraint);

  template <const bool t_doMaster, const bool t_doSlave>
  bool ComputeNodeFacetRateResponse(NodeFacetConstraintItem &r_constraint, const float uNexts[], const float uCurrs[]);
  template <const bool t_doMaster, const bool t_doSlave>
  bool ComputeEdgeEdgeRateResponse(EdgeEdgeConstraintItem &r_constraint, const float uNexts[], const float uCurrs[]);

  /** 
   * \brief Computes node-facet contact frictions.
   *
   * 
   * Does rewrite the constraint forces to the corresponding friction force. Also, requires consolidation of normal forces, hence must be called <i>after</i> ApplyContactForces.
   */
  template <const bool t_doMaster, const bool t_doSlave>
  void ComputeNodeFacetFrictionResponse(std::vector<NodeFacetConstraintItem> &r_constraints, const float uNexts[], const float uCurrs[], const float frictionCoefficient);

  /** 
   * \brief Computes edge-edge contact frictions.
   *
   * 
   * Does rewrite the constraint forces to the corresponding friction force. Also, requires consolidation of normal forces, hence must be called <i>after</i> ApplyContactForces.
   */
  template <const bool t_doMaster, const bool t_doSlave>
  void ComputeEdgeEdgeFrictionResponse(std::vector<EdgeEdgeConstraintItem> &r_constraints, const float uNexts[], const float uCurrs[], const float frictionCoefficient);

  /** Performs the contact search, returns true if any contacts (node-facet or edge-edge) were found */
  virtual bool FindCollisions(const float uNexts[], const float uCurrs[]);

  /** Computes the actual forces arising from the found contacts, must return true if there were any non-zero forces were computed. */
  virtual bool ComputeContactForces(float *p_f, const float uNexts[], const float uCurrs[]) = 0;

  /** Applies the contact force to the slave node of a node-facet contact and marks it as active. */
  void ApplyForceToSlaveNode(NodeFacetConstraintItem &r_constraint, const float df[]);

  /** Applies the contact force to a master-facet node of a node-facet contact and marks it as active. */
  void ApplyForceToMasterFacetNode(NodeFacetConstraintItem &r_constraint, const float df[], const int vertexIndex);

  /** Applies the contact force to a master-edge node of an edge-edge contact and marks it as active. */
  void ApplyForceToMasterEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex);

  /** Applies the contact force to a slave-edge node of an edge-edge contact and marks it as active. */
  void ApplyForceToSlaveEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex);

  /** Applies the friction force to the slave node of a node-facet contact. */
  void ApplyFrictionToSlaveNode(NodeFacetConstraintItem &r_constraint, const float df[]);

  /** Applies the friction force to a master-facet node of a node-facet contact. */
  void ApplyFrictionToMasterFacetNode(NodeFacetConstraintItem &r_constraint, const float df[], const int vertexIndex);

  /** Applies the friction force to a master-edge node of an edge-edge contact. */
  void ApplyFrictionToMasterEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex);

  /** Applies the friction force to a slave-edge node of an edge-edge contact. */
  void ApplyFrictionToSlaveEdgeNode(EdgeEdgeConstraintItem &r_constraint, const float df[], const int nodeIndex);

  /** Consolidates all normal contact constraints and applies the consolidated forces to the global internal forces vector, does not reset the forces. */
  template <const bool t_doMaster, const bool t_doSlave>
  void ApplyContactForces(float *p_f, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints);

  /** Consolidates all friction constraints and applies the forces to the global internal forces vector, does not reset the forces. */
  template <const bool t_doMaster, const bool t_doSlave>
  void ApplyFrictionForces(float *p_f, const std::vector<NodeFacetConstraintItem> &nodeFacetConstraints, const std::vector<EdgeEdgeConstraintItem> &edgeEdgeConstraints);

public:
  virtual bool ComputeContactResponses(float *p_f, const float uNexts[], const float uCurrs[]);
  /** @} */

  /**
   * \name Contact Search
   * @{
   */
private:
  tledBVHTraverserCPU *mp_BVHTraverser;
  
protected:
  /** Must instantiate but not initialise the BVH traverser */
  virtual tledBVHTraverserCPU* InstantiateBVHTraverser(void) = 0;
  virtual void InitBVHTraverser(void);

  tledBVHTraverserCPU& GetBVHTraverser(void) { return *mp_BVHTraverser; }
  const tledBVHTraverserCPU& GetBVHTraverser(void) const { return *mp_BVHTraverser; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(void);

  tledContactSolverImplCPU(tledUnstructuredContactManager &r_manager) : Superclass(r_manager), mp_BVHTraverser(NULL) {}
  virtual ~tledContactSolverImplCPU(void);
  /** @} */
};

#include "tledContactSolverCPU.tpp"
#endif

