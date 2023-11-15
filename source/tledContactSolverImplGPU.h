// =========================================================================
// File:       tledContactSolverImplGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSolverImplGPU_H
#define tledContactSolverImplGPU_H

#include "tledContactSolverGPU.h"
#include "tledBVHTraverserGPU.h"

/**
 * GPU contact handler implementation base class
 * \ingroup contact
 */
template <class TContactMesh, class TAPI>
class tledContactSolverImplGPU : public tledContactSolverImpl<TContactMesh, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSolverImpl<TContactMesh, TAPI> Superclass;
  typedef TContactMesh ContactMesh;
  typedef typename TContactMesh::GPUSurface GPUSurface;
  typedef typename ContactMesh::Facet Facet;
  typedef tledBVHTraverserGPU::NodeFacetNarrowPhaseResult<ContactMesh::Facet::NumberOfVertices> NodeProjection;
  typedef tledBVHTraverserGPU::EdgeEdgeNarrowPhaseResult EdgeProjection;
  /** @} */

  /**
   * \name Mesh
   * @{
   */
protected:
  GPUSurface& GetHostGPUSurface(void) { return static_cast<GPUSurface&>(this->GetMesh().GetHostGPUSurface()); }
  const GPUSurface& GetHostGPUSurface(void) const { return static_cast<const GPUSurface&>(this->GetMesh().GetHostGPUSurface()); }
  /** @} */

  /**
   * \name Contact Search
   * @{
   */
private:
  tledBVHTraverserGPU *mp_BVHTraverser;

protected:
  /** This class takes care of instantiation through this pure-virtual member function, but not initialisation! */
  virtual tledBVHTraverserGPU* InstantiateBVHTraverser(void) = 0;
  virtual void InitBVHTraverser(void);

  /** Executes the collision detection pipeline. Returns true if the latter finds any contacts. */
  virtual bool FindCollisions(void);

  tledBVHTraverserGPU& GetBVHTraverser(void) { return *mp_BVHTraverser; }
  const tledBVHTraverserGPU& GetBVHTraverser(void) const { return *mp_BVHTraverser; }
  /** @} */  

  /**
   * \name Contact Response Calculation
   * @{
   */
public:
  typedef tledContactSolverGPU::ContactResponse ContactResponse;
  class ProjectionToResponseTransform;
  class NodeProjectionToResponseTransform;
  class EdgeProjectionToResponseTransform;

private:
  tledCUDADeviceMemoryBlock *mp_NodeVelocityBuffer;

protected:
  /** Function for node-facet response computation */
  virtual void ComputeNodeFacetResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t nodeFacetStream) = 0;
  /** Hook called after initial node-facet responses are computed and stream is synchronised, descendant must release all intermediate result memory in it. */
  virtual void ReleaseNodeFacetResponseComputationMemory(void) = 0;
  /** Function for edge-edge response computation */
  virtual void ComputeEdgeEdgeResponses(ContactResponse *dp_dst, const float4 *dpc_uNexts, const float4 *dpc_uCurrs, const cudaStream_t edgeEdgeStream) = 0;
  /** Hook called after initial edge-edge responses are computed and stream is synchronised, descendant must release all intermediate result memory in it. */
  virtual void ReleaseEdgeEdgeResponseComputationMemory(void) = 0;
  /** Sub-classes must provide this function that returns the number of un-consolidated node-facet contact responses for the current BVH traverser state */
  virtual int GetNumberOfNodeFacetResponses(void) const = 0;
  /** Sub-classes must provide this function that returns the number of un-consolidated edge-edge contact responses for the current BVH traverser state */
  virtual int GetNumberOfEdgeEdgeResponses(void) const = 0;

  /** Access to node velocities for rate response computation (device ptr, no division by dt). */
  const float3* GetNodeVelocities(void) const { return mp_NodeVelocityBuffer->GetBuffer<float3>(); }
  
  /** Releases the node-velocity buffer (normally no call by descendant required) */
  virtual void ReleaseNodeVelocities(void);

  /** Computes all (deformable) mesh-node velocities (normally no call by descendant required). Use GetNodeVelocities for subsequent queries. */
  virtual void ComputeNodeVelocities(const float4 *dpc_uNexts, const float4 *dpc_uCurrs);

  static void ConsolidateResponses(tledCUDADeviceMemoryBlock* &rp_consolidatedResponsesn, int &r_numConsolidated, tledCUDADeviceMemoryBlock &r_responses, const int numResponses);  

  /** Consolidates per-constraint responses and adds them to the global external force vector */
  void ConsolidateResponses(float4 *p_f, tledCUDADeviceMemoryBlock &r_responses, const int numResponses) const;

public:
  virtual bool ComputeContactResponses(float4 *dp_R, const float4 *dpc_uNexts, const float4 *dpc_uCurrs);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(void);

  tledContactSolverImplGPU(tledUnstructuredContactManager &r_manager) : Superclass(r_manager), mp_BVHTraverser(NULL), mp_NodeVelocityBuffer(NULL) {}
  virtual ~tledContactSolverImplGPU(void);
  /** @} */
};

/** 
 * \brief Base class for projection to response force transform.
 * \ingroup contact
 */
template <class TContactMesh, class TAPI>
class tledContactSolverImplGPU<TContactMesh, TAPI>::ProjectionToResponseTransform {
  /**
   * \name Simulation Parameters
   * @{
   */
private:
  const float mc_Dt, mc_SafetyOffset, mc_RateDist;

protected:
  __device__ float GetTimeStep(void) const { return mc_Dt; }
  __device__ float GetSafetyOffset(void) const { return mc_SafetyOffset; }
  __device__ float GetRateResponseMaxDistance(void) const { return mc_RateDist; }
  /** @} */
  
  /**
   * \name Current-Configuration, Mesh Queries
   * @{
   */
private:
  const float3 *mdpc_Velocities;
  const float *mdpc_NodeMasses;

protected:
  __device__ float GetNodeMass(const int nInd) const { return mdpc_NodeMasses[nInd]; }
  __device__ const float3& GetNodeVelocity(const int nInd) const { return mdpc_Velocities[nInd]; }
  /** @} */

  /**
   * \name Response Computation
   * @{
   */
protected:
  __device__ float ComputeRateResponseWeight(const float effectiveDepth) const;
  __device__ float ComputeBaseLambdaFromDepth(const float effectiveDepth) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  ProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : mc_Dt(dt), mc_SafetyOffset(off), mc_RateDist(rateDelta), mdpc_NodeMasses(dpc_masses), mdpc_Velocities(dpc_v) {}
  virtual ~ProjectionToResponseTransform(void) {}
  /** @} */
};

/** 
 * \brief Base class for node-facet projection to response force transform.
 * \ingroup contact
 */
template <class TContactMesh, class TAPI>
class tledContactSolverImplGPU<TContactMesh, TAPI>::NodeProjectionToResponseTransform : public ProjectionToResponseTransform {

  /**
   * \name Response-Force Computation
   * @{
   */
protected:
  __device__ float ComputeShapeSquareSum(const NodeProjection &proj) const;
  
  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeSlaveBeta(const NodeProjection &proj) const; 

  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeStoppingLambda(const float vMag, const NodeProjection &proj) const;

  /** Unlike ComputeStoppingLambda: Returns 0 if vMag > 0. */
  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeStoppingLambdaClamped(const float vMag, const NodeProjection &proj) const;

  __device__ float ComputeEffectiveDepth(const NodeProjection &proj) const { return proj.GapValue - this->GetSafetyOffset(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  NodeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : ProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta) {}
  virtual ~NodeProjectionToResponseTransform(void) {}
  /** @} */
};

/** 
 * \brief Base class for edge-edge projection to response force transform.
 * \ingroup contact
 */
template <class TContactMesh, class TAPI>
class tledContactSolverImplGPU<TContactMesh, TAPI>::EdgeProjectionToResponseTransform : public ProjectionToResponseTransform {
  /**
   * \name Response-Force Computation
   * @{
   */
protected:
  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float4 ComputeShapeValues(const EdgeProjection &proj) const;

  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeStoppingLambda(const float vMag, const float4 &shapeVals, const EdgeProjection &proj) const;

  /** Unlike ComputeStoppingLambda: Returns 0 if vMag > 0. */
  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeStoppingLambdaClamped(const float vMag, const float4 &shapeVals, const EdgeProjection &proj) const;

  template <const bool t_doSlave, const bool t_doMaster>
  __device__ float ComputeShapeSquareSum(const float4 &shapeVals) const;

  __device__ float ComputeEffectiveDepth(const EdgeProjection &proj) const { return proj.Xi.x - this->GetSafetyOffset(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  EdgeProjectionToResponseTransform(const float *dpc_masses, const float3 *dpc_v, const float dt, const float off, const float rateDelta) : ProjectionToResponseTransform(dpc_masses, dpc_v, dt, off, rateDelta) {}
  /** @} */
};

#ifdef __CUDACC__
#include "tledContactSolverImplGPU.tpp"
#endif

#endif
