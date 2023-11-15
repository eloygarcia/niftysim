// =========================================================================
// File:       tledGeometricSelfCollisionBVHUpdaterDevelGPU.h
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
#ifndef tledGeometricSelfCollisionBVHUpdaterDevelGPU_H
#define tledGeometricSelfCollisionBVHUpdaterDevelGPU_H

#include "tledGeometricSelfCollisionBVHUpdaterCPU.h"
#include "tledCUDAHelpers.h"

template <class TBVH>
class tledGeometricSelfCollisionBVHUpdaterDevelGPU : public tledGeometricSelfCollisionBVHUpdaterCPU<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef tledGeometricSelfCollisionBVHUpdaterCPU<TBVH> Superclass;
  typedef typename Superclass::UpdateNodeInfo UpdateNodeInfo;
  /** @} */

  /**
   * \name Update nodes
   * @{
   */
public:
  struct GPUUpdateNodeInfo;  

private:
  int *mdp_UpdateNodeNodeIndices;
  GPUUpdateNodeInfo *mdp_UpdateNodeData;

protected:
  virtual void ConvertToGPU(GPUUpdateNodeInfo &r_gpuUN, const typename Superclass::UpdateNodeInfo &hostUN) const;
  virtual void ConvertFromGPU(typename Superclass::UpdateNodeInfo &r_hostUN, const GPUUpdateNodeInfo &gpuUN) const;

  virtual void CopyUpdateNodesToGPU(void);
  virtual void CopyUpdateNodesFromGPU(void);    

  GPUUpdateNodeInfo*& GetOnDeviceUpdateNodes(void) { return mdp_UpdateNodeData; }
  const GPUUpdateNodeInfo* GetOnDeviceUpdateNodes(void) const { return mdp_UpdateNodeData; }
  const int* GetOnDeviceUpdateNodeIndices(void) const { return mdp_UpdateNodeNodeIndices; }

  void ComputeMaxNonTranslationalDisplacement(void);
  void ComputeMaxNonRigidDisplacement(void);
  void SetOnDeviceUpdateStatuses(void);

  void ComputeUpdateNodeCentroids(void);
  /** @} */

  /**
   * \name Development Members
   * @{
   */
private:
  void _SetUpdateStatusesCPU(void);
  /** @} */

  /**
   * \name Update API
   * @{
   */
public:
  virtual void ComputeNodeUpdateStatuses(void);
  virtual void UpdateBVH(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  virtual void AllocateOnDeviceUpdateNodes(void);

public:
  virtual void Init(void);

  tledGeometricSelfCollisionBVHUpdaterDevelGPU(BVH &r_bvh) : Superclass(r_bvh), mdp_UpdateNodeNodeIndices(NULL) {}
  tledGeometricSelfCollisionBVHUpdaterDevelGPU(void) : mdp_UpdateNodeNodeIndices(NULL) {}
  virtual ~tledGeometricSelfCollisionBVHUpdaterDevelGPU(void);
  /** @} */
};

template <class TBVH>
struct tledGeometricSelfCollisionBVHUpdaterDevelGPU<TBVH>::GPUUpdateNodeInfo {
  float3 Rotation[3], Translation;
  float3 CurrentCentroid, LastUpdateCentroid;
  float LastUpdateConeAngle, MaxNonRigidDisplacement, MaxNonTranslational;
  int BVIndex;
  int NodeStartIndex, NodeEndIndex;
  int Status, RigidUpdateCounter;
};

#endif
