// =========================================================================
// File:       tledDeformableContactSurfaceGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableContactSurfaceGPU_H
#define tledDeformableContactSurfaceGPU_H

#ifdef _GPU_
#include "tledDeformableContactSurface.h"
#include "tledDynamicContactSurfaceGPU.h"
#include "tledMesh.h"
#include "tledCUDAHelpers.h"

/**
 * \brief Deformable contact surface (GPU) interface.
 * \ingroup contact
 */
class tledDeformableContactSurfaceGPU : public tledDeformableContactSurface, public tledContactSurfaceGPU {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurface Superclass;
  /** @} */

  /**
   * \name Geometry Update Routines
   * @{
   */
protected:
  static void UpdateNodePositions(float3 *dp_x, const float3 *dpc_x0, const float4 *dpc_uNext, const int *dpc_surfaceToVolumeNodeIndexMap, const int numNodes);

public:
  /** Updates the surface current configuration to accurately reflect existing contacts/mesh intersections. */
  virtual void Update(const float4 *dpc_U) { tledDeformableContactSurface::IncUpdateCounter(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDeformableContactSurfaceGPU(void) {}

public:
  static tledDeformableContactSurfaceGPU* CreateSurface(const tledMesh &mesh);
  static tledDeformableContactSurfaceGPU* CreateSurface(const XMLNode xml);
  static tledDeformableContactSurfaceGPU* CreateSurface(const std::string &type);

  virtual ~tledDeformableContactSurfaceGPU(void) {}
  /** @} */
};

/**
 * \brief GPU deformable contact-surface implementation.
 * \ingroup contact
 */
template <class TBaseSurface>
class tledDeformableContactSurfaceImplGPU : public tledDynamicContactSurfaceGPU<tledContactSurfaceImplGPU<TBaseSurface> > {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDynamicContactSurfaceGPU<tledContactSurfaceImplGPU<TBaseSurface> > Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name GPU
   * @{
   */
public:
  struct GPUSurface : public Superclass::GPUSurface {
    int2 *NodeFacetIndexRanges;
    int *NodeFacetList;
    float *NodeMasses;
    int *VolumeToSurfaceNodeIndexMap;
    int *SurfaceToVolumeNodeIndexMap;
  };

private:
  float3 *mdp_OldCoordinateBuffer;

private:
  GPUSurface& _GetHostSurface(void) { return static_cast<GPUSurface&>(this->GetHostGPUSurface()); }
  GPUSurface* _GetDeviceSurface(void) { return static_cast<GPUSurface*>(this->GetDeviceGPUSurface()); }

protected:
  virtual void AllocateDeviceSurface(void);
  virtual void CopySurfaceToDevice(void);
  virtual void InitNodes(void);  
  virtual void InitFacets(void);
  virtual void CopyNodes(void);  
  virtual int GetCoordinateHistoryBufferMultiplier(void) const { return this->GetCoordinateHistorySize(); }
  /** @} */

  /**
   * \name CPU Compatibility
   * @{
   */
private:
  std::vector<float> m_FacetNormals;

public:
  /** Only an alias for GetAllNodeCoordinates, needed for BVH construction. */
  const float* GetAllOldNodeCoordinates(void) const { return this->GetAllNodeCoordinates(); }

  /** No actual caching, always recomputes the normal. Used in BVH construction */
  const float* GetFacetNormalCached(const int facetIndex) { return this->ComputeFacetNormal(&m_FacetNormals[3*facetIndex], facetIndex); }

  /** Used in BVH construction */
  float* ComputeNormalisedFacetNormalCached(float *p_dst, const int facetIndex) { return tledVectorArithmetic::ScalarDiv(p_dst, tledVectorArithmetic::Norm(this->ComputeFacetNormal(p_dst, facetIndex))); }
  /** @} */

  /**
   * \name Pre-Computation/Geometry Update Routines
   * @{
   */
public:  
  /** Update the vertex normals of a list of nodes (assumes the list is already unique) */
  virtual void ComputeNodeNormals(const int *dpc_nodeIndices, const int numNodes);

  /** Updates all unnormalised facet normals and normalised node normals. */
  virtual void UpdateAllNormals(void);

  virtual void InitNodeMasses(const float globalMasses[]);
  virtual void Init(void);
  virtual void Update(const float4 *dp_us);  
  virtual void Save(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void LoadFromXMLPostloadHook(void);

  tledDeformableContactSurfaceImplGPU(const tledMesh &mesh) { this->ConstructFromSolidMesh(mesh); }
  tledDeformableContactSurfaceImplGPU(void) {}

  virtual ~tledDeformableContactSurfaceImplGPU(void);
  /** @} */
};

typedef tledDeformableContactSurfaceImplGPU<tledDeformableContactSurfaceImpl<3, tledDeformableContactSurfaceGPU>  > tledDeformableContactSurfaceT3GPU;
typedef tledDeformableContactSurfaceImplGPU<tledDeformableContactSurfaceImpl<4, tledDeformableContactSurfaceGPU>  > tledDeformableContactSurfaceQ4GPU;

#ifdef __CUDACC__
#include "tledDeformableContactSurfaceGPU.tpp"
#endif

#endif
#endif
