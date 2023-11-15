// =========================================================================
// File:       tledContactSurfaceGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    December 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSurfaceGPU_H
#define tledContactSurfaceGPU_H

#include "tledContactSurface.h"
#include "tledCUDAHelpers.h"

class tledContactSurfaceGPU {
  /**
   * \name On-Device Representation
   * @{
   */
public:
  struct GPUSurface {
    float3 *NodeCoordinates;
    int2 *Edges;
    float3 *NodeNormals;
    int NumberOfFacets, NumberOfNodes, NumberOfEdges;
  };

private:
  float3 *mdp_NodeCoordinates;
  int *mdp_FacetVertexIndices, *mdp_FacetEdgeIndices;
  GPUSurface *mhp_Surface, *mdp_Surface;

protected:
  void SetOnDeviceNodeCoordinateBuffer(float3 *dp_X) { mdp_NodeCoordinates = dp_X; }
  float3* GetAllOnDeviceNodeCoordinates(void) { return mdp_NodeCoordinates; }

  void InitNodes(const int numNodes);

  virtual void AllocateDeviceSurface(void) = 0;
  virtual void InitNodes(void) = 0;
  virtual void InitFacets(void) = 0;
  virtual void CopyNodes(void) = 0;
  virtual void CopySurfaceToDevice(void) = 0;
  virtual void InitDeviceSurface(void);
  virtual void InitEdges(void) = 0;  

  GPUSurface*& GetHostGPUSurfacePointer(void) { return mhp_Surface; }
  void SetOnDeviceFacetVertexIndices(const int *inds, const int numInds);
  void SetOnDeviceFacetEdgeIndices(const int *inds, const int numInds);

public:
  GPUSurface& GetHostGPUSurface(void) { return *mhp_Surface; }
  const GPUSurface& GetHostGPUSurface(void) const { return *mhp_Surface; }

  const GPUSurface* GetDeviceGPUSurface(void) const { return mdp_Surface; }  
  GPUSurface*& GetDeviceGPUSurface(void) { return mdp_Surface; }

  /** Access to an array containing only the facet vertex indices */
  const int* GetOnDeviceFacetVertexIndices(void) const { return mdp_FacetVertexIndices; }

  /** Access to an array containing only the facet edge indices */
  const int* GetOnDeviceFacetEdgeIndices(void) const { return mdp_FacetEdgeIndices; }

  const float3* GetAllOnDeviceNodeCoordinates(void) const { return mdp_NodeCoordinates; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactSurfaceGPU(void) : mdp_NodeCoordinates(NULL), mdp_FacetVertexIndices(NULL), mdp_FacetEdgeIndices(NULL), mhp_Surface(NULL), mdp_Surface(NULL) {}
  virtual ~tledContactSurfaceGPU(void);
  /** @} */
};

template <class TBaseSurface>
class tledContactSurfaceImplGPU : public TBaseSurface {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseSurface Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Device Representation
   * @{
   */
public:
  typedef tledContactFacet<Facet::NumberOfVertices> GPUFacet;

  struct GPUSurface : public Superclass::GPUSurface {
    static const int NumberOfFacetVertices = Facet::NumberOfVertices;

    GPUFacet *Facets;
  };

private:
  GPUSurface& _GetHostGPUSurface(void) {
    return static_cast<GPUSurface&>(this->GetHostGPUSurface());
  }

protected:
  /** Copies a standard 3-component node vector to the device. Overwrites NodeCopyBuffer */
  void CopyNodeVectorToDevice(float3 *dp_dst, const float hostV[]);
  void CopyNodeVectorToHost(float *p_dst, const float3 *dpc_gpuV);

  virtual void AllocateDeviceSurface(void);
  virtual void InitNodes(void);
  virtual void InitEdges(void);
  virtual void CopyNodes(void);
  virtual void InitFacets(void);
  virtual void CopySurfaceToDevice(void);
  /** @} */  

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactSurfaceImplGPU(void) {}
  virtual ~tledContactSurfaceImplGPU(void);
  /** @} */
};

#include "tledContactSurfaceGPU.tpp"
#endif
