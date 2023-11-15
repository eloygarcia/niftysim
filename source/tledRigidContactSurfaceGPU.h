// =========================================================================
// File:       tledRigidContactSurfaceGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledRigidContactSurfaceGPU_H
#define tledRigidContactSurfaceGPU_H

#include "tledRigidContactSurface.h"
#include "tledContactSurfaceGPU.h"
#include "tledCUDAMemoryBlock.h"

class tledRigidContactSurfaceGPU : public tledRigidContactSurface, public tledContactSurfaceGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledRigidContactSurfaceGPU* CreateSurface(const XMLNode &meshSpec);
  virtual ~tledRigidContactSurfaceGPU(void) {}
  /** @} */
};

template <class TBaseSurface>
class tledRigidContactSurfaceImplGPU : public tledContactSurfaceImplGPU<TBaseSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceImplGPU<TBaseSurface> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name GPU
   * @{
   */
public:
  struct GPUSurface : public Superclass::GPUSurface {
    float3 *FacetNormals;
    float4 *FacetProjectionOperators;
  };

private:
  GPUSurface& _GetHostSurface(void) { return static_cast<GPUSurface&>(this->GetHostGPUSurface()); }

protected:
  virtual void AllocateDeviceSurface(void);
  virtual void CopySurfaceToDevice(void);
  virtual void InitNodes(void);  
  virtual void InitFacets(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  virtual void InitNormals(void) {}

public:
  virtual void Init(void);

  /** Preferably instantiated through tledRigidContactSurface::CreateSurface */
  tledRigidContactSurfaceImplGPU(void) {}

  virtual ~tledRigidContactSurfaceImplGPU(void) {}
  /** @} */
};

/**
 * \brief Shorthand for rigid triangle surface meshes.
 * \ingroup contact
 */
typedef tledRigidContactSurfaceImplGPU<tledRigidContactSurfaceImpl<3, tledRigidContactSurfaceGPU> > tledRigidContactSurfaceT3GPU;

/**
 * \brief Shorthand for rigid quad surface meshes.
 * \ingroup contact
 */
typedef tledRigidContactSurfaceImplGPU<tledRigidContactSurfaceImpl<4, tledRigidContactSurfaceGPU> > tledRigidContactSurfaceQ4GPU;

#include "tledRigidContactSurfaceGPU.tpp"
#endif
