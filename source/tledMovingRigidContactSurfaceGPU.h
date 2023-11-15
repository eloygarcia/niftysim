// =========================================================================
// File:       tledMovingRigidContactSurfaceGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMovingRigidContactSurfaceGPU_H
#define tledMovingRigidContactSurfaceGPU_H

#include "tledMovingRigidContactSurface.h"
#include "tledRigidContactSurfaceGPU.h"
#include "tledDynamicContactSurfaceGPU.h"

/**
 * \brief Moving rigid contact surface GPU interface
 * \ingroup contact
 */
class tledMovingRigidContactSurfaceGPU : public tledMovingRigidContactSurface, public tledRigidContactSurfaceGPU {
  /**
   * \name Transformation Utilities
   * @{
   */
protected:
  static void TranslateNodes(float3 *dp_cds, float3 *dp_oldCds, const float3 *dpc_x0, const int numNodes, const float currentT[], const float oldT[]);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledMovingRigidContactSurfaceGPU(void) {}

public:
  static tledMovingRigidContactSurfaceGPU* CreateSurface(const std::string &type);
  virtual ~tledMovingRigidContactSurfaceGPU(void) {}
  /** @} */
};

/**
 * \brief Moving rigid surface GPU implementation.
 * \ingroup contact
 */
template <class TBaseSurface> 
class tledMovingRigidContactSurfaceImplGPU : public tledDynamicContactSurfaceGPU<tledRigidContactSurfaceImplGPU<TBaseSurface> > {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledRigidContactSurfaceImplGPU<TBaseSurface> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Device Memory Management
   */
protected:
  virtual int GetCoordinateHistoryBufferMultiplier(void) const { return 1; }
  /** @} */

  /**
   * \name Transform
   * @{
   */
protected:
  /** Updates current and historical node coordinates. */
  void TranslateNodes(void);

public:
  virtual void Update(void);
  /** @} */

  /**
   * \name CPU Compatibility Interface (mainly for BVH construction)
   * @{
   */
public:
  const float* GetAllOldNodeCoordinates(void) const { return this->GetAllNodeCoordinates(); }  
  const float* GetAllOldNodeCoordinates0(void) const { return this->GetAllNodeCoordinates(); }  
  /** @} */

  /**
   * \name Construction, Destruction
   */
public:
  /** Should be instantiated through tledRigidContactSurface::CreateSurface */
  tledMovingRigidContactSurfaceImplGPU(void) {}

  virtual ~tledMovingRigidContactSurfaceImplGPU(void) {}
  /** @} */
};

/**
 * \brief Shorthand for moving rigid triangle surface meshes.
 * \ingroup contact
 */
typedef tledMovingRigidContactSurfaceImplGPU<tledMovingRigidContactSurfaceImpl<3, tledMovingRigidContactSurfaceGPU> > tledMovingRigidContactSurfaceT3GPU;

/**
 * \brief Shorthand for moving rigid quad surface meshes.
 * \ingroup contact
 */
typedef tledMovingRigidContactSurfaceImplGPU<tledMovingRigidContactSurfaceImpl<4, tledMovingRigidContactSurfaceGPU> > tledMovingRigidContactSurfaceQ4GPU;

#include "tledMovingRigidContactSurfaceGPU.tpp"
#endif
