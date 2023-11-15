// =========================================================================
// File:       tledDynamicBVHGPU.h
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
#ifndef tledDynamicBVHGPU_H
#define tledDynamicBVHGPU_H
#ifdef _GPU_
#include "tledDynamicBVH.h"
#include "tledBVHGPU.h"

template <class TBV>
class tledDynamicGPUBV : public TBV {
public:
  int UpdateCounter;


};

/**
 * \name Interface of GPU dynamic BVHs (mainly for moving rigid bodies)
 * \ingroup contact
 */
class tledDynamicBVHGPU : public tledDynamicBVH, public tledBVHGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDynamicBVHGPU* CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin);

  /** XML import */
  static tledDynamicBVHGPU* CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root);

  virtual ~tledDynamicBVHGPU(void) {}
  /** @} */  
};

template <class TSurface, class TBV, class TAPI = tledDynamicBVHGPU>
class tledDynamicBVHImplGPU : public tledBVHImplGPU<tledDynamicBVHImpl<TSurface, TBV, TAPI> > {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHImplGPU<tledDynamicBVHImpl<TSurface, TBV, TAPI> > Superclass;
  typedef TBV BoundingVolume;
  typedef TSurface ContactMesh;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDynamicBVHImplGPU(const ContactMesh &mesh) : Superclass(const_cast<ContactMesh&>(mesh)) {}
  virtual ~tledDynamicBVHImplGPU(void) {}
  /** @} */
};

#endif
#endif
