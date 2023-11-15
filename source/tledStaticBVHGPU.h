// =========================================================================
// File:       tledStaticBVHGPU.h
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
#ifndef tledStaticBVHGPU_H
#define tledStaticBVHGPU_H

#include "tledBVHGPU.h"
#include "tledStaticBVH.h"
#include "tledRigidContactSurfaceGPU.h"
#include "xmlParser.h"

class tledStaticBVHGPU : public tledBVHGPU, public tledStaticBVH {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Only suitable for instantiation of rigid surface BVHs */
  static tledStaticBVHGPU* CreateBVH(const tledRigidContactSurfaceGPU &mesh, const float margin);

  /** Only suitable for instantiation of rigid surface BVHs from previous XML export */
  static tledStaticBVHGPU* CreateBVH(const tledRigidContactSurfaceGPU &mesh, const XMLNode root);

  virtual ~tledStaticBVHGPU(void) {}
  /** @} */
};

#endif
