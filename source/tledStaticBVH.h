// =========================================================================
// File:       tledStaticBVH.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledStaticBVH_H
#define tledStaticBVH_H

#include "tledBVH.h"
#include "tledRigidContactSurface.h"
#include "xmlParser.h"

class tledStaticBVH : public tledBVH {
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Factory for rigid contact-surface BVHs */
  static tledStaticBVH* CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin, const bool useGPU);

  /** XML import factory */
  static tledStaticBVH* CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root, const bool useGPU);

  virtual ~tledStaticBVH(void) {}
  /** @} */
};

#endif
