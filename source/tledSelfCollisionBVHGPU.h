// =========================================================================
// File:       tledSelfCollisionBVHGPU.h
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
#ifndef tledSelfCollisionBVHGPU_H
#define tledSelfCollisionBVHGPU_H

#include "tledSelfCollisionBVH.h"
#include "tledBVHGPU.h"
#include "xmlParser.h"

class tledSelfCollisionBVHGPU : public tledSelfCollisionBVH, public tledBVHGPU {
  /**
   * \name Construction, Destruction
   * @{
   */
private:
  template <class TBV>
  static tledSelfCollisionBVHGPU* _CreateBVHForBVType(tledDeformableContactSurface &r_mesh, const float margin, const float maxDisplacement);

public:
  static tledSelfCollisionBVHGPU* CreateBVH(tledDeformableContactSurface &r_mesh, const std::string &bvType, const float margin, const float maxDisplacement);

  /** Rebuild from previous export */
  static tledSelfCollisionBVHGPU* CreateBVH(tledDeformableContactSurface &r_mesh, const XMLNode root);

  virtual ~tledSelfCollisionBVHGPU(void) {}
  /** @} */
};

template <class TSurface, class TBV, class TAPI = tledSelfCollisionBVH>
class tledSelfCollisionBVHImplGPU : public tledDynamicBVHImpl<TSurface, tledSelfCollisionBV<TBV>, TAPI> {
  
};

#endif
