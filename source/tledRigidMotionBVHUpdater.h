// =========================================================================
// File:       tledRigidMotionBVHUpdater.h
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
#ifndef tledRigidMotionBVHUpdater_H
#define tledRigidMotionBVHUpdater_H

#include "tledDynamicBVHUpdater.h"
#include "tledVectorArithmetic.h"
#include "tledMovingRigidContactSurface.h"

#include <algorithm>
#include <cassert>

/**
 * \brief Rigid-body motion BVH updater, only compatible with tledMovingRigidContactSurface
 * \ingroup contact
 */
template <class TBVH>
class tledRigidMotionBVHUpdater : public tledDynamicBVHUpdaterImpl<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename BVH::ContactMesh ContactMesh;
  typedef tledDynamicBVHUpdaterImpl<TBVH> Superclass;
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
private:
  float m_RootTranslation[3];

protected:
  const float* GetRootTranslation(void) const { return m_RootTranslation; }
  void SetRootTranslation(const float t[]) { std::copy(t, t + 3, m_RootTranslation); }

public:
  virtual void Init(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledRigidMotionBVHUpdater(BVH &r_bvh) : Superclass(r_bvh) { std::fill(m_RootTranslation, m_RootTranslation + 3, 0.0f); }
  tledRigidMotionBVHUpdater(void) { std::fill(m_RootTranslation, m_RootTranslation + 3, 0.0f); }
  virtual ~tledRigidMotionBVHUpdater(void) {}
  /** @} */
};

#include "tledRigidMotionBVHUpdater.tpp"
#endif
