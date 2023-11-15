// =========================================================================
// File:       tledContactSurfaceCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactSurfaceCreator_H
#define tledContactSurfaceCreator_H

#include "tledBasicSurfaceCreator.h"
#include "tledSurfaceTopology.h"
#include "tledHelper.h"

#include <algorithm>

template <class TSurface>
class tledContactSurfaceCreator : public tledBasicSurfaceCreator<TSurface> {
  /** 
   * \name Types
   * @{
   */
public:
  typedef tledBasicSurfaceCreator<TSurface> Superclass;
  /** @} */

  /**
   * \name Main 
   * @{
   */
private:
  tledBasicSurfaceCreator<TSurface> *mp_BaseCreator;
  
protected:
  /** Adds edges to the output surface */
  void InitEdges(void);

public:
  virtual void Create();
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactSurfaceCreator(tledBasicSurfaceCreator<TSurface> &r_baseCreator) : mp_BaseCreator(&r_baseCreator) {}
  virtual ~tledContactSurfaceCreator(void) {}
  /** @} */
};

#include "tledContactSurfaceCreator.tpp"
#endif
