// =========================================================================
// File:       tledBasicSurfaceCreator.h
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
#ifndef tledBasicSurfaceCreator_H
#define tledBasicSurfaceCreator_H

#include "tledBasicMeshCreator.h"

/**
 * \brief Creates a surface through some means other than reading from a file.
 * \ingroup surface
 */
template <class TSurface>
class tledBasicSurfaceCreator : public tledBasicMeshCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBasicMeshCreator<TSurface> Superclass;
  /** @} */

  /**
   * \name Main 
   * @{
   */
public:
  virtual void Create(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBasicSurfaceCreator(void) {}
  /** @} */
};
#endif
