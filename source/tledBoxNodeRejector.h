// =========================================================================
// File:       tledBoxNodeRejector.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBoxNodeRejector_H
#define tledBoxNodeRejector_H

#include "tledNodeRejector.h"

/**
 * Axis-aligned box boundary node rejector.
 * \ingroup model
 */
class tledBoxNodeRejector : public tledNodeRejector {
  /**
   * \name Geometry
   * @{
   */
private:
  float m_Bounds[3][2];
  
public:
  void SetMinBound(const int axis, const float bound) { m_Bounds[axis][0] = bound; }
  void SetMaxBound(const int axis, const float bound) { m_Bounds[axis][1] = bound; }

  float GetMinBound(const int axis) const { return m_Bounds[axis][0]; }
  float GetMaxBound(const int axis) const { return m_Bounds[axis][1]; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
public:
  virtual bool DoReject(const int nodeIndex) const;
  /** @} */

  /**
   * \name Construction, Init, Destruction
   * @{
   */
public:
  virtual void InitialiseFromXMLSpec(const XMLNode rootNode);
  virtual ~tledBoxNodeRejector(void) {}
  /** @} */
};

#endif
