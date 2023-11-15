// =========================================================================
// File:       tledSphereNodeRejector.h
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
#ifndef tledSphereNodeRejector_H
#define tledSphereNodeRejector_H

#include "tledNodeRejector.h"

/**
 * Sphere boundary node rejector.
 * \ingroup model
 */
class tledSphereNodeRejector : public tledNodeRejector {
  /**
   * \name Geometry
   * @{
   */
private:
  float m_Centre[3], m_Radius;
  
public:
  void SetCentre(const float cnt[]) { std::copy(cnt, cnt + 3, m_Centre); }
  void SetRadius(const float r) { m_Radius = r; }

  const float* GetCentre(void) const { return m_Centre; }
  float GetRadius(void) const { return m_Radius; }
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
  virtual ~tledSphereNodeRejector(void) {}
  /** @} */
};

#endif
