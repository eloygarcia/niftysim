// =========================================================================
// File:       tledSelfCollisionBVHXMLExporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHXMLExporter_H
#define tledSelfCollisionBVHXMLExporter_H

#include "tledBVHXMLExporter.h"
#include "tledSelfCollisionBVXMLExporter.h"
#include "tledAABBXMLExporter.h"
#include "tledOBBXMLExporter.h"

#include <iostream>
#include <cstdlib>
#include <typeinfo>

/**
 * \brief Self-collision BVH XML exporters
 * \ingroup contact
 */
template <class TBVH>
class tledSelfCollisionBVHXMLExporter : public tledBVHXMLExporter<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHXMLExporter<TBVH> Superclass;
  typedef TBVH BVH;
  typedef typename BVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  /** Writes information concerning primitive adjacency */
  virtual void WriteAdjacencyData(void);

  /** Writes self-collision candidate lists */
  virtual void WriteCollisionCandidates(void);

  virtual void WriteBody(void);

public:
  virtual const char* GetRootElementName(void) const { return "SelfCollisionBVH"; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledSelfCollisionBVHXMLExporter(void) {}
  /** @} */
};

#include "tledSelfCollisionBVHXMLExporter.tpp"
#endif
