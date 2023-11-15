// =========================================================================
// File:       tledSelfCollisionBVHXMLImporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSelfCollisionBVHXMLImporter_H
#define tledSelfCollisionBVHXMLImporter_H

#include "tledBVHXMLImporter.h"
#include "tledSelfCollisionBVXMLImporter.h"
#include "tledSelfCollisionBVHXMLExporter.h"
#include "tledAABBXMLImporter.h"

/**
 * \brief Parser for XML created with tledSelfCollisionBVHXMLExporter
 * \ingroup contact
 */
template <class TBVH>
class tledSelfCollisionBVHXMLImporter : public tledBVHXMLImporter<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHXMLImporter<TBVH> Superclass;
  typedef TBVH BVH;
  typedef typename BVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  /** Reads information concerning primitive adjacency */
  virtual void ParseAdjacencyData(void);

  /** Reads self-collision candidate lists */
  virtual void ParseCollisionCandidates(void);

public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSelfCollisionBVHXMLImporter(void) : Superclass() {}
  virtual ~tledSelfCollisionBVHXMLImporter(void) {}
  /** @} */
};

#include "tledSelfCollisionBVHXMLImporter.tpp"
#endif
