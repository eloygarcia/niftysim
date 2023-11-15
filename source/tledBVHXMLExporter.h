// =========================================================================
// File:       tledBVHXMLExporter.h
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
#ifndef tledBVHXMLExporter_H
#define tledBVHXMLExporter_H

#include "tledXMLExporter.h"
#include "tledBVXMLExporter.h"
#include "tledAABB.h"
#include "tledAABBXMLExporter.h"
#include "tledOBB.h"
#include "tledOBBXMLExporter.h"

#include <iostream>
#include <cstdlib>
#include <typeinfo>

/**
 * \brief Base class for BVH XML exporters
 * \ingroup contact
 */
template <class TBVH>
class tledBVHXMLExporter : public tledXMLExporter<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLExporter<TBVH> Superclass;
  typedef TBVH BVH;
  typedef typename BVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  /** Writes global BVH stats (margins, etc) to the XML tree. */
  virtual void WriteBVHBaseStats(void);
  
  /** Writes the data associated with leaf BVs, primitives */
  virtual void WriteLeafData(void);

  /** Writes the actual BVs */
  virtual void WriteBVs(void);

  virtual void WriteBody(void);

public:
  virtual const char* GetRootElementName(void) const { return "BVH"; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHXMLExporter(void) {}
  /** @} */
};

#include "tledBVHXMLExporter.tpp"
#endif
