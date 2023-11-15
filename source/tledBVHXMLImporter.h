// =========================================================================
// File:       tledBVHXMLImporter.h
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
#ifndef tledBVHXMLImporter_H
#define tledBVHXMLImporter_H

#include "tledXMLImporter.h"
#include "tledAABBXMLImporter.h"
#include "tledAABBXMLExporter.h"
#include "tledOBBXMLImporter.h"
#include "tledOBBXMLExporter.h"

template <class TBVH>
class tledBVHXMLImporter : public tledXMLImporter<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLImporter<TBVH> Superclass;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  /** Reads information pertaining to BVH leafs */
  virtual void ParseLeafData(void);

  /** Reads the global parameters of the BVH */
  virtual void ParseBaseStatsData(void);
  
  /** Reads the BV data */
  virtual void ParseBVData(void);

public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBVHXMLImporter(void) : Superclass() {}
  virtual ~tledBVHXMLImporter(void) {}
  /** @} */
};

#include "tledBVHXMLImporter.tpp"
#endif
