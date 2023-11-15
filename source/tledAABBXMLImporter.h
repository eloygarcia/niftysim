// =========================================================================
// File:       tledAABBXMLImporter.h
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
#ifndef tledAABBXMLImporter_H
#define tledAABBXMLImporter_H

#include "tledBVXMLImporter.h"

#include <algorithm>

/**
 * \brief XML reader for AABBs
 * \ingroup contact
 * \ingroup fileio
 */
template <class TBV>
class tledAABBXMLImporter : public tledBVXMLImporter<TBV> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVXMLImporter<TBV> Superclass;
  typedef TBV BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledAABBXMLImporter(void) : Superclass() {}
  virtual ~tledAABBXMLImporter(void) {}
  /** @} */
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledAABB<2> > : public tledAABBXMLImporter<tledAABB<2> > {
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledAABB<4> > : public tledAABBXMLImporter<tledAABB<4> > {
};

#include "tledAABBXMLImporter.tpp"
#endif
