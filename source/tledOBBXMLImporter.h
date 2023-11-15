// =========================================================================
// File:       tledOBBXMLImporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledOBBXMLImporter_H
#define tledOBBXMLImporter_H
#include "tledBVXMLImporter.h"

#include <algorithm>

/**
 * \brief XML reader for OBBs
 * \ingroup contact
 * \ingroup fileio
 */
template <class TBV>
class tledOBBXMLImporter : public tledBVXMLImporter<TBV> {
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
  tledOBBXMLImporter(void) : Superclass() {}
  virtual ~tledOBBXMLImporter(void) {}
  /** @} */
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledOBB<2> > : public tledOBBXMLImporter<tledOBB<2> > {
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledOBB<4> > : public tledOBBXMLImporter<tledOBB<4> > {
};

#include "tledOBBXMLImporter.tpp"
#endif
