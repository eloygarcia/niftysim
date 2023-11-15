// =========================================================================
// File:       tledSelfCollisionBVXMLImporter.h
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
#ifndef tledSelfCollisionBVXMLImporter_H
#define tledSelfCollisionBVXMLImporter_H

#include "tledBVXMLImporter.h"
#include "tledSelfCollisionBV.h"
#include "tledAABBXMLImporter.h"
#include "tledOBBXMLImporter.h"

template <class TBaseImporter>
class tledSelfCollisionBVXMLImporter : public TBaseImporter {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseImporter Superclass;
  typedef typename Superclass::BoundingVolume BoundingVolume;
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
  tledSelfCollisionBVXMLImporter(void) : Superclass() {}
  virtual ~tledSelfCollisionBVXMLImporter(void) {}
  /** @} */
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledSelfCollisionBV<tledAABB<2> > > : public tledSelfCollisionBVXMLImporter<tledAABBXMLImporter<tledSelfCollisionBV<tledAABB<2> > > > {
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledSelfCollisionBV<tledAABB<4> > > : public tledSelfCollisionBVXMLImporter<tledAABBXMLImporter<tledSelfCollisionBV<tledAABB<4> > > > {
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledSelfCollisionBV<tledOBB<2> > > : public tledSelfCollisionBVXMLImporter<tledOBBXMLImporter<tledSelfCollisionBV<tledOBB<2> > > > {
};

template <>
class tledBVHImporterBVXMLImporterAdapter<tledSelfCollisionBV<tledOBB<4> > > : public tledSelfCollisionBVXMLImporter<tledOBBXMLImporter<tledSelfCollisionBV<tledOBB<4> > > > {
};

#include "tledSelfCollisionBVXMLImporter.tpp"
#endif
