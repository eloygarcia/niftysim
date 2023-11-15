// =========================================================================
// File:       tledSelfCollisionBVXMLExporter.h
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
#ifndef tledSelfCollisionBVXMLExporter_H
#define tledSelfCollisionBVXMLExporter_H

#include "tledBVXMLExporter.h"
#include "tledSelfCollisionBV.h"
#include "tledAABBXMLExporter.h"

/**
 * \brief Self-collision BV XML exporter
 * \ingroup contact
 */
template <class TBaseExporter>
class tledSelfCollisionBVXMLExporter : public TBaseExporter {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseExporter Superclass;
  typedef typename Superclass::BoundingVolume BoundingVolume;  
  /** @} */

  /** 
   * \name Processing
   * @{
   */
protected:
  virtual void WriteBody(void);

public:
  virtual const char* GetRootElementName(void) const { return "SelfCollisionBV"; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledSelfCollisionBVXMLExporter(void) {}
  /** @} */
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledSelfCollisionBV<tledAABB<2> > > : public tledSelfCollisionBVXMLExporter<tledAABBXMLExporter<tledSelfCollisionBV<tledAABB<2> > > > {
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledSelfCollisionBV<tledAABB<4> > > : public tledSelfCollisionBVXMLExporter<tledAABBXMLExporter<tledSelfCollisionBV<tledAABB<4> > > > {
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledSelfCollisionBV<tledOBB<2> > > : public tledSelfCollisionBVXMLExporter<tledOBBXMLExporter<tledSelfCollisionBV<tledOBB<2> > > > {
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledSelfCollisionBV<tledOBB<4> > > : public tledSelfCollisionBVXMLExporter<tledOBBXMLExporter<tledSelfCollisionBV<tledOBB<4> > > > {
};

#include "tledSelfCollisionBVXMLExporter.tpp"
#endif
