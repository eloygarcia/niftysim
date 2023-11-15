// =========================================================================
// File:       tledAABBXMLExporter.h
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
#ifndef tledAABBXMLExporter_H
#define tledAABBXMLExporter_H

#include "tledBVXMLExporter.h"
#include "tledAABB.h"

/**
 * \brief AABB XML exporter
 * \ingroup contact
 */
template <class TBV>
class tledAABBXMLExporter : public tledBVXMLExporter<TBV> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVXMLExporter<TBV> Superclass;
  typedef TBV BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  virtual void WriteBody(void);

public:
  virtual const char* GetRootElementName(void) const { return "AABB"; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledAABBXMLExporter(void) : Superclass() {}
  virtual ~tledAABBXMLExporter(void) {}
  /** @} */
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledAABB<2> > : public tledAABBXMLExporter<tledAABB<2> > {
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledAABB<4> > : public tledAABBXMLExporter<tledAABB<4> > {
};

#include "tledAABBXMLExporter.tpp"
#endif
