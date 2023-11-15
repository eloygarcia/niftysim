// =========================================================================
// File:       tledOBBXMLExporter.h
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
#ifndef tledOBBXMLExporter_H
#define tledOBBXMLExporter_H

#include "tledBVXMLExporter.h"
#include "tledOBB.h"

/**
 * \brief OBB XML exporter
 * \ingroup contact
 */
template <class TBV>
class tledOBBXMLExporter : public tledBVXMLExporter<TBV> {
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
  virtual const char* GetRootElementName(void) const { return "OBB"; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledOBBXMLExporter(void) : Superclass() {}
  virtual ~tledOBBXMLExporter(void) {}
  /** @} */
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledOBB<2> > : public tledOBBXMLExporter<tledOBB<2> > {
};

template <>
class tledBVHExporterBVXMLExporterAdapter<tledOBB<4> > : public tledOBBXMLExporter<tledOBB<4> > {
};

#include "tledOBBXMLExporter.tpp"
#endif
