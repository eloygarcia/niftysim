// =========================================================================
// File:       tledBVXMLExporter.h
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
#ifndef tledBVXMLExporter_H
#define tledBVXMLExporter_H

#include "tledXMLExporter.h"
#include "tledBV.h"

/**
 * \brief Base class for bounding volume XML exporters
 * \ingroup contact
 */
template <class TBV>
class tledBVXMLExporter : public tledXMLExporter<TBV> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLExporter<TBV> Superclass;
  typedef TBV BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  virtual void WriteBody(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVXMLExporter(void) {}
  /** @} */
};

/**
 * \brief Dummy class used in BVH exporters to get the right exporter for a given BV type.
 */
template <class TBV>
class tledBVHExporterBVXMLExporterAdapter : public tledBVXMLExporter<TBV> {
  /**
   * \name Processing
   * @{
   */
protected:
  virtual void WriteBody(void) { std::abort(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHExporterBVXMLExporterAdapter(void) {}
  /** @} */  
};

#include "tledBVXMLExporter.tpp"
#endif
