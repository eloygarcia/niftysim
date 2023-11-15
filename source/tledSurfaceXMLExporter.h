// =========================================================================
// File:       tledSurfaceXMLExporter.h
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
#ifndef tledSurfaceXMLExporter_H
#define tledSurfaceXMLExporter_H
#include "tledXMLExporter.h"
#include "tledSurface.h"

/**
 * \brief Converts a surface object to XML
 * \ingroup model
 */
template <class TSurface>
class tledSurfaceXMLExporter : public tledXMLExporter<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef TSurface Surface;
  typedef tledXMLExporter<TSurface> Superclass;
  /** @} */

  /** 
   * \name Processing
   * @{
   */
protected:
  /** Writes node-related data */
  virtual void WriteNodes(void);

  /** Writes facet-related data */
  virtual void WriteFacets(void);

  virtual void WriteBody(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledSurfaceXMLExporter(void) {}
  /** @} */
};

#include "tledSurfaceXMLExporter.tpp"
#endif
