// =========================================================================
// File:       tledContactSurfaceXMLExporter.h
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
#ifndef tledContactSurfaceXMLExporter_H
#define tledContactSurfaceXMLExporter_H

#include "tledSurfaceXMLExporter.h"

/**
 * \brief Converts a tledContactSurface to XML, not to be confused with the file writers derived from tledBasicMeshFileWriter
 */
template <class TContactSurface>
class tledContactSurfaceXMLExporter : public tledSurfaceXMLExporter<TContactSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSurfaceXMLExporter<TContactSurface> Superclass;
  typedef TContactSurface Surface;
  /** @} */

  /** 
   * \name Processing
   * @{
   */
protected:
  /** Writes edge-related data */
  virtual void WriteEdges(void);

  virtual void WriteFacets(void);
  virtual void WriteBody(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledContactSurfaceXMLExporter(void) {}
  /** @} */
};

#include "tledContactSurfaceXMLExporter.tpp"
#endif
