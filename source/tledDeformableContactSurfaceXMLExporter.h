// =========================================================================
// File:       tledDeformableContactSurfaceXMLExporter.h
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
#ifndef tledDeformableContactSurfaceXMLExporter_H
#define tledDeformableContactSurfaceXMLExporter_H

#include "tledContactSurfaceXMLExporter.h"

/**
 * \name Exports a deformable contact surface. Assumes that the surface has never been updated or saved. 
 * \ingroup contact
 */
template <class TContactSurface>
class tledDeformableContactSurfaceXMLExporter : public tledContactSurfaceXMLExporter<TContactSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceXMLExporter<TContactSurface> Superclass;
  typedef TContactSurface Surface;
  /** @} */

  /** 
   * \name Processing
   * @{
   */
protected:
  /** Writes edge-related data */
  virtual void WriteMaps(void);

  virtual void WriteNodes(void);
  virtual void WriteBody(void);

public:
  static const char* GetDeformableContactSurfaceXMLTag(void) { return "DeformableContactSurface"; }
  virtual const char* GetRootElementName(void) const { return GetDeformableContactSurfaceXMLTag(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledDeformableContactSurfaceXMLExporter(void) {}
  /** @} */
};

#include "tledDeformableContactSurfaceXMLExporter.tpp"
#endif
