// =========================================================================
// File:       tledDeformableMembraneContactSurfaceXMLExporter.h
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
#ifndef tledDeformableMembraneContactSurfaceXMLExporter_H
#define tledDeformableMembraneContactSurfaceXMLExporter_H
#include "tledDeformableContactSurfaceXMLExporter.h"

/**
 * \brief Exports a deformable contact surface comprising a membrane. Assumes that the surface has never been updated or saved. 
 * \ingroup contact
 */
template <class TContactSurface>
class tledDeformableMembraneContactSurfaceXMLExporter : public tledDeformableContactSurfaceXMLExporter<TContactSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableContactSurfaceXMLExporter<TContactSurface> Superclass;
  /** @} */

  /** 
   * \name Processing
   * @{
   */
protected:
  /** Writes membrane base indices (see tledDeformableMembraneContactSurface) */
  virtual void WriteIndexBounds(void);
  virtual void WriteNodes0(void);
  virtual void WriteBody(void);

public:
  static const char* GetDeformableMembraneContactSurfaceXMLTag(void) { return "DeformableMembraneContactSurface"; }
  virtual const char* GetRootElementName(void) const { return GetDeformableMembraneContactSurfaceXMLTag(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledDeformableMembraneContactSurfaceXMLExporter(void) {}
  /** @} */
};

#include "tledDeformableMembraneContactSurfaceXMLExporter.tpp"
#endif
