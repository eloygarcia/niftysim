// =========================================================================
// File:       tledXMLDeformableMembraneContactSurfaceCreator.h
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
#ifndef tledXMLDeformableMembraneContactSurfaceCreator_H
#define tledXMLDeformableMembraneContactSurfaceCreator_H

#include "tledXMLDeformableContactSurfaceCreator.h"

/** 
 * \brief Parser for files created with tledDeformableMembraneContactSurfaceXMLExporter 
 * \ingroup contact
 */
template <class TSurface>
class tledXMLDeformableMembraneContactSurfaceCreator : public tledXMLDeformableContactSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLDeformableContactSurfaceCreator<TSurface> Superclass;
  /** @} */

  /**
   * \name Parsing
   * @{
   */
protected:
  class XMLImporter;

protected:
  virtual typename tledXMLSurfaceCreator<TSurface>::XMLImporter* InstantiateParser(void) const { return new XMLImporter(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledXMLDeformableMembraneContactSurfaceCreator(void) : Superclass() {}
  virtual ~tledXMLDeformableMembraneContactSurfaceCreator(void) {}
  /** @} */
};

/** 
 * \brief XML parsing back-end 
 */
template <class TSurface>
class tledXMLDeformableMembraneContactSurfaceCreator<TSurface>::XMLImporter : public Superclass::XMLImporter {
  /**
   * \name Import
   * @{
   */
protected:
  virtual void ParseIndexBounds(void);
  virtual void ParseNodeData(void);
  virtual void ParseCoordinates0(void);

public:
  virtual void Import(void);
  /** @} */
};

#include "tledXMLDeformableMembraneContactSurfaceCreator.tpp"
#endif
