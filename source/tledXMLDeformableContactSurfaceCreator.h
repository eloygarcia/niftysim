// =========================================================================
// File:       tledXMLDeformableContactSurfaceCreator.h
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
#ifndef tledXMLDeformableContactSurfaceCreator_H
#define tledXMLDeformableContactSurfaceCreator_H

#include "tledXMLContactSurfaceCreator.h"

/** 
 * \brief Parser for files created with tledDeformableContactSurfaceXMLExporter 
 * \ingroup contact
 */
template <class TSurface>
class tledXMLDeformableContactSurfaceCreator : public tledXMLContactSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLContactSurfaceCreator<TSurface> Superclass;
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
  tledXMLDeformableContactSurfaceCreator(void) : Superclass() {}
  virtual ~tledXMLDeformableContactSurfaceCreator(void) {}
  /** @} */
};

/** \brief XML parsing back-end */
template <class TSurface>
class tledXMLDeformableContactSurfaceCreator<TSurface>::XMLImporter : public Superclass::XMLImporter {
  /**
   * \name Import
   * @{
   */
protected:
  virtual void ParseMapsData(void);
  virtual void ParseNodeData(void);

public:
  virtual void Import(void);
  /** @} */
};

#include "tledXMLDeformableContactSurfaceCreator.tpp"
#endif
