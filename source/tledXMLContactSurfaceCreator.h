// =========================================================================
// File:       tledXMLContactSurfaceCreator.h
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
#ifndef tledXMLContactSurfaceCreator_H
#define tledXMLContactSurfaceCreator_H

#include "tledXMLSurfaceCreator.h"

#include <algorithm>

/**
 * \brief XML importer class for parsing of XML representations of contact surfaces generated with NiftySim exporters (not user definitions!)
 * \ingroup contact
 */
template <class TSurface>
class tledXMLContactSurfaceCreator : public tledXMLSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLSurfaceCreator<TSurface> Superclass;
  /** @} */

  /**
   * \name Parsing
   * @{
   */
protected:
  class XMLImporter;

protected:
  virtual typename Superclass::XMLImporter* InstantiateParser(void) const { return new XMLImporter(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledXMLContactSurfaceCreator(void) : Superclass() {}
  virtual ~tledXMLContactSurfaceCreator(void) {}
  /** @} */
};

/** \brief XML parsing back-end */
template <class TSurface>
class tledXMLContactSurfaceCreator<TSurface>::XMLImporter : public Superclass::XMLImporter {
  /**
   * \name Processing
   * @{
   */
protected:
  virtual void ParseFacetData(void);
  virtual void ParseEdgeData(void);

public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~XMLImporter(void) {}
  /** @} */
};

#include "tledXMLContactSurfaceCreator.tpp"
#endif
