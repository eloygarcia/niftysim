// =========================================================================
// File:       tledXMLSurfaceCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledXMLSurfaceCreator_H
#define tledXMLSurfaceCreator_H

#include "xmlParser.h"
#include "tledBasicSurfaceCreator.h"
#include "tledXMLImporter.h"
#include "tledModel.h"

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <sstream>

/**
 * \brief Reads a surface from XML
 * \ingroup import
 *
 *
 * Format follows that of solid mesh spec:
\verbatim
<Nodes NumNodes="XXX">
   Node1_X Node1_Y Node1_Z
   ...
   NodeXXX_X NodeXXX_Y NodeXXX_Z
</Nodes>
<Elements NumEls="YYY" Type="T3|Q4">
  Element1_V0 Element1_V1 Element1_V2 ...
  ....
  ElementYYY_V0 ElementYYY_V1 ElementYYY_V2 ...
</Elements>
\endverbatim
 * Supported element types are:
 * <ul>
 * <li>T3: Triangles</li>
 * <li>Q4: Quadrilaterals</li>
 * </ul>
 */
template <class TSurface>
class tledXMLSurfaceCreator : public tledBasicSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBasicSurfaceCreator<TSurface> Superclass;  
  /** @} */

  /**
   * \name Main 
   * @{
   */
protected:
  class XMLImporter;

private:
  const XMLNode *pc_RootNode;

protected:
  virtual XMLImporter* InstantiateParser(void) const { return new XMLImporter(); }

public:
  /** Expects the root node of XML tree defining the surface */
  void SetXMLRoot(const XMLNode &surfaceRootNode) { pc_RootNode = &surfaceRootNode; }

  virtual void Create(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledXMLSurfaceCreator(void) : pc_RootNode(NULL) {}
  virtual ~tledXMLSurfaceCreator(void) {}
  /** @} */
};

/** \brief XML parsing back-end */
template <class TSurface>
class tledXMLSurfaceCreator<TSurface>::XMLImporter : public tledXMLImporter<TSurface> {
  /**
   * \name Processing
   * @{
   */
protected:
  virtual void ParseNodeData(void);
  virtual void ParseFacetData(void);

public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  XMLImporter(void) : tledXMLImporter<TSurface>() {}
  virtual ~XMLImporter(void) {}
  /** @} */
};

#include "tledXMLSurfaceCreator.tpp"
#endif
