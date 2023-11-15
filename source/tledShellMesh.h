// =========================================================================
// File:       tledShellMesh.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellMesh_H
#define tledShellMesh_H

#include "tledModel.h"
#include "tledMeshSurface.h"

/**
 * \defgroup shell Shell/Membrane Modelling
 */

/**
 * \brief Shell Element Mesh
 * \ingroup shell
 */
template <const int t_numFacetVertices> 
class tledShellMesh : public tledMeshSurface<t_numFacetVertices> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledMeshSurface<t_numFacetVertices> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Surface Nodes
   * @{
   */
public:  
  /** Setter for global node vector */
  void SetNodeVector(const float nodes[], const int numNodes) { Superclass::SetNodeVector(nodes, numNodes); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledShellMesh(void) {}

  tledShellMesh(void) {}

  /** Builds mesh by parsing XML specs. */
  tledShellMesh(XMLNode xModel);

  /** Extracts surface from a solid mesh. */
  tledShellMesh(const tledMesh &solidMesh) : Superclass(solidMesh) {}
  /** @} */
};
#include "tledShellMesh.tpp"
#endif
