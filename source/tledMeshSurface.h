// =========================================================================
// File:       tledMeshSurface.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMeshSurface_H
#define tledMeshSurface_H

#include "tledMesh.h"
#include "tledSurface.h"

/**
 * \brief Holds a surface extracted from the sim's solid mesh.
 * \ingroup surface
 */
template <const int t_numFacetVertices>
class tledMeshSurface : public tledSurfaceImpl<tledBasicSurfaceFacet<t_numFacetVertices>, tledSurface> {
  /**
   * \name Types
   * @{
   */
private:
  class _Extractor;

public:
  typedef tledSurfaceImpl<tledBasicSurfaceFacet<t_numFacetVertices>, tledSurface> Superclass;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledMeshSurface(void) {}

public:
  tledMeshSurface(const tledMesh &mesh);

  virtual ~tledMeshSurface(void) {}
  /** @} */
};

template <>
tledMeshSurface<3>::tledMeshSurface(const tledMesh &mesh);

template <>
tledMeshSurface<4>::tledMeshSurface(const tledMesh &mesh);
#endif
