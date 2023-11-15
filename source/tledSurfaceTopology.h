// =========================================================================
// File:       tledSurfaceTopology.h
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
#ifndef tledSurfaceTopology_H
#define tledSurfaceTopology_H

#include "tledSurface.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include <cassert>

/**
 * \brief Extended surface-mesh topology information
 * \ingroup surface
 */
template <class TSurface>
class tledSurfaceTopology {
  /**
   * \name Types
   * @{
   */
public:
  typedef TSurface Surface;
  typedef typename Surface::Facet Facet;
  /** @} */

  /**
   * \name Surface Access
   * @{
   */
private:
  const Surface *mpc_Surface;

public:
  const Surface& GetSurface(void) const { return *mpc_Surface; }
  /** @} */

  /**
   * \name Edges
   * @{
   */
private:
  class _EdgeComparator;

private:
  std::vector<std::pair<int, int> > m_Edges;
  std::vector<std::pair<int, int> > m_EdgeNeighbours;
  std::vector<tledArray<int, Facet::NumberOfVertices> > m_FacetEdges;

public:
  const std::vector<std::pair<int, int> >& GetEdges(void) const { return m_Edges; }

  /** 
   * Returns the max. two facets adjacent to an edge (meshes with holes leading to configurations with more than two facets adjacent to an edge are not supported.).<br>
   * Boundary edges can be identified by the -1 in the second component of the neighbour tuple.
   */
  const std::vector<std::pair<int, int> >& GetEdgeNeighbours(void) const { return m_EdgeNeighbours; }

  /** List of global indices of facet edges */ 
  const std::vector<tledArray<int, Facet::NumberOfVertices> >& GetFacetEdges(void) const { return m_FacetEdges; }

  /** Computes edge and adjacent-facet lists */
  void ComputeEdges(void);
  /** @} */
  
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSurfaceTopology(const Surface &surface) : mpc_Surface(&surface) {}
  /** @} */
};

#include "tledSurfaceTopology.tpp"
#endif
