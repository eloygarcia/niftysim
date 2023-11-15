// =========================================================================
// File:       tledSurface.h
// Purpose:    General surface mesh representation
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef tledSurface_H
#define tledSurface_H

#include "tledVectorArithmetic.h"
#include "tledHelper.h"

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>

/**
 * \defgroup surface Surface representations
 * \ingroup mesh
 * \brief Classes for representing 2D elements embedded in 3D space.
 */

/**
 * \brief 3-space 2D element
 * \ingroup surface
 * 
 * Node indices refer to the surface the facet belongs to not the mesh from which the surface was derived (if applicable).
 */
template <const int t_numFacetVertices>
struct tledBasicSurfaceFacet {
  int NodeIndices[t_numFacetVertices];  

  static const int NumberOfVertices = t_numFacetVertices;
}; 

/**
 * \brief Basic, untemplated interface for 2D surfaces in 3-space
 * \ingroup surface
 *
 *
 * For performance-critical operations, a template-based approach is preferable over using this interface.
 */
class tledSurface {
  /**
   * \name Facets
   * @{
   */
public:
  virtual void SetNumberOfFacets(const int numFacets) {}
  virtual int GetNumberOfFacets(void) const = 0;
  virtual const int* GetFacetNodeIndices(const int facetIndex) const = 0;
  virtual int GetNumberOfFacetVertices(void) const = 0;
  /** @} */

  /**
   * \name Nodes
   * @{
   */
public:
  virtual int GetNumberOfNodes(void) const = 0;
  virtual const float* GetNodeCoordinates(const int nodeInd) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledSurface(void) {}
  /** @} */
};

/**
 * \brief Basic 2D surfaces in 3-space implementation.
 */
template <class TFacet, class TInterface>
class tledSurfaceImpl : public TInterface {
  /**
   * \name Facets
   * @{
   */
public:
  typedef TFacet Facet;
  
private:
  std::vector<Facet> m_Facets;

public:
  /** Allocates facet related memory */
  virtual void SetNumberOfFacets(const int numFacets);
  
  /** R/W access to facet data structure */
  Facet& GetFacet(const int facetInd) { return m_Facets[facetInd]; }

  virtual int GetNumberOfFacets(void) const { return m_Facets.size(); }
  virtual int GetNumberOfFacetVertices(void) const { return Facet::NumberOfVertices; }  
  const std::vector<Facet>& GetAllFacets(void) const { return m_Facets; }
  const Facet& GetFacet(const int facetInd) const { return m_Facets[facetInd]; }

  /**
   * \brief Computes a (piece-wise constant) covariant basis for a facet given by a node array and vertex indices.
   *
   * The basis is given by 
   * \f$\frac{dx(\xi, \eta)}{d\xi},\;\frac{dx(\xi, \eta)}{d\eta},\;\frac{dx(\xi, \eta)}{d\xi}\times\;\frac{dx(\xi, \eta)}{d\eta}\f$
   * The facet normal is not normalised, for triangular facets its magnitude is 2 x the facet area.
   * For quadrilateral facets the basis at the facet mid-point is returned (\f$\xi = \eta = 0.5\f$).
   */
  static void ComputeFacetBasis(float *p_dxDXi, float *p_dxDEta, float *p_n, const int facetVertexIndices[], const float nodes[]);

  /**
   * \brief Computes an unnormalised facet normal for an explicitly given facet. See ComputeFacetBasis for details.
   */
  static float* ComputeFacetNormal(float *p_normalDst, const int facetVertexIndices[], const float nodes[]);

  /**
   * \brief Computes an unnormalised facet normal for a facet given through an index.
   */
  float* ComputeFacetNormal(float *p_normalDst, const int facetInd) const { return ComputeFacetNormal(p_normalDst, GetFacet(facetInd).NodeIndices, GetAllNodeCoordinates()); }

  /**
   * \brief Computes the full covariant facet basis for a given facet.
   */
  void ComputeFacetBasis(float *p_dxDXi, float *p_dxDEta, float *p_n, const int facetInd) { return ComputeFacetBasis(p_dxDXi, p_dxDEta, p_n, GetFacet(facetInd).NodeIndices, GetAllNodeCoordinates()); }

  /** \brief Facet area, if possible avoid using this and compute area from a precomputed facet normal. */
  float ComputeFacetArea(const int facetInd) const;

  /**
   * \brief Computes a facet's centroid
   */
  float* ComputeCentroid(float *p_centroid, const int facetInd) const;

  /**
   * \brief Computes normalised facet normal
   */
  float* ComputeNormalisedFacetNormal(float *p_normalDst, const int facetInd) const;  

  /** Returns the list of nodes comprised in a list facets. */
  template <typename TConstIterator>
  std::vector<int> CompileNodeListFromFacetList(const TConstIterator ic_pIndsBegin, const TConstIterator ic_pIndsEnd) const;

  virtual const int* GetFacetNodeIndices(const int facetIndex) const { return GetFacet(facetIndex).NodeIndices; }
  /** @} */

  /**
   * \name Surface Nodes
   * @{
   */
private:
  const float *mpc_NodeCoordinates;
  int m_NumNodes;

protected:
  void SetNodeVector(const float nodes[], const int numNodes);

public:
  virtual int GetNumberOfNodes(void) const { return m_NumNodes; }
  virtual const float* GetNodeCoordinates(const int nodeInd) const { return &mpc_NodeCoordinates[3*nodeInd]; }
  const float* GetAllNodeCoordinates(void) const { return mpc_NodeCoordinates; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledSurfaceImpl(void) {}
  /** @} */
};

#include "tledSurface.tpp"

#endif
