// =========================================================================
// File:       tledContactSurface.h
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
#ifndef tledContactSurface_H
#define tledContactSurface_H

#include "tledSurface.h"
#include "tledVectorArithmetic.h"
#include "tledHelper.h"

#include <limits>
#include <utility>
#include <algorithm>
#include <cmath>

class tledUnstructuredContactManager;

/**
 * \brief Facets comprising edge information (only used in contact modelling)
 * \ingroup contact 	 
 */
template <const int t_numFacetVertices>
struct tledContactFacet : public tledBasicSurfaceFacet<t_numFacetVertices> {
  int EdgeIndices[t_numFacetVertices];  
};

/**
 * \ingroup Contact surface untemplated interface
 * \ingroup contact 	 
 */
class tledContactSurface : public tledSurface {
  /**
   * \name Edges
   * @{
   */
private:
  std::vector<std::pair<int, int> > m_Edges;
  
public:
  int GetNumberOfEdges(void) const { return m_Edges.size(); }

  /** Returns a reference to the vector of all edges in the surface (R/O) */
  const std::vector<std::pair<int, int> >& GetAllEdges(void) const { return m_Edges; }

  /** Returns a reference to the vector of all edges in the surface (R/W) */
  std::vector<std::pair<int, int> >& GetAllEdges(void) { return m_Edges; }

  const std::pair<int, int>& GetEdge(const int edgeIndex) const { return GetAllEdges()[edgeIndex]; }
  /** @} */

  /**
   * \brief Nodes
   * @{
   */
public:
  /** Allocates node memory */
  virtual void SetNumberOfNodes(const int numNodes) {}
  /** @} */

  /**
   * \name Misc. Geometry Queries
   * @{
   */
protected:
  float m_MaxH, m_MinH;
  
public:
  /** Computes max./min. facet diameters (only possible after initialisation) */
  virtual void ComputeDiameters(void) = 0;
  
  /** Maximum surface element thickness, defaults to safety margin value of manager */
  float GetMaxDiameter(void) const { return m_MaxH; }
  /** Minimum surface element thickness, defaults to safety margin value of manager */
  float GetMinDiameter(void) const { return m_MinH; }
  /** @} */

  /**
   * \name Friction
   * @{
   */
private:
  float m_CoulombMu;

public:
  void SetFrictionCoefficient(const float mu) { m_CoulombMu = mu; }
  float GetFrictionCoefficient(void) const { return m_CoulombMu; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactSurface(void) : m_MaxH(std::numeric_limits<float>::quiet_NaN()), m_MinH(std::numeric_limits<float>::quiet_NaN()), m_CoulombMu(std::numeric_limits<float>::quiet_NaN()) {}
  virtual ~tledContactSurface(void) {}
  /** @} */
};

/**
 * \brief Contact surface base class implementation
 * \ingroup contact 	 
 */
template <const int t_numFacetVertices, class TInterface>
class tledContactSurfaceImpl : public tledSurfaceImpl<tledContactFacet<t_numFacetVertices>, TInterface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactFacet<t_numFacetVertices> Facet;
  typedef tledSurfaceImpl<tledContactFacet<t_numFacetVertices>, TInterface> Superclass;
  /** @} */

  /**
   * \name Surface Nodes
   * @{
   */
 private:
  std::vector<float> m_NodeCoordinates;

 public:
  /** Allocates node memory */
  virtual void SetNumberOfNodes(const int numNodes);

  /** R/W access to nodes */
  float* GetAllNodeCoordinates(void) { return &m_NodeCoordinates.front(); }

  /** R/O access to nodes */
  const float* GetAllNodeCoordinates(void) const { return this->Superclass::GetAllNodeCoordinates(); }
  /** @} */

  /**
   * \name Surface Projection
   * @{
   */
 public:
  /**
   * \brief Computes a reduced matrix operator for projection of a coordinate onto a facet, specified by facetInd   
   *
   * Based on Moller-Trumbore algorithm.<br>
   * Use with ProjectOntoFacet.
   */
  float* ComputeFacetProjectionOperator(float *p_projOp, const int facetInd) const;
  /** @} */

  /**
   * \name Misc. Geometry Queries
   * @{
   */
public:
  virtual void ComputeDiameters(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledContactSurfaceImpl(void) {}

public:
  virtual ~tledContactSurfaceImpl(void) {}  
  /** @} */
};

#include "tledContactSurface.tpp"
#endif
