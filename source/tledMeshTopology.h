// =========================================================================
// File:       tledMeshTopology.h
// Purpose:    Mesh topology class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#ifndef tledMeshTopology_H
#define tledMeshTopology_H
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <utility>
#include <vector>

#include "tledMesh.h"
#include "tledHelper.h"
#include "tledVectorArithmetic.h"

/**
 * \brief Provides extended topological information for tledMesh objects
 * \ingroup mesh
 */
template <const int t_numElNodes = 4>
class tledMeshTopology {
public:
  tledMeshTopology(const tledMesh &mesh) : mc_Mesh(mesh) {}
  virtual ~tledMeshTopology(void) {}

  /**
   * \name Mesh 
   * @{
   */
private:
  const tledMesh &mc_Mesh;

public:
  const tledMesh& GetMesh(void) const { return mc_Mesh; }
  /** @} */

  /**
   * \name Facets
   * @{
   */
public:
  /**
   * \brief Facet representation
   */
  class Facet : public tledArray<int, t_numElNodes == 4? 3 : 4> {
  public:
    static const int NumFacetVertices = t_numElNodes == 4? 3 : 4;

  public:
    typedef tledArray<int, NumFacetVertices> Superclass;

  protected:
    Facet(void) {}

  public:
    Facet(const int v0, const int v1, const int v2, const int v3 = -1);
  };

  /**
   * \brief Fixed length arrays of element facets
   */
  typedef tledArray<int, t_numElNodes == 4? 4 : 6> ElementFacetList;  

public:
  const Facet& GetFacet(const int facetInd) const { return m_Facets[facetInd]; }

  int GetNumFacets(void) const { return m_Facets.size(); }

  /**
   * Compiles a list of all facets that are on the surface of the mesh.
   */
  std::vector<int> GetSurfaceFacets(void) const;
  /** @} */
  
  /**
   * \name Computation routines
   * @{
   */
public:
  /**
   * \brief Computes for all nodes lists of adjacent elements.
   */
  void ComputeNodeAdjacency(void);

  /**
   * Computes lists of facets and adjacent elements.
   */
  void ComputeFacets(void);

  /**
   * \brief Computes min./max. element diameters accesible through GetElementMinH() and GetElementMaxH(), respectively.   
   */
  void ComputeElementDiameters(void);
  /** @} */

  /**
   * \name Getters
   * @{
   */
public:
  /**
   * Returns the elements adjacent to a facet. If a facet is part of the mesh surface, the second component is 
   * set to -1
   */
  const std::pair<int,int>& GetFacetElements(const int facetInd) const { return m_FacetEls[facetInd]; }
  
  const ElementFacetList& GetElementFacets(const int elInd) const { return m_ElFacets[elInd]; }

  /**
   * \brief Returns a list of facets adjacent to a node.
   */
  const std::vector<int>& GetNodeFacets(const int nodeInd) const { return m_NodeFacets[nodeInd]; }

  /**
   * \brief Returns a list of elements adjacent to a node.
   */
  const std::vector<int>& GetNodeElements(const int nodeInd) const { return m_NodeEls[nodeInd]; }
  /** @} */

  /**
   * \name Mesh Stats
   * @{
   */
private:
  float m_MinH, m_MaxH;

public:
  /**
   * \brief Returns the shortest element <i>edge</i> (not true element diameter)
   */
  float GetElementMinH(void) const { return m_MinH; }

  /**
   * \brief Returns the longest element <i>edge</i> (not true element diameter)
   */
  float GetElementMaxH(void) const { return m_MaxH; }
  /** @} */
  
private:
  std::vector<std::pair<int,int> > m_FacetEls;
  std::vector<Facet> m_Facets;
  std::vector<std::vector<int> > m_NodeEls;
  std::vector<ElementFacetList> m_ElFacets;
  std::vector<std::vector<int> > m_NodeFacets; 
};


typedef tledMeshTopology<4> tledTetrahedronMeshTopology;
typedef tledMeshTopology<8> tledHexedronMeshTopology;

#include "tledMeshTopology.tpp"
#endif
