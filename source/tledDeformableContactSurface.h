// =========================================================================
// File:       tledDeformableContactSurface.h
// Purpose:    
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
#ifndef tledDeformableContactSurface_H
#define tledDeformableContactSurface_H

#include "tledVectorArithmetic.h"
#include "tledMesh.h"
#include "tledContactSurface.h"
#include "tledSurfaceTopology.h"
#include "tledContactVolumeSurfaceExtractor.h"
#include "tledDeformableContactSurfaceXMLExporter.h"

#include <limits>
#include <algorithm>
#include <cassert>
#include <string>

/**
 * \brief Deformable contact surface interface.
 * \ingroup contact
 */
class tledDeformableContactSurface : public tledContactSurface {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurface Superclass;
  /** @} */

  /**
   * \name Masses
   * @{
   */
private:
  std::vector<float> m_SurfaceNodeMasses;

public:
  /** Returns the vector of surface node masses */
  const std::vector<float>& GetAllSurfaceNodeMasses(void) const { return m_SurfaceNodeMasses; }
  /** Returns the vector of surface node masses */
  float GetSurfaceNodeMass(const int nodeIndex) const { return m_SurfaceNodeMasses[nodeIndex]; }
  /** Initialise the nodal masses used in contact force calculation */
  virtual void InitNodeMasses(const float globalMasses[]);
  /** @} */

  /**
   * \name Pre-Computation/Geometry Update Routines
   * @{
   */
private:
  int m_UpdateCounter;
  int m_SaveCounter;

protected:
  /** Increments the update counter used in lazy evaluation of various quanities. */
  void IncUpdateCounter(void) { m_UpdateCounter += 1; }

public:
  /** Number of updates since initialisation */
  int GetUpdateCount(void) const { return m_UpdateCounter; }

  /** Number of saves since initialisation */
  int GetSaveCount(void) const { return m_SaveCounter; }

  void ResetSaveCount(void) { m_SaveCounter = 0; }
  void ResetUpdateCount(void) { m_UpdateCounter = 0; }

  /**
   * \brief Stores the surface after it has been ensured there are no intersections.
   */
  virtual void Save(void) { m_SaveCounter += 1; }
  /** @} */

  /**
   * \name Surface - Volume Mappings
   * @{
   */
private:
  std::vector<int> m_Volume2SurfaceNodeMap;
  std::vector<int> m_Surface2VolumeNodeMap;
  int m_CoordinateHistorySize;
  
public:
  void SetVolume2SurfaceNodeMap(const std::vector<int> &map) { m_Volume2SurfaceNodeMap = map; }
  void SetSurface2VolumeNodeMap(const std::vector<int> &map) { m_Surface2VolumeNodeMap = map; }

  const std::vector<int>& GetVolume2SurfaceNodeMap(void) const { return m_Volume2SurfaceNodeMap; }
  const std::vector<int>& GetSurface2VolumeNodeMap(void) const { return m_Surface2VolumeNodeMap; }

  /** @{ */
  /**
   * \brief Controls how many time steps to go back for "old" node positions (default 10).
   */  
  int GetCoordinateHistorySize(void) const { return m_CoordinateHistorySize; }
  void SetCoordinateHistorySize(const int numTSteps) { m_CoordinateHistorySize = numTSteps; }
  /** @} */

  /**
   * \brief Maps a surface node index to the corresponding volume mesh node index.
   */
  int MapSurface2VolumeNode(const int surfNodeInd) const { return m_Surface2VolumeNodeMap[surfNodeInd]; }

  /**
   * \brief Maps a volume mesh node index to a surface node one, provided the node is contained in the surface mesh.
   *
   * Returns -1 if there is no corresponding surface node.
   */
  int MapVolume2SurfaceNode(const int volumeNodeInd) const { return m_Volume2SurfaceNodeMap[volumeNodeInd]; }
  /** @} */

  /**
   * \name XML Export
   * @{
   */
public:
  virtual XMLNode ExportToXML(void) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Called at the end of loading from an XML file. */
  virtual void LoadFromXMLPostloadHook(void) {}

  /** Finalises the surface definition after all facets and nodes have been defined */
  virtual void Init(void) { m_SaveCounter = m_UpdateCounter = 0; }

  static tledDeformableContactSurface* CreateSurface(const tledMesh &mesh, const bool useGPU);
  static tledDeformableContactSurface* CreateSurface(const XMLNode xmlRep, const bool useGPU);

  /** Creates an empty surface of the type represented by the type string (T3, Q4, etc.), mainly used in testing */
  static tledDeformableContactSurface* CreateSurface(const std::string &type, const bool useGPU);

  tledDeformableContactSurface(void) : m_UpdateCounter(0), m_SaveCounter(0), m_CoordinateHistorySize(10) {}
  virtual ~tledDeformableContactSurface(void) {}
  /** @} */
};

/**
 * \brief Deformable contact surface implementation
 * \ingroup contact
 * 
 *
 * Caches various quantities, e.g. normals. Keeps track of node positions via the *OldNodeCoordinates API, Save.
 */
template <const int t_numFacetVertices, class TAPI = tledDeformableContactSurface>
class tledDeformableContactSurfaceImpl : public tledContactSurfaceImpl<t_numFacetVertices, TAPI> {
  /**
   * \name Types
   * @{
   */
protected:
  struct CachedQuantity {
    int UpdateCounter;

    CachedQuantity(void) : UpdateCounter(-1) {}
  };

public:
  typedef tledContactSurfaceImpl<t_numFacetVertices, TAPI> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Normals
   * @{
   */
public:
  struct CachedNormal : public CachedQuantity {
    float Normal[3];
  };

protected:
  struct NodeNormalData {
    CachedNormal Normal;
    int NodeFacetsStartIndex, NodeFacetsEndIndex;
  };

private:
  std::vector<NodeNormalData> m_NodeNormals;

protected:
  CachedNormal& GetNodeNormalCache(const int nodeIndex); 
  const CachedNormal& GetNodeNormalCache(const int nodeIndex) const;

  NodeNormalData& GetNodeNormalData(const int nodeIndex) { return m_NodeNormals[nodeIndex]; }
  const NodeNormalData& GetNodeNormalData(const int nodeIndex) const { return m_NodeNormals[nodeIndex]; }
  /** @} */

  /**
   * \name Nodes
   * @{
   */
private:
  std::vector<int> m_NodeFacetInds;
  std::vector<float> m_NodeCoordinates0;

protected:
  const int* GetAllNodeFacetIndices(void) const { return &m_NodeFacetInds.front(); }
  int GetNumberOfAllNodeFacetIndices(void) const { return m_NodeFacetInds.size(); }
  std::vector<int>& GetAllNodeFacetIndices(void) { return m_NodeFacetInds; }

public:
  virtual void SetNumberOfNodes(const int numNodes);

  /** All initial configuration node positions */
  const float* GetAllNodeCoordinates0(void) const { return &m_NodeCoordinates0.front(); }

  /** Initial configuration node positions */
  const float* GetNodeCoordinates0(const int nInd) const { return this->GetAllNodeCoordinates0() + 3*nInd; }

  /** Initial configuration node positions (R/W) */
  std::vector<float>& GetAllNodeCoordinates0(void) { return m_NodeCoordinates0; }
  
  /** Number of facets adjacent to node */
  int GetNumberOfNodeFacets(const int nInd) const { return m_NodeNormals[nInd].NodeFacetsEndIndex - m_NodeNormals[nInd].NodeFacetsStartIndex; }

  /** Start of list of facets adjacent to node (needed for node normal computation) */
  const int* GetNodeFacetIndices(const int nInd) const { return &m_NodeFacetInds.front() + m_NodeNormals[nInd].NodeFacetsStartIndex; }  

  /** Setter for node facet-fan lists */
  void SetNodeFacetIndices(const std::vector<int> &allNodeFacetIndices, const std::vector<int> &numNodeFacets);
  /** @} */

  /**
   * \name Pre-Computation/Geometry Update Routines
   * @{
   */
public:
  virtual void Init(void);
  /** @} */

  /**
   * \name XML Export
   * @{
   */
public:
  virtual XMLNode ExportToXML(void) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  void ConstructFromSolidMesh(const tledMesh &volMesh);

public:
  tledDeformableContactSurfaceImpl(void) {}
  virtual ~tledDeformableContactSurfaceImpl(void) {}
  /** @} */
};

#include "tledDeformableContactSurface.tpp"
#endif
